"""
_scoring.py — H-DNA triplex thermodynamic stability scoring.

Uses an O(n²) constrained maximum-subarray search over all (L, R) arm
windows to find the globally optimal arm boundaries; for each L the window
sums are extended one base at a time as R grows. The constraint
2*(al-R)+ll <= (R-L)*v enforces that the new spacer (right tail x 2 +
original spacer) does not exceed v times the selected arm length.
If no window satisfies that constraint, score the original detected arms.

No external dependencies (no Biopython, no attrs).
"""

from __future__ import annotations

_REVCOMP = str.maketrans("ACGTacgt", "TGCAtgca")
_G_BITS = str.maketrans("GACT", "1000")
_A_BITS = str.maketrans("GACT", "0100")


def _bits(seq: str, table: dict) -> int:
    """Bitmask of the positions of one base in an ACGT string."""
    return int(seq.translate(table), 2) if seq else 0


def _popcount(x: int) -> int:
    return bin(x).count("1")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_REVCOMP)[::-1]


class _Scorer:
    """Stateless scoring engine — instantiated once as a module singleton."""

    def __init__(self):
        self.scoring = {"G": 7.6, "A": 3.83}
        self.stacking_score = 5
        self.nuc = {"A", "G", "C", "T"}
        self.d = 2
        self.dv = 2
        self.min_al = 8
        self.v = 1.0  # max spacer-to-arm ratio

    def _stacking_step(self, prev: str, cur: str, tail: int) -> tuple[float, int]:
        """One term of :meth:`_calc_stacking` and the updated mismatch tail."""
        if prev == cur == "1":
            return float(self.stacking_score), -1
        return -self._discount(max(tail, 0)), tail + 1

    def _discount(self, n: int) -> float:
        n = min(n, 12)
        p = self.stacking_score ** (n + 1)
        return p if n + 1 < 2 else round(p / self.d, self.dv)

    def _calc_stacking(self, scoring_array: str) -> list[float]:
        scores = [0.0]
        incorrect_tail = -1
        for j in range(1, len(scoring_array)):
            if scoring_array[j - 1] == scoring_array[j] == "1":
                scores.append(float(self.stacking_score))
                incorrect_tail = -1
            else:
                scores.append(-self._discount(max(incorrect_tail, 0)))
                incorrect_tail += 1
        return scores

    def _pair_score(self, n1: str, n2: str) -> float:
        if (n1 == n2 == "G") or (n1 == n2 == "A"):
            return self.scoring[n1]
        return -(self.scoring["A"] + self.scoring["G"]) / 2.0

    def score(self, full_sequence: str, arm_length: int, *,
              include_stacking: bool = True,
              optimize_boundaries: bool = True) -> dict | None:
        """Score an H-DNA hit and return adjusted triplex information.

        Parameters
        ----------
        full_sequence : str
            The complete motif: left_arm + spacer + right_arm.
        arm_length : int
            Length of each arm in the original hit.

        Returns
        -------
        dict or None
            ``{stacking_score, pairing_score, total_score, putative_triplex}``
            or ``None`` if scoring could not be applied.
        """
        try:
            return self._score_impl(full_sequence, arm_length,
                                    include_stacking=include_stacking,
                                    optimize_boundaries=optimize_boundaries)
        except (ValueError, IndexError):
            return None

    def _arms(self, s: str, al: int) -> tuple[str, str, str, str]:
        """Orient the motif on its purine strand and split it into arms."""
        s = s.upper()
        s1 = s[:al]
        s2_rev = s[-al:]
        s3 = s[al:-al]

        # If left arm is CT-rich, use the reverse-complement strand
        if s1.count("C") + s1.count("T") > s1.count("G") + s1.count("A"):
            s = _reverse_complement(s)
            s1 = s[:al]
            s2_rev = s[-al:]
            s3 = s[al:-al]

        # s2 read 3'->5' so position j of s1 pairs with position j of s2
        s2 = s2_rev[::-1]
        ll = len(s3)

        if not self.nuc.issuperset(s1 + s2):
            raise ValueError("non-ACGT base in arms")
        return s1, s2, s3, ll

    def upper_bound(self, s: str, al: int) -> float | None:
        """An upper bound on ``total_score`` for this motif, in O(al).

        No window's pairing sum exceeds the sum of the positive per-position
        pairing scores, and its stacking sum gains +stacking_score only for
        adjacent matched positions (every other term is a penalty). Returns
        None when the motif cannot be scored.
        """
        try:
            s1, s2, _, _ = self._arms(s, al)
        except (ValueError, IndexError):
            return None
        # Bit j marks a G:G or A:A pair at arm position j; only those pairs
        # score positively, and only adjacent ones add a stacking bonus.
        g_pairs = _bits(s1, _G_BITS) & _bits(s2, _G_BITS)
        a_pairs = _bits(s1, _A_BITS) & _bits(s2, _A_BITS)
        matched = g_pairs | a_pairs
        pairs = (_popcount(g_pairs) * self.scoring["G"]
                 + _popcount(a_pairs) * self.scoring["A"])
        stacks = _popcount(matched & (matched >> 1))
        # Scores are round(x, 3), and round is monotone, so rounding the bound
        # the same way keeps it valid while letting a window that attains it
        # tie exactly. The 1e-6 nudge absorbs float summation-order error.
        return max(round(pairs + stacks * self.stacking_score + 1e-6, 3), 0)

    def _score_impl(self, s: str, al: int, *,
                    include_stacking: bool = True,
                    optimize_boundaries: bool = True) -> dict:
        s1, s2, s3, ll = self._arms(s, al)

        # Per-position pairing scores and match array for full arm
        scoring_array = ""
        pairing_scores: list[float] = []
        for j in range(al):
            val = (s1[j] == s2[j] == "G") or (s1[j] == s2[j] == "A")
            scoring_array += "1" if val else "0"
            pairing_scores.append(self._pair_score(s1[j], s2[j]))

        # O(n^2) search: maximize pairing+stacking score over window [L, R)
        # subject to: new_spacer_length = 2*(al-R)+ll <= (R-L)*v.
        # Never shrink an arm below min_al, unless the detected arm was
        # already shorter. The old min(L + min_al, al) lower bound let late
        # windows shrink to just one base.
        best_score = -float("inf")
        best_pair: tuple[int, int] | None = None
        if optimize_boundaries:
            # For each L, extend [L, R) one base at a time. The running sums
            # add the same terms in the same order as sum(pairing_scores[L:R])
            # and sum(self._calc_stacking(scoring_array[L:R])). Python >= 3.12
            # sum() compensates rounding, so the raw sums can differ in the
            # last bits, but windows are compared after round(..., 3) and the
            # reported scores below are still computed with sum().
            min_window = min(self.min_al, al)
            for L in range(al):
                pair_sum = 0
                stack_sum = 0.0 if include_stacking else 0
                tail = -1
                for R in range(L + 1, al + 1):
                    pair_sum += pairing_scores[R - 1]
                    if include_stacking and R - 1 > L:
                        term, tail = self._stacking_step(
                            scoring_array[R - 2], scoring_array[R - 1], tail)
                        stack_sum += term
                    if R < L + min_window:
                        continue
                    if 2 * (al - R) + ll <= (R - L) * self.v:
                        cur = round(pair_sum + stack_sum, 3)
                        if cur > best_score:
                            best_score = cur
                            best_pair = (L, R)
        else:
            best_pair = (0, al)

        if best_pair is None:
            # The constraint governs boundary adjustment, not whether a
            # detected repeat can be scored at its original boundaries.
            best_pair = (0, al)

        L, R = best_pair
        sub_scoring = scoring_array[L:R]
        stacking_scores = self._calc_stacking(sub_scoring) if include_stacking else []
        stacking_score = sum(stacking_scores)
        pairing_score = sum(pairing_scores[L:R])
        total_score = round(stacking_score + pairing_score, 3)

        new_spacer = s1[R:] + s3 + s2[R:][::-1]
        # The detector returns lowercase arms and full_sequence. Keep the
        # derived motif in the same case in every public scoring path.
        putative_triplex = (
            s1[L:R] + "[" + new_spacer + "]" + s2[L:R][::-1]
        ).lower()

        return {
            "stacking_score": stacking_score,
            "pairing_score": pairing_score,
            "total_score": max(total_score, 0),
            "putative_triplex": putative_triplex,
        }


# Module-level singleton
_scorer = _Scorer()


def score_hit(left_arm: str, spacer: str, right_arm: str, arm_length: int) -> dict | None:
    """Public entry point for scoring a single H-DNA hit.

    Parameters
    ----------
    left_arm : str
        Left arm sequence (5'->3').
    spacer : str
        Spacer/loop sequence.  "." (C-extension placeholder for zero
        spacer) is treated as an empty string.
    right_arm : str
        Right arm sequence (5'->3').
    arm_length : int
        Length of each arm in bases.

    Returns
    -------
    dict or None
        ``{stacking_score, pairing_score, total_score, putative_triplex}``
        or ``None`` if the sequence is not scorable. ``putative_triplex``
        uses lowercase DNA bases, matching detector output.
    """
    sp = "" if spacer == "." else spacer
    full_seq = left_arm + sp + right_arm
    return _scorer.score(full_seq, arm_length)


def score_hit_components(
    left_arm: str, spacer: str, right_arm: str, arm_length: int, *,
    include_stacking: bool = True, optimize_boundaries: bool = True,
) -> dict | None:
    """Score one hit with optional stacking and boundary optimization.

    Intended for component ablations. With both options enabled, this is
    equivalent to :func:`score_hit`. With boundary optimization disabled,
    the entire detected arm is scored, including any mismatch penalties.
    """
    sp = "" if spacer == "." else spacer
    return _scorer.score(left_arm + sp + right_arm, arm_length,
                         include_stacking=include_stacking,
                         optimize_boundaries=optimize_boundaries)
