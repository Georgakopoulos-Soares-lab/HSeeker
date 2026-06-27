"""
_scoring.py — H-DNA triplex thermodynamic stability scoring.

Uses an O(n²) constrained maximum-subarray search over all (L, R) arm
windows to find the globally optimal arm boundaries. The constraint
2*(al-R)+ll ≤ (R-L)*v enforces that the new spacer (right tail × 2 +
original spacer) does not exceed v times the selected arm length.

No external dependencies (no Biopython, no attrs).
"""

from __future__ import annotations

_REVCOMP = str.maketrans("ACGTacgt", "TGCAtgca")


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

    def _discount(self, n: int) -> float:
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

    def score(self, full_sequence: str, arm_length: int) -> dict | None:
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
            return self._score_impl(full_sequence, arm_length)
        except (ValueError, IndexError):
            return None

    def _score_impl(self, s: str, al: int) -> dict:
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

        # s2 read 3'→5' so position j of s1 pairs with position j of s2
        s2 = s2_rev[::-1]
        ll = len(s3)

        if any(n not in self.nuc for n in s1 + s2):
            raise ValueError("non-ACGT base in arms")

        # Per-position pairing scores and match array for full arm
        scoring_array = ""
        pairing_scores: list[float] = []
        for j in range(al):
            val = (s1[j] == s2[j] == "G") or (s1[j] == s2[j] == "A")
            scoring_array += "1" if val else "0"
            pairing_scores.append(self._pair_score(s1[j], s2[j]))

        # O(n²) search: maximize pairing+stacking score over window [L, R)
        # subject to: new_spacer_length = 2*(al-R)+ll ≤ (R-L)*v
        best_score = -float("inf")
        best_pair: tuple[int, int] | None = None
        for L in range(al):
            for R in range(min(L + self.min_al, al), al + 1):
                if 2 * (al - R) + ll <= (R - L) * self.v:
                    sub = scoring_array[L:R]
                    stacked = self._calc_stacking(sub)
                    cur = round(sum(pairing_scores[L:R]) + sum(stacked), 3)
                    if cur > best_score:
                        best_score = cur
                        best_pair = (L, R)

        if best_pair is None:
            raise ValueError("No valid arm window satisfies the spacer constraint")

        L, R = best_pair
        sub_scoring = scoring_array[L:R]
        stacking_scores = self._calc_stacking(sub_scoring)
        stacking_score = sum(stacking_scores)
        pairing_score = sum(pairing_scores[L:R])
        total_score = round(stacking_score + pairing_score, 3)

        new_spacer = s1[R:] + s3 + s2[R:][::-1]
        putative_triplex = s1[L:R] + "[" + new_spacer + "]" + s2[L:R][::-1]

        return {
            "stacking_score": stacking_score,
            "pairing_score": pairing_score,
            "total_score": total_score,
            "putative_triplex": putative_triplex,
        }


# Module-level singleton
_scorer = _Scorer()


def score_hit(left_arm: str, spacer: str, right_arm: str, arm_length: int) -> dict | None:
    """Public entry point for scoring a single H-DNA hit.

    Parameters
    ----------
    left_arm : str
        Left arm sequence (5'→3').
    spacer : str
        Spacer/loop sequence.  ``"."`` (C-extension placeholder for zero
        spacer) is treated as an empty string.
    right_arm : str
        Right arm sequence (5'→3').
    arm_length : int
        Length of each arm in bases.

    Returns
    -------
    dict or None
        ``{stacking_score, pairing_score, total_score, putative_triplex}``
        or ``None`` if the sequence is not scorable.
    """
    sp = "" if spacer == "." else spacer
    full_seq = left_arm + sp + right_arm
    return _scorer.score(full_seq, arm_length)
