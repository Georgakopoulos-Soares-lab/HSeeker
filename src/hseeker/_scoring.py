"""
_scoring.py — H-DNA triplex thermodynamic stability scoring.

Based on the Farmer scoring algorithm.  Evaluates stacking and pairing
energetics to produce an adjusted triplex with optimal arm/spacer boundaries.

No external dependencies (no Biopython, no attrs).
"""

_REVCOMP = str.maketrans("ACGTacgt", "TGCAtgca")


def _reverse_complement(seq: str) -> str:
    """Reverse-complement a DNA string (no Biopython needed)."""
    return seq.translate(_REVCOMP)[::-1]


class _Scorer:
    """Stateless scoring engine — instantiated once as a module singleton."""

    def __init__(self):
        self.scoring = {"G": 7.6, "A": 3.83}
        self.stacking_score = 5
        self.nuc = {"A", "G", "C", "T"}
        self.RT = 0.6
        # Default parameters (same as Farmer.defaults)
        self.d = 2
        self.dv = 2
        self.MS = 50
        self.m = 2
        self.min_al = 8
        self.min_spacer = 4

    # ------------------------------------------------------------------
    #  Core helpers (identical logic to Farmer)
    # ------------------------------------------------------------------

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

    def _trim_ends(self, s1: str, s2: str) -> tuple[str, str]:
        cutoff = 0
        for i in range(len(s1)):
            if (s1[i] == s2[i] == "G") or (s1[i] == s2[i] == "A"):
                break
            cutoff += 1
        return s1[cutoff:], s2[cutoff:]

    # ------------------------------------------------------------------
    #  Main entry point
    # ------------------------------------------------------------------

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
            or ``None`` if scoring could not be applied (e.g. non-ACGT bases).
        """
        try:
            return self._score_impl(full_sequence, arm_length)
        except (ValueError, IndexError):
            return None

    def _score_impl(self, s: str, al: int) -> dict:
        s = s.upper()
        s1 = s[:al]          # left arm
        s2_rev = s[-al:]     # right arm (still in 5'→3' orientation)
        s3 = s[al:-al]       # spacer

        # If left arm is CT-rich, reverse-complement the whole sequence
        if s1.count("C") + s1.count("T") > s1.count("G") + s1.count("A"):
            s = _reverse_complement(s)
            s1 = s[:al]
            s2_rev = s[-al:]
            s3 = s[al:-al]

        # s2 is the right arm read 3'→5' (reverse of genomic orientation)
        # This is how the original algorithm compares arm positions
        s2 = s2_rev[::-1]
        ll = len(s3)  # original spacer length

        if any(n not in self.nuc for n in s1 + s2):
            raise ValueError("non-ACGT base in arms")

        # --- Trim ends ---
        s1, s2 = self._trim_ends(s1, s2)
        al = len(s1)

        # --- Step 1: compute scores, find minimum prefix, cut ---
        scoring_array = ""
        pairing_scores: list[float] = []
        for j in range(al):
            val = (s1[j] == s2[j] == "G") or (s1[j] == s2[j] == "A")
            scoring_array += "1" if val else "0"
            pairing_scores.append(self._pair_score(s1[j], s2[j]))

        stacking_scores = self._calc_stacking(scoring_array)
        score_list = [pairing_scores[i] + stacking_scores[i] for i in range(al)]

        prefix_sum = 0.0
        min_prefix = 0.0
        excl = 0
        for i in range(al):
            prefix_sum += score_list[i]
            if prefix_sum < min_prefix:
                min_prefix = prefix_sum
                excl = i + 1

        cut = min(excl, max(al - max(self.min_al, ll), 0))
        while cut < al and ((s1[cut] != s2[cut]) or s1[cut] in ("C", "T")):
            cut += 1
        if cut >= al:
            cut = al

        s1 = s1[cut:]
        s2 = s2[cut:]
        al = len(s1)
        putative_triplex = s1 + "[" + s3 + "]" + s2[::-1]

        stacking_score_final = sum(stacking_scores)
        pairing_score_final = sum(pairing_scores)
        total_score = round(stacking_score_final + pairing_score_final, 3)

        if ll > al:
            return {
                "stacking_score": stacking_score_final,
                "pairing_score": pairing_score_final,
                "total_score": total_score,
                "putative_triplex": putative_triplex,
            }

        # --- Step 2: spacer avalanche ---
        pairing_scores = pairing_scores[cut:]
        scoring_array = scoring_array[cut:]
        stacking_scores = self._calc_stacking(scoring_array)
        score_list = [
            pairing_scores[i] + stacking_scores[i]
            for i in range(len(pairing_scores))
        ]

        avalance = max(0, min(al // self.m, self.MS) - ll)
        prefix_sum = 0.0
        for i in range(al - avalance):
            prefix_sum += score_list[i]
        max_prefix = prefix_sum
        excl = al - avalance
        for i in range(al - avalance, al):
            prefix_sum += score_list[i]
            if prefix_sum > max_prefix:
                max_prefix = prefix_sum
                excl = i + 1

        s_l_adj = s1[: max(excl, self.min_al)]
        s_r_adj = s2[: max(excl, self.min_al)]
        al_adj = len(s_l_adj)
        sp_adj = s1[al_adj:] + s3 + s2[al_adj:][::-1]

        stacking_score_final = sum(stacking_scores)
        pairing_score_final = sum(pairing_scores)
        total_score = round(stacking_score_final + pairing_score_final, 3)
        putative_triplex = s_l_adj + "[" + sp_adj + "]" + s_r_adj[::-1]

        if al_adj < len(sp_adj):
            return {
                "stacking_score": stacking_score_final,
                "pairing_score": pairing_score_final,
                "total_score": total_score,
                "putative_triplex": putative_triplex,
            }

        # --- Step 3: handle zero-spacer case ---
        pairing_scores = pairing_scores[: max(excl, self.min_al)]
        scoring_array = scoring_array[: max(excl, self.min_al)]
        stacking_scores = self._calc_stacking(scoring_array)

        if len(sp_adj) == 0:
            cutoff = min(self.min_spacer, len(s_l_adj)) // 2
            while (
                cutoff > 1
                and (
                    (s_l_adj[-cutoff] == s_r_adj[-cutoff] == "G")
                    or (s_l_adj[-cutoff] == s_r_adj[-cutoff] == "A")
                )
            ):
                cutoff -= 1
            sp_adj = s_l_adj[-cutoff:] + s_r_adj[-cutoff:][::-1]
            s_l_adj = s_l_adj[:-cutoff]
            s_r_adj = s_r_adj[:-cutoff]

        stacking_score_final = sum(stacking_scores)
        pairing_score_final = sum(pairing_scores)
        total_score = round(stacking_score_final + pairing_score_final, 3)
        putative_triplex = s_l_adj + "[" + sp_adj + "]" + s_r_adj[::-1]

        return {
            "stacking_score": stacking_score_final,
            "pairing_score": pairing_score_final,
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
