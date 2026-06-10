from typing import Optional
from Bio.Seq import Seq 
import math
from attr import field 
import attr 

@attr.s(slots=True, frozen=True)
class HDNA:
    s: str = field()
    al: int = field()
    s1: str = field()
    s2: str = field()
    s3: str = field()
    s_l_adj: str = field()
    s_r_adj: str = field()
    sp_adj: str = field()
    putative_triplex: str = field()
    stacking_scores: list[float] = field()
    pairing_scores: list[float] = field()
    stacking_score: float = field()
    pairing_score: float = field()
    total_score: float = field()

    @classmethod
    def from_dict(cls, data: dict) -> "HDNA":
        return cls(
            s=data['s'],
            al=data['al'],
            s1=data['s1'],
            s2=data['s2'],
            s3=data['s3'],
            s_l_adj=data['s_l_adj'],
            s_r_adj=data['s_r_adj'],
            sp_adj=data['sp_adj'],
            putative_triplex=data['putative_triplex'],
            stacking_scores=data['stacking_scores'],
            pairing_scores=data['pairing_scores'],
            stacking_score=data['stacking_score'],
            pairing_score=data['pairing_score'],
            total_score=data['total_score']
        )

class Farmer:
    defaults = {
                'd': 2,
                'dv': 2,
                'MS': 50,
                'm': 3,
                'min_al': 8
                }
    def __init__(self, d = None, dv = None, min_al = None, MS = None, m = None):
        self.scoring = {
                        'G': 7.6,
                        'A': 3.83
                        }
        self.stacking_score = 5
        self.nuc = {'A', 'G', 'C', 'T'}
        self.RT = 0.6
        self.MS = MS
        self.dv = dv
        self.d = d
        self.m = m
        self.min_al = min_al
        self.min_spacer = 4
        if self.min_al is None:
            self.min_al = Farmer.defaults['min_al']
        if self.m is None:
            self.m = Farmer.defaults['m']
        if self.d is None:
            self.d = Farmer.defaults['d']
        if self.dv is None:
            self.dv = Farmer.defaults['dv']
        if self.MS is None:
            self.MS = Farmer.defaults['MS']
    def calc_stacking_scores(self, scoring_array: str) -> list:
        scores = [0]
        incorrect_tail = -1
        for j in range(1, len(scoring_array)):
            if scoring_array[j-1] == scoring_array[j] == '1':
                scores.append(self.stacking_score)
                incorrect_tail = -1
            else:
                scores.append(-self.discount(max(incorrect_tail, 0)))
                incorrect_tail += 1
        return scores
    def discount(self, n: int) -> float:
        p = self.stacking_score ** (n+1) 
        return p if n+1<2 else round(p / self.d, self.dv) 
    def score(self, nuc1: str, nuc2: str) -> float:
        if (nuc1 == nuc2 == 'G') or (nuc1 == nuc2 == 'A'):
            return self.scoring[nuc1]
        return - (self.scoring['A'] + self.scoring['G']) / 2
    def decouple(self, s: str) -> tuple[str, str, str]:
        s = s.upper(); a = 0; l = len(s)
        for i in range(l):
            if s[i] != s[l-i-1] or i >= l -i - 1:
                break
            a += 1
        b = l - 2 * a
        return (s[:a], s[-a:], s[a:a+b])
    def trim_ends(self, s1: str, s2: str) -> tuple[str, str]:
        cutoff = 0
        for i in range(len(s1)):
            if (s1[i] == s2[i] == 'G') or (s1[i] == s2[i] == 'A'):
                break
            cutoff += 1
        return s1[cutoff:], s2[cutoff:]
    def print(self, data: dict) -> None:
        print()
        print(f"H-DNA sequence: {data['s']} - {data['al']=}")
        print(f"H-DNA sequence: {data['s1']}-{data['s3']}-{data['s2'][::-1]}\n" + "*" * 20)
        print(f"H-DNA sequence: {data['s_l_adj']}-{data['sp_adj']}-{data['s_r_adj']}\n" + "*" * 20)
        print(f"Stacking array: {data['stacking_scores']}")
        print(f"Total stacking score: {data['stacking_score']}.")
        print("PAIRING")
        print("*" * 20)
        print(f"Pairing array: {data['pairing_scores']}")
        print(f"Total pairing score: {data['pairing_score']}.")
        print(f"Total score: {data['total_score']}.")
    def stability(self, 
                  s: str, 
                  al: Optional[int] = None, 
                  loop : Optional[str] = None, 
                  trim: bool = True) -> dict:
        rev = lambda seq: str(Seq(seq).reverse_complement())
        s = s.upper()
        if al:
            s1, s2, s3 = s[:al], s[-al:], s[al:-al]
        else:
            s1, s2, s3 = self.decouple(s)
            al = len(s1)
        if s.count("C") + s.count("T") > s.count("G") + s.count("A"):
            s = rev(s)
        if al:
            s1, s2, s3 = s[:al], s[-al:], s[al:-al]
        else:
            s1, s2, s3 = self.decouple(s)
            al = len(s1)
        ll = len(s3)
        s2 = s2[::-1]
        if any(n not in self.nuc for n in s1+s2):
            raise ValueError()
        # TRIM - pre-STEP - omit
        if trim:
            s1, s2 = self.trim_ends(s1, s2)
            s = s1 + s3 + s2[::-1]
            al = len(s1)
        # Step 1
        pairing_scores = []
        stacking_scores = [0]
        scoring_array = ''
        for j in range(len(s1)):
            val = (s1[j] == s2[j] == 'G') or (s1[j] == s2[j] == 'A')
            scoring_array += str(int(val))
            pairing_scores.append(self.score(s1[j], s2[j]))
        al_adj = al
        stacking_scores = self.calc_stacking_scores(scoring_array)
        stacking_score = sum(stacking_scores)
        pairing_score = sum(pairing_scores)
        total_score = round(stacking_score + pairing_score, 3)
        score_list = [pairing_scores[i] + stacking_scores[i] for i in range(len(pairing_scores))]
        prefix_sum = 0
        min_prefix = 0
        excl = 0
        for i in range(al):
            prefix_sum += score_list[i]
            if prefix_sum < min_prefix:
                min_prefix = prefix_sum
                excl = i + 1
        cut = min(excl, max(al - max(self.min_al, ll), 0))
        while (s1[cut] != s2[cut]) or s1[cut] == 'C' or s1[cut] == 'T':
            cut += 1
        s1 = s1[cut:]
        s2 = s2[cut:]
        al = len(s1)
        putative_triplex = s1 + "[" + s3 + "]" + s2[::-1]
        if ll > al:
            return HDNA.from_dict({
                    's': s,
                    'al': al,
                    's1': s1,
                    's2': s2[::-1],
                    's3': s3,
                    's_l_adj': s1,
                    's_r_adj': s2[::-1],
                    'sp_adj': s3,
                    'putative_triplex': putative_triplex,
                    'stacking_scores': stacking_scores,
                    'pairing_scores': pairing_scores,
                    'stacking_score': stacking_score,
                    'pairing_score': pairing_score,
                    'total_score': total_score})
        # Step 2
        pairing_scores = pairing_scores[cut:]
        scoring_array = scoring_array[cut:]
        stacking_scores = self.calc_stacking_scores(scoring_array)
        score_list = [pairing_scores[i] + stacking_scores[i] for i in range(len(pairing_scores))]
        avalance = max(0, min(al // self.m, self.MS) - ll)
        prefix_sum = 0
        print(f"Avail: {avalance}. {al}")
        for i in range(al-avalance):
            prefix_sum += score_list[i]
        max_prefix = prefix_sum
        excl = al
        for i in range(al-avalance, al, 1):
            prefix_sum += score_list[i]
            if prefix_sum > max_prefix:
                max_prefix = prefix_sum
                excl = i + 1
        s_l_adj = s1[:max(excl, self.min_al)]
        s_r_adj = s2[:max(excl, self.min_al)]
        al_adj = len(s_l_adj)
        sp_adj = s1[al_adj:] + s3 + s2[al_adj:][::-1]
        s = s_l_adj + sp_adj + s_r_adj[::-1]
        stacking_score = sum(stacking_scores)
        pairing_score = sum(pairing_scores)
        total_score = round(stacking_score + pairing_score, 3)
        putative_triplex = s_l_adj+ "[" + sp_adj + "]" + s_r_adj[::-1]
        if al_adj < len(sp_adj):
            return HDNA.from_dict({
                    's': s, 
                    'al': al, 
                    's1': s1, 
                    's2': s2[::-1], 
                    's3': s3,
                    's_l_adj': s_l_adj, 
                    's_r_adj': s_r_adj, 
                    'sp_adj': sp_adj,
                    'putative_triplex': putative_triplex,
                    'stacking_scores': stacking_scores, 
                    'pairing_scores': pairing_scores,
                    'stacking_score': stacking_score,
                    'pairing_score': pairing_score,
                    'total_score': total_score
                    })
        # Step 3
        pairing_scores = pairing_scores[:max(excl, self.min_al)]
        scoring_array = scoring_array[:max(excl, self.min_al)]
        stacking_scores = self.calc_stacking_scores(scoring_array)
        score_list = [pairing_scores[i] + stacking_scores[i] for i in range(len(pairing_scores))]
        if len(sp_adj) == 0:
            cutoff = min(self.min_spacer, len(s_l_adj)) // 2
            while ((s_l_adj[-cutoff] == s_r_adj[-cutoff] == 'G') \
                  or (s_l_adj[-cutoff] == s_r_adj[-cutoff] == 'A')) \
                  and cutoff > 1:
                cutoff -= 1
            sp_adj = s_l_adj[-cutoff:] + s_r_adj[-cutoff:][::-1]
            s_l_adj = s_l_adj[:-cutoff]
            s_r_adj = s_r_adj[:-cutoff]
            pairing_scores = pairing_scores[:-cutoff]
            scoring_array = scoring_array[:-cutoff]
            stacking_scores = self.calc_stacking_scores(scoring_array)
        stacking_score = sum(stacking_scores)
        pairing_score = sum(pairing_scores)
        total_score = round(stacking_score + pairing_score, 3)
        if loop is not None:
            loop = loop.upper()
            if any(n not in self.nuc for n in loop):
                raise ValueError()
            n = len(loop)
            DG_loop = math.log(n) * 1.5 * self.RT
            total_score -= DG_loop
        putative_triplex = s_l_adj + "[" + sp_adj + "]" + s_r_adj[::-1]
        data = {
                's': s,
                'al': al,
                's1': s1,
                's2': s2[::-1],
                's3': s3,
                's_l_adj': s_l_adj,
                's_r_adj': s_r_adj,
                'sp_adj': sp_adj,
                'putative_triplex': putative_triplex,
                'stacking_scores': stacking_scores,
                'pairing_scores': pairing_scores,
                'stacking_score': stacking_score,
                'pairing_score': pairing_score,
                'total_score': total_score
                }
        return HDNA.from_dict(data)
    
            # if avalance > 0:
        #     j = 0
        #     while j < avalance:
        #         if score_list[-avalance+j] > 0:
        #             j += 1
        #             continue

        #         for i in range(1, avalance-j+1, 1):
        #             end = -avalance + j + i
        #             ssum = sum(score_list[-avalance+j:end if end != 0 else None])  # O(n) — use prefix sums to make O(1)
        #             if ssum > 0:
        #                 j += i
        #                 break
        #         else:
        #             excl = - avalance + j
        #             break
        #     if excl != 0:
        #         stacking_scores = stacking_scores[:excl]
        #         pairing_scores = pairing_scores[:excl]
        #         al_adj -= abs(excl)
        #         s_l_adj = s[:al_adj]
        #         s_r_adj = s[-al_adj:]
        #         sp_adj = s[al_adj:len(s)-al_adj]
        # if excl == 0:
        #     s_l_adj, s_r_adj, sp_adj = s1, s2[::-1], s3
        
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Calculate H-DNA stability score.")
    # parser.add_argument("sequence", type=str, help="The DNA sequence to analyze.")
    # parser.add_argument("--al", type=int, default=None, help="The length of the arms. If not provided, the program will attempt to determine it automatically.")
    # parser.add_argument("--loop", type=str, default=None, help="The sequence of the loop. If not provided, no loop penalty will be applied.")
    # parser.add_argument("--no-trim", action="store_true", help="Do not trim the ends of the arms.")
    # args = parser.parse_args()
    # data = Farmer().stability(args.sequence, 
    #                           args.al, 
    #                           args.loop, 
    #                           not args.no_trim)

    # TEST SUITE
    seq = [
        ("AAGGGAGAAGGGG", "GGGGAAGAGGGAA"),
        ("AAGGGAGAAAGGG", "GGGGAAGAGGGAA"),
        ("GGGAGGGG", "GGGGAGGG"),
        ("GGGGGAGG", "GGAGGGGG"),
        ("GGGAGGGG", "GGGAGGGG"),
        ("GGAGGGG", "GGAGGGTG"),


        ("AAGGGAGAAGGGGCC", "CCGGGGAAGAGGGAA"),
        ("AAGGGAGAAGGGGC", "CGGGGAAGAGGGAA"),
        ("AAAAAAAGGGAGAAGGCGGCCCGG", "GGCCCCGGCGGAAGAGGGAAAAAAA")
        ]
    seq = [
        ("AAGGGAGAAGCGGTATAGGGCGAAGAGGGAA", None),
        ("AAGGGAGAAGCCGGTATAGGGCCGAAGAGGGAA", None),
        ("AAGGGAGAAGCCGGTATAGGGCCGAAGAGGGAA", None),


        ("AAGGGAGAAGGGGTATAGGGGAAGAGGGAA", None), #GG32
        ("AAGGGAGAAAGGGTATAGGGGAAGAGGGAA", 13), #AG32
        ("GGGAGGGGCGCTTATGGGGAGGG", None), # MYC
        ("GGGGGAGGCGCTTATGGAGGGGG", None), # MYCAA
        ("GGGAGGGGCGCTTATGGAGGGGG", 8), # MYC-AG NO
        ("GGAGGGGCGCTTATGGAGGGGG", 7), # Mycm1 NO
        ("GGAGGGGCGCTTATGGAGGGTG", 7), # MycM2 NO
        ("AAGGGAGAAGGGGCCTCGCCGGGGAAGAGGGAA", 15),
        ("AAGGGAGAAGGGGCTCGCGGGGAAGAGGGAA", 14),
        ("AAAAAAAGGGAGACAGGCGGCCCGGTCGGGCCCGGCGGACAGAGGGAAAAAAA", 24),
        ("AAAAAAAGGGAGACAGGCGGCCCGGTCGGGCCCGGCGGATAGAGGGAAAAAAA", 24),
        ("AACAAAAGGGAGACAGGCGGCCCGGTCGGGCCCGGCGGATAGAGGCAAAAAAA", 24),

        ]
    seq = [
        
           ("CCCCCCCCCCCGGGTATGGGCCCCCCCCCCC", None),
           
           ("GCGGGTTGGGCG", None),

           ("CGCCGGGCCCGGGGGAAGGGGGCCCGGGCCGC", None),

           ("CCCGCGGGGGGGGGGGGGCCCGGGGGGGGGGGGGGGGCCCGGGGGGGGGGGGGCGCCC", None),

            ("CCCGCGAGAGAGAGAGAGCCCGGGGGGGGGGGGGGGGCCCGAGAGAGAGAGAGCGCCC", None)]

    for s1, al in seq:
        print()
        data = Farmer().stability(s1, al)
        print(data)
        # breakpoint()