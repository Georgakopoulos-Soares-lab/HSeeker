# WS3-A1 / WS3-A2 — comparator builds and target-class characterization (reviewer R3.6)

**Status:** **PARTIAL — builds and interface characterization done; benchmark runs
still blocked on WS0-A1 (benchmark freeze).** Not to be marked complete.
**Date:** 2026-09-23
**Reviewer item R3.6 (extract):** "The manuscript also discusses NeSSie and nBMST as
relevant tools but does not include them in the main experimental benchmark. A more
comprehensive comparison with algorithms designed specifically to identify mirror
repeats/H-DNA would strengthen the manuscript considerably. If technical limitations
prevent inclusion of particular tools, these should be explicitly documented."

## Why this was done now

Both actions are graph-blocked on WS0-A1, but only for the **run** — the frozen
benchmark decides what sequences go in. Building the tools and establishing how
their output maps to an H-DNA call is independent of benchmark content, so it was
done ahead of the freeze to de-risk the comparison.

---

## WS3-A1 — NeSSie

**Source:** `github.com/B3rse/nessie`, commit
`dbe6cb2d13050d9b9877b3ffd317bb746143f21a` (2018-06-26)
**Build:** `make` → "Successfully built nessie!" (g++, warnings only, no errors)
**Publication:** Berselli, Lavezzo, Toppo, *Bioinformatics* 2018, doi:10.1093/bioinformatics/bty142

### Relevant modes

- `-M/--mirror` — mirror-symmetric motifs
- `-T/--triplex` — "DNA-triplex forming potential", with `-p` (percent non-purine
  bases permitted)
- degeneracy via `-m` (mismatch %), `-g` (gap %), `-k`/`-K` (min/max motif length)

### Target-class finding — NeSSie does not search the same object as HSeeker

Measured, on benchmark sequences:

| Input | NeSSie `-M -k 8` | HSeeker |
|---|---|---|
| HDNA0002 `pGG32` (forming, canonical H-DNA) | **no hit** | hit, arm 15, spacer 2 |
| HDNA0053 poly-A20 (non-forming) | hit: `AAAAAAAA` | hit, arm 10, spacer 0 |
| HDNA0063 `FXN_GAA66` | `-T`: `AAGAAGAA` ×64; `-M -k 10`: `GAAGAAGAAG` ×63 | — |

NeSSie's `-M` reports **k-mers that are themselves mirror-symmetric**, and `-T`
reports **purine-rich k-mers**. Neither models an *arm–spacer–arm* motif, which is
what H-DNA is and what HSeeker detects. That is why a canonical H-DNA positive
(pGG32, arms separated by a `TATA` spacer) produces no `-M` hit, while a
homopolymer does.

Output is `motif string / @counts / @indexes` per sequence — **no per-locus score
and no arm/spacer decomposition**.

**Consequence for the response letter:** NeSSie cannot be scored as a drop-in H-DNA
classifier. Any comparison must either (a) state an explicit hit→call mapping rule,
or (b) present NeSSie as a symmetry-search tool with a different target class. This
is the documented "technical limitation" R3.6 asks for, and it is a substantive
point in HSeeker's favour rather than an evasion.

---

## WS3-A2 — nBMST / non-B_gfa

**Source:** `github.com/abcsFrederick/non-B_gfa`, commit
`a891b59a4f262c83f990739d5802863e6c621175` (2026-06-22)
**Build:** `make` → `gfa` binary (gcc -O2, clean)

### Closest comparator of the two

`gfa` has a dedicated Mirror_Repeat module whose parameters map directly onto
HSeeker's:

| gfa | HSeeker |
|---|---|
| `-minMRrep` (min length of half a mirror repeat, default 10) | `minrep` |
| `-maxMRspacer` (max distance between halves, default 100) | `maxspacer` |

### Operational defect found — fixed-size path buffer

```
./gfa -seq /long/absolute/path/smoke.fa -out /long/absolute/path/nbmst ...
  *** buffer overflow detected *** : terminated      (no output produced)

./gfa -seq s.fa -out o ...                            (short relative paths)
  -> o_MR.tsv, o_MR.gff, plus APR/DR/GQ/IR/STR/Z files
```

`gfa` aborts on long path arguments. **Workaround for the TACC run: invoke it from
inside the working directory with short relative paths.** This must be recorded in
Methods, or the run will fail silently on a cluster with deep scratch paths.

### Output and a first comparison

`gfa -seq s.fa -out o -minMRrep 8 -maxMRspacer 10`:

```
Sequence_name           Type            Start Stop Length Score Strand Repeat Spacer Composition
HDNA0002_forming        Mirror_Repeat   1     32   32     NA    +      14     4      5A/0C/9G/0T
HDNA0053_non_forming    Mirror_Repeat   1     20   20     NA    +      10     0      10A/0C/0G/0T
```

against HSeeker on the same four sequences:

| Record | Label | nBMST MR | HSeeker |
|---|---|---|---|
| HDNA0002 `pGG32` | forming | hit (repeat 14, spacer 4) | hit (arm 15, spacer 2, score 152.55) |
| HDNA0030 `pcMyc` | forming | **no hit** | hit (arm 8, spacer 7, score 92.03) |
| HDNA0053 poly-A20 | non-forming | hit | hit (score 83.30) |
| SYN0001 (synthetic negative) | non-forming | no hit | **no hit** |

Three things this establishes:

1. **`Score` is literally `NA`.** nBMST emits motif geometry, not a propensity
   score, so a hit→call mapping rule is mandatory — confirming the plan's concern
   empirically rather than by assumption.
2. **nBMST missed `pcMyc`**, a canonical experimentally-verified H-DNA former that
   HSeeker detects, consistent with a stricter/exact mirror criterion.
3. **Both tools call poly-A20 a mirror repeat.** HSeeker's production CLI drops it
   via the homopolymer filter; the library default used here does not. This is a
   concrete illustration of why the composition filters exist (R2.6) — and a
   reminder that the two arms must be run under one pinned configuration (WS2-A1).

---

## What remains before these actions can be marked complete

- Freeze the benchmark (**WS0-A1**, blocked on DEC-2/DEC-3).
- Run both tools on the frozen set and deposit `nessie_benchmark.csv`,
  `nbmst_benchmark.csv` with exact commands.
- State the hit→call mapping rule for each in Methods.
- **Hard rule carried over from the plan:** the RTR's existing sentence "We run
  nBMST with the sequences and show performance lower than ours" must not be
  submitted until those outputs exist. Nothing measured here yet supports a
  performance claim — the 4-sequence smoke test is an interface check, not a
  benchmark.

## Reproduce

```
git clone --depth 1 https://github.com/B3rse/nessie.git && cd nessie && make
./nessie -I in.fa -O out.txt -T -k 8 -K 50 -p 10 -m 10

git clone --depth 1 https://github.com/abcsFrederick/non-B_gfa.git && cd non-B_gfa && make
./gfa -seq s.fa -out o -minMRrep 8 -maxMRspacer 10     # short relative paths only
```
