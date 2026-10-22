# Plan: Benchmark curation, 1:1 balancing, and dependency-graph model of the revision plan

## User decisions (locked)
- **Removal scope**: keep every record that has paper provenance; remove only what cannot be traced; justify every decision.
- **Balance**: 1:1 forming:non-forming.
- **Synthetic negatives**: a subset of the reviewer's five classes, with the perturbation-focused class as the core.
- **Graph**: interactive HTML.
- **Immutability**: reviewer comments are quoted verbatim; existing deliverables (`HSeeker_Reviewer_Answer_Plan_v2.md`, `hseeker_revision_ledger.yaml`, `hdna_benchmark_provenance_audit.csv`) are NOT modified — all outputs are new files.

## Step 1 — Curate the experimental benchmark (v3)

Rule: a record stays iff a primary-source paper was traced in the provenance audit. Applying it:

**Remove — no paper provenance (3 records):**
| Record | Label | Justification |
|---|---|---|
| HDNA0035 GG32-long | forming | 42-nt oligo = pGG32 core + EcoRI/BamHI overhangs; no primary source for this exact construct; near-duplicate of HDNA0002 |
| HDNA0036 Mycm1 | non-forming | pMycAG-derived 22-nt mutant; not present in Belotserkovskii 2007 (constructs were WT/TT/Q/AG/AA); no primary source found |
| HDNA0037 MycM2 | non-forming | Mycm1 + G21T; no primary source found |

**Exclude — provenance exists but the record is invalid as an ACGT benchmark entry (2 records, recommended default, user-overridable):**
| Record | Label | Justification |
|---|---|---|
| HDNA0026 R2-ino | non-forming | The experimental label belongs to an inosine-substituted oligo (Del Mundo 2017); the FASTA silently converts inosine→G, manufacturing a perfect mirror repeat that was never tested experimentally — it produced 1 of 2 false positives |
| HDNA0080 G10TTAA_AG5 | forming | Label-scope error: 'forming' refers to an intermolecular TFO complex; the tool predicts intramolecular H-DNA only — it manufactured a false negative |

**Keep with flag (1 record):** HDNA0033 AG32-30nt — provenance verified (Rooney & Moore 1995, exact 30-nt match) but the paper reports the H-DNA transition as strongly impaired vs the CSV's `forming` label. Kept with `label_conflict=true`; WS2 sensitivity analysis reports metrics with both labelings.

**Result: 75 records = 69 forming / 6 non-forming.** Output: `hdna_benchmark_experimental_v3.csv` with added columns `curation_decision` (kept/removed/excluded), `curation_justification`, `label_conflict`.

## Step 2 — Generate 63 synthetic negatives → 1:1 balanced set (138 records)

Seeded generator script `generate_synthetic_negatives.py` (Python, seed fixed, all parameters inside). Mix:

| Class | n | Method | Grounded rationale |
|---|---|---|---|
| Mirror-disrupted point mutants | 20 | 1–3 purine-preserving substitutions (G→A) in ONE arm of real forming positives | Directly answers reviewer R2.1 "artificially perturb the positive cases to disrupt H-DNA"; hardest negatives (purity kept, mirror broken) |
| Dinucleotide-preserving shuffles | 12 | Shuffle real positives preserving dinucleotide composition; re-draw until mirror score destroyed | Composition-matched controls; destroys long-range mirror symmetry |
| Random GC/length-matched | 12 | Random sequences matched to positives' length and GC distribution | Reviewer's "random sequences" |
| Other-structure | 12 | G4 (c-MYC Pu27, telomeric, thrombin aptamer variants), Z-DNA (CG)n, mixed B-DNA | Reviewer's "sequences with other known structures, like B-DNA, or DNA quadruplexes" |
| poly-A/poly-T homopolymers | 4 | Lengths 15–40 | Reviewer's "including poly-A"; A20 experimentally non-forming (Hanvey 1988) |
| Perfect-mirror, balanced-composition | 3 | Mirror-symmetric but ~50% GC, mixed purine/pyrimidine strands | Reviewer's "test mirror sequences" — mirror without purine richness must not be called H-DNA |

**Grounded guardrails:**
- NO poly-G/poly-C negatives: poly-dG·dC genuinely forms H-DNA (Kohwi 1988; our own HDNA0054 is forming poly-dC30).
- Validation is structural only: mirror score below threshold, purity check, exact + reverse-complement dedup against all 75 experimental sequences and within the synthetic set, length/GC match verified. Synthetic records are NEVER filtered by HSeeker's own predictions (avoids circular benchmarking).
- Every synthetic record: `record_type=synthetic`, `generation_method`, `source_record`, `seed`, `expected_label_rationale`.

Output: `hdna_benchmark_balanced_v3.csv` (138 rows = 75 experimental + 63 synthetic; 69:69). Dual-reporting policy unchanged: metrics reported on experimental-only AND balanced sets.

## Step 3 — Dependency analysis + interactive graph model

Build `hseeker_revision_graph.json` (nodes/edges source of truth), then render `hseeker_revision_graph.html` (pyvis interactive; install if needed) + `hseeker_revision_graph_preview.png` (Graphviz static snapshot for QA/pasting).

- **Nodes**: 26 reviewer items (R1, R2.1–R2.11, R3.1–R3.12; verbatim quotes, unaltered), 5 author decisions (DEC-1..5), ~34 actions (existing WS0–WS6 actions from the ledger + new CUR-1 curate, SYN-1 generate, SYN-2 validate, SYN-3 freeze), key artifacts (CSVs, figures, response letter), 2 gates (benchmark-v3 frozen; numbers verified).
- **Edges**: `requires` (action←artifact/decision), `addresses` (action→reviewer item), `informs` (decision→action), `invalidates` (curation/balancing→every downstream metric, comparator, and manuscript number that must be re-derived), `cites` (response→artifact).
- **Highlights**: critical path (curation → freeze → re-benchmark → stats → verified numbers → response letter); keystone gate "benchmark v3 frozen"; the full cascade of what a benchmark change forces to re-run (WS2 metrics, WS3 comparators, R3.1/R3.5/R2.1 responses, affected manuscript consistency-register rows).
- Visual encoding: color by node type, border/shape by status (verified/pending), critical path emphasized, details in hover tooltips.

## Step 4 — Verification
- CSV checks: 75/138 rows, 69:69 balance, no duplicate sequences (exact + revcomp), all synthetic rows have method+seed+rationale, all removals have justifications.
- Graph checks: JSON parses; every ledger action + all 26 reviewer items present; every edge endpoint resolves; HTML file well-formed; PNG preview passes media check.

## Deliverables (all new files; nothing existing modified)
1. `hdna_benchmark_experimental_v3.csv`
2. `hdna_benchmark_balanced_v3.csv`
3. `generate_synthetic_negatives.py`
4. `hseeker_revision_graph.html`
5. `hseeker_revision_graph.json`
6. `hseeker_revision_graph_preview.png`

## Out of scope (unchanged)
Re-running HSeeker/comparators on v3, genome-scale validation, response-letter drafting — remain ledger-tracked future work (now represented as pending nodes in the graph).
