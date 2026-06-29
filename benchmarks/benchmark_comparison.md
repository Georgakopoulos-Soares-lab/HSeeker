# Computational Benchmark: HSeeker, triplex, and Triplexator on Human Chromosome 1

## Benchmark Configuration

All three tools were evaluated on the same input: the GRCh38 primary assembly of human chromosome 1
(`chr1.fa`, 248,956,422 bp, 252.1 MB), retrieved from UCSC.
Experiments were run on a single node of the TACC Lonestar6 HPC cluster
(AMD EPYC 7763, 32 GB RAM, Rocky Linux 8.10).
Each tool was allocated 16 CPU cores.
Wall-clock time and peak resident set size (RSS) were measured using `/usr/bin/time -v`.

---

## Tools

### HSeeker (v0.1.0)

HSeeker detects H-DNA mirror repeat motifs — intramolecular triplex-forming sequences
characterised by a purine/pyrimidine-rich arm pair flanking a central spacer.
The core scan is implemented as a compiled C extension that releases the GIL, enabling
true CPU parallelism across a Python `ThreadPoolExecutor` thread pool.
Chromosome-scale records are split into overlapping 1 Mbp chunks; each chunk is dispatched
as an independent task so all 16 threads remain fully utilised throughout the scan.

For this benchmark, HSeeker was invoked with the RMQ (Range Maximum Query) purity
prefilter enabled (`-purity-rmq`), which uses sparse-table range-maximum queries over
prefix-sum arrays to discard (center, spacer) pairs whose right arm cannot satisfy the
purity threshold before the extension loop begins.
Thermodynamic stability scoring (stacking and base-pairing scores) was applied to every
reported hit.

**Parameters:** minimum arm length 10 bp, purity threshold 90%, mismatch tolerance 10%,
maximum spacer 10 bp, overlap removal enabled, scoring enabled, purity-RMQ prefilter enabled.

**Parallelism model:** intra-record chunk parallelism via Python `ThreadPoolExecutor`;
16 worker threads; GIL released in C extension during each chunk scan.

---

### triplex (Bioconductor, v1.50.0)

The `triplex` Bioconductor R package detects intramolecular triplex-forming sequences
in double-stranded DNA by scoring all eight canonical triplex types against a statistical
model parameterised for eukaryotic sequences.
Internally, `triplex.search` is a single-threaded C routine with no native parallelism
parameter; the function signature exposes no `ncore` or equivalent argument.

To utilise 16 cores, the benchmark script (`benchmarks/triplex_benchmark.R`) manually
splits the chromosome into 249 non-overlapping 1 Mbp windows with 1,000 bp of boundary
padding, then dispatches each window as an independent `triplex.search` call via
`parallel::mclapply(mc.cores = 16)`.
Hits falling within the padding region of non-final chunks are discarded during result
merging to prevent double-counting across boundaries.
Wall time reported here is the R-internal elapsed time (i.e., excludes Rscript and
Bioconductor package load time).

**Parameters:** all defaults — `min_score = 15`, `p_value = 0.05`, `min_len = 6`,
`max_len = 25`, triplex types 0–7, eukaryotic scoring model.

**Parallelism model:** manual chunk-level parallelism via `parallel::mclapply`; 16 forked
R worker processes; no shared memory between workers.

---

### Triplexator (v1.3.2)

Triplexator detects triplex target sites (TTS) in duplex DNA and/or triplex-forming
oligonucleotides (TFO) in single-stranded sequences using canonical Hoogsteen and
reverse-Hoogsteen bond rules.
Unlike the two tools above, Triplexator has native OpenMP parallelism built into its C++
core; the `-p NUM` flag sets the thread count and the `-rm` flag selects the parallelism
strategy.

For chromosome-scale single-sequence input, `-rm 1` (parallelize TTSs) is the appropriate
mode: each putative TTS locus within the duplex sequence is processed by an independent
OpenMP thread, keeping all 16 hardware threads active with minimal coordination overhead.
Triplexator was built from source (GitHub commit on the `master` branch, CMake 4.1.1,
GCC 9.4.0, OpenMP 4.5) after patching two C++11 compatibility issues in the bundled SeqAn
library headers.

**Parameters:** all defaults — minimum length 16 bp, maximum length 30 bp, error rate 5%,
minimum guanine content 10%, repeat/low-complexity filter enabled, all motifs
(R/Y/M/P/A), duplicate detection off.

**Parallelism model:** native OpenMP; `-rm 1 -p 16`; shared-memory threading within a
single process.

---

## Results at 16 Cores

| Tool | Version | Wall time (s) | Peak RSS (MB) | Throughput (Mbp/s) | Hits | Parallelism |
|---|---|---|---|---|---|---|
| HSeeker | 0.1.0 | **25.2** | 1,965 | **10.0** | 96,729 | ThreadPoolExecutor, 16 threads |
| triplex | 1.50.0 | 57.7 | **578** | 4.4 | 81,906 | mclapply, 16 forked processes |
| Triplexator | 1.3.2 | 67.5 | 1,709 | 3.7 | 275,476 | OpenMP -rm 1, 16 threads |

Input: chr1 GRCh38, 248,956,422 bp, 252.1 MB. Node: AMD EPYC 7763, 32 GB RAM, 16 vCPUs.
HSeeker RSS includes the full-chromosome string held in memory during chunking.
triplex RSS reflects 16 forked processes each holding a 1 Mbp subsequence.
Wall time for triplex is R-internal elapsed time (excludes Rscript and package startup).

---

## Parallel Scaling Sweep (1, 2, 4, 8, 16 workers)

### HSeeker

| Workers | Wall (s) | RSS (MB) | Speedup | Efficiency |
|---|---|---|---|---|
| 1 | 187.7 | 1,612 | 1.00× | 100% |
| 2 | 105.3 | 1,614 | 1.78× | 89% |
| 4 | 59.5 | 1,599 | 3.15× | 79% |
| 8 | 35.9 | 1,719 | 5.23× | 65% |
| 16 | 25.2 | 1,965 | 7.45× | 47% |

### triplex R v1.50.0

Wall time is R-internal (excludes startup); RSS from `/usr/bin/time -v` on the full process.
The 16-core point was measured in a standalone run at lower system memory commit; all other
points are from a single continuous sweep.

| Workers | Wall (s) | RSS (MB) | Speedup | Efficiency |
|---|---|---|---|---|
| 1 | 816.3 | 644 | 1.00× | 100% |
| 2 | 407.7 | 606 | 2.00× | 100% |
| 4 | 209.9 | 599 | 3.89× | 97% |
| 8 | 112.4 | 593 | 7.26× | 91% |
| 16 | 57.7 | 578 | 14.15× | 88% |

### Triplexator v1.3.2 (native OpenMP, -rm 1)

Triplexator's `-rm 1` mode parallelises over TTS loci. On a single long chromosome the
locus discovery phase itself is single-threaded, so parallel efficiency is negligible
regardless of thread count.

| Workers | Wall (s) | RSS (MB) | Speedup | Efficiency |
|---|---|---|---|---|
| 1 | 68.3 | 1,709 | 1.00× | 100% |
| 2 | 67.8 | 1,709 | 1.01× | 50% |
| 4 | 66.5 | 1,709 | 1.03× | 26% |
| 8 | 67.8 | 1,710 | 1.01× | 13% |
| 16 | 67.5 | 1,709 | 1.01× | 6% |

---

## Discussion

**Speed.** HSeeker is the fastest of the three tools at 16 cores, completing the chr1 scan
in 25.2 s at 10.0 Mbp/s — 2.3× faster than triplex (57.7 s) and 2.7× faster than
Triplexator (67.5 s).
Although triplex achieves near-ideal parallel scaling (14.2× at 16 cores, 88% efficiency),
it starts from a much higher single-core baseline (816 s), so HSeeker remains faster in
absolute wall time at every worker count.
HSeeker's more modest relative speedup (7.4× at 16 cores, 47% efficiency) reflects
overheads from the Python chunk dispatcher and the O(n²) scoring pass, which are not
fully parallelised.

**Memory.** triplex has by far the lowest peak RSS (578 MB at 16 cores) because each
forked worker process holds only a 1 Mbp subsequence at any one time; memory stays
flat or slightly decreases as more workers split the I/O load.
HSeeker (1,965 MB) and Triplexator (1,709 MB) both load the full chromosome into a
single process address space before parallel work begins; their RSS scales with
chromosome length and, in HSeeker's case, grows modestly with worker count due to
thread-local buffers.

**Parallel scaling comparison.** HSeeker and triplex R both exhibit meaningful speedup
with increasing worker count; triplex R achieves higher relative efficiency because each
chunk is purely serial C code with no Python dispatcher overhead.
Triplexator shows essentially no scaling under `-rm 1` on a single long chromosome: the
bottleneck is sequential locus enumeration in the duplex pass, not the per-locus work
that OpenMP parallelises.

**Hit counts and comparability.** The three tools report substantially different hit counts
(81,906–275,476) because they implement distinct detection models with different default
thresholds and motif definitions.
HSeeker targets H-DNA mirror repeats with a minimum arm length of 10 bp and a 90%
purine/pyrimidine purity constraint, producing 96,729 non-overlapping hits.
triplex scores all eight intramolecular triplex types against a statistical eukaryotic
model with a minimum score threshold of 15 and reports 81,906 hits.
Triplexator searches for triplex target sites using Hoogsteen bond rules with a minimum
length of 16 bp, a 5% error rate, and an active repeat filter, reporting 275,476 TTS loci.
Direct numerical comparison of hit counts across tools is therefore not meaningful without
re-parameterisation to a common detection criterion.
