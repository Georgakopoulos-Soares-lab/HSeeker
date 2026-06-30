#!/usr/bin/env Rscript
# benchmarks/triplex_benchmark.R
#
# Benchmark the triplex Bioconductor package on a FASTA file using
# chunk-level parallelism via mclapply (fork).
#
# Chr1 (~249 Mbp) is split into CHUNK_BP-sized windows with OVERLAP bp
# of padding at each boundary to avoid truncating hits that straddle
# chunk edges.  The master pre-extracts each batch of chunk strings before
# each mclapply call so workers receive data via CoW (no FASTA re-reads),
# while the fork-time commit per worker stays at package overhead + ~1 Mbp
# of strings rather than the full chr1 sequence.
#
# Note on parallelism limits: under Linux overcommit_memory=2 with no swap,
# fork() can fail when Committed_AS + ncores×VmData > CommitLimit.  On the
# TACC LS6 node used here, up to 14 workers are stable; the 16-core data
# point in triplex_sweep.json was collected in a standalone run under lower
# system load.
#
# Peak RSS and wall time tracked externally via:
#   /usr/bin/time -v Rscript benchmarks/triplex_benchmark.R ...
#
# Usage:
#   Rscript benchmarks/triplex_benchmark.R [fasta] [output.csv] [ncores]
#
# Defaults:
#   fasta   = benchmarks/data/chr1.fa
#   output  = benchmarks/data/triplex_chr1_bench.csv
#   ncores  = 1

suppressPackageStartupMessages({
  library(triplex)
  library(Biostrings)
  library(parallel)
})

# Worker receives pre-extracted chunk string via mclapply fork.
# Fork inherits loaded packages — no per-worker dyn.load from disk.
.triplex_chunk <- function(task) {
  tryCatch({
    sub <- DNAString(task$seq_str)
    res <- triplex.search(sub)
    if (length(res) == 0L) return(NULL)
    df  <- as.data.frame(res)

    if (!is.na(task$owned_end_rel)) {
      df <- df[df$start <= task$owned_end_rel, , drop = FALSE]
    }
    if (nrow(df) == 0L) return(NULL)

    df$start  <- df$start + task$cs - 1L
    df$end    <- df$end   + task$cs - 1L
    df$seq_id <- task$seq_id
    df
  }, error = function(e) {
    message(sprintf("  chunk %d error: %s", task$ci, conditionMessage(e)))
    NULL
  })
}

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
argv    <- commandArgs(trailingOnly = TRUE)
fa_path <- if (length(argv) >= 1) argv[[1]] else "benchmarks/data/chr1.fa"
out_csv <- if (length(argv) >= 2) argv[[2]] else "benchmarks/data/triplex_chr1_bench.csv"
ncores  <- if (length(argv) >= 3) as.integer(argv[[3]]) else 1L

CHUNK_BP <- 1000000L   # 1 Mbp per chunk
OVERLAP  <- 1000L      # bp padding to catch hits that span chunk boundaries

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
cat(sprintf("%-20s %s\n", "triplex version:",   as.character(packageVersion("triplex"))))
cat(sprintf("%-20s %s\n", "Biostrings version:", as.character(packageVersion("Biostrings"))))
cat(sprintf("%-20s %s\n", "FASTA:",              fa_path))
cat(sprintf("%-20s %d\n", "Cores:",              ncores))
cat(sprintf("%-20s %s bp\n", "Chunk size:",      format(CHUNK_BP, big.mark = ",")))
cat(sprintf("%-20s %s bp\n", "Overlap:",         format(OVERLAP, big.mark = ",")))
cat(sprintf("%-20s %s\n", "Started:",            format(Sys.time())))
cat(strrep("-", 60), "\n")
flush.console()

# ---------------------------------------------------------------------------
# Load FASTA
# ---------------------------------------------------------------------------
t_load <- proc.time()[["elapsed"]]
seqs   <- readDNAStringSet(fa_path)
cat(sprintf("Loaded %d record(s), %s bp total (%.1fs)\n",
            length(seqs),
            format(sum(width(seqs)), big.mark = ","),
            proc.time()[["elapsed"]] - t_load))
flush.console()

# ---------------------------------------------------------------------------
# Chunk builder
# ---------------------------------------------------------------------------
build_chunks <- function(seq_len, chunk_bp, overlap) {
  starts <- seq(1L, seq_len, by = chunk_bp)
  ends   <- pmin(starts + chunk_bp - 1L + overlap, seq_len)
  data.frame(chunk_id  = seq_along(starts),
             abs_start = starts,
             abs_end   = ends,
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Main scan
# ---------------------------------------------------------------------------
t0       <- proc.time()[["elapsed"]]
all_hits <- vector("list", length(seqs))

for (rec_idx in seq_along(seqs)) {
  seq_id  <- names(seqs)[[rec_idx]]
  seqlen  <- length(seqs[[rec_idx]])
  chunks  <- build_chunks(seqlen, CHUNK_BP, OVERLAP)
  n_ch    <- nrow(chunks)

  cat(sprintf("Scanning '%s' (%s bp) → %d chunks @ %d cores ...\n",
              seq_id, format(seqlen, big.mark = ","), n_ch, ncores))
  flush.console()

  dna_seq     <- seqs[[rec_idx]]
  chunk_hits  <- vector("list", n_ch)
  batch_size  <- max(ncores, 1L)

  # Process in batches of ncores: extract batch strings, fork, collect.
  # Keeps pre-fork RSS = packages + chr1 + batch_size×1MB instead of
  # packages + chr1 + all_chunks, reducing fork commit charge.
  for (b_start in seq(1L, n_ch, by = batch_size)) {
    b_end   <- min(b_start + batch_size - 1L, n_ch)
    b_tasks <- lapply(b_start:b_end, function(ci) {
      cs <- chunks$abs_start[[ci]]
      ce <- chunks$abs_end[[ci]]
      owned_end_rel <- if (ci < n_ch)
        chunks$abs_start[[ci + 1L]] - cs
      else
        NA_integer_
      list(
        seq_str       = as.character(subseq(dna_seq, start = cs, end = ce)),
        cs            = cs,
        ci            = ci,
        owned_end_rel = owned_end_rel,
        seq_id        = seq_id
      )
    })

    if (ncores == 1L) {
      b_hits <- lapply(b_tasks, .triplex_chunk)
    } else {
      b_hits <- mclapply(b_tasks, .triplex_chunk, mc.cores = length(b_tasks))
    }
    chunk_hits[b_start:b_end] <- b_hits
  }
  rm(dna_seq); gc(verbose = FALSE)

  valid <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0L,
                  chunk_hits)
  if (length(valid) > 0L)
    all_hits[[rec_idx]] <- do.call(rbind, valid)
}

wall_sec <- proc.time()[["elapsed"]] - t0

# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------
cat(strrep("-", 60), "\n")

if (all(vapply(all_hits, is.null, logical(1L)))) {
  out <- data.frame(seq_id = character(),
                    start  = integer(),
                    end    = integer(),
                    width  = integer(),
                    score  = numeric(),
                    stringsAsFactors = FALSE)
} else {
  out <- do.call(rbind, Filter(Negate(is.null), all_hits))
  rownames(out) <- NULL
}

cat(sprintf("Total hits  : %s\n",              format(nrow(out), big.mark = ",")))
cat(sprintf("Wall time   : %.1f s  (%.2f min)\n", wall_sec, wall_sec / 60))
cat(sprintf("Throughput  : %.2f Mbp/s\n",      sum(width(seqs)) / 1e6 / wall_sec))
cat(sprintf("Output CSV  : %s\n",              out_csv))
cat(sprintf("Finished    : %s\n",              format(Sys.time())))

write.csv(out, out_csv, row.names = FALSE)

# Machine-readable summary alongside the CSV
summary_txt <- paste0(out_csv, ".summary.txt")
writeLines(c(
  sprintf("triplex_version=%s",   as.character(packageVersion("triplex"))),
  sprintf("fasta=%s",             fa_path),
  sprintf("ncores=%d",            ncores),
  sprintf("total_hits=%d",        nrow(out)),
  sprintf("wall_sec=%.3f",        wall_sec),
  sprintf("wall_min=%.3f",        wall_sec / 60),
  sprintf("throughput_mbps=%.4f", sum(width(seqs)) / 1e6 / wall_sec)
), summary_txt)

cat(sprintf("Summary     : %s\n", summary_txt))
