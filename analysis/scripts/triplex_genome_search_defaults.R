#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    stop("Usage: triplex_genome_search_defaults.R <input_fasta> <output_csv>", call. = FALSE)
  }
  .libPaths(c(".r-lib", .libPaths()))
  library(triplex)
  library(Biostrings)
})

input_fasta <- args[[1]]
output_csv <- args[[2]]

seqs <- readDNAStringSet(input_fasta)
rows <- list()
row_idx <- 1L
t0 <- proc.time()[["elapsed"]]

for (record_idx in seq_along(seqs)) {
  seq_id <- names(seqs)[[record_idx]]
  dna <- seqs[[record_idx]]
  res <- triplex.search(dna)
  hit_count <- length(res)
  if (hit_count > 0) {
    scores <- score(res)
    pvalues <- pvalue(res)
    types <- type(res)
    starts <- start(res)
    ends <- end(res)
    widths <- width(res)
    for (i in seq_len(hit_count)) {
      rows[[row_idx]] <- data.frame(
        seq_id = seq_id,
        start = starts[[i]],
        end = ends[[i]],
        width = widths[[i]],
        triplex_score = scores[[i]],
        triplex_pvalue = pvalues[[i]],
        triplex_type = types[[i]],
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1L
    }
  }
}

runtime_sec <- proc.time()[["elapsed"]] - t0

if (length(rows) == 0) {
  out <- data.frame(
    seq_id = character(),
    start = integer(),
    end = integer(),
    width = integer(),
    triplex_score = numeric(),
    triplex_pvalue = numeric(),
    triplex_type = integer(),
    stringsAsFactors = FALSE
  )
} else {
  out <- do.call(rbind, rows)
}

attr(out, "runtime_sec") <- runtime_sec
write.csv(out, output_csv, row.names = FALSE)
writeLines(sprintf("triplex_runtime_sec,%.6f", runtime_sec), paste0(output_csv, ".runtime.txt"))
