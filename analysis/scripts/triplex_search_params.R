#!/usr/bin/env Rscript
# triplex.search on every record of a FASTA file with explicit parameters (tool comparison, 2026-10-11).
# Usage: triplex_search_params.R <in.fasta> <out.csv> <min_score> <p_value> <min_len> <max_len> <seq_type>
suppressPackageStartupMessages({
  library(triplex)
  library(Biostrings)
})
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 7) stop("usage: in.fasta out.csv min_score p_value min_len max_len seq_type", call. = FALSE)
seqs <- readDNAStringSet(a[[1]])
min_score <- as.numeric(a[[3]]); p_value <- as.numeric(a[[4]])
min_len <- as.integer(a[[5]]); max_len <- as.integer(a[[6]]); seq_type <- a[[7]]
rows <- list()
for (i in seq_along(seqs)) {
  res <- triplex.search(seqs[[i]], min_score = min_score, p_value = p_value,
                        min_len = min_len, max_len = max_len, seq_type = seq_type)
  if (length(res) > 0) {
    rows[[length(rows) + 1]] <- data.frame(record_id = names(seqs)[[i]], start = start(res), end = end(res),
                                           score = score(res), pvalue = pvalue(res), type = type(res))
  }
}
out <- if (length(rows)) do.call(rbind, rows) else
  data.frame(record_id = character(), start = integer(), end = integer(), score = numeric(),
             pvalue = numeric(), type = integer())
write.csv(out, a[[2]], row.names = FALSE)
