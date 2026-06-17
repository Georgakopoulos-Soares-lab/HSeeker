#!/usr/bin/env Rscript

# Run the Bioconductor Triplex package on the curated validation CSV.
#
# This is a helper for scripts/compare_triplex.py. It writes the raw Triplex
# presence/absence calls into results/fast_gated/triplex_predictions_raw.csv.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[startsWith(args, file_arg)][1])
if (is.na(script_path)) {
  script_path <- file.path("scripts", "run_triplex_benchmark.R")
}
script_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
local_lib <- file.path(root, ".r-lib")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(triplex)
  library(Biostrings)
})

input_csv <- file.path(root, "hdna_experimental_sequences_final.csv")
output_csv <- file.path(root, "results", "fast_gated", "triplex_predictions_raw.csv")
dir.create(dirname(output_csv), showWarnings = FALSE, recursive = TRUE)

df <- read.csv(input_csv, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", check.names = FALSE)

get_col <- function(df, candidates) {
  normalized <- tolower(sub("^\ufeff", "", names(df)))
  for (candidate in candidates) {
    idx <- match(tolower(candidate), normalized)
    if (!is.na(idx)) return(df[[idx]])
  }
  stop(sprintf("Missing required column. Tried: %s", paste(candidates, collapse = ", ")))
}

clean_sequence <- function(x) {
  gsub("[^A-Za-z]", "", toupper(x))
}

record_id <- get_col(df, c("record_id", "sequence_id", "id"))
sequence_name <- get_col(df, c("sequence_name", "name"))
label <- get_col(df, c("label", "forming_label", "experimental_label"))
sequence <- clean_sequence(get_col(df, c("sequence_5to3", "sequence", "dna_sequence")))

run_one <- function(seq) {
  dna <- DNAString(seq)
  elapsed <- system.time({
    invisible(capture.output({
      hits <- triplex.search(dna)
    }))
  })[["elapsed"]]
  n <- length(hits)
  if (n == 0) {
    return(data.frame(
      triplex_n_hits = 0,
      triplex_pred = "non-forming",
      triplex_elapsed_seconds = elapsed,
      triplex_best_score = NA_real_,
      triplex_best_pvalue = NA_real_,
      triplex_best_type = NA_integer_,
      triplex_best_strand = NA_character_,
      triplex_best_start = NA_integer_,
      triplex_best_end = NA_integer_,
      triplex_best_width = NA_integer_,
      triplex_best_loop_start = NA_integer_,
      triplex_best_loop_end = NA_integer_,
      triplex_best_loop_width = NA_integer_,
      triplex_best_insertions = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }
  scores <- score(hits)
  best <- which.max(scores)
  data.frame(
    triplex_n_hits = n,
    triplex_pred = "forming",
    triplex_elapsed_seconds = elapsed,
    triplex_best_score = scores[best],
    triplex_best_pvalue = pvalue(hits)[best],
    triplex_best_type = type(hits)[best],
    triplex_best_strand = as.character(strand(hits)[best]),
    triplex_best_start = start(hits)[best],
    triplex_best_end = end(hits)[best],
    triplex_best_width = width(hits)[best],
    triplex_best_loop_start = lstart(hits)[best],
    triplex_best_loop_end = lend(hits)[best],
    triplex_best_loop_width = lwidth(hits)[best],
    triplex_best_insertions = ins(hits)[best],
    stringsAsFactors = FALSE
  )
}

rows <- vector("list", length(sequence))
for (i in seq_along(sequence)) {
  cat(sprintf("Triplex %d/%d: %s\n", i, length(sequence), record_id[i]))
  rows[[i]] <- run_one(sequence[i])
}

out <- cbind(
  data.frame(
    sequence_id = record_id,
    sequence_name = sequence_name,
    label = label,
    sequence_length = nchar(sequence),
    stringsAsFactors = FALSE
  ),
  do.call(rbind, rows)
)

write.csv(out, output_csv, row.names = FALSE, na = "")
cat(sprintf("Wrote %s\n", output_csv))
