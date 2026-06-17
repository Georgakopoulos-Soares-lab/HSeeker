#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) {
    stop("Usage: triplex_validate.R <input_csv> <output_csv>", call. = FALSE)
  }
  .libPaths(c(".r-lib", .libPaths()))
  library(triplex)
  library(Biostrings)
})

input_csv <- args[[1]]
output_csv <- args[[2]]

df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("sequence_id", "label", "sequence")
missing <- setdiff(required, names(df))
if (length(missing) > 0) {
  stop(paste("Missing columns:", paste(missing, collapse = ", ")), call. = FALSE)
}

rows <- vector("list", nrow(df))
for (i in seq_len(nrow(df))) {
  sid <- df$sequence_id[[i]]
  seq <- toupper(gsub("[^ACGT]", "N", df$sequence[[i]]))
  triplex_count <- 0L
  best_score <- NA_real_
  best_pvalue <- NA_real_
  best_type <- NA_integer_
  best_start <- NA_integer_
  best_end <- NA_integer_
  err <- ""

  if (nchar(seq) > 0) {
    tryCatch({
      res <- triplex.search(
        DNAString(seq),
        type = 0:7,
        min_score = 0,
        p_value = 1,
        min_len = 6,
        max_len = min(25, max(6, nchar(seq))),
        min_loop = 1,
        max_loop = 50
      )
      triplex_count <- length(res)
      if (triplex_count > 0) {
        scores <- score(res)
        idx <- which.max(scores)
        best_score <- scores[[idx]]
        best_pvalue <- pvalue(res)[[idx]]
        best_type <- type(res)[[idx]]
        best_start <- start(res)[[idx]]
        best_end <- end(res)[[idx]]
      }
    }, error = function(e) {
      err <<- conditionMessage(e)
    })
  }

  rows[[i]] <- data.frame(
    sequence_id = sid,
    label = df$label[[i]],
    triplex_pred = triplex_count > 0,
    triplex_count = triplex_count,
    triplex_best_score = best_score,
    triplex_best_pvalue = best_pvalue,
    triplex_best_type = best_type,
    triplex_best_start = best_start,
    triplex_best_end = best_end,
    triplex_error = err,
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, rows)
write.csv(out, output_csv, row.names = FALSE)
