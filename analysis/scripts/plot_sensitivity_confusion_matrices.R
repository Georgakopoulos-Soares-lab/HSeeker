#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(grid)
})

root <- normalizePath(getwd(), mustWork = TRUE)

read_metrics <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(df) <- df$method
  df
}

as_metrics <- function(metrics, method_key) {
  row <- metrics[method_key, , drop = FALSE]
  list(
    TP = as.integer(row$TP),
    FN = as.integer(row$FN),
    TN = as.integer(row$TN),
    FP = as.integer(row$FP),
    sensitivity = as.numeric(row$sensitivity),
    specificity = as.numeric(row$specificity),
    f1 = as.numeric(row$F1),
    mcc = as.numeric(row$MCC)
  )
}

cell_fill <- function(correct, count, max_count) {
  frac <- if (max_count > 0) count / max_count else 0
  if (correct) {
    # Green scale for correct classifications.
    rgb(
      red = 0.88 - 0.42 * frac,
      green = 0.96 - 0.20 * frac,
      blue = 0.88 - 0.42 * frac
    )
  } else {
    # Light red scale for classification errors.
    rgb(
      red = 1.00,
      green = 0.91 - 0.25 * frac,
      blue = 0.91 - 0.25 * frac
    )
  }
}

draw_label <- function(label, x, y, size = 10, fontface = "plain", rot = 0,
                       just = "centre", col = "#1A1A1A") {
  grid.text(
    label,
    x = unit(x, "in"),
    y = unit(y, "in"),
    gp = gpar(fontsize = size, fontface = fontface, col = col, fontfamily = "Helvetica"),
    rot = rot,
    just = just
  )
}

draw_cell <- function(x, y, size, count, total, correct, max_count) {
  grid.rect(
    x = unit(x + size / 2, "in"),
    y = unit(y + size / 2, "in"),
    width = unit(size * 0.96, "in"),
    height = unit(size * 0.96, "in"),
    gp = gpar(fill = cell_fill(correct, count, max_count), col = "white", lwd = 1.4)
  )
  pct <- if (total > 0) 100 * count / total else 0
  draw_label(as.character(count), x + size / 2, y + size * 0.60, size = 24, fontface = "bold")
  draw_label(sprintf("(%.1f%%)", pct), x + size / 2, y + size * 0.36, size = 11, col = "#333333")
}

draw_matrix <- function(x0, y0, title, metrics, show_y_title = FALSE) {
  cell <- 0.92
  total <- metrics$TP + metrics$FN + metrics$TN + metrics$FP
  max_count <- max(metrics$TP, metrics$FN, metrics$TN, metrics$FP)

  # Panel title and compact summary.
  draw_label(title, x0 + cell, y0 + 2.62, size = 15, fontface = "bold")
  draw_label(
    sprintf("Sensitivity %.3f   Specificity %.3f", metrics$sensitivity, metrics$specificity),
    x0 + cell, y0 + 2.39, size = 9.6, fontface = "bold", col = "#2A2A2A"
  )
  draw_label(
    sprintf("F1 %.3f   MCC %.3f", metrics$f1, metrics$mcc),
    x0 + cell, y0 + 2.20, size = 9.6, fontface = "bold", col = "#2A2A2A"
  )

  # Column labels: left/right = non-forming / H-DNA forming.
  draw_label("Non-forming", x0 + cell * 0.50, y0 + 1.93, size = 8.9)
  draw_label("H-DNA forming", x0 + cell * 1.50, y0 + 1.93, size = 8.9)

  # Row labels: top/bottom = non-forming / H-DNA forming.
  draw_label("Non-forming", x0 - 0.16, y0 + cell * 1.50, size = 8.9, just = "right")
  draw_label("H-DNA forming", x0 - 0.16, y0 + cell * 0.50, size = 8.9, just = "right")

  # Optional y-axis title for the left matrix.
  if (show_y_title) {
    draw_label("Experimental class", x0 - 1.08, y0 + cell, size = 11.5, fontface = "bold", rot = 90)
  }

  # Top row: actual non-forming. Bottom row: actual H-DNA forming.
  # Left column: predicted non-forming. Right column: predicted H-DNA forming.
  draw_cell(x0, y0 + cell, cell, metrics$TN, total, TRUE, max_count)
  draw_cell(x0 + cell, y0 + cell, cell, metrics$FP, total, FALSE, max_count)
  draw_cell(x0, y0, cell, metrics$FN, total, FALSE, max_count)
  draw_cell(x0 + cell, y0, cell, metrics$TP, total, TRUE, max_count)

  # Border and axis title for each matrix.
  grid.rect(
    x = unit(x0 + cell, "in"),
    y = unit(y0 + cell, "in"),
    width = unit(cell * 1.92, "in"),
    height = unit(cell * 1.92, "in"),
    gp = gpar(fill = NA, col = "#333333", lwd = 0.6)
  )
  draw_label("Predicted class", x0 + cell, y0 - 0.36, size = 12, fontface = "bold")
}

draw_legend <- function(x0, y0) {
  sw <- 0.22
  grid.rect(unit(x0, "in"), unit(y0, "in"), unit(sw, "in"), unit(sw, "in"),
            gp = gpar(fill = "#6CC070", col = NA))
  draw_label("Correct", x0 + 0.22, y0, size = 9.2, just = "left")
  grid.rect(unit(x0 + 0.95, "in"), unit(y0, "in"), unit(sw, "in"), unit(sw, "in"),
            gp = gpar(fill = "#F5B6B6", col = NA))
  draw_label("Incorrect", x0 + 1.17, y0, size = 9.2, just = "left")
}

plot_one <- function(metrics, hseeker_key, triplex_key, title, subtitle, out_pdf) {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  cairo_pdf(out_pdf, width = 8.7, height = 5.15, family = "Helvetica")
  grid.newpage()

  draw_label(title, 0.52, 4.82, size = 18, fontface = "bold", just = "left")
  draw_label(subtitle, 0.52, 4.55, size = 10.5, just = "left", col = "#3A3A3A")
  draw_legend(6.35, 4.62)

  draw_matrix(1.42, 1.20, "HSeeker", as_metrics(metrics, hseeker_key), show_y_title = TRUE)
  draw_matrix(5.12, 1.20, "Triplex", as_metrics(metrics, triplex_key), show_y_title = FALSE)

  dev.off()
}

direct_dir <- file.path(root, "analysis", "results", "sensitivity_direct")
injected_dir <- file.path(root, "analysis", "results", "sensitivity_injected")

plot_one(
  read_metrics(file.path(direct_dir, "summary_metrics.csv")),
  "HSeeker direct hit",
  "Triplex direct hit",
  "Direct sequence validation",
  "HSeeker and Triplex run directly on the 80 curated experimental sequences.",
  file.path(direct_dir, "plots", "confusion_matrices_hseeker_triplex.pdf")
)

plot_one(
  read_metrics(file.path(injected_dir, "summary_metrics.csv")),
  "HSeeker overlap",
  "Triplex overlap",
  "Injected E. coli validation",
  "Hits matched to inserted intervals using the 80% overlap rule.",
  file.path(injected_dir, "plots", "confusion_matrices_hseeker_triplex.pdf")
)

cat(file.path(direct_dir, "plots", "confusion_matrices_hseeker_triplex.pdf"), "\n")
cat(file.path(injected_dir, "plots", "confusion_matrices_hseeker_triplex.pdf"), "\n")
