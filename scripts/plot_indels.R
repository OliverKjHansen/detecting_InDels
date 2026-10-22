#!/usr/bin/env Rscript

# Plot InDel quality distributions, filter impacts, and spectra
# Usage: Rscript plot_indels.R <indels.tsv> <output.pdf>
# The plots are of the calls only: the opportunity tables are too big to load.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript plot_indels.R <indels.tsv> <output.pdf>")
}

indel_file <- args[1]
out_pdf <- args[2]

# Load data
indels <- read.delim(indel_file, sep = "\t", header = TRUE, comment.char = "#")

pdf(out_pdf, width = 12, height = 8)

# ---------- Plot 1: kmerQ distribution for indels ----------
if ("kmerQ" %in% colnames(indels)) {
  p1 <- ggplot(indels, aes(x = kmerQ, fill = indel_type)) +
    geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7) +
    labs(title = "kmerQ Distribution of Detected InDels",
         x = "kmerQ Score", y = "Count", fill = "InDel Type") +
    theme_minimal() +
    geom_vline(xintercept = 93, linetype = "dashed", color = "red") +
    annotate("text", x = 93, y = Inf, label = "Q93 threshold",
             vjust = 2, hjust = -0.1, color = "red", size = 3)
  print(p1)
} else {
  message("Skipping kmerQ plot: 'kmerQ' column not found")
}

# ---------- Plot 2: InDel spectrum (context counts) ----------
if ("indel_string" %in% colnames(indels) && "indel_type" %in% colnames(indels)) {
  spectrum <- indels %>%
    count(indel_type, indel_string) %>%
    arrange(desc(n))

  # Show top 20 most common motifs
  top_motifs <- spectrum %>% slice_head(n = 20)

  p2 <- ggplot(top_motifs, aes(x = reorder(indel_string, -n), y = n, fill = indel_type)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    labs(title = "InDel Spectrum (Top 20 Motifs)",
         x = "InDel Motif", y = "Count", fill = "InDel Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p2)
} else {
  message("Skipping InDel spectrum plot: 'indel_string' or 'indel_type' column not found")
}

# ---------- Plot 3: Distance to read end ----------
if ("ToreadEnd" %in% colnames(indels)) {
  p3 <- ggplot(indels, aes(x = ToreadEnd, fill = indel_type)) +
    geom_histogram(binwidth = 100, position = "dodge", alpha = 0.7) +
    labs(title = "Distance from InDel to Read End",
         x = "Distance to Read End (bp)", y = "Count", fill = "InDel Type") +
    theme_minimal() +
    geom_vline(xintercept = 100, linetype = "dashed", color = "red") +
    annotate("text", x = 100, y = Inf, label = "100bp threshold",
             vjust = 2, hjust = -0.1, color = "red", size = 3)
  print(p3)
} else {
  message("Skipping distance-to-read-end plot: 'ToreadEnd' column not found")
}

# ---------- Plot 4: Unique vs non-unique indels ----------
if ("unique" %in% colnames(indels)) {
  p4 <- ggplot(indels, aes(x = unique, fill = indel_type)) +
    geom_bar(position = "dodge", alpha = 0.8) +
    labs(title = "Unique vs Non-Unique InDels",
         x = "Uniqueness", y = "Count", fill = "InDel Type") +
    theme_minimal()
  print(p4)
} else {
  message("Skipping uniqueness plot: 'unique' column not found")
}

# ---------- Plot 5: Haploid check results ----------
if ("haploid_check" %in% colnames(indels)) {
  p5 <- ggplot(indels, aes(x = as.factor(haploid_check), fill = indel_type)) +
    geom_bar(position = "dodge", alpha = 0.8) +
    labs(title = "Haploid Concordance Check",
         x = "Haploid Check (True = concordant)", y = "Count", fill = "InDel Type") +
    theme_minimal()
  print(p5)
} else {
  message("Skipping haploid check plot: 'haploid_check' column not found")
}

# ---------- Plot 6: Tissue distribution ----------
if ("tissue" %in% colnames(indels)) {
  tissue_counts <- indels %>%
    count(tissue, indel_type) %>%
    filter(tissue != "NA")

  if (nrow(tissue_counts) > 0) {
    p6 <- ggplot(tissue_counts, aes(x = reorder(tissue, -n), y = n, fill = indel_type)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
      labs(title = "InDels by Tissue Type",
           x = "Tissue", y = "Count", fill = "InDel Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p6)
  }
} else {
  message("Skipping tissue plot: 'tissue' column not found")
}

# ---------- Plot 7: kmerQ vs haploid distance (if available) ----------
if ("kmerQ" %in% colnames(indels) && "distance_to_indel_in_hap" %in% colnames(indels)) {
  p7 <- ggplot(indels, aes(x = kmerQ, y = distance_to_indel_in_hap, color = indel_type)) +
    geom_point(alpha = 0.3, size = 0.5) +
    labs(title = "kmerQ vs Distance to InDel in Haploid Mapping",
         x = "kmerQ Score", y = "Distance to InDel in Haploid",
         color = "InDel Type") +
    theme_minimal() +
    geom_vline(xintercept = 93, linetype = "dashed", color = "red")
  print(p7)
} else {
  message("Skipping kmerQ vs haploid distance plot: required columns not found")
}

# ---------- Plot 8: Subread coverage (sa tag) distribution ----------
if ("sa" %in% colnames(indels)) {
  sa_data <- indels %>% filter(!is.na(sa) & sa != "NA")
  if (nrow(sa_data) > 0) {
    sa_data$sa <- as.numeric(sa_data$sa)
    p8 <- ggplot(sa_data, aes(x = sa, fill = indel_type)) +
      geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7) +
      labs(title = "Subread Coverage (sa tag) at InDel Positions",
           x = "Subread Coverage", y = "Count", fill = "InDel Type") +
      theme_minimal()
    print(p8)
  }
} else {
  message("Skipping subread coverage plot: 'sa' column not found")
}

dev.off()
message("Plots saved to: ", out_pdf)
