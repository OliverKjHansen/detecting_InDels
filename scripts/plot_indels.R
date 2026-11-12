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

# ---------- Workflow 3's own columns (ticket 49) ----------
# read.delim turns True/False into TRUE/FALSE: a check passes if it is TRUE, or
# NA (it couldn't be made).
passes_or_na <- function(x) is.na(x) | toupper(as.character(x)) == "TRUE"
has <- function(...) all(c(...) %in% colnames(indels))
if (has("indel_string")) indels$size <- factor(nchar(indels$indel_string))

# ---------- Plot 9: kmerQ by position ----------
if (has("inner_flank_q", "outer_flank_q")) {
  levels_q <- c("<60", "60-93", ">=93")
  by_flank <- indels %>%
    mutate(inner_flank_q = factor(inner_flank_q, levels_q), outer_flank_q = factor(outer_flank_q, levels_q)) %>%
    count(indel_type, inner_flank_q, outer_flank_q)
  p9 <- ggplot(by_flank, aes(x = outer_flank_q, y = inner_flank_q, fill = n)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = n), size = 3) +
    scale_fill_gradient(low = "#cde2fb", high = "#184f95", trans = "log10") +
    facet_wrap(~ indel_type) +
    labs(title = "kmerQ by position: the lowest quality of L1, R1 (inner) and of L2, R2 (outer)",
         subtitle = "The proposed rule keeps inner >=93 and outer 60-93 or >=93 (kmerq_results.md)",
         x = "Outer flanks (L2, R2)", y = "Inner flanks (L1, R1)", fill = "Calls") +
    theme_minimal()
  print(p9)
}

# ---------- Plot 10: an insertion's inserted bases ----------
if (has("inserted_min_quality", "size")) {
  inserted <- indels %>% filter(indel_type == "insertion", !is.na(inserted_min_quality))
  if (nrow(inserted) > 0) {
    p10 <- ggplot(inserted, aes(x = as.numeric(inserted_min_quality))) +
      geom_histogram(binwidth = 1, fill = "#2a78d6") +
      facet_wrap(~ size, scales = "free_y", labeller = label_both) +
      geom_vline(xintercept = 93, linetype = "dashed", color = "red") +
      labs(title = "Insertion calls: the lowest base quality among the inserted bases",
           subtitle = "Final filtering keeps the calls at 93 (calls only)",
           x = "Lowest inserted base quality", y = "Calls") +
      theme_minimal()
    print(p10)
  }
}

# ---------- Plot 11: subreads disagreeing at the junction ----------
if (has("junction_disagreement", "size")) {
  p11 <- ggplot(indels, aes(x = size, fill = addNA(factor(junction_disagreement)))) +
    geom_bar(position = "fill") +
    facet_wrap(~ indel_type) +
    labs(title = "A subread that disagrees with the read at L1 or R1",
         x = "Size (bp)", y = "Share of calls", fill = "Disagreeing") +
    theme_minimal()
  print(p11)
}

# ---------- Plot 12: the distance to the mask ----------
if (has("mask_distance", "size")) {
  distance_bins <- c("0", "1-5", "6-10", "11-20", "21-50", "51-100", "101-1000", "none within 1000")
  by_distance <- indels %>%
    mutate(d = suppressWarnings(as.numeric(mask_distance)),
           bin = factor(case_when(is.na(d) ~ "none within 1000", d == 0 ~ "0", d <= 5 ~ "1-5", d <= 10 ~ "6-10",
                                  d <= 20 ~ "11-20", d <= 50 ~ "21-50", d <= 100 ~ "51-100", TRUE ~ "101-1000"),
                        distance_bins))
  p12 <- ggplot(by_distance, aes(x = bin, fill = size)) +
    geom_bar(position = "dodge") +
    facet_wrap(~ indel_type, scales = "free_y") +
    labs(title = "Distance from the event (anchor to base after) to the nearest masked base",
         subtitle = "0: the base after a deletion is masked; final filtering drops those calls and their opportunities",
         x = "Distance (bp)", y = "Calls", fill = "Size (bp)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p12)
}

# ---------- Plot 13: deletions whose bases repeat next to them ----------
if (has("repeats_next", "size")) {
  deletions <- indels %>% filter(indel_type == "deletion")
  if (nrow(deletions) > 0) {
    p13 <- ggplot(deletions, aes(x = size, fill = addNA(factor(repeats_next)))) +
      geom_bar(position = "fill") +
      labs(title = "Deletion calls whose deleted bases repeat right before or after them (ACAC -> AC)",
           subtitle = "false_positive_deletions.md: 80% of 2 bp calls, against 42% of real 2 bp deletions",
           x = "Size (bp)", y = "Share of calls", fill = "Repeat next") +
      theme_minimal()
    print(p13)
  }
}

# ---------- Plot 14: the calls through the filters (calls only) ----------
# Final filtering's steps for workflow 3 (config.yaml final_filters), then the
# candidates of false_positive_deletions.md. Counted on the calls only: final
# filtering applies the keys to the opportunities too, and its filter report
# counts both.
steps <- list(
  list("all calls", "", function(d) rep(TRUE, nrow(d))),
  list("+ unique", "unique", function(d) d$unique == "unique"),
  list("+ kmerQ by position", c("inner_flank_q", "outer_flank_q"),
       function(d) d$inner_flank_q == ">=93" & d$outer_flank_q %in% c("60-93", ">=93")),
  list("+ base after not masked", "base_after_masked", function(d) d$base_after_masked == "no"),
  list("+ inserted bases at Q93", "inserted_bases_pass",
       function(d) is.na(d$inserted_bases_pass) | d$inserted_bases_pass == "pass"),
  list("+ haploid check", "haploid_check", function(d) passes_or_na(d$haploid_check)),
  list("+ hap1/hap2 test", c("hap1_check", "hap2_check"),
       function(d) passes_or_na(d$hap1_check) & passes_or_na(d$hap2_check)),
  list("candidate: no disagreeing subread at L1, R1", "junction_disagreement",
       function(d) is.na(d$junction_disagreement) | d$junction_disagreement == "no"),
  list("candidate: read >= 8 kb", "read_length", function(d) as.numeric(d$read_length) >= 8000),
  list("candidate: 2 bp deletions don't repeat next", c("repeats_next", "size"),
       function(d) !(d$indel_type == "deletion" & d$size == "2" & d$repeats_next %in% "yes"))
)
if (has("indel_type", "size")) {
  keep <- rep(TRUE, nrow(indels))
  funnel <- list()
  for (step in steps) {
    if (!all(step[[2]] %in% c("", colnames(indels)))) {
      message("Skipping funnel step '", step[[1]], "': columns not found")
      next
    }
    keep <- keep & step[[3]](indels)
    keep[is.na(keep)] <- FALSE
    funnel[[length(funnel) + 1]] <- indels[keep, ] %>% count(indel_type, size) %>% mutate(step = step[[1]])
  }
  funnel <- bind_rows(funnel) %>% mutate(step = factor(step, unique(step)))
  p14 <- ggplot(funnel, aes(x = step, y = n, colour = size, group = size)) +
    geom_line() +
    geom_point() +
    scale_y_log10() +
    facet_wrap(~ indel_type) +
    labs(title = "Calls left after each filter (calls only)",
         subtitle = "Final filtering's steps, then the candidates of false_positive_deletions.md",
         x = NULL, y = "Calls (log scale)", colour = "Size (bp)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p14)
}

dev.off()
message("Plots saved to: ", out_pdf)
