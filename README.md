# detecting_InDels

Pipeline for detecting small InDels (1-3 bp insertions and deletions) from PacBio HiFi long reads at a single-read level.

## Overview

This GWF workflow takes a de novo assembly `.fasta` and a `.bam` file with reads mapped back (self-alignment) and identifies small InDels using PySAM. 

## Pipeline Steps

1. **Build filter BED** — Merge microsatellite intervals and unresolved haplotype regions into a single filter BED file to mask unreliable regions.
2. **Detect small InDels** — Scan each contig for 1-3 bp insertions/deletions using PySAM pileup, applying coverage and homopolymer filters.
3. **Extract mutation FASTQ** — Extract reads supporting detected InDels into FASTQ format for haploid remapping.
4. **Count opportunities** — Count the number of non-indel positions (opportunities) per context, tissue, and quality score.
5. **Haploid remapping** — Remap mutation-supporting reads to the haploid assembly using pbmm2.
6. **Haploid concordance filter** — Check whether InDels detected in the diploid alignment are also present in the haploid mapping.
7. **Merge & plot** — Combine per-contig results and generate quality diagnostic plots using ggplot2 in R.
8. **Final filtering** — Apply quality filters (kmerQ, haploid concordance, distance-to-read-end) with a step-by-step report of rows removed at each stage.

## Filtering Approaches

- **Coverage filter**: Positions must have coverage between 0.5× and 2× the sample average.
- **Microsatellite/unresolved region masking**: InDels in tandem repeat or unresolved haplotype regions are excluded.
- **Homopolymer filter**: InDels in homopolymer runs (≤1 break in base identity) are excluded.
- **kmerQ quality filter**: Minimum base quality (Phred) at the 4-mer surrounding the InDel.
- **Haploid concordance**: InDels must be confirmed in the haploid mapping (matching type and size).
- **Distance to read end**: InDels too close to read ends are unreliable and can be filtered.(not implemented)
- **Uniqueness**: Whether the InDel is seen on one or multiple reads at that position.

## Outputs

- `output/final/<sample>_combined_filtered.tsv` — Filtered InDels
- `output/final/<sample>_opportunities_combined_filtered.tsv` — Filtered opportunities
- `output/final/<sample>_filter_report.tsv` — Step-by-step filter report
- `output/plots/<sample>_indel_plots.pdf` — Quality and spectrum diagnostic plots

## Requirements

- Python 3 with `pysam`, `pandas`, `intervaltree`
- R with `ggplot2`, `dplyr`, `tidyr`
- `samtools`, `bedtools`, `pbmm2`
- [GWF](https://gwf.app/) workflow system
