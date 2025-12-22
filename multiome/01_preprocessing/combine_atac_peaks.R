#!/usr/bin/env Rscript
# ==============================================================================
# Combine ATAC peaks across all samples
#
# This script combines MACS2-called peaks from all individual samples into
# a unified peak set for downstream analysis. Peaks are filtered to remove
# blacklist regions and non-standard chromosomes.
#
# This should be run after all individual samples have been processed with
# process_atac_sample.R
#
# Usage: Rscript combine_atac_peaks.R <input_path>
#   input_path: Path to directory containing processed samples (with peak_macs2.csv files)
# ==============================================================================

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript combine_atac_peaks.R <input_path>")
}

inPath <- args[1]

# ==============================================================================
# Step 1: Load peaks from all samples
# ==============================================================================
cat("Loading peaks from all samples...\n")
samples <- list.files(inPath, "peak_macs2.csv", recursive = TRUE)
peak_list <- c()

for (sample_i in samples) {
  cat("Processing:", sample_i, "\n")
  # Convert to genomic ranges
  peak_gr <- makeGRangesFromDataFrame(read.csv(file.path(inPath, sample_i)))
  peak_gr <- keepStandardChromosomes(peak_gr, pruning.mode = "coarse")
  peak_gr <- subsetByOverlaps(x = peak_gr, ranges = blacklist_hg38_unified, invert = TRUE)
  peak_list <- c(peak_gr, peak_list)
  cat("  Peaks:", length(peak_gr), "\n")
}

# ==============================================================================
# Step 2: Create unified peak set
# ==============================================================================
cat("Creating unified peak set...\n")
peaks_combined <- GenomicRanges::reduce(x = peak_list)
peakwidths <- width(peaks_combined)

cat("Total combined peaks:", length(peaks_combined), "\n")
cat("Peak width summary:\n")
print(summary(peakwidths))

# ==============================================================================
# Step 3: Save combined peaks
# ==============================================================================
write.csv(peaks_combined, paste0(inPath, "/peak_macs2_combined_across_samples.csv"))

pdf(paste0(inPath, "/histogram_peak_macs2_combined_across_samples.pdf"), width = 8, height = 11)
hist(peakwidths, main = "Distribution of Combined Peak Widths", xlab = "Peak Width (bp)")
dev.off()

cat("Combined peaks saved to:", paste0(inPath, "/peak_macs2_combined_across_samples.csv"), "\n")

