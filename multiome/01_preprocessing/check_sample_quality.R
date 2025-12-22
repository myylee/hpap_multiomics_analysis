#!/usr/bin/env Rscript
# ==============================================================================
# Sample Quality Control Summary
#
# This script aggregates quality control metrics from all processed samples
# and creates summary tables for:
# 1. GEX QC metrics per sample
# 2. ATAC QC metrics per sample
# 3. Combined multiome metrics (cells with both GEX and ATAC data)
#
# This summary is used to:
# - Evaluate sample quality across all donors
# - Identify high-quality samples suitable for use as reference samples
#   in reference-based integration (RPCA reference-based GEX and RLSI ATAC)
# - Make decisions about which samples to include in downstream analysis
#
# The output file (sampleQC_summary.csv) contains per-sample metrics that
# can be used to select reference samples (typically high-quality control samples
# with good QC metrics across both GEX and ATAC modalities).
#
# Usage: Rscript check_sample_quality.R <input_path>
# ==============================================================================

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript check_sample_quality.R <input_path>")
}

inPath <- args[1]

# ==============================================================================
# Step 1: Load metadata from all samples
# ==============================================================================
cat("Loading metadata from all samples...\n")
gex_md_paths <- list.files(inPath, "metadata_sctNormGEX.csv", recursive = TRUE)
atac_md_paths <- list.files(inPath, "metadata_peakMacs2ATAC.csv", recursive = TRUE)

gex_md_list <- lapply(lapply(gex_md_paths, fread), as.data.frame)
atac_md_list <- lapply(lapply(atac_md_paths, fread), as.data.frame)

# ==============================================================================
# Step 2: Combine GEX metadata
# ==============================================================================
cat("Combining GEX metadata...\n")
gex_common_columns <- Reduce(intersect, lapply(gex_md_list, colnames))
gex_md_combined <- do.call(
  rbind,
  lapply(gex_md_list, function(md) { md[, gex_common_columns] })
)
rownames(gex_md_combined) <- gex_md_combined$V1

# ==============================================================================
# Step 3: Combine ATAC metadata
# ==============================================================================
cat("Combining ATAC metadata...\n")
atac_common_columns <- Reduce(intersect, lapply(atac_md_list, colnames))
atac_md_combined <- do.call(
  rbind,
  lapply(atac_md_list, function(md) { md[, atac_common_columns] })
)
rownames(atac_md_combined) <- atac_md_combined$V1

# ==============================================================================
# Step 4: Calculate per-sample GEX QC metrics
# ==============================================================================
cat("Calculating per-sample QC metrics...\n")
sample_ncells <- gex_md_combined %>%
  group_by(donor) %>%
  summarize(
    n_gex = n(),
    median_nCountRNA = median(nCount_RNA),
    median_nFeatureRNA = median(nFeature_RNA),
    median_percentMt = median(percent.mt)
  ) %>%
  full_join(
    atac_md_combined %>%
      group_by(donor) %>%
      summarize(
        n_atac = n(),
        median_fragCount = median(frequency_count),
        median_nCountPeakMacs2 = median(nCount_peak_macs2),
        median_tssEnrichment = median(TSS.enrichment),
        median_nucleosomeSignal = median(nucleosome_signal)
      )
  ) %>%
  arrange(donor)

# ==============================================================================
# Step 5: Combine multiome metadata and calculate metrics
# ==============================================================================
cat("Identifying multiome cells...\n")
bc_common <- intersect(rownames(gex_md_combined), rownames(atac_md_combined))
duplicated_col <- intersect(colnames(gex_md_combined), colnames(atac_md_combined))

multiome_md_combined <- cbind(
  gex_md_combined[bc_common, ],
  atac_md_combined[bc_common, setdiff(colnames(atac_md_combined), duplicated_col)]
)

# Add multiome metrics to summary
sample_ncells <- full_join(
  sample_ncells,
  multiome_md_combined %>%
    group_by(donor) %>%
    summarize(
      n_multiome = n(),
      multiome_median_nCountRNA = median(nCount_RNA),
      multiome_median_nFeatureRNA = median(nFeature_RNA),
      multiome_median_percentMt = median(percent.mt),
      multiome_median_fragCount = median(frequency_count),
      multiome_median_nCountPeakMacs2 = median(nCount_peak_macs2),
      multiome_median_tssEnrichment = median(TSS.enrichment),
      multiome_median_nucleosomeSignal = median(nucleosome_signal)
    )
)

# ==============================================================================
# Step 6: Save results
# ==============================================================================
cat("Saving QC summary...\n")
fwrite(sample_ncells, file = paste0(inPath, '/sampleQC_summary.csv'), row.names = TRUE, quote = FALSE)
fwrite(multiome_md_combined, file = paste0(inPath, '/multiome_md_combined.csv'), row.names = TRUE, quote = FALSE)

cat("QC summary saved to:", paste0(inPath, '/sampleQC_summary.csv'), "\n")
cat("Combined metadata saved to:", paste0(inPath, '/multiome_md_combined.csv'), "\n")

cat("\nSummary:\n")
cat("  Total samples:", nrow(sample_ncells), "\n")
cat("  Total GEX cells:", sum(sample_ncells$n_gex, na.rm = TRUE), "\n")
cat("  Total ATAC cells:", sum(sample_ncells$n_atac, na.rm = TRUE), "\n")
cat("  Total multiome cells:", sum(sample_ncells$n_multiome, na.rm = TRUE), "\n")

