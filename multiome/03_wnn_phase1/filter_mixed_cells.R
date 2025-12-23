#!/usr/bin/env Rscript
# ==============================================================================
# Phase 1: Filter Mixed Cells for Phase 2 Re-integration
#
# This script performs the filtering step at the end of Phase 1:
# 1. Manual annotation of WNN5 clusters
# 2. Three-step filtering of mixed/contaminated clusters using:
#    - WNN5 annotations (wnn5Annot)
#    - ATAC RLSI reference integration clusters (atacRlsiRef)
#    - ATAC Harmony integration clusters (atacHarmony)
#
# IMPORTANT FOR REPRODUCIBILITY:
# - Use seed.use = 1234 for all UMAP runs to ensure identical layouts
# - Apply annotations in the exact order of sorted cluster IDs
# - Filter using the exact cluster annotations specified below
#
# WORKFLOW SUMMARY:
# After Phase 1 WNN integration, we evaluated multiple WNN results (wnn2-wnn6).
# Each WNN result was manually annotated by examining:
# - Marker gene expression patterns
# - Comparison with per-modality annotations (gexHarmony, gexRpcaRef, etc.)
# - Cluster quality metrics
# - UMAP visualizations
#
# After evaluation, WNN5 (RPCA Ref GEX + Harmony ATAC) was selected as the
# best integration based on cluster separation and biological coherence.
#
# The filtering uses a multi-modal approach:
# - Step 1: Filter based on WNN5 annotations (acinar_alpha_beta, alpha_beta,
#           acinar_ductal, lowQual)
# - Step 2: Additional filtering based on ATAC RLSI clusters (acinar_alpha_beta)
# - Step 3: Additional filtering based on ATAC Harmony clusters
#           (acinar_alpha_beta, alpha_beta)
#
# This ensures that cells showing mixed signals in any modality are removed.
# The filtered dataset then undergoes Phase 2 re-integration (see integrate_gex_phase2.R
# and integrate_atac_phase2.R), followed by subclustering and final annotation
# (see subcluster_and_annotate_final.R).
#
# Usage: Rscript filter_mixed_cells.R <input_path> <output_path> [n_workers]
# ==============================================================================

library(Seurat)
library(Signac)
library(tidyverse)
library(future)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript filter_mixed_cells.R <input_path> <output_path> [n_workers]")
}

inPath <- args[1]
outPath <- args[2]
n_workers <- ifelse(length(args) >= 3, as.numeric(args[3]), 4)

plan("multisession", workers = n_workers)
options(future.globals.maxSize = 200000 * 1024^2)  # 200GB

# ==============================================================================
# Step 1: Load object with per-modality annotations
# ==============================================================================
cat("Loading object with per-modality annotations...\n")
seurat <- readRDS(paste0(inPath, '/seurat_gexAndAtac.rds'))

# ==============================================================================
# Step 2: Manual Annotation of WNN5 Clusters
# ==============================================================================
# WNN5 uses: RPCA Ref GEX + Harmony ATAC reductions
#
# Annotation was performed by examining:
# - Marker gene expression (DotPlots, FeaturePlots)
# - Cluster quality metrics (VlnPlots of QC metrics)
# - Comparison with per-modality annotations
# - UMAP visualization
#
# IMPORTANT: Annotation vector must match the sorted numeric order of cluster IDs.
# Cluster IDs are sorted as 0, 1, 2, ..., n before applying annotations.

cat("Annotating WNN5 clusters...\n")
# WNN5 clusters (36 clusters: 0-35)
wnn5Annot <- c(
  "alpha", "acinar", "alpha", "ductal", "beta",
  "beta", "ductal", "activated.stellate", "acinar", "acinar_alpha_beta",
  "delta", "endothelial", "alpha", "beta", "ductal",
  "pp.gamma", "c11", "quiescent.stellate", "acinar_alpha_beta", "acinar_alpha_beta",
  "acinar_ductal", "alpha", "alpha", "acinar_alpha_beta", "immune",
  "ductal", "immune", "alpha_beta", "acinar_alpha_beta", "immune",
  "lowQual", "acinar_alpha_beta", "acinar_alpha_beta", "acinar_alpha_beta",
  "beta", "alpha"
)
seurat@meta.data$wnn5Annot <- seurat$clusters_wnn5
levels(seurat$wnn5Annot) <- wnn5Annot

cat("WNN5 cluster distribution:\n")
print(table(seurat$wnn5Annot))

# ==============================================================================
# Step 3: Filter Mixed/Contaminated Clusters (Three-Step Filtering)
# ==============================================================================
# Clusters identified as mixed or contaminated are filtered out before Phase 2.
# These clusters show mixed signals that may represent:
# - Doublets/multiplets
# - Technical artifacts
# - Transitional cell states that cannot be clearly assigned
#
# The filtering uses a three-step approach combining:
# 1. WNN5 annotations (wnn5Annot)
# 2. ATAC RLSI reference integration clusters (atacRlsiRef)
# 3. ATAC Harmony integration clusters (atacHarmony)
#
# This multi-modal filtering ensures that cells showing mixed signals in any
# modality are removed.

cat("Filtering mixed/contaminated clusters using three-step approach...\n")

# Step 1: Filter based on WNN5 annotations
clusters_to_filter_wnn5 <- c("acinar_alpha_beta", "alpha_beta", "acinar_ductal", "lowQual")
cat("Step 1: Filtering based on wnn5Annot...\n")
seurat_wnn5Filtered <- subset(
  seurat,
  subset = wnn5Annot %in% clusters_to_filter_wnn5,
  invert = TRUE
)
cat("  Cells after Step 1:", ncol(seurat_wnn5Filtered), "\n")

# Step 2: Additional filtering based on ATAC RLSI reference clusters
# Remove cells annotated as "acinar_alpha_beta" in ATAC RLSI integration
cat("Step 2: Additional filtering based on atacRlsiRef clusters...\n")
if ("atacRlsiRef" %in% colnames(seurat_wnn5Filtered@meta.data)) {
  atac_rlsi_clusters_to_remove <- "acinar_alpha_beta"
  atac_rlsi_clusters_to_keep <- setdiff(unique(seurat_wnn5Filtered$atacRlsiRef), atac_rlsi_clusters_to_remove)
  seurat_wnn5Filtered <- subset(
    seurat_wnn5Filtered,
    subset = atacRlsiRef %in% atac_rlsi_clusters_to_keep
  )
  cat("  Cells after Step 2:", ncol(seurat_wnn5Filtered), "\n")
} else {
  warning("atacRlsiRef not found in metadata. Skipping Step 2.\n")
}

# Step 3: Additional filtering based on ATAC Harmony clusters
# Remove cells annotated as "acinar_alpha_beta" or "alpha_beta" in ATAC Harmony integration
cat("Step 3: Additional filtering based on atacHarmony clusters...\n")
if ("atacHarmony" %in% colnames(seurat_wnn5Filtered@meta.data)) {
  atac_harmony_clusters_to_remove <- c("acinar_alpha_beta", "alpha_beta")
  atac_harmony_clusters_to_keep <- setdiff(unique(seurat_wnn5Filtered$atacHarmony), atac_harmony_clusters_to_remove)
  seurat_wnn5Filtered <- subset(
    seurat_wnn5Filtered,
    subset = atacHarmony %in% atac_harmony_clusters_to_keep
  )
  cat("  Cells after Step 3:", ncol(seurat_wnn5Filtered), "\n")
} else {
  warning("atacHarmony not found in metadata. Skipping Step 3.\n")
}

cat("\nFinal filtering summary:\n")
cat("  Cells before filtering:", ncol(seurat), "\n")
cat("  Cells after filtering:", ncol(seurat_wnn5Filtered), "\n")
cat("  Cells filtered:", ncol(seurat) - ncol(seurat_wnn5Filtered), "\n")
cat("  Percentage retained:", round(100 * ncol(seurat_wnn5Filtered) / ncol(seurat), 2), "%\n")

# ==============================================================================
# Step 4: Save filtered object
# ==============================================================================
cat("Saving filtered object...\n")
seurat <- seurat_wnn5Filtered
saveRDS(seurat, paste0(outPath, '/seurat_wnn5Filtered.rds'))

cat("\nWNN5 annotation and filtering complete!\n")
cat("Filtered object saved to:", paste0(outPath, '/seurat_wnn5Filtered.rds'), "\n")
cat("\nNext steps for Phase 2 re-integration:\n")
cat("1. Run integrate_gex_phase2.R on the filtered object\n")
cat("2. Run integrate_atac_phase2.R on the filtered object\n")
cat("3. Run integrate_wnn_phase2.R to combine re-integrated modalities\n")
cat("\nIMPORTANT: Ensure all UMAP runs use seed.use = 1234 for reproducibility!\n")

