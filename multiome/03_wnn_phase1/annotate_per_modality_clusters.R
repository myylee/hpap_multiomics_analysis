#!/usr/bin/env Rscript
# ==============================================================================
# Phase 1: Manual Annotation of Per-Modality Clusters
#
# This script documents the manual annotation process for clusters from
# per-modality integrations (GEX and ATAC). These annotations are used to
# guide the evaluation of WNN integration results.
#
# IMPORTANT: Manual annotation is performed interactively by examining:
# 1. Marker gene expression (DotPlots and FeaturePlots)
# 2. Cluster quality metrics (VlnPlots)
# 3. UMAP visualizations
# 4. Comparison with per-sample scSorter annotations
#
# The annotation vectors provided here are the final annotations used in the
# published analysis. For reproducibility, ensure you:
# - Use the same seed (1234) for UMAP and clustering in integration steps
# - Examine the same cluster quality plots
# - Apply annotations in the same order as cluster IDs
#
# Usage: Rscript annotate_per_modality_clusters.R <input_path> <output_path>
# ==============================================================================

library(Seurat)
library(tidyverse)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript annotate_per_modality_clusters.R <input_path> <output_path>")
}

inPath <- args[1]
outPath <- args[2]

# ==============================================================================
# Step 1: Load combined object
# ==============================================================================
cat("Loading combined GEX and ATAC object...\n")
seurat <- readRDS(paste0(inPath, '/seurat_gexAndAtac.rds'))

# ==============================================================================
# Step 2: Manual Annotation Vectors
# ==============================================================================
# These annotation vectors map cluster IDs (in sorted numeric order) to
# cell type annotations. The annotations were determined by manual inspection
# of marker gene expression, cluster quality metrics, and UMAP visualizations.
#
# Annotation key:
# - Pure cell types: alpha, beta, delta, gamma/pp.gamma, epsilon, acinar, ductal,
#                    endothelial, immune, activated.stellate, quiescent.stellate
# - Mixed populations: acinar_alpha_beta, acinar_ductal, alpha_beta, etc.
# - Unknown/low quality: unknown_c##, lowQual, c11 (specific cell type)
#
# CRITICAL: These annotations must be applied in the exact order of sorted
# cluster IDs. Cluster IDs are sorted numerically (0, 1, 2, ...) and the
# annotation vector index corresponds to this order.

# GEX Harmony integration clusters
# Cluster IDs: 0, 1, 2, ..., 22 (23 clusters)
harmonyAnnot <- c(
  "alpha", "ductal", "acinar", "acinar_alpha_beta", "beta",
  "activated.stellate", "delta", "endothelial", "c11", "beta",
  "acinar", "ductal", "pp.gamma", "quiescent.stellate", "acinar_ductal",
  "immune", "alpha", "acinar", "alpha", "unknown_c19",
  "beta", "activated.stellate", "alpha"
)
seurat@meta.data$gexHarmony <- seurat$combined_sct_harmony_clusters
levels(seurat$gexHarmony) <- harmonyAnnot

# GEX RPCA Reference integration clusters
# Cluster IDs: 0, 1, 2, ..., 18 (19 clusters)
rpcaRefAnnot <- c(
  "alpha", "acinar", "ductal", "alpha", "beta",
  "acinar", "activated.stellate", "beta", "delta", "c11",
  "endothelial", "beta", "acinar_ductal", "ductal", "pp.gamma",
  "acinar_alpha_beta", "quiescent.stellate", "alpha_c17", "immune"
)
seurat@meta.data$gexRpcaRef <- seurat$combined_sct_rpcaRefIntegrated_clusters
levels(seurat$gexRpcaRef) <- rpcaRefAnnot

# GEX RPCA All integration clusters
# Cluster IDs: 0, 1, 2, ..., 21 (22 clusters)
rpcaAllAnnot <- c(
  "alpha", "ductal", "acinar", "alpha", "acinar",
  "beta", "alpha", "beta", "activated.stellate", "delta",
  "c11", "endothelial", "beta", "acinar_alpha_beta", "ductal",
  "pp.gamma", "acinar_ductal", "alpha", "quiescent.stellate",
  "immune", "unknown_c20", "unknown_c21"
)
seurat@meta.data$gexRpcaAll <- seurat$combined_sct_rpcaAllIntegrated_clusters
levels(seurat$gexRpcaAll) <- rpcaAllAnnot

# ATAC Harmony integration clusters
# Cluster IDs: 0, 1, 2, ..., 38 (39 clusters)
atacHarmonyAnnot <- c(
  "alpha_beta", "alpha", "ductal", "acinar", "beta",
  "alpha", "acinar", "alpha", "alpha", "beta",
  "alpha", "beta", "ductal", "delta", "endothelial",
  "beta", "ductal", "acinar", "activated.stellate", "ductal",
  "ductal", "pp.gamma", "activated.stellate", "acinar_alpha_beta", "acinar",
  "c11", "acinar", "quiescent.stellate", "c13", "acinar",
  "acinar", "activated.stellate", "acinar", "immune", "ductal",
  "delta", "ductal_c36", "acinar_alpha_beta", "acinar_alpha_beta", "acinar_alpha_beta"
)
seurat@meta.data$atacHarmony <- seurat$macs2Combined_harmony_clusters
levels(seurat$atacHarmony) <- atacHarmonyAnnot

# ATAC RLSI Reference integration clusters
# Cluster IDs: 0, 1, 2, ..., 15 (16 clusters)
rlsiRefAnnot <- c(
  "alpha", "alpha", "ductal", "beta", "acinar",
  "acinar", "beta", "delta", "activated.stellate", "endothelial",
  "c11", "ductal", "ductal", "acinar_alpha_beta", "quiescent.stellate",
  "immune"
)
seurat@meta.data$atacRlsiRef <- seurat$macs2Combined_rlsiRefIntegrated_clusters
levels(seurat$atacRlsiRef) <- rlsiRefAnnot

# ==============================================================================
# Step 3: Save annotated object
# ==============================================================================
cat("Saving annotated object...\n")
saveRDS(seurat, paste0(outPath, '/seurat_gexAndAtac_annotated.rds'))

# Also update the main object
saveRDS(seurat, paste0(inPath, '/seurat_gexAndAtac.rds'))

cat("Per-modality cluster annotation complete!\n")
cat("Annotated columns: gexHarmony, gexRpcaRef, gexRpcaAll, atacHarmony, atacRlsiRef\n")
cat("Next step: Evaluate WNN integration results using these annotations as reference.\n")

