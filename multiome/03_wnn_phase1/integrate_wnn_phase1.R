#!/usr/bin/env Rscript
# ==============================================================================
# Phase 1: WNN Integration
#
# This script performs Weighted Nearest Neighbor (WNN) integration combining
# GEX and ATAC modalities using different integration method combinations:
# 1. Harmony GEX + Harmony ATAC (wnn1)
# 2. Harmony GEX + RLSI ATAC (wnn2)
# 3. RPCA All GEX + Harmony ATAC (wnn3)
# 4. RPCA All GEX + RLSI ATAC (wnn4)
# 5. RPCA Ref GEX + Harmony ATAC (wnn5) - SELECTED FOR FINAL ANALYSIS
# 6. RPCA Ref GEX + RLSI ATAC (wnn6)
#
# WNN finds the optimal weighting between modalities to create a unified
# representation of the data.
#
# IMPORTANT FOR REPRODUCIBILITY:
# - All UMAP runs use seed.use = 1234 to ensure identical layouts
# - Clustering uses algorithm = 3 (Louvain) for consistency
# - Cluster IDs are sorted numerically before any annotation
#
# WORKFLOW:
# After running this script, you need to:
# 1. Manually annotate each WNN result by examining cluster quality plots
# 2. Compare WNN results to select the best integration (we selected wnn5)
# 3. Filter mixed/contaminated clusters from the selected WNN result
# 4. Proceed to Phase 2 re-integration on filtered cells
#
# See annotate_per_modality_clusters.R for per-modality cluster annotations.
# See annotate_wnn_filter_phase2.R for WNN5 annotation and filtering workflow.
#
# Usage: Rscript integrate_wnn_phase1.R <input_path> [n_workers]
# ==============================================================================

library(Seurat)
library(stringr)
library(sctransform)
library(harmony)
library(tidyverse)
library(data.table)
library(Signac)
library(future)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript integrate_wnn_phase1.R <input_path> [n_workers]")
}

inPath <- args[1]
n_workers <- ifelse(length(args) >= 2, as.numeric(args[2]), 4)

plan("multisession", workers = n_workers)
options(future.globals.maxSize = 200000 * 1024^2)  # 200GB

# ==============================================================================
# Step 1: Load combined GEX and ATAC object
# ==============================================================================
cat("Loading combined GEX and ATAC object...\n")
seurat <- readRDS(paste0(inPath, '/seurat_gexAndAtac.rds'))

# ==============================================================================
# Step 2: Define WNN integration combinations
# ==============================================================================
# Different combinations of GEX and ATAC integration methods
wnn_integration_list <- list(
  "1" = c("combinedSctNormGEX.harmony", "macs2Combined.harmony"),
  "2" = list("combinedSctNormGEX.harmony", "macs2Combined.rlsiRefIntegrated"),
  "3" = list("combinedSctNormGEX.rpcaAllIntegrated", "macs2Combined.harmony"),
  "4" = list("combinedSctNormGEX.rpcaAllIntegrated", "macs2Combined.rlsiRefIntegrated"),
  "5" = list("combinedSctNormGEX.rpcaRefIntegrated", "macs2Combined.harmony"),
  "6" = list("combinedSctNormGEX.rpcaRefIntegrated", "macs2Combined.rlsiRefIntegrated")
)

# ==============================================================================
# Step 3: Perform WNN integration for each combination
# ==============================================================================
cat("Performing WNN integration...\n")
for (i in 2:length(wnn_integration_list)) {
  rd_list <- wnn_integration_list[[i]]
  key_i <- names(wnn_integration_list)[i]
  
  cat("WNN integration", key_i, ":", 
      paste(rd_list, collapse = " + "), "\n")
  
  # Find multimodal neighbors
  seurat <- FindMultiModalNeighbors(
    seurat,
    reduction.list = rd_list,
    dims.list = list(1:30, 2:30),
    knn.graph.name = paste0("wknn", key_i),
    snn.graph.name = paste0("wsnn", key_i),
    weighted.nn.name = paste0("weighted.nn", key_i)
  )
  
  # Run UMAP
  # CRITICAL: Use seed.use = 1234 for reproducibility
  # Changing this seed will result in different UMAP layouts and may affect
  # cluster assignments, making it difficult to reproduce published results.
  seurat <- RunUMAP(
    seurat,
    nn.name = paste0("weighted.nn", key_i),
    reduction.name = paste0("wnn", key_i, ".umap"),
    reduction.key = paste0("wnn", key_i, "UMAP_"),
    seed.use = 1234  # Fixed seed for reproducibility
  )
  
  # Find clusters
  seurat <- FindClusters(
    seurat,
    graph.name = paste0("wsnn", key_i),
    algorithm = 3,
    verbose = FALSE,
    cluster.name = paste0("clusters_wnn", key_i)
  )
  
  # Ensure cluster names are properly formatted
  clusterCol <- paste0("clusters_wnn", key_i)
  seurat@meta.data[, clusterCol] <- factor(
    as.character(seurat@meta.data[, clusterCol]),
    levels = as.character(sort(as.integer(levels(seurat@meta.data[, clusterCol]))))
  )
  
  # Visualize
  pdf(paste0(inPath, '/wnn', key_i, '_umap.pdf'), width = 12, height = 8)
  p1 <- DimPlot(
    seurat,
    group.by = "donor",
    reduction = paste0("wnn", key_i, ".umap"),
    label = FALSE
  ) + NoLegend()
  p2 <- DimPlot(
    seurat,
    group.by = clusterCol,
    reduction = paste0("wnn", key_i, ".umap"),
    label = TRUE
  ) + NoLegend()
  gridExtra::grid.arrange(grobs = list(p1, p2), ncol = 2)
  dev.off()
  
  cat("  Completed WNN", key_i, "\n")
}

# ==============================================================================
# Step 4: Save object
# ==============================================================================
cat("Saving object with WNN integrations...\n")
saveRDS(seurat, paste0(inPath, '/seurat_gexAndAtac.rds'))

cat("\nPhase 1 WNN integration complete!\n")
cat("WNN integrations performed:", paste(names(wnn_integration_list)[-1], collapse = ", "), "\n")
cat("\nNext steps:\n")
cat("1. Generate cluster quality plots for each WNN result\n")
cat("2. Manually annotate WNN clusters by examining marker expression\n")
cat("3. Compare WNN results and select the best integration\n")
cat("4. Filter mixed/contaminated clusters from selected WNN result\n")
cat("5. Run annotate_wnn_filter_phase2.R to annotate and filter wnn5\n")
cat("6. Proceed to Phase 2 re-integration on filtered cells\n")
cat("\nIMPORTANT: For reproducibility, ensure all downstream steps also use seed.use = 1234!\n")

