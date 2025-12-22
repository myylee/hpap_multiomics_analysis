#!/usr/bin/env Rscript
# ==============================================================================
# Phase 2: GEX Re-integration After Filtering
#
# After filtering low-quality cells based on Phase 1 WNN results, this script
# re-runs GEX integration using RPCA with reference samples.
#
# IMPORTANT FOR REPRODUCIBILITY:
# - UMAP uses seed.use = 1234 to ensure identical layouts
# - Clustering uses method = "igraph" (Louvain) for consistency
# - Reference sample indices must match the original analysis
#
# Usage: Rscript integrate_gex_phase2.R <input_path> <filtered_object_path> [npc] [n_workers] [reference_indices]
# ==============================================================================

library(Seurat)
library(stringr)
library(sctransform)
library(harmony)
library(tidyverse)
library(data.table)
library(future)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript integrate_gex_phase2.R <input_path> <filtered_object_path> [npc] [n_workers] [reference_indices]")
}

inPath <- args[1]
filteredObjectPath <- args[2]
npc <- ifelse(length(args) >= 3, as.numeric(args[3]), 30)
n_workers <- ifelse(length(args) >= 4, as.numeric(args[4]), 4)

# Parse reference indices
if (length(args) >= 5) {
  reference <- as.numeric(strsplit(args[5], ",")[[1]])
} else {
  # Default: high-quality control samples
  reference <- c(5, 6, 9, 16, 21, 27, 29, 31, 32)
}

plan("multisession", workers = n_workers)
options(future.globals.maxSize = 200000 * 1024^2)  # 200GB

# ==============================================================================
# Step 1: Load filtered object
# ==============================================================================
cat("Loading filtered object...\n")
seurat <- readRDS(filteredObjectPath)

# ==============================================================================
# Step 2: RPCA re-integration with reference
# ==============================================================================
cat("Running RPCA re-integration (reference-based)...\n")
cat("Using reference samples (indices):", reference, "\n")

seurat <- IntegrateLayers(
  object = seurat,
  method = RPCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "combinedSctNormGEX.pca",
  new.reduction = "combinedSctNormGEX.rpcaRefIntegrated.run2",
  reference = reference,
  verbose = TRUE
)

# ==============================================================================
# Step 3: UMAP, clustering, and visualization
# ==============================================================================
cat("Running UMAP and clustering...\n")
# CRITICAL: Use seed.use = 1234 for reproducibility
seurat <- RunUMAP(
  seurat,
  reduction = "combinedSctNormGEX.rpcaRefIntegrated.run2",
  reduction.key = 'RUN2RPCREFASCTNORMGEXUMAP_',
  dims = 1:npc,
  verbose = FALSE,
  reduction.name = 'rpcaRefIntegratedSctNormGEXRun2.umap',
  seed.use = 1234  # Fixed seed for reproducibility
)

seurat <- FindNeighbors(
  seurat,
  reduction = "combinedSctNormGEX.rpcaRefIntegrated.run2",
  dims = 1:npc,
  verbose = FALSE,
  graph.name = c('rpcarefrun2.nn', 'rpcarefrun2.snn')
)

seurat <- FindClusters(
  seurat,
  resolution = 0.5,
  method = "igraph",
  verbose = FALSE,
  graph.name = 'rpcarefrun2.snn',
  cluster.name = "combined_sct_rpcaRefIntegratedRun2_clusters"
)

# Visualize
pdf(paste0(inPath, '/GEX_rpcaRefIntegrated_run2_umap.pdf'), width = 11, height = 11)
p1 <- DimPlot(seurat, group.by = "donor", reduction = "rpcaRefIntegratedSctNormGEXRun2.umap", label = FALSE)
p2 <- DimPlot(seurat, group.by = "combined_sct_rpcaRefIntegratedRun2_clusters", 
              reduction = "rpcaRefIntegratedSctNormGEXRun2.umap", label = TRUE) + NoLegend()
gridExtra::grid.arrange(grobs = list(p1, p2), nrow = 1)
dev.off()

# ==============================================================================
# Step 4: Save metadata and reductions
# ==============================================================================
cat("Saving metadata and reductions...\n")
saveRDS(seurat@reductions, paste0(inPath, '/reduction_sctNormGEX_wnn5Filtered_harmonyIntegrated_run2.rds'))
saveRDS(seurat@meta.data, paste0(inPath, '/metadata_sctNormGEX_wnn5Filtered_harmonyIntegrated_run2.rds'))

# Save updated object
saveRDS(seurat, filteredObjectPath)

cat("Phase 2 GEX re-integration complete!\n")

