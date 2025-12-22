#!/usr/bin/env Rscript
# ==============================================================================
# Phase 1: GEX RPCA Integration (Reference-Based)
#
# Performs RPCA integration using selected reference samples.
# Reference samples are high-quality control samples.
# This is run after integrate_gex_phase1.R
#
# Usage: Rscript integrate_gex_rpca_ref.R <input_path> [npc] [n_workers] [reference_indices]
#   reference_indices: Comma-separated list of sample indices to use as reference
#                      Default: 5,6,9,16,21,27,29,31,32 (high-quality control samples)
# ==============================================================================

library(Seurat)
library(stringr)
library(sctransform)
library(harmony)
library(cowplot)
library(tidyverse)
library(data.table)
library(future)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript integrate_gex_rpca_ref.R <input_path> [npc] [n_workers] [reference_indices]")
}

inPath <- args[1]
npc <- ifelse(length(args) >= 2, as.numeric(args[2]), 30)
n_workers <- ifelse(length(args) >= 3, as.numeric(args[3]), 4)

# Parse reference indices
if (length(args) >= 4) {
  reference <- as.numeric(strsplit(args[4], ",")[[1]])
} else {
  # Default: high-quality control samples
  reference <- c(5, 6, 9, 16, 21, 27, 29, 31, 32)
}

plan("multisession", workers = n_workers)
options(future.globals.maxSize = 200000 * 1024^2)


# Load unintegrated object
cat("Loading unintegrated object...\n")
seurat_combined_v5 <- readRDS(paste0(inPath, '/seurat_sctNormGEX_unintegrated.rds'))

cat("Using reference samples:", reference, "\n")

# ==============================================================================
# RPCA Integration (Reference-Based)
# ==============================================================================
cat("Running RPCA integration (reference-based)...\n")
seurat_combined_v5 <- IntegrateLayers(
  object = seurat_combined_v5,
  method = RPCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "combinedSctNormGEX.pca",
  new.reduction = "combinedSctNormGEX.rpcaRefIntegrated",
  reference = reference,
  verbose = TRUE
)

seurat_combined_v5 <- RunUMAP(
  seurat_combined_v5, 
  reduction = "combinedSctNormGEX.rpcaRefIntegrated", 
  reduction.key = 'RPCREFASCTNORMGEXUMAP_',
  dims = 1:npc, 
  verbose = FALSE, 
  reduction.name = 'rpcaRefIntegratedSctNormGEX.umap'
)

seurat_combined_v5 <- FindNeighbors(
  seurat_combined_v5, 
  reduction = "combinedSctNormGEX.rpcaRefIntegrated", 
  dims = 1:npc, 
  verbose = FALSE, 
  graph.name = c('rpcaref.nn', 'rpcaref.snn')
)

seurat_combined_v5 <- FindClusters(
  seurat_combined_v5, 
  resolution = 0.5, 
  method = "igraph", 
  verbose = FALSE, 
  graph.name = 'rpcaref.snn',
  cluster.name = "combined_sct_rpcaRefIntegrated_clusters"
)

# Visualize
pdf(paste0(inPath, '/GEX_sctNormGEX_rpcaRefIntegrated_umap.pdf'), width = 11, height = 11)
p1 <- DimPlot(seurat_combined_v5, group.by = "donor", reduction = "rpcaRefIntegratedSctNormGEX.umap", label = FALSE)
p2 <- DimPlot(seurat_combined_v5, group.by = "combined_sct_rpcaRefIntegrated_clusters", 
              reduction = "rpcaRefIntegratedSctNormGEX.umap", label = TRUE) + NoLegend()
gridExtra::grid.arrange(grobs = list(p1, p2), nrow = 1)
dev.off()

# Save metadata and reductions
saveRDS(seurat_combined_v5@meta.data, paste0(inPath, '/metadata_sctNormGEX_rpcaRefIntegrated.rds'))
saveRDS(seurat_combined_v5@reductions, paste0(inPath, '/reduction_sctNormGEX_rpcaRefIntegrated.rds'))

cat("RPCA integration (reference-based) complete.\n")

