#!/usr/bin/env Rscript
# ==============================================================================
# Phase 1: GEX RPCA Integration (All Samples)
#
# Performs RPCA integration using all samples (no reference).
# This is run after integrate_gex_phase1.R
#
# Usage: Rscript integrate_gex_rpca_all.R <input_path> [npc] [n_workers]
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
  stop("Usage: Rscript integrate_gex_rpca_all.R <input_path> [npc] [n_workers]")
}

inPath <- args[1]
npc <- ifelse(length(args) >= 2, as.numeric(args[2]), 30)
n_workers <- ifelse(length(args) >= 3, as.numeric(args[3]), 4)

plan("multisession", workers = n_workers)
options(future.globals.maxSize = 200000 * 1024^2)


# Load unintegrated object
cat("Loading unintegrated object...\n")
seurat_combined_v5 <- readRDS(paste0(inPath, '/seurat_sctNormGEX_unintegrated.rds'))

# ==============================================================================
# RPCA Integration (All Samples)
# ==============================================================================
cat("Running RPCA integration (all samples)...\n")
seurat_combined_v5 <- IntegrateLayers(
  object = seurat_combined_v5,
  method = RPCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "combinedSctNormGEX.pca",
  new.reduction = "combinedSctNormGEX.rpcaAllIntegrated",
  verbose = TRUE
)

seurat_combined_v5 <- RunUMAP(
  seurat_combined_v5, 
  reduction = "combinedSctNormGEX.rpcaAllIntegrated", 
  reduction.key = 'RPCALLASCTNORMGEXUMAP_',
  dims = 1:npc, 
  verbose = FALSE, 
  reduction.name = 'rpcaAllIntegratedSctNormGEX.umap'
)

seurat_combined_v5 <- FindNeighbors(
  seurat_combined_v5, 
  reduction = "combinedSctNormGEX.rpcaAllIntegrated", 
  dims = 1:npc, 
  verbose = FALSE, 
  graph.name = c('rpcaall.nn', 'rpcaall.snn')
)

seurat_combined_v5 <- FindClusters(
  seurat_combined_v5, 
  resolution = 0.5, 
  method = "igraph", 
  verbose = FALSE, 
  graph.name = 'rpcaall.snn',
  cluster.name = "combined_sct_rpcaAllIntegrated_clusters"
)

# Visualize
pdf(paste0(inPath, '/GEX_sctNormGEX_rpcaAllIntegrated_umap.pdf'), width = 11, height = 11)
p1 <- DimPlot(seurat_combined_v5, group.by = "donor", reduction = "rpcaAllIntegratedSctNormGEX.umap", label = FALSE)
p2 <- DimPlot(seurat_combined_v5, group.by = "combined_sct_rpcaAllIntegrated_clusters", 
              reduction = "rpcaAllIntegratedSctNormGEX.umap", label = TRUE) + NoLegend()
gridExtra::grid.arrange(grobs = list(p1, p2), nrow = 1)
dev.off()

# Save metadata and reductions
saveRDS(seurat_combined_v5@meta.data, paste0(inPath, '/metadata_sctNormGEX_rpcaAllIntegrated.rds'))
saveRDS(seurat_combined_v5@reductions, paste0(inPath, '/reduction_sctNormGEX_rpcaAllIntegrated.rds'))

cat("RPCA integration (all samples) complete.\n")

