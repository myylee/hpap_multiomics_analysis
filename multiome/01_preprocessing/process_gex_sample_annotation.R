#!/usr/bin/env Rscript
# ==============================================================================
# Complete GEX sample processing
#
# This script continues from process_gex_sample.R:
# 1. Runs UMAP and clustering
# 2. Adds donor metadata
# 3. Saves final processed object
#
# Usage: Rscript process_gex_sample_annotation.R <path> [npc]
# ==============================================================================

library(tidyverse)
library(Seurat)
library(sctransform)
library(fs)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript process_gex_sample_annotation.R <path> [npc]")
}

path <- args[1]
npc <- ifelse(length(args) >= 2, as.numeric(args[2]), 30)

# Load intermediate object
seurat <- readRDS(paste0(path, '/seurat_GEX.tmp.rds'))

# ==============================================================================
# Step 1: UMAP and clustering
# ==============================================================================
cat("Running UMAP and clustering...\n")
seurat <- RunUMAP(
  seurat, 
  reduction = "sctNormGEX.pca", 
  reduction.key = 'SCTNORMGEXUMAP_',
  dims = 1:npc, 
  verbose = FALSE, 
  reduction.name = 'sctNormGEX.umap'
)

seurat <- FindNeighbors(
  seurat, 
  reduction = "sctNormGEX.pca", 
  dims = 1:npc, 
  verbose = FALSE, 
  graph.name = c('sct.nn', 'sct.snn')
)

seurat <- FindClusters(
  seurat, 
  resolution = 0.5, 
  method = "igraph", 
  verbose = FALSE, 
  graph.name = 'sct.snn'
)
seurat$sct.clusters <- seurat$seurat_clusters

# ==============================================================================
# Step 2: Visualization
# ==============================================================================
cat("Creating visualization plots...\n")
pdf(paste0(path, '/GEX_sctNormGEX_umap.pdf'))
d <- DimPlot(
  seurat, 
  group.by = "sct.clusters", 
  reduction = "sctNormGEX.umap", 
  label = TRUE
)
print(d)
dev.off()

cat("Cluster distribution:\n")
print(table(seurat$sct.clusters))

# ==============================================================================
# Step 3: Add metadata and rename cells
# ==============================================================================
donorFull <- path_file(path)
arr <- str_split(donorFull, "-")[[1]]
donor <- arr[1]

seurat@meta.data$donor <- donor
seurat@meta.data$modality <- 'multiomeRNA'
seurat <- RenameCells(seurat, new.names = paste0(donor, "_", colnames(seurat)))

# ==============================================================================
# Step 4: Save final object and metadata
# ==============================================================================
saveRDS(seurat, file = paste0(path, '/seurat_sctNormGEX.rds'))
data.table::fwrite(
  seurat@meta.data, 
  file = paste0(path, '/metadata_sctNormGEX.csv'),
  row.names = TRUE,
  quote = FALSE
)

cat("GEX sample processing complete!\n")

