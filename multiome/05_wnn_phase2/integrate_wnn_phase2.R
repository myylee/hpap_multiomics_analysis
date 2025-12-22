#!/usr/bin/env Rscript
# ==============================================================================
# Phase 2: Final WNN Integration
#
# After Phase 2 re-integration of GEX and ATAC, this script performs the
# final WNN integration using the re-integrated reductions.
#
# Usage: Rscript integrate_wnn_phase2.R <input_path> <filtered_object_path> [n_workers]
# ==============================================================================

library(Seurat)
library(stringr)
library(sctransform)
library(harmony)
library(tidyverse)
library(Signac)
library(future)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript integrate_wnn_phase2.R <input_path> <filtered_object_path> [n_workers]")
}

inPath <- args[1]
filteredObjectPath <- args[2]
n_workers <- ifelse(length(args) >= 3, as.numeric(args[3]), 4)

plan("multisession", workers = n_workers)
options(future.globals.maxSize = 200000 * 1024^2)  # 200GB

# ==============================================================================
# Step 1: Load filtered and re-integrated object
# ==============================================================================
cat("Loading filtered and re-integrated object...\n")
seurat <- readRDS(filteredObjectPath)

# ==============================================================================
# Step 2: Final WNN integration
# ==============================================================================
cat("Performing final WNN integration...\n")
rd_list <- list("combinedSctNormGEX.rpcaRefIntegrated.run2", "macs2Combined.harmony2")
key_i <- "5filteredRun2"

cat("Using reductions:", paste(rd_list, collapse = " + "), "\n")

# Find multimodal neighbors
seurat <- FindMultiModalNeighbors(
  seurat,
  reduction.list = rd_list,
  dims.list = list(1:30, 2:30),
  knn.graph.name = paste0("wknn", key_i),
  snn.graph.name = paste0("wsnn", key_i),
  weighted.nn.name = paste0("weighted.nn", key_i)
)

# Run UMAP (with model for potential projection)
seurat <- RunUMAP(
  seurat,
  nn.name = "weighted.nn5filteredRun2",
  reduction.name = "wnn5filteredRun2.umap",
  reduction.key = "wnn5filteredRun2Run2UMAP_",
  return.model = TRUE,
  seed.use = 1234
)

# ==============================================================================
# Step 3: Final clustering
# ==============================================================================
cat("Performing final clustering...\n")
seurat <- FindClusters(
  seurat,
  graph.name = paste0("wsnn", key_i),
  algorithm = 3,
  verbose = FALSE,
  cluster.name = paste0("clusters_wnn", key_i)
)

clusterCol <- paste0("clusters_wnn", key_i)
seurat@meta.data[, clusterCol] <- factor(
  as.character(seurat@meta.data[, clusterCol]),
  levels = as.character(sort(as.integer(levels(seurat@meta.data[, clusterCol]))))
)

# ==============================================================================
# Step 4: Visualization
# ==============================================================================
cat("Creating visualizations...\n")
pdf(paste0(inPath, '/wnn5filteredRun2_umap.pdf'), width = 12, height = 8)
p1 <- DimPlot(seurat, group.by = "donor", reduction = "wnn5filteredRun2.umap", label = FALSE)
p2 <- DimPlot(seurat, group.by = clusterCol, reduction = "wnn5filteredRun2.umap", label = TRUE) + NoLegend()
gridExtra::grid.arrange(grobs = list(p1, p2), ncol = 2)
dev.off()

# ==============================================================================
# Step 5: Save final object
# ==============================================================================
cat("Saving final integrated object...\n")
saveRDS(seurat, filteredObjectPath)

# Also save with a descriptive name
finalOutputPath <- paste0(inPath, '/multiome40DonorsRef-v2.rds')
saveRDS(seurat, finalOutputPath)

cat("Phase 2 final WNN integration complete!\n")
cat("Final object saved to:", finalOutputPath, "\n")
cat("Final UMAP reduction: wnn5filteredRun2.umap\n")

