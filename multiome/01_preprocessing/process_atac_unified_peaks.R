#!/usr/bin/env Rscript
# ==============================================================================
# Create Unified Peak Set Assay for ATAC Sample
#
# This script creates the macs2_combined assay by quantifying fragments
# against the unified peak set (from combine_atac_peaks.R). This is essential
# for integration because all samples must use the same peak set.
#
# Workflow:
# 1. Sample-specific peaks (peak_macs2) - created in process_atac_sample.R
# 2. Unified peaks (peak_macs2_combined_across_samples.csv) - created in combine_atac_peaks.R
# 3. This script: Quantify each sample against unified peaks → macs2_combined assay
#
# Usage: Rscript process_atac_unified_peaks.R <path> <npc> <unified_peaks_path>
#   path: Path to sample directory (contains seurat_ATAC.rds)
#   npc: Number of principal components for LSI
#   unified_peaks_path: Path to unified peaks CSV file (from combine_atac_peaks.R)
# ==============================================================================

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(fs)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript process_atac_unified_peaks.R <path> <npc> <unified_peaks_path>")
}

path <- args[1]
npc <- as.numeric(args[2])
peaks <- args[3]

# ==============================================================================
# Step 1: Load sample object and unified peaks
# ==============================================================================
cat("Loading sample object...\n")
seurat <- readRDS(paste0(path, '/seurat_ATAC.rds'))

cat("Loading unified peaks from:", peaks, "\n")
peak_gr <- makeGRangesFromDataFrame(read.csv(peaks))

# ==============================================================================
# Step 2: Quantify fragments against unified peaks
# ==============================================================================
cat("Quantifying fragments against unified peaks...\n")
DefaultAssay(seurat) <- 'ATAC'

macs2Combined_counts <- FeatureMatrix(
  fragments = Fragments(seurat),
  features = peak_gr,
  cells = rownames(seurat@meta.data)
)

# ==============================================================================
# Step 3: Create macs2_combined assay
# ==============================================================================
cat("Creating macs2_combined assay...\n")
chrom_assay <- CreateChromatinAssay(
  counts = macs2Combined_counts,
  sep = c("-", "-"),
  fragments = Fragments(seurat),
  annotation = Annotation(seurat)
)

seurat[['macs2_combined']] <- chrom_assay

# ==============================================================================
# Step 4: LSI dimensionality reduction on unified peaks
# ==============================================================================
cat("Running LSI on unified peaks...\n")
DefaultAssay(seurat) <- "macs2_combined"
seurat <- seurat[rowSums(seurat) > 0, ]
seurat <- RunTFIDF(seurat)
seurat <- FindTopFeatures(seurat, min.cutoff = 5)
seurat <- RunSVD(
  seurat,
  assay = 'macs2_combined',
  reduction.name = 'macs2Combined.lsi',
  verbose = FALSE,
  reduction.key = 'MACS2COMBINEDLSI_'
)

# ==============================================================================
# Step 5: Visualize LSI
# ==============================================================================
pdf(paste0(path, '/ATAC_lsi_elbow_macs2Combined.pdf'))
e <- ElbowPlot(seurat, ndims = 50, reduction = 'macs2Combined.lsi')
print(e)
dev.off()

pdf(paste0(path, '/ATAC_lsi_depthCor_macs2Combined.pdf'))
DepthCor(seurat, reduction = 'macs2Combined.lsi')
dev.off()

# ==============================================================================
# Step 6: UMAP and clustering
# ==============================================================================
cat("Running UMAP and clustering...\n")
seurat <- RunUMAP(
  seurat,
  reduction = "macs2Combined.lsi",
  reduction.key = 'MACS2COMBINEDUMAP_',
  dims = 2:npc,
  verbose = FALSE,
  reduction.name = 'macs2Combined.umap'
)

seurat <- FindNeighbors(
  seurat,
  reduction = "macs2Combined.lsi",
  dims = 1:npc,
  verbose = FALSE,
  graph.name = c('macs2Combined.nn', 'macs2Combined.snn')
)

seurat <- FindClusters(
  seurat,
  resolution = 0.5,
  method = "igraph",
  verbose = FALSE,
  graph.name = 'macs2Combined.snn'
)
seurat$macs2Combined.clusters <- seurat$seurat_clusters

pdf(paste0(path, '/ATAC_macs2Combined_umap.pdf'))
d <- DimPlot(
  seurat,
  group.by = "macs2Combined.clusters",
  reduction = "macs2Combined.umap",
  label = TRUE
)
print(d)
dev.off()

# ==============================================================================
# Step 7: Save updated object
# ==============================================================================
cat("Saving updated object with macs2_combined assay...\n")
saveRDS(seurat, paste0(path, '/seurat_ATAC.rds'))

cat("Unified peak set processing complete!\n")
cat("macs2_combined assay created with", nrow(seurat[['macs2_combined']]), "unified peaks\n")
cat("Next step: Run integration scripts using macs2_combined assay.\n")

