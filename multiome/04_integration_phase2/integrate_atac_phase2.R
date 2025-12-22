#!/usr/bin/env Rscript
# ==============================================================================
# Phase 2: ATAC Re-integration After Filtering
#
# After filtering low-quality cells based on Phase 1 WNN results, this script
# re-runs ATAC integration using Harmony.
#
# IMPORTANT FOR REPRODUCIBILITY:
# - UMAP uses seed.use = 1234 to ensure identical layouts
# - Clustering resolution = 2 matches the original analysis
# - LSI dimensions 2:30 are used (dim 1 typically captures sequencing depth)
#
# Usage: Rscript integrate_atac_phase2.R <input_path> <filtered_object_path> [n_workers]
# ==============================================================================

library(Seurat)
library(Signac)
library(stringr)
library(harmony)
library(tidyverse)
library(future)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript integrate_atac_phase2.R <input_path> <filtered_object_path> [n_workers]")
}

inPath <- args[1]
filteredObjectPath <- args[2]
n_workers <- ifelse(length(args) >= 3, as.numeric(args[3]), 4)

plan("multisession", workers = n_workers)
options(future.globals.maxSize = 800000 * 1024^2)  # 800GB

# ==============================================================================
# Step 1: Load filtered object
# ==============================================================================
cat("Loading filtered object...\n")
seurat <- readRDS(filteredObjectPath)

# ==============================================================================
# Step 2: Prepare ATAC assay for re-integration
# ==============================================================================
cat("Preparing ATAC assay for re-integration...\n")
DefaultAssay(seurat) <- "macs2_combined"

seurat <- FindTopFeatures(seurat, min.cutoff = 5)
seurat <- RunTFIDF(seurat)
seurat <- RunSVD(seurat)

# Split by donor for integration
seurat[["macs2_combined"]] <- split(seurat[["macs2_combined"]], f = seurat$donor)
VariableFeatures(object = seurat, assay = 'macs2_combined') <- rownames(seurat)[1:2000]

# ==============================================================================
# Step 3: Harmony re-integration
# ==============================================================================
cat("Running Harmony re-integration...\n")
seurat <- IntegrateLayers(
  object = seurat,
  method = HarmonyIntegration,
  assay = 'macs2_combined',
  orig.reduction = "lsi",
  new.reduction = "macs2Combined.harmony2",
  verbose = FALSE
)

# ==============================================================================
# Step 4: Clustering and visualization
# ==============================================================================
cat("Clustering and creating visualizations...\n")
seurat <- FindNeighbors(seurat, reduction = "macs2Combined.harmony2", dims = 2:30)
seurat <- FindClusters(
  seurat,
  resolution = 2,
  cluster.name = "macs2Combined_wnn5Filtered_harmony_clusters"
)

# CRITICAL: Use seed.use = 1234 for reproducibility
seurat <- RunUMAP(
  seurat,
  reduction = "macs2Combined.harmony2",
  dims = 2:30,
  reduction.name = "macs2Combined.harmony2.umap",
  reduction.key = 'HARMONY2UMAP_',
  seed.use = 1234  # Fixed seed for reproducibility
)

# Visualize
pdf(paste0(inPath, '/ATAC_harmony2_umap.pdf'), width = 12, height = 8)
p1 <- DimPlot(seurat, group.by = "donor", reduction = "macs2Combined.harmony2.umap", label = FALSE)
p2 <- DimPlot(seurat, group.by = "macs2Combined_wnn5Filtered_harmony_clusters", 
              reduction = "macs2Combined.harmony2.umap", label = TRUE) + NoLegend()
gridExtra::grid.arrange(grobs = list(p1, p2), ncol = 2)
dev.off()

# ==============================================================================
# Step 5: Save metadata and reductions
# ==============================================================================
cat("Saving metadata and reductions...\n")
saveRDS(seurat@reductions, paste0(inPath, '/reduction_macs2Combined_wnn5Filtered_harmonyIntegrated_run2.rds'))
saveRDS(seurat@meta.data, paste0(inPath, '/metadata_macs2Combined_wnn5Filtered_harmonyIntegrated_run2.rds'))

# Save updated object
saveRDS(seurat, filteredObjectPath)

cat("Phase 2 ATAC re-integration complete!\n")

