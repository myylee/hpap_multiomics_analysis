#!/usr/bin/env Rscript
# ==============================================================================
# Phase 1: ATAC Harmony Integration
#
# This script integrates ATAC-seq data from all multiome samples using Harmony:
# 1. Loads individual ATAC sample objects
# 2. Combines into multi-layer Seurat object (Seurat v5)
# 3. Runs LSI dimensionality reduction
# 4. Performs Harmony integration
# 5. Clusters and visualizes results
#
# Usage: Rscript integrate_atac_harmony.R <input_path> <metadata_path> [npc] [n_workers]
# ==============================================================================

library(Signac)
library(Seurat)
library(tidyverse)
library(fs)
library(stringr)
library(data.table)
library(future)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript integrate_atac_harmony.R <input_path> <metadata_path> [npc] [n_workers]")
}

inPath <- args[1]
mdPath <- args[2]
npc <- ifelse(length(args) >= 3, as.numeric(args[3]), 30)
n_workers <- ifelse(length(args) >= 4, as.numeric(args[4]), 4)

# Set up parallel processing
plan("multisession", workers = n_workers)
options(future.globals.maxSize = 200000 * 1024^2)  # 200GB

# ==============================================================================
# Step 1: Load individual sample objects
# ==============================================================================
cat("Loading individual ATAC sample objects...\n")
multiome_md_combined <- fread(mdPath) %>% 
  as.data.frame() %>% 
  column_to_rownames("V1")

paths <- list.files(inPath, "seurat_ATAC.rds", recursive = TRUE)
obj_list <- list()
sample_names <- gsub("/islets.*", "", paths)

for (i in seq_along(paths)) {
  path_i <- paths[i]
  cat("Loading:", path_i, "\n")
  donor_i <- sample_names[i]
  obj_i <- readRDS(file.path(inPath, path_i))
  
  # Subset to cells in metadata
  bc_sel <- multiome_md_combined %>% 
    filter(donor == donor_i) %>% 
    pull(V1)
  
  if (length(bc_sel) == 0) {
    next
  }
  
  obj_i <- subset(obj_i, cells = bc_sel)
  obj_list[[i]] <- obj_i
}

names(obj_list) <- sample_names
# Filter objects with no cells
obj_list <- obj_list[which(!unlist(lapply(obj_list, is.null)))]

# ==============================================================================
# Step 2: Extract counts and metadata, create combined object
# ==============================================================================
cat("Creating combined Seurat object...\n")
counts_list <- list()
md_list <- list()

for (i in seq_along(obj_list)) {
  obj_i <- obj_list[[i]]
  counts_list[[i]] <- obj_i@assays$macs2_combined@counts
  md_list[[i]] <- obj_i@meta.data
}

names(counts_list) <- names(obj_list)

common_columns <- Reduce(intersect, lapply(md_list, colnames))
md_combined <- do.call(
  rbind,
  lapply(md_list, function(md) { md[, common_columns] })
)

# Create Seurat object with multiple layers (Seurat v5)
seurat_combined <- CreateSeuratObject(
  counts = counts_list,
  meta.data = md_combined,
  assay = 'macs2_combined'
)

# ==============================================================================
# Step 3: LSI dimensionality reduction (unintegrated)
# ==============================================================================
cat("Running LSI dimensionality reduction...\n")
merged <- JoinLayers(seurat_combined)
merged <- FindTopFeatures(merged, min.cutoff = 5)
merged <- RunTFIDF(merged)
merged <- RunSVD(merged)
merged <- RunUMAP(
  merged, 
  reduction = "lsi", 
  dims = 2:30,
  reduction.name = "macs2CombinedUnintegrated.umap",
  reduction.key = 'MACS2COMBINEDUNINTEGRATEDUMAP_'
)

# ==============================================================================
# Step 4: Harmony integration
# ==============================================================================
cat("Running Harmony integration...\n")
merged[["macs2_combined"]] <- split(merged[["macs2_combined"]], f = merged$donor)
VariableFeatures(object = merged, assay = 'macs2_combined') <- rownames(merged)[1:2000]

merged <- IntegrateLayers(
  object = merged,
  method = HarmonyIntegration,
  assay = 'macs2_combined',
  orig.reduction = "lsi",
  new.reduction = "macs2Combined.harmony",
  verbose = FALSE
)

# ==============================================================================
# Step 5: Clustering and visualization
# ==============================================================================
cat("Clustering and creating visualizations...\n")
merged <- FindNeighbors(merged, reduction = "macs2Combined.harmony", dims = 2:30)
merged <- FindClusters(
  merged, 
  resolution = 2, 
  cluster.name = "macs2Combined_harmony_clusters"
)

merged <- RunUMAP(
  merged, 
  reduction = "macs2Combined.harmony", 
  dims = 2:30, 
  reduction.name = "macs2CombinedHarmony.umap",
  reduction.key = 'MACS2COMBINEDHARMONYUMAP_'
)

# Visualize
p1 <- DimPlot(
  merged,
  reduction = "macs2CombinedUnintegrated.umap",
  group.by = c("macs2Combined_harmony_clusters", "donor")
)

p2 <- DimPlot(
  merged,
  reduction = "macs2CombinedHarmony.umap",
  group.by = c("macs2Combined_harmony_clusters", "donor")
)

pdf(paste0(inPath, '/ATAC_macs2Combined_wHarmony_umap.pdf'), width = 20, height = 20)
print(p1 / p2)
dev.off()

# Save object
saveRDS(merged, paste0(inPath, '/seurat_macs2Combined_wHarmony.rds'))

cat("ATAC Harmony integration complete!\n")
cat("Saved object:", paste0(inPath, '/seurat_macs2Combined_wHarmony.rds'), "\n")

