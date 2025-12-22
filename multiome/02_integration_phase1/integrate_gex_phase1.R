#!/usr/bin/env Rscript
# ==============================================================================
# Phase 1: GEX (Gene Expression) Integration
#
# This script integrates GEX data from all multiome samples using multiple
# integration methods:
# 1. Creates combined Seurat object with layers (Seurat v5)
# 2. Runs SCTransform normalization
# 3. Performs Harmony integration
# 4. Performs RPCA integration (all samples and reference-based)
#
# The output includes multiple integration approaches for comparison.
#
# Usage: Rscript integrate_gex_phase1.R <input_path> <metadata_path> [npc] [n_workers]
# ==============================================================================

library(Seurat)
library(stringr)
library(sctransform)
library(harmony)
library(cowplot)
library(tidyverse)
library(data.table)
library(future)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript integrate_gex_phase1.R <input_path> <metadata_path> [npc] [n_workers]")
}

inPath <- args[1]
mdPath <- args[2]
npc <- ifelse(length(args) >= 3, as.numeric(args[3]), 30)
n_workers <- ifelse(length(args) >= 4, as.numeric(args[4]), 4)

# Set up parallel processing
plan("multisession", workers = n_workers)
options(future.globals.maxSize = 200000 * 1024^2)  # 200GB


# ==============================================================================
# Step 1: Load individual sample objects and combine
# ==============================================================================
cat("Loading individual sample objects...\n")
multiome_md_combined <- fread(mdPath) %>% 
  as.data.frame() %>% 
  column_to_rownames("V1")

paths <- list.files(inPath, "seurat_sctNormGEX.rds", recursive = TRUE)
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
  obj_i <- subset(obj_i, cells = bc_sel)
  obj_list[[i]] <- obj_i
}

# ==============================================================================
# Step 2: Extract counts and metadata, create combined object
# ==============================================================================
cat("Creating combined Seurat object...\n")
counts_list <- list()
md_list <- list()

for (i in seq_along(obj_list)) {
  obj_i <- obj_list[[i]]
  counts_list[[i]] <- obj_i[["RNA"]]$counts
  md_list[[i]] <- obj_i@meta.data
}

names(counts_list) <- sample_names

common_columns <- Reduce(intersect, lapply(md_list, colnames))
md_combined <- do.call(
  rbind,
  lapply(md_list, function(md) { md[, common_columns] })
)

# Create Seurat object with multiple layers (Seurat v5)
seurat_combined_v5 <- CreateSeuratObject(
  counts = counts_list,
  meta.data = md_combined,
  assay = 'RNA'
)

# ==============================================================================
# Step 3: SCTransform normalization
# ==============================================================================
cat("Running SCTransform normalization...\n")
seurat_combined_v5 <- SCTransform(
  seurat_combined_v5, 
  vars.to.regress = "percent.mt"
)

seurat_combined_v5 <- RunPCA(
  seurat_combined_v5, 
  assay = 'SCT', 
  reduction.name = 'combinedSctNormGEX.pca', 
  verbose = FALSE, 
  reduction.key = 'COMBINEDSCTNORMGEXPC_'
)

# Unintegrated UMAP and clustering
seurat_combined_v5 <- RunUMAP(
  seurat_combined_v5, 
  reduction = "combinedSctNormGEX.pca", 
  reduction.key = 'COMBINEDSCTNORMGEXUMAP_',
  dims = 1:npc, 
  verbose = FALSE, 
  reduction.name = 'combinedSctNormGEX.umap'
)

seurat_combined_v5 <- FindNeighbors(
  seurat_combined_v5, 
  reduction = "combinedSctNormGEX.pca", 
  dims = 1:npc, 
  verbose = FALSE, 
  graph.name = c('combinedsct.nn', 'combinedsct.snn')
)

seurat_combined_v5 <- FindClusters(
  seurat_combined_v5, 
  resolution = 0.5, 
  method = "igraph", 
  verbose = FALSE, 
  graph.name = 'combinedsct.snn',
  cluster.name = "combined_sct_unintegrated_clusters"
)

# Save unintegrated object
saveRDS(seurat_combined_v5, paste0(inPath, '/seurat_sctNormGEX_unintegrated.rds'))
cat("Unintegrated object saved.\n")

# ==============================================================================
# Step 4: Harmony integration
# ==============================================================================
cat("Running Harmony integration...\n")
set.seed(1234)
seurat_combined_v5 <- IntegrateLayers(
  object = seurat_combined_v5,
  method = HarmonyIntegration,
  orig.reduction = "combinedSctNormGEX.pca",
  new.reduction = "combinedSctNormGEX.harmony",
  verbose = FALSE
)

seurat_combined_v5 <- FindNeighbors(
  seurat_combined_v5, 
  reduction = "combinedSctNormGEX.harmony", 
  dims = 1:npc, 
  verbose = FALSE, 
  graph.name = c('harmony.nn', 'harmony.snn')
)

seurat_combined_v5 <- FindClusters(
  seurat_combined_v5, 
  resolution = 0.5, 
  method = "igraph", 
  verbose = FALSE, 
  graph.name = 'harmony.snn',
  cluster.name = "combined_sct_harmony_clusters"
)

seurat_combined_v5 <- RunUMAP(
  seurat_combined_v5, 
  reduction = "combinedSctNormGEX.harmony", 
  reduction.key = 'HARMONYSCTNORMGEXUMAP_',
  dims = 1:npc, 
  verbose = FALSE, 
  reduction.name = 'harmonySctNormGEX.umap'
)

# Visualize Harmony integration
pdf(paste0(inPath, '/GEX_sctNormGEX_harmonyIntegrated_umap.pdf'), width = 11, height = 11)
p1 <- DimPlot(seurat_combined_v5, group.by = "donor", reduction = "harmonySctNormGEX.umap", label = FALSE)
p2 <- DimPlot(seurat_combined_v5, group.by = "combined_sct_harmony_clusters", 
              reduction = "harmonySctNormGEX.umap", label = TRUE) + NoLegend()
gridExtra::grid.arrange(grobs = list(p1, p2), nrow = 1)
dev.off()

# Save metadata and reductions for later use
saveRDS(seurat_combined_v5@meta.data, paste0(inPath, '/metadata_sctNormGEX_harmonyIntegrated.rds'))
saveRDS(seurat_combined_v5@reductions, paste0(inPath, '/reduction_sctNormGEX_harmonyIntegrated.rds'))

cat("Harmony integration complete.\n")
cat("Next: Run integrate_gex_rpca_all.R and integrate_gex_rpca_ref.R for RPCA integration.\n")

