#!/usr/bin/env Rscript
# ==============================================================================
# Phase 2: Subcluster and Create Final Annotation (celltype2)
#
# After Phase 2 WNN integration, this script performs subclustering on specific
# clusters to identify rare cell types (e.g., epsilon cells) and creates the
# final cell type annotation (celltype2) used in downstream analysis.
#
# WORKFLOW:
# 1. Subcluster specific clusters from Phase 2 WNN results (clusters 11, 12, 13, 14, 15)
# 2. Combine subcluster annotations into a unified annotation
# 3. Create celltype2 annotation with refined cell type assignments (including epsilon)
# 4. Filter out "mixed" clusters to create final cleaned object
# 5. Re-run final WNN integration on cleaned object
#
# IMPORTANT FOR REPRODUCIBILITY:
# - All subclustering uses algorithm = 3 (Louvain) for consistency
# - Final WNN uses seed.use = 1234 for identical UMAP layouts
# - Cluster IDs are sorted numerically before annotation mapping
#
# Usage: Rscript subcluster_and_annotate_final.R <input_path> <filtered_object_path> [n_workers]
# ==============================================================================

library(Seurat)
library(Signac)
library(tidyverse)
library(future)
library(gtools)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript subcluster_and_annotate_final.R <input_path> <filtered_object_path> [n_workers]")
}

inPath <- args[1]
filteredObjectPath <- args[2]
n_workers <- ifelse(length(args) >= 3, as.numeric(args[3]), 4)

plan("multisession", workers = n_workers)
options(future.globals.maxSize = 200000 * 1024^2)  # 200GB

# ==============================================================================
# Step 1: Load Phase 2 WNN integrated object
# ==============================================================================
cat("Loading Phase 2 WNN integrated object...\n")
cat("Expected file:", filteredObjectPath, "\n")
seurat <- readRDS(filteredObjectPath)

# Check that we have the Phase 2 WNN clusters
expected_cluster_col <- "clusters_wnn5filteredRun2"
if (!expected_cluster_col %in% colnames(seurat@meta.data)) {
  stop("Expected cluster column '", expected_cluster_col, "' not found in metadata.\n",
       "Please ensure you've run integrate_wnn_phase2.R first.\n")
}

cat("Object loaded:", ncol(seurat), "cells\n")
# Use $ notation for consistency with original code
cat("Phase 2 WNN clusters:", length(unique(seurat[[expected_cluster_col]])), "clusters\n")

# ==============================================================================
# Step 2: Subcluster Specific Clusters
# ==============================================================================
# Subclustering is performed on specific clusters to identify rare cell types
# (e.g., epsilon cells from cluster 15) and refine annotations
#
# Clusters to subcluster: 11, 12, 13, 14, 15
# These clusters were identified as potentially containing rare or mixed populations

cat("\nSubclustering specific clusters to refine cell type annotations...\n")
graph_name <- "wsnn5filteredRun2"

# Subcluster cluster 13
cat("Subclustering cluster 13...\n")
key_sub13 <- "5filteredRun2_wSubC13"
cluster_col_sub13 <- paste0("clusters_wnn", key_sub13)
seurat <- FindSubCluster(
  seurat,
  cluster = 13,
  graph.name = graph_name,
  algorithm = 3,
  subcluster.name = cluster_col_sub13
)
seurat@meta.data[, cluster_col_sub13] <- factor(
  seurat@meta.data[, cluster_col_sub13],
  levels = gtools::mixedsort(unique(seurat@meta.data[, cluster_col_sub13]))
)

# Subcluster cluster 11
cat("Subclustering cluster 11...\n")
key_sub11 <- "5filteredRun2_wSubC11"
cluster_col_sub11 <- paste0("clusters_wnn", key_sub11)
seurat <- FindSubCluster(
  seurat,
  cluster = 11,
  graph.name = graph_name,
  algorithm = 3,
  subcluster.name = cluster_col_sub11
)
seurat@meta.data[, cluster_col_sub11] <- factor(
  seurat@meta.data[, cluster_col_sub11],
  levels = gtools::mixedsort(unique(seurat@meta.data[, cluster_col_sub11]))
)

# Subcluster cluster 12
cat("Subclustering cluster 12...\n")
key_sub12 <- "5filteredRun2_wSubC12"
cluster_col_sub12 <- paste0("clusters_wnn", key_sub12)
seurat <- FindSubCluster(
  seurat,
  cluster = 12,
  graph.name = graph_name,
  algorithm = 3,
  subcluster.name = cluster_col_sub12
)
seurat@meta.data[, cluster_col_sub12] <- factor(
  seurat@meta.data[, cluster_col_sub12],
  levels = gtools::mixedsort(unique(seurat@meta.data[, cluster_col_sub12]))
)

# Subcluster cluster 14
cat("Subclustering cluster 14...\n")
key_sub14 <- "5filteredRun2_wSubC14"
cluster_col_sub14 <- paste0("clusters_wnn", key_sub14)
seurat <- FindSubCluster(
  seurat,
  cluster = 14,
  graph.name = graph_name,
  algorithm = 3,
  subcluster.name = cluster_col_sub14
)
seurat@meta.data[, cluster_col_sub14] <- factor(
  seurat@meta.data[, cluster_col_sub14],
  levels = gtools::mixedsort(unique(seurat@meta.data[, cluster_col_sub14]))
)

# Subcluster cluster 15 (epsilon cells identified here)
cat("Subclustering cluster 15 (epsilon cell identification)...\n")
key_sub15 <- "5filteredRun2_wSubC15"
cluster_col_sub15 <- paste0("clusters_wnn", key_sub15)
seurat <- FindSubCluster(
  seurat,
  cluster = 15,
  graph.name = graph_name,
  algorithm = 3,
  subcluster.name = cluster_col_sub15
)
seurat@meta.data[, cluster_col_sub15] <- factor(
  seurat@meta.data[, cluster_col_sub15],
  levels = gtools::mixedsort(unique(seurat@meta.data[, cluster_col_sub15]))
)

cat("Subclustering complete!\n")

# ==============================================================================
# Step 3: Combine Subcluster Annotations
# ==============================================================================
cat("\nCombining subcluster annotations...\n")
key_combined <- "5filteredRun2_wSubCombined"
cluster_col_combined <- paste0("clusters_wnn", key_combined)

# Start with original Phase 2 WNN clusters
seurat@meta.data[[cluster_col_combined]] <- as.character(seurat[[expected_cluster_col]])
cell_barcodes <- rownames(seurat@meta.data)

# Replace subclustered clusters with their subcluster assignments
# Use $ notation for consistency with original code
seurat@meta.data[
  cell_barcodes[which(seurat[[expected_cluster_col]] == "13")],
  cluster_col_combined
] <- as.character(seurat@meta.data[
  cell_barcodes[which(seurat[[expected_cluster_col]] == "13")],
  cluster_col_sub13
])

seurat@meta.data[
  cell_barcodes[which(seurat[[expected_cluster_col]] == "11")],
  cluster_col_combined
] <- as.character(seurat@meta.data[
  cell_barcodes[which(seurat[[expected_cluster_col]] == "11")],
  cluster_col_sub11
])

seurat@meta.data[
  cell_barcodes[which(seurat[[expected_cluster_col]] == "14")],
  cluster_col_combined
] <- as.character(seurat@meta.data[
  cell_barcodes[which(seurat[[expected_cluster_col]] == "14")],
  cluster_col_sub14
])

seurat@meta.data[
  cell_barcodes[which(seurat[[expected_cluster_col]] == "15")],
  cluster_col_combined
] <- as.character(seurat@meta.data[
  cell_barcodes[which(seurat[[expected_cluster_col]] == "15")],
  cluster_col_sub15
])

# Convert to factor with sorted levels
seurat@meta.data[[cluster_col_combined]] <- factor(
  seurat@meta.data[[cluster_col_combined]],
  levels = gtools::mixedsort(unique(seurat@meta.data[[cluster_col_combined]]))
)

cat("Combined annotation created:", length(levels(seurat@meta.data[[cluster_col_combined]])), "clusters\n")

# ==============================================================================
# Step 4: Create Final Annotation (celltype2)
# ==============================================================================
cat("\nCreating final celltype2 annotation...\n")

# First create detailed annotation (wnn5SubCombinedAnnot)
# This annotation includes all subcluster details
wnn5SubCombinedAnnot <- c(
  "alpha", "acinar", "alpha", "ductal", "beta", "beta", "ductal",
  "activated.stellate", "acinar", "delta",
  "endothelial", rep("beta", 3), "alpha", rep("beta", 2), "ductal",
  rep("pp.gamma", 3), "epsilon", rep("pp.gamma", 2), rep("c11", 3), "c11_exocrine",
  rep("c13", 3), "c13_beta", "c13", "c13_acinar", "c13",
  "quiescent.stellate", "alpha", "acinar", "immune", "ductal",
  "alpha_c21", "immune_endocrine", "beta_c23", "alpha_beta_c24"
)

seurat@meta.data$wnn5SubCombinedAnnot <- seurat@meta.data[[cluster_col_combined]]
levels(seurat@meta.data$wnn5SubCombinedAnnot) <- wnn5SubCombinedAnnot

# Create celltype2 annotation (simplified, with "mixed" for ambiguous clusters)
seurat@meta.data$celltype2 <- seurat@meta.data[[cluster_col_combined]]

celltype2Annt <- c(
  "alpha", "acinar", "alpha", "ductal", "beta", "beta", "ductal",
  "activated.stellate", "acinar", "delta",
  "endothelial", rep("beta", 3), "alpha", rep("beta", 2), "ductal",
  rep("pp.gamma", 3), "epsilon", rep("pp.gamma", 2), rep("c11", 3), "mixed",
  rep("c13", 3), "mixed", "c13", "mixed", "c13",
  "quiescent.stellate", "alpha", "acinar", "immune", "ductal",
  "mixed", "mixed", "mixed", "mixed"
)

levels(seurat@meta.data$celltype2) <- celltype2Annt

cat("Final annotation created:\n")
cat("  Total clusters:", length(unique(seurat$celltype2)), "\n")
cat("  Cell type distribution:\n")
print(table(seurat$celltype2))

# ==============================================================================
# Step 5: Filter Mixed Clusters and Create Final Object
# ==============================================================================
cat("\nFiltering out 'mixed' clusters...\n")
seurat_final <- subset(
  seurat,
  subset = celltype2 %in% setdiff(unique(seurat$celltype2), "mixed")
)

cat("Cells before filtering:", ncol(seurat), "\n")
cat("Cells after filtering:", ncol(seurat_final), "\n")
cat("Cells filtered:", ncol(seurat) - ncol(seurat_final), "\n")
cat("Percentage retained:", round(100 * ncol(seurat_final) / ncol(seurat), 2), "%\n")

# ==============================================================================
# Step 6: Final WNN Integration on Cleaned Object
# ==============================================================================
cat("\nPerforming final WNN integration on cleaned object...\n")
rd_list <- list("combinedSctNormGEX.rpcaRefIntegrated.run2", "macs2Combined.harmony2")
key_final <- "5filteredRun2"

cat("Using reductions:", paste(rd_list, collapse = " + "), "\n")

# Find multimodal neighbors
seurat_final <- FindMultiModalNeighbors(
  seurat_final,
  reduction.list = rd_list,
  dims.list = list(1:30, 2:30),
  knn.graph.name = paste0("wknn", key_final),
  snn.graph.name = paste0("wsnn", key_final),
  weighted.nn.name = paste0("weighted.nn", key_final)
)

# Run UMAP (with model for potential projection)
# CRITICAL: Use seed.use = 1234 for reproducibility
seurat_final <- RunUMAP(
  seurat_final,
  nn.name = paste0("weighted.nn", key_final),
  reduction.name = "wnn5filteredRun2.umap",
  reduction.key = "wnn5filteredRun2Run2UMAP_",
  return.model = TRUE,  # Allow projection of query data
  seed.use = 1234  # Fixed seed for reproducibility
)

# ==============================================================================
# Step 7: Save Final Object
# ==============================================================================
cat("\nSaving final object...\n")
final_output_path <- paste0(inPath, '/multiome40DonorsRef-v2.rds')
saveRDS(seurat_final, final_output_path)

# Also save the object with subclusters (before final filtering) for reference
saveRDS(seurat, paste0(inPath, '/seurat_wnnCombined_run5FilteredClean2Run2.rds'))

cat("\nFinal annotation and filtering complete!\n")
cat("Final object saved to:", final_output_path, "\n")
cat("Final annotation column: celltype2\n")
cat("Final UMAP reduction: wnn5filteredRun2.umap\n")
cat("Epsilon cells identified through subclustering of cluster 15\n")
cat("\nThe final object (multiome40DonorsRef-v2.rds) is ready for downstream analysis!\n")

