#!/usr/bin/env Rscript
# ==============================================================================
# Phase 1: ATAC RLSI Integration (Reference-Based)
#
# This script integrates ATAC-seq data using Reciprocal LSI (RLSI) with
# reference-based integration:
# 1. Loads individual ATAC sample objects
# 2. Prepares objects for integration (LSI reduction)
# 3. Merges objects
# 4. Performs RLSI integration using reference samples
# 5. Clusters and visualizes results
#
# Usage: Rscript integrate_atac_rlsi.R <input_path> <metadata_path> [npc] [n_workers] [reference_indices]
#   reference_indices: Comma-separated list of sample indices to use as reference
#                      Default: 5,6,9,16,21,27,29,31,32 (high-quality control samples)
# ==============================================================================

library(Seurat)
library(Signac)
library(stringr)
library(tidyverse)
library(data.table)
library(future)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript integrate_atac_rlsi.R <input_path> <metadata_path> [npc] [n_workers] [reference_indices]")
}

inPath <- args[1]
mdPath <- args[2]
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
options(future.globals.maxSize = 800000 * 1024^2)  # 800GB
options(Seurat.object.assay.version = "v3")

set.seed(1234)

# ==============================================================================
# Step 1: Load and prepare individual sample objects
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
  
  # Prepare for integration: filter peaks, run LSI
  peak_keep <- which(rowSums(obj_i[['macs2_combined']]@counts) > 0)
  obj_i <- CreateSeuratObject(
    obj_i[['macs2_combined']]@counts[peak_keep, ],
    meta.data = obj_i@meta.data,
    assay = 'macs2_combined'
  ) %>%
    FindTopFeatures(min.cutoff = 5) %>%
    RunTFIDF() %>%
    RunSVD()
  
  obj_list[[i]] <- obj_i
}

names(obj_list) <- sample_names
# Filter objects with no cells
obj_list <- obj_list[which(!unlist(lapply(obj_list, is.null)))]

cat("Reference samples (indices):", reference, "\n")
cat("Reference sample names:", names(obj_list)[reference], "\n")

# ==============================================================================
# Step 2: Merge objects and compute LSI
# ==============================================================================
cat("Merging objects and computing LSI...\n")
atac.combined <- obj_list[[1]]
for (i in 2:length(obj_list)) {
  atac.combined <- merge(atac.combined, obj_list[[i]])
}

atac.combined <- FindTopFeatures(atac.combined, min.cutoff = 10)
atac.combined <- RunTFIDF(atac.combined)
atac.combined <- RunSVD(atac.combined)

features <- atac.combined@assays$macs2_combined@var.features

# ==============================================================================
# Step 3: Find integration anchors
# ==============================================================================
cat("Finding integration anchors...\n")
integration.anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  reduction = "rlsi",
  dims = 2:30,
  reference = reference,
  anchor.features = features
)

saveRDS(integration.anchors, file = paste0(inPath, '/macs2Combined_anchors_v3.tmp.rds'))

# ==============================================================================
# Step 4: Integrate embeddings
# ==============================================================================
cat("Integrating LSI embeddings...\n")
atac.integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = atac.combined[["lsi"]],
  new.reduction.name = "macs2Combined.rlsiRefIntegrated",
  dims.to.integrate = 1:30
)

# ==============================================================================
# Step 5: UMAP, clustering, and visualization
# ==============================================================================
cat("Running UMAP and clustering...\n")
atac.integrated <- RunUMAP(
  atac.integrated,
  reduction = 'macs2Combined.rlsiRefIntegrated',
  dims = 2:30,
  reduction.name = "rlsiRefIntegrated.umap",
  reduction.key = "RLSIREFINTEGRATEDUMAP",
  seed.use = 1234
)

atac.integrated <- FindNeighbors(
  atac.integrated, 
  reduction = "macs2Combined.rlsiRefIntegrated", 
  dims = 1:npc, 
  verbose = FALSE, 
  graph.name = c('rlsirefintegrated.nn', 'rlsirefintegrated.snn')
)

atac.integrated <- FindClusters(
  atac.integrated, 
  resolution = 0.5, 
  method = "igraph", 
  verbose = FALSE, 
  graph.name = 'rlsirefintegrated.snn',
  cluster.name = "macs2Combined_rlsiRefIntegrated_clusters"
)

# Visualize
pdf(paste0(inPath, '/ATAC_macs2Combined_rlsiRefIntegrated_umap.pdf'), width = 12, height = 8)
p1 <- DimPlot(
  atac.integrated, 
  group.by = "donor", 
  reduction = "rlsiRefIntegrated.umap", 
  label = FALSE
)
p2 <- DimPlot(
  atac.integrated, 
  group.by = "macs2Combined_rlsiRefIntegrated_clusters", 
  reduction = "rlsiRefIntegrated.umap", 
  label = TRUE
) + NoLegend()
gridExtra::grid.arrange(grobs = list(p1, p2), ncol = 2)
dev.off()

# Save object
saveRDS(atac.integrated, file = paste0(inPath, '/ATAC_macs2Combined_rlsiRefIntegrated.rds'))

cat("ATAC RLSI integration complete!\n")
cat("Saved object:", paste0(inPath, '/ATAC_macs2Combined_rlsiRefIntegrated.rds'), "\n")

