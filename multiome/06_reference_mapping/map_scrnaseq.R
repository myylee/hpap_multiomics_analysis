#!/usr/bin/env Rscript
# ==============================================================================
# Map scRNA-seq Query Data to Multiome Reference
#
# This script maps scRNA-seq query data to the multiome reference using
# Seurat's reference mapping approach. The script assumes the query data
# contains cells from multiple donors merged together. It will:
# 1. Create a clean query object and split by donor
# 2. Process each donor sample separately
# 3. Map each sample to the reference
# 4. Combine results
#
# NOTE: If analyzing a single sample (one donor), you can skip the splitting
# and batch processing steps (Step 2) and directly process the query object.
# Simply comment out the splitting code and use the query object directly in
# the mapping loop (or replace the loop with a single mapping call).
#
# Usage: Rscript map_scrnaseq.R <reference_path> <query_path> <output_dir> [reference_umap]
#   reference_path: Path to multiome reference object
#   query_path: Path to query scRNA-seq Seurat object (may contain multiple donors)
#   output_dir: Directory to save mapping results
#   reference_umap: Name of UMAP reduction in reference (default: wnn5filteredRun2.umap)
# ==============================================================================

library(Seurat)
library(stringr)
library(sctransform)
library(harmony)
library(cowplot)
library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript map_scrnaseq.R <reference_path> <query_path> <output_dir> [reference_umap]")
}

referencePath <- args[1]
queryPath <- args[2]
out_dir <- args[3]
reference_umap <- ifelse(length(args) >= 4, args[4], "wnn5filteredRun2.umap")

dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# Step 1: Prepare reference
# ==============================================================================
cat("Loading and preparing reference...\n")
reference <- readRDS(referencePath)

DefaultAssay(reference) <- 'RNA'
reference <- NormalizeData(reference) %>%
  FindVariableFeatures() %>%
  ScaleData()

# Use supervised PCA based on WNN graph
reference <- RunSPCA(reference, graph = 'wsnn5filteredRun2')

# Set cell type annotation
if ("celltype2" %in% colnames(reference@meta.data)) {
  reference@meta.data$celltype <- reference$celltype2
} else if ("celltype" %in% colnames(reference@meta.data)) {
  # Use existing celltype if celltype2 not available
} else {
  stop("No cell type annotation found in reference object")
}

Idents(reference) <- 'celltype'

# Build reference neighbor graph for fast mapping
reference <- FindNeighbors(
  object = reference,
  reduction = "spca",
  dims = 1:30,
  graph.name = "spca.annoy.neighbors",
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

# ==============================================================================
# Step 2: Prepare query data
# ==============================================================================
cat("Loading and preparing query data...\n")
query <- readRDS(file = queryPath)

# Keep only RNA assay to keep dataset small
query_clean <- CreateSeuratObject(
  counts = query@assays$RNA@counts,
  meta.data = query@meta.data
)
rm(query)
gc()

query <- query_clean
rm(query_clean)

# Split query by donor for processing multiple samples
# NOTE: If you have a single sample (one donor), you can skip this splitting
# step and directly normalize and process the query object. See comments in
# Step 3 for single-sample processing.
cat("Splitting query data by donor...\n")
query@meta.data$donor <- gsub("-.*", "", query$orig.ident)
query.batches <- SplitObject(query, split.by = "donor")
query.batches <- lapply(X = query.batches, FUN = NormalizeData, verbose = FALSE)

cat("Found", length(query.batches), "donor(s) in query data.\n")

# ==============================================================================
# Step 3: Map each query batch to reference
# ==============================================================================
# NOTE: If you skipped the splitting step (single sample analysis), you can
# replace this loop with a single mapping call. First normalize the query:
#   query <- NormalizeData(query, verbose = FALSE)
# Then use 'query' directly instead of 'query.batches[[i]]' in the mapping code below.

cat("Mapping query data to reference...\n")
format_mapping_res <- function(query_seurat) {
  require(dplyr)
  df <- cbind(
    query_seurat@meta.data %>%
      select(predicted.celltype, predicted.celltype.score),
    query_seurat@reductions$ref.umap@cell.embeddings
  )
  return(df)
}

mapped_res <- c()
for (i in seq_along(query.batches)) {
  cat("Mapping batch", i, "of", length(query.batches), "\n")
  query <- query.batches[[i]]
  
  # Find transfer anchors
  anchors_i <- FindTransferAnchors(
    reference = reference,
    query = query,
    reference.reduction = "spca",
    reference.neighbors = "spca.annoy.neighbors",
    dims = 1:30,
    k.filter = NA
  )
  
  # Map query to reference
  query <- MapQuery(
    anchorset = anchors_i,
    query = query,
    reference = reference,
    refdata = list(celltype = "celltype"),
    reference.reduction = "spca",
    reduction.model = reference_umap,
    transferdata.args = list(k.weight = 20)
  )
  
  Idents(query) <- 'predicted.celltype'
  
  # Visualize mapping
  pdf(file.path(out_dir, paste0('mapping_', unique(query.batches[[i]]$donor), '.pdf')), 
      width = 16, height = 12)
  p1 <- DimPlot(query, reduction = "ref.umap", label = TRUE, label.size = 4, 
                repel = TRUE, shuffle = TRUE) + 
    NoLegend() + 
    ggtitle(paste("Mapped labels:", unique(query.batches[[i]]$donor)))
  p2 <- FeaturePlot(query, features = "predicted.celltype.score", combine = TRUE, 
                    pt.size = 0.1, reduction = "ref.umap") &
    scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
  p3 <- DimPlot(reference, reduction = reference_umap, label = TRUE, label.size = 4,
                group.by = 'celltype', repel = TRUE, shuffle = TRUE) + 
    NoLegend() + 
    ggtitle("Reference annotations")
  gridExtra::grid.arrange(grobs = list(p1, p2, p3), nrow = 2)
  dev.off()
  
  # Store results
  res_i <- format_mapping_res(query)
  mapped_res <- rbind(mapped_res, res_i)
  query.batches[[i]] <- query
}

# ==============================================================================
# Step 4: Save results
# ==============================================================================
cat("Saving mapping results...\n")
write.csv(mapped_res, file = file.path(out_dir, 'mapping_results.csv'), 
          quote = FALSE, row.names = TRUE)

# Merge all query batches
# NOTE: If you have only one sample, skip the merge and use query.batches[[1]] directly
if (length(query.batches) > 1) {
  pancreas.query <- merge(query.batches[[1]], query.batches[2:length(query.batches)], 
                          merge.dr = "ref.umap")
} else {
  pancreas.query <- query.batches[[1]]
}

saveRDS(pancreas.query@reductions, 
        file.path(out_dir, 'reduction_scRNAseq_mapped.rds'))
write.csv(pancreas.query@meta.data, 
          file.path(out_dir, 'metadata_scRNAseq_mapped.csv'), 
          quote = FALSE, row.names = TRUE)

cat("Mapping complete!\n")
cat("Results saved to:", out_dir, "\n")

