#!/usr/bin/env Rscript
# ==============================================================================
# Map snATAC-seq Query Data to Multiome Reference
#
# This script maps snATAC-seq query data to the multiome reference:
# 1. Prepares reference (multiome object) for mapping
# 2. Processes query snATAC-seq data
# 3. Finds transfer anchors using LSI reduction
# 4. Transfers cell type labels and projects onto reference UMAP
#
# Usage: Rscript map_snatacseq.R <reference_path> <query_path> <output_dir> [reference_umap]
#   reference_umap: Name of UMAP reduction in reference (default: wnn5filteredRun2.umap)
# ==============================================================================

library(Seurat)
library(Signac)
library(stringr)
library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript map_snatacseq.R <reference_path> <query_path> <output_dir> [reference_umap]")
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

DefaultAssay(reference) <- "macs2_combined"

# Compute LSI if not present
if (!"lsi" %in% names(reference@reductions)) {
  cat("Computing LSI for reference...\n")
  reference <- FindTopFeatures(reference, min.cutoff = 5)
  reference <- RunTFIDF(reference)
  reference <- RunSVD(reference)
}

# Set cell type annotation
if ("celltype2" %in% colnames(reference@meta.data)) {
  reference@meta.data$celltype <- reference$celltype2
} else if ("celltype" %in% colnames(reference@meta.data)) {
  # Use existing celltype if celltype2 not available
} else {
  stop("No cell type annotation found in reference object")
}

Idents(reference) <- 'celltype'

# Build reference neighbor graph
reference <- FindNeighbors(
  object = reference,
  reduction = "lsi",
  dims = 2:30,
  graph.name = "lsi.annoy.neighbors",
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

DefaultAssay(query) <- "peaks"  # or appropriate ATAC assay name

# Compute LSI for query
cat("Computing LSI for query...\n")
query <- FindTopFeatures(query, min.cutoff = 5)
query <- RunTFIDF(query)
query <- RunSVD(query)

# ==============================================================================
# Step 3: Find anchors and map query to reference
# ==============================================================================
cat("Finding transfer anchors...\n")
transfer.anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  reference.reduction = "lsi",
  reference.neighbors = "lsi.annoy.neighbors",
  dims = 2:30,
  k.filter = NA
)

cat("Mapping query to reference...\n")
query <- MapQuery(
  anchorset = transfer.anchors,
  query = query,
  reference = reference,
  refdata = list(celltype = "celltype"),
  reference.reduction = "lsi",
  reduction.model = reference_umap,
  transferdata.args = list(k.weight = 20)
)

Idents(query) <- 'predicted.celltype'

# ==============================================================================
# Step 4: Visualize and save results
# ==============================================================================
cat("Creating visualizations...\n")
pdf(file.path(out_dir, 'mapping_snATACseq.pdf'), width = 16, height = 12)
p1 <- DimPlot(query, reduction = "ref.umap", label = TRUE, label.size = 4,
              repel = TRUE, shuffle = TRUE) +
  NoLegend() +
  ggtitle("Mapped snATAC-seq labels")
p2 <- FeaturePlot(query, features = "predicted.celltype.score", combine = TRUE,
                  pt.size = 0.1, reduction = "ref.umap") &
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
p3 <- DimPlot(reference, reduction = reference_umap, label = TRUE, label.size = 4,
              group.by = 'celltype', repel = TRUE, shuffle = TRUE) +
  NoLegend() +
  ggtitle("Reference annotations")
gridExtra::grid.arrange(grobs = list(p1, p2, p3), nrow = 2)
dev.off()

# Save results
format_mapping_res <- function(query_seurat) {
  require(dplyr)
  df <- cbind(
    query_seurat@meta.data %>%
      select(predicted.celltype, predicted.celltype.score),
    query_seurat@reductions$ref.umap@cell.embeddings
  )
  return(df)
}

mapped_res <- format_mapping_res(query)
write.csv(mapped_res, file = file.path(out_dir, 'mapping_results.csv'),
          quote = FALSE, row.names = TRUE)

saveRDS(query@reductions, file.path(out_dir, 'reduction_snATACseq_mapped.rds'))
write.csv(query@meta.data, file.path(out_dir, 'metadata_snATACseq_mapped.csv'),
          quote = FALSE, row.names = TRUE)

cat("Mapping complete!\n")
cat("Results saved to:", out_dir, "\n")

