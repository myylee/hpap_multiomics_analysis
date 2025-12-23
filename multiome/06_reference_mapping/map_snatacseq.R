#!/usr/bin/env Rscript
# ==============================================================================
# Map snATAC-seq Query Data to Downsampled Multiome Reference
#
# This script maps snATAC-seq query data to a prepared downsampled multiome
# reference (1000CellsPerDonor). The reference should be prepared using
# prepare_reference_for_mapping.R, which downsamples the multiome reference
# to 1000 cells per donor, re-runs WNN integration, and prepares it with
# Supervised LSI (SLSI) for mapping.
#
# Steps:
# 1. Loads prepared 1000CellsPerDonor reference (created by prepare_reference_for_mapping.R)
# 2. Processes query snATAC-seq data
# 3. Finds transfer anchors using SLSI reduction
# 4. Transfers cell type labels and projects onto reference UMAP
# 5. Visualizes and saves mapping results
#
# Usage: Rscript map_snatacseq.R <reference_path> <query_path> <output_dir> [query_assay] [reference_umap]
#   reference_path: Path to prepared 1000CellsPerDonor reference (from prepare_reference_for_mapping.R)
#   query_path: Path to query snATAC-seq Seurat object
#   output_dir: Directory to save mapping results
#   query_assay: Name of ATAC assay in query object (default: peaks)
#   reference_umap: Name of UMAP reduction in reference (default: wnn5filteredRun2.umap)
# ==============================================================================

library(Seurat)
library(Signac)
library(stringr)
library(tidyverse)
library(data.table)
library(gridExtra)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript map_snatacseq.R <reference_path> <query_path> <output_dir> [query_assay] [reference_umap]")
}

referencePath <- args[1]
queryPath <- args[2]
out_dir <- args[3]
query_assay <- ifelse(length(args) >= 4, args[4], "peaks")
reference_umap <- ifelse(length(args) >= 5, args[5], "wnn5filteredRun2.umap")

dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# Step 1: Load prepared reference
# ==============================================================================
cat("Loading prepared downsampled reference...\n")
reference <- readRDS(referencePath)

# Verify reference has required components
required_reductions <- c("slsi", reference_umap)
missing_reductions <- setdiff(required_reductions, names(reference@reductions))
if (length(missing_reductions) > 0) {
  stop("Reference missing required reductions: ", paste(missing_reductions, collapse = ", "),
       "\nPlease run prepare_reference_for_mapping.R first.")
}

if (!"celltype" %in% colnames(reference@meta.data)) {
  stop("Reference missing 'celltype' annotation. Please check reference object.")
}

cat("Reference loaded:", ncol(reference), "cells\n")
cat("Reference reductions:", paste(names(reference@reductions), collapse = ", "), "\n")

# ==============================================================================
# Step 2: Prepare query data
# ==============================================================================
cat("Loading and preparing query data...\n")
query <- readRDS(file = queryPath)

# Set assay (should match the ATAC assay name in query object)
DefaultAssay(query) <- query_assay

# Filter peaks with zero counts
cat("Filtering peaks with zero counts...\n")
query <- query[rowSums(query[[query_assay]]@counts) > 0, ]

cat("Query data:", ncol(query), "cells,", nrow(query), "peaks\n")

# Compute LSI for query
cat("Computing LSI for query...\n")
query <- RunTFIDF(query)
query <- FindTopFeatures(query, min.cutoff = 5)
query <- RunSVD(query)

# ==============================================================================
# Step 3: Find anchors and map query to reference
# ==============================================================================
cat("Finding transfer anchors using SLSI...\n")
transfer.anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  reference.assay = "macs2_combined",
  query.assay = query_assay,
  reference.reduction = "slsi",
  reduction = "lsiproject",
  dims = 2:30
)

cat("Mapping query to reference...\n")
query <- MapQuery(
  anchorset = transfer.anchors,
  query = query,
  reference = reference,
  refdata = list(ct = "celltype"),
  reference.reduction = "slsi",
  reduction.model = reference_umap,
  transferdata.args = list(k.weight = 20)
)

Idents(query) <- 'predicted.ct'

# ==============================================================================
# Step 4: Visualize and save results
# ==============================================================================
cat("Creating visualizations...\n")

# Identify high-confidence predictions (score > 0.8)
cells_high_confidence <- colnames(query)[query$predicted.ct.score > 0.8]
cat("High-confidence predictions (score > 0.8):", length(cells_high_confidence), "cells\n")

pdf(file.path(out_dir, 'mapping_snATACseq.pdf'), width = 11, height = 11)
p1 <- DimPlot(query, reduction = "ref.umap", label = TRUE, label.size = 4,
              repel = TRUE, shuffle = TRUE, group.by = "predicted.ct") +
  NoLegend() +
  ggtitle("Predicted annotation")
p2 <- FeaturePlot(query, features = "predicted.ct.score", reduction = "ref.umap",
                  pt.size = 0.1) +
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle("Prediction score")
p3 <- DimPlot(query, reduction = "ref.umap", label = TRUE, label.size = 4,
              repel = TRUE, shuffle = TRUE, group.by = "predicted.ct",
              cells = cells_high_confidence) +
  NoLegend() +
  ggtitle("Predicted annotation (score > 0.8)")
p4 <- DimPlot(reference, reduction = reference_umap, label = TRUE,
              label.size = 4, group.by = 'celltype', repel = TRUE,
              shuffle = TRUE) +
  NoLegend() +
  ggtitle("Reference annotations")
gridExtra::grid.arrange(grobs = list(p1, p2, p3, p4), nrow = 2)
dev.off()

# Summary of prediction scores
cat("\nPrediction score summary:\n")
print(summary(query$predicted.ct.score))
cat("\nHigh-confidence predictions (score > 0.8):", 
    sum(query$predicted.ct.score > 0.8, na.rm = TRUE), "cells\n")

# Save results
format_mapping_res <- function(query_seurat) {
  require(dplyr)
  df <- cbind(
    query_seurat@meta.data %>%
      select(predicted.ct, predicted.ct.score),
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

cat("\nMapping complete!\n")
cat("Results saved to:", out_dir, "\n")
