#!/usr/bin/env Rscript
# ==============================================================================
# Prepare Downsampled Reference for snATAC-seq Mapping
#
# This script prepares a downsampled multiome reference for mapping snATAC-seq
# query data:
# 1. Loads reference and downsamples to 1000 cells per donor
# 2. Re-runs WNN integration on downsampled reference
# 3. Prepares reference using Supervised LSI (SLSI)
# 4. Saves the prepared reference object
#
# The downsampled reference reduces computational requirements while maintaining
# mapping accuracy. This reference can be reused for multiple query samples.
#
# Usage: Rscript prepare_reference_for_mapping.R <reference_path> <output_path>
# ==============================================================================

library(Seurat)
library(Signac)
library(stringr)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript prepare_reference_for_mapping.R <reference_path> <output_path>")
}

referencePath <- args[1]
outputPath <- args[2]

# ==============================================================================
# Step 1: Load and downsample reference
# ==============================================================================
cat("Loading reference...\n")
reference <- readRDS(referencePath)

# Set cell type annotation
if ("celltype2" %in% colnames(reference@meta.data)) {
  reference@meta.data$celltype <- reference$celltype2
} else if ("celltype" %in% colnames(reference@meta.data)) {
  # Use existing celltype if celltype2 not available
} else {
  stop("No cell type annotation found in reference object")
}

# Downsample to 1000 cells per donor for mapping
cat("Downsampling reference to 1000 cells per donor...\n")
Idents(reference) <- "donor"
reference_1000perDonor <- subset(x = reference, downsample = 1000)

cat("Downsampled reference:", ncol(reference_1000perDonor), "cells (from", ncol(reference), ")\n")
cat("Number of donors:", length(unique(reference_1000perDonor$donor)), "\n")

# ==============================================================================
# Step 2: Re-run WNN integration on downsampled reference
# ==============================================================================
cat("Re-running WNN integration on downsampled reference...\n")
# Use Phase 2 reductions for WNN
rd_list <- list("combinedSctNormGEX.rpcaRefIntegrated.run2", "macs2Combined.harmony2")
key_i <- "5filteredRun2"

reference_1000perDonor <- FindMultiModalNeighbors(
  reference_1000perDonor,
  reduction.list = rd_list,
  dims.list = list(1:30, 2:30),
  knn.graph.name = paste0("wknn", key_i),
  snn.graph.name = paste0("wsnn", key_i),
  weighted.nn.name = paste0("weighted.nn", key_i)
)

reference_1000perDonor <- RunUMAP(
  reference_1000perDonor,
  nn.name = "weighted.nn5filteredRun2",
  reduction.name = "wnn5filteredRun2.umap",
  reduction.key = "wnn5filteredRun2Run2UMAP_",
  return.model = TRUE,
  seed.use = 1234
)

cat("WNN integration complete.\n")

# ==============================================================================
# Step 3: Prepare reference using Supervised LSI (SLSI)
# ==============================================================================
cat("Preparing reference for mapping using Supervised LSI...\n")
DefaultAssay(reference_1000perDonor) <- "macs2_combined"

# Join layers if reference was saved with layers
reference_1000perDonor <- JoinLayers(reference_1000perDonor)

# Find top features and convert assay class
reference_1000perDonor <- FindTopFeatures(reference_1000perDonor, min.cutoff = 'q50')
reference_1000perDonor[["macs2_combined"]] <- as(
  object = reference_1000perDonor[["macs2_combined"]],
  Class = "Assay"
)

# Run Supervised LSI using WNN graph
reference_1000perDonor <- RunSLSI(
  reference_1000perDonor,
  assay = 'macs2_combined',
  graph = 'wsnn5filteredRun2',
  n = 30,
  features = VariableFeatures(reference_1000perDonor)
)

Idents(reference_1000perDonor) <- 'celltype'

cat("SLSI reduction complete.\n")

# ==============================================================================
# Step 4: Save prepared reference
# ==============================================================================
cat("Saving prepared reference to:", outputPath, "\n")
saveRDS(reference_1000perDonor, outputPath)

cat("\nReference preparation complete!\n")
cat("Downsampled reference saved with:\n")
cat("  -", ncol(reference_1000perDonor), "cells\n")
cat("  - WNN integration (wnn5filteredRun2)\n")
cat("  - SLSI reduction for mapping\n")
cat("  - Cell type annotations (celltype column)\n")
cat("\nThis reference can now be used with map_snatacseq.R for mapping query samples.\n")

