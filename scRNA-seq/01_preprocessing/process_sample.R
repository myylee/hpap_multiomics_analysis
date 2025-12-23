#!/usr/bin/env Rscript
# ==============================================================================
# Process scRNA-seq Sample Data
#
# This script processes raw CellRanger output for a single scRNA-seq sample:
# 1. Loads 10X data using SoupX for ambient RNA correction
# 2. Removes doublets using scDblFinder
# 3. Filters cells based on quality metrics
# 4. Normalizes data using LogNormalize
# 5. Runs PCA
#
# Usage: Rscript process_sample.R <input_path> <output_path> [npc]
#   input_path: Path to CellRanger output directory
#   output_path: Path to save processed Seurat object
#   npc: Number of principal components for downstream analysis (default: 40)
# ==============================================================================

library(tidyverse)
library(fs)
library(stringr)
library(SoupX)
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript process_sample.R <input_path> <output_path> [npc]")
}

inPath <- args[1]
outPath <- args[2]
npc <- ifelse(length(args) >= 3, as.numeric(args[3]), 40)

options(Seurat.object.assay.version = "v3")

# Extract donor information from path
donorFull <- path_file(inPath)
arr <- str_split(donorFull, "-")[[1]]
donor <- arr[1]
source <- arr[2]
sample <- arr[3]

# ==============================================================================
# Step 1: Load and correct for ambient RNA using SoupX
# ==============================================================================
cat("Loading 10X data and correcting for ambient RNA...\n")
sc <- load10X(paste0(inPath, '/outs'))

pdf(paste0(outPath, '/autoSoupX.pdf'))
sc <- autoEstCont(sc, forceAccept = TRUE)
dev.off()

out <- adjustCounts(sc, roundToInt = TRUE)
rm(sc)
gc()

# ==============================================================================
# Step 2: Create Seurat object and calculate QC metrics
# ==============================================================================
cat("Creating Seurat object...\n")
seurat.tmp <- CreateSeuratObject(
  out,
  min.cells = 3,
  min.features = 200,
  project = donorFull
)
rm(out)
gc()

# Calculate mitochondrial percentage
seurat.tmp[["percent.mt"]] <- PercentageFeatureSet(seurat.tmp, pattern = "^MT-")

# ==============================================================================
# Step 3: Remove doublets
# ==============================================================================
cat("Removing doublets...\n")
sce <- as.SingleCellExperiment(seurat.tmp)
set.seed(123)
sce <- scDblFinder(sce)
stopifnot(all(colnames(sce) == colnames(seurat.tmp)))

seurat.tmp$scDblFinder.class <- sce$scDblFinder.class
seurat <- subset(seurat.tmp, subset = scDblFinder.class == "singlet")
rm(sce, seurat.tmp)
gc()

# ==============================================================================
# Step 4: Quality filtering
# ==============================================================================
cat("Applying quality filters...\n")
pdf(paste0(outPath, '/vlnPlotBefore.pdf'))
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Filter cells based on quality metrics
# Filters: 200 < nFeature_RNA < 10000, percent.mt < 25, 500 < nCount_RNA < 100000
seurat <- subset(
  seurat,
  subset = nFeature_RNA > 200 & nFeature_RNA < 10000 &
    percent.mt < 25 & nCount_RNA > 500 & nCount_RNA < 100000
)

pdf(paste0(outPath, '/vlnPlotAfter.pdf'))
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# ==============================================================================
# Step 5: LogNormalize and PCA
# ==============================================================================
cat("Running LogNormalize and PCA...\n")
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 3000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(
  seurat,
  assay = 'RNA',
  reduction.name = 'rna.pca',
  reduction.key = 'RNAPC_'
)

pdf(paste0(outPath, '/rna_elbow.pdf'))
e <- ElbowPlot(seurat, ndims = npc, reduction = 'rna.pca')
print(e)
dev.off()

# ==============================================================================
# Step 6: Save intermediate object
# ==============================================================================
cat("Saving intermediate object...\n")
saveRDS(seurat, paste0(outPath, '/seurat.tmp.rds'))

cat("Sample processing complete. Intermediate object saved.\n")
cat("Next step: Run cluster_and_annotate.R for clustering and annotation.\n")

