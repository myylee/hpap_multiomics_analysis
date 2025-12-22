#!/usr/bin/env Rscript
# ==============================================================================
# Process GEX (Gene Expression) data for a single multiome sample
# 
# This script processes raw CellRanger Arc output for RNA-seq data:
# 1. Loads 10X data using SoupX for ambient RNA correction
# 2. Removes doublets using scDblFinder
# 3. Filters cells based on quality metrics
# 4. Normalizes using SCTransform
# 5. Runs PCA and basic clustering
# 6. Annotates cells using scSorter
#
# Usage: Rscript process_gex_sample.R <input_path> <output_path> [npc]
#   npc: Number of principal components for downstream analysis (default: 30)
#   input_path: Path to CellRanger Arc output directory
#   output_path: Path to save processed Seurat object
#   npc: Number of principal components for downstream analysis (default: 30)
# ==============================================================================

# Load required libraries
library(tidyverse)
library(fs)
library(stringr)
library(SoupX)
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(sctransform)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript process_gex_sample.R <input_path> <output_path> [npc]")
}

inPath <- args[1]
outPath <- args[2]
npc <- ifelse(length(args) >= 3, as.numeric(args[3]), 30)

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

pdf(paste0(outPath, '/GEX_autoSoupX.pdf'))
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
pdf(paste0(outPath, '/GEX_vlnPlotBefore.pdf'))
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Filter cells based on quality metrics
# Changed percent.mt to 30 from 25, and nCount_RNA 1000 from 500
seurat <- subset(
  seurat, 
  subset = nFeature_RNA > 200 & nFeature_RNA < 10000 &
    percent.mt < 30 & nCount_RNA > 1000 & nCount_RNA < 100000
)

pdf(paste0(outPath, '/GEX_vlnPlotAfter.pdf'))
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# ==============================================================================
# Step 5: SCTransform normalization
# ==============================================================================
cat("Running SCTransform...\n")
seurat <- SCTransform(
  seurat, 
  vst.flavor = "v2", 
  vars.to.regress = "percent.mt", 
  verbose = FALSE
)

# ==============================================================================
# Step 6: Dimensionality reduction and clustering
# ==============================================================================
cat("Running PCA...\n")
seurat <- RunPCA(
  seurat, 
  assay = 'SCT', 
  reduction.name = 'sctNormGEX.pca', 
  verbose = FALSE, 
  reduction.key = 'SCTNORMGEXPC_'
)

pdf(paste0(outPath, '/GEX_sct_elbow.pdf'))
e <- ElbowPlot(seurat, ndims = 50, reduction = 'sctNormGEX.pca')
print(e)
dev.off()

# Save intermediate object
saveRDS(seurat, paste0(outPath, '/seurat_GEX.tmp.rds'))

cat("GEX preprocessing complete. Intermediate object saved.\n")
cat("Next step: Run process_gex_sample_annotation.R for cell type annotation.\n")

