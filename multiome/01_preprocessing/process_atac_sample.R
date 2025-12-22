#!/usr/bin/env Rscript
# ==============================================================================
# Process ATAC-seq data for a single multiome sample
#
# This script processes raw CellRanger Arc output for ATAC-seq data:
# 1. Loads ATAC count matrix and fragments
# 2. Calls peaks using MACS2
# 3. Creates chromatin assay with called peaks
# 4. Removes doublets using scDblFinder and Amulet
# 5. Calculates QC metrics (TSS enrichment, nucleosome signal, etc.)
#
# Usage: Rscript process_atac_sample.R <input_path> <output_path> <macs2_path> <annotation_path>
#   input_path: Path to CellRanger Arc output directory
#   output_path: Path to save processed Seurat object
#   macs2_path: Path to MACS2 executable
#   annotation_path: Path to GRanges annotation object
# ==============================================================================

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(fs)
library(stringr)
library(SingleCellExperiment)
library(scDblFinder)
library(sctransform)
library(harmony)
library(cowplot)
library(GenomicRanges)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript process_atac_sample.R <input_path> <output_path> <macs2_path> <annotation_path>")
}

inPath <- args[1]
outPath <- args[2]
macs2Path <- args[3]
annotationPath <- args[4]

set.seed(123)

# ==============================================================================
# Step 1: Load 10X ATAC data
# ==============================================================================
cat("Loading 10X ATAC data...\n")
counts <- Read10X_h5(file.path(inPath, "/outs/filtered_feature_bc_matrix.h5"))
multiome_fragpath <- file.path(inPath, "/outs/atac_fragments.tsv.gz")

# Load annotation
# If annotation not present, create using:
# annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# ucsc.levels <- str_replace(string = paste("chr", seqlevels(annotation), sep = ""), 
#                            pattern = "chrMT", replacement = "chrM")
# seqlevels(annotation) <- ucsc.levels
# saveRDS(annotation, annotationPath)
annotation <- readRDS(annotationPath)

# ==============================================================================
# Step 2: Create chromatin assay and call peaks
# ==============================================================================
cat("Creating chromatin assay and calling peaks...\n")
chrom_assay <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = multiome_fragpath,
  annotation = annotation
)

seurat <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC"
)

# Call peaks using all cells with MACS2
DefaultAssay(seurat) <- "ATAC"
peaks <- CallPeaks(
  object = seurat,
  macs2.path = macs2Path
)
peaks <- seurat@assays$ATAC@ranges
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# ==============================================================================
# Step 3: Create MACS2 peak-based assay
# ==============================================================================
cat("Creating MACS2 peak-based assay...\n")
macs2_counts <- FeatureMatrix(
  fragments = CreateFragmentObject(multiome_fragpath),
  features = peaks,
  cells = rownames(seurat@meta.data)
)

seurat[["peak_macs2"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  sep = c("-", "-"),
  fragments = multiome_fragpath,
  annotation = annotation
)

write.csv(peaks, paste0(outPath, "/peak_macs2.csv"))

DefaultAssay(seurat) <- "peak_macs2"
saveRDS(seurat, paste0(outPath, '/seurat_ATAC.tmp.rds'))

# ==============================================================================
# Step 4: Remove doublets
# ==============================================================================
cat("Removing doublets...\n")
sce <- as.SingleCellExperiment(seurat)
set.seed(123)

# Create mock doublets for scDblFinder
sce <- mockDoubletSCE(ngenes = 300)
sce <- scDblFinder(
  sce, 
  artificialDoublets = 1, 
  aggregateFeatures = TRUE, 
  nfeatures = 25, 
  processing = "normFeatures"
)

# Use Amulet for additional doublet detection
fragfile <- multiome_fragpath
blacklist_to_remove <- blacklist_hg38_unified
otherChroms <- GRanges(
  c("M", "chrM", "MT", "X", "Y", "chrX", "chrY"),
  IRanges(1L, width = 10^8)
)
toExclude <- suppressWarnings(c(blacklist_to_remove, otherChroms))
res <- amulet(fragfile, regionsToExclude = toExclude)

# Combine doublet scores
res$scDblFinder.p <- 1 - colData(sce)[row.names(res), "scDblFinder.score"]
res$combined <- apply(
  res[, c("scDblFinder.p", "p.value")], 
  1, 
  FUN = function(x) {
    x[x < 0.001] <- 0.001  # Prevent too much skew from very small or 0 p-values
    suppressWarnings(aggregation::fisher(x))
  }
)

seurat$scDblFinder.combined.p <- res$combined
saveRDS(seurat, paste0(outPath, '/seurat_ATAC.tmp.rds'))

cat("ATAC preprocessing (part 1) complete. Intermediate object saved.\n")
cat("Next step: Run process_atac_sample_qc.R for quality control and filtering.\n")

