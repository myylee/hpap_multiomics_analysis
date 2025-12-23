#!/usr/bin/env Rscript
# ==============================================================================
# Cluster and Annotate scRNA-seq Sample
#
# This script performs clustering and cell type annotation for a processed
# scRNA-seq sample:
# 1. Runs clustering and UMAP using RNA (LogNormalize) assay
# 2. Performs cell type annotation using scSorter with marker genes
#    (Note: scSorter annotation is a first pass to get a general idea of cell types)
# 3. Visualizes annotations
#
# Usage: Rscript cluster_and_annotate.R <path> <npc>
#   path: Path to directory containing seurat.tmp.rds (from process_sample.R)
#   npc: Number of principal components to use for clustering/UMAP
# ==============================================================================

library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(scSorter)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript cluster_and_annotate.R <path> <npc>")
}

path <- args[1]
npc <- as.numeric(args[2])

options(Seurat.object.assay.version = "v3")

# Marker genes for broad cell type annotation (from Golnaz marker list)
markers <- data.frame(
  Type = c(
    "acinar", "acinar", "acinar", "acinar",
    "alpha", "alpha", "alpha",
    "beta", "beta", "beta",
    "ductal", "ductal", "ductal", "ductal",
    "stellates.mesenchymal", "stellates.mesenchymal", "stellates.mesenchymal", "stellates.mesenchymal",
    "endothelial", "endothelial",
    "immune", "immune",
    "pp.gamma", "pp.gamma", "pp.gamma", "pp.gamma",
    "delta", "delta", "delta",
    "epsilon"
  ),
  Marker = c(
    "PRSS1", "REG1A", "CPA1", "CPA2",
    "GCG", "GC", "TTR",
    "INS", "IGF2", "IAPP",
    "KRT19", "CFTR", "SFRP5", "MMP7",
    "COL1A1", "PDGFRB", "RGS10", "THY1",
    "VWF", "CD93",
    "NCF2", "PTPRC",
    "PPY", "CARTPT", "PCDH10", "PLAC8",
    "SST", "LEPR", "PRG4",
    "GHRL"
  ),
  stringsAsFactors = FALSE
)

# Define colors for each cell type
cell_types <- c(as.character(levels(as.factor(markers$Type))), 'unknown')
colors <- data.frame(
  Type = cell_types,
  col = c("#e30800", "#f56505", "#dec400", "#006630", "#0223c7", "#5b02c7",
          "#00b0e6", "#c40080", "#02f00a", "#7d3301", "#000000")
)

# Load processed object
cat("Loading processed object...\n")
seurat <- readRDS(paste0(path, '/seurat.tmp.rds'))

# ==============================================================================
# Step 1: RNA (LogNormalize) clustering and UMAP
# ==============================================================================
cat("Running RNA (LogNormalize) clustering and UMAP...\n")
seurat <- FindNeighbors(
  seurat,
  reduction = "rna.pca",
  dims = 1:npc,
  verbose = FALSE,
  graph.name = c('rna.nn', 'rna.snn')
)
seurat <- FindClusters(
  seurat,
  resolution = 0.5,
  method = "igraph",
  verbose = FALSE,
  graph.name = 'rna.snn'
)
seurat$rna.clusters <- seurat$seurat_clusters
seurat <- RunUMAP(
  seurat,
  reduction = "rna.pca",
  reduction.key = 'RNAUMAP_',
  dims = 1:npc,
  verbose = FALSE,
  reduction.name = 'rna.umap'
)

# ==============================================================================
# Step 2: Cell type annotation using scSorter
# Note: scSorter annotation is a first pass to get a general idea of cell types
# ==============================================================================
cat("Annotating cell types using scSorter (first pass annotation)...\n")
rna.expr <- seurat[['RNA']]$data
rna.markers <- markers[markers$Marker %in% rownames(rna.expr), ]
DefaultAssay(seurat) <- 'RNA'
rna.topgenes <- VariableFeatures(seurat)
rna.topgene.filter <- rowSums(rna.expr[rna.topgenes, ] != 0) > ncol(rna.expr) * 0.1
rna.topgenes <- rna.topgenes[rna.topgene.filter]
rna.picked.genes <- unique(c(rna.markers$Marker, rna.topgenes))
rna.expr <- rna.expr[rownames(rna.expr) %in% rna.picked.genes, ]
rna.rts <- scSorter(rna.expr, rna.markers)
rna.preds <- rna.rts$Pred_Type
rna.preds[rna.preds == 'Unknown'] <- 'unknown'
seurat$donorScSorter.rna <- rna.preds
colors.new <- colors[colors$Type %in% as.character(levels(as.factor(seurat$donorScSorter.rna))), ]$col

pdf(paste0(path, '/rna_umap.pdf'))
d <- DimPlot(
  seurat,
  group.by = "donorScSorter.rna",
  reduction = "rna.umap",
  label = FALSE,
  cols = colors.new
)
print(d)
dev.off()

# ==============================================================================
# Step 3: Save final object
# ==============================================================================
DefaultAssay(seurat) <- 'RNA'
saveRDS(seurat, file = paste0(path, '/seurat.rds'))

cat("\nClustering and annotation complete!\n")
cat("Final object saved to:", paste0(path, '/seurat.rds'), "\n")

