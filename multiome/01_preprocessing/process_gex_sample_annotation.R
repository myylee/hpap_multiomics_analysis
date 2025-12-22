#!/usr/bin/env Rscript
# ==============================================================================
# Complete GEX sample processing
#
# This script continues from process_gex_sample.R:
# 1. Runs UMAP and clustering
# 2. Annotates cells using scSorter with marker genes
# 3. Adds donor metadata
# 4. Saves final processed object
#
# Usage: Rscript process_gex_sample_annotation.R <path> [npc]
# ==============================================================================

library(tidyverse)
library(Seurat)
library(sctransform)
library(SingleCellExperiment)
library(scSorter)
library(fs)
library(stringr)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript process_gex_sample_annotation.R <path> [npc]")
}

path <- args[1]
npc <- ifelse(length(args) >= 2, as.numeric(args[2]), 30)

# Marker genes for broad cell type annotation (from Golnaz marker list, excluding c11)
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

# Load intermediate object
seurat <- readRDS(paste0(path, '/seurat_GEX.tmp.rds'))

# ==============================================================================
# Step 1: UMAP and clustering
# ==============================================================================
cat("Running UMAP and clustering...\n")
seurat <- RunUMAP(
  seurat, 
  reduction = "sctNormGEX.pca", 
  reduction.key = 'SCTNORMGEXUMAP_',
  dims = 1:npc, 
  verbose = FALSE, 
  reduction.name = 'sctNormGEX.umap'
)

seurat <- FindNeighbors(
  seurat, 
  reduction = "sctNormGEX.pca", 
  dims = 1:npc, 
  verbose = FALSE, 
  graph.name = c('sct.nn', 'sct.snn')
)

seurat <- FindClusters(
  seurat, 
  resolution = 0.5, 
  method = "igraph", 
  verbose = FALSE, 
  graph.name = 'sct.snn'
)
seurat$sct.clusters <- seurat$seurat_clusters

# ==============================================================================
# Step 2: Cell type annotation using scSorter
# ==============================================================================
cat("Annotating cells using scSorter...\n")
sct.expr <- seurat[["SCT"]]$data
sct.markers <- markers[markers$Marker %in% rownames(sct.expr), ]
DefaultAssay(seurat) <- 'SCT'
sct.topgenes <- VariableFeatures(seurat)
sct.topgene.filter <- rowSums(sct.expr[sct.topgenes, ] != 0) > ncol(sct.expr) * 0.1
sct.topgenes <- sct.topgenes[sct.topgene.filter]
sct.picked.genes <- unique(c(sct.markers$Marker, sct.topgenes))
sct.expr <- sct.expr[rownames(sct.expr) %in% sct.picked.genes, ]
sct.rts <- scSorter(sct.expr, sct.markers)
sct.preds <- sct.rts$Pred_Type
sct.preds[sct.preds == 'Unknown'] <- 'unknown'
seurat$donorScSorter.sct <- sct.preds

cat("Cell type distribution:\n")
print(table(seurat$donorScSorter.sct))

# ==============================================================================
# Step 3: Visualization
# ==============================================================================
cat("Creating visualization plots...\n")
pdf(paste0(path, '/GEX_sctNormGEX_umap.pdf'))
d1 <- DimPlot(
  seurat, 
  group.by = "sct.clusters", 
  reduction = "sctNormGEX.umap", 
  label = TRUE
)
d2 <- DimPlot(
  seurat, 
  group.by = "donorScSorter.sct", 
  reduction = "sctNormGEX.umap", 
  label = FALSE
)
gridExtra::grid.arrange(grobs = list(d1, d2), ncol = 2)
dev.off()

# ==============================================================================
# Step 4: Add metadata and rename cells
# ==============================================================================
donorFull <- path_file(path)
arr <- str_split(donorFull, "-")[[1]]
donor <- arr[1]

seurat@meta.data$donor <- donor
seurat@meta.data$modality <- 'multiomeRNA'
seurat <- RenameCells(seurat, new.names = paste0(donor, "_", colnames(seurat)))

# ==============================================================================
# Step 5: Save final object and metadata
# ==============================================================================
saveRDS(seurat, file = paste0(path, '/seurat_sctNormGEX.rds'))
data.table::fwrite(
  seurat@meta.data, 
  file = paste0(path, '/metadata_sctNormGEX.csv'),
  row.names = TRUE,
  quote = FALSE
)

cat("GEX sample processing complete!\n")

