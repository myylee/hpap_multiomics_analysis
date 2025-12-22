#!/usr/bin/env Rscript
# ==============================================================================
# Quality control and filtering for ATAC-seq sample
#
# This script continues from process_atac_sample.R:
# 1. Filters doublets based on combined p-value
# 2. Calculates QC metrics (TSS enrichment, nucleosome signal, etc.)
# 3. Filters cells based on quality thresholds
# 4. Runs LSI dimensionality reduction
# 5. Adds donor metadata and saves final object
#
# Usage: Rscript process_atac_sample_qc.R <path> <npc>
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
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript process_atac_sample_qc.R <path> [npc]")
}

path <- args[1]
npc <- ifelse(length(args) >= 2, as.numeric(args[2]), 30)

# ==============================================================================
# Step 1: Load and filter doublets
# ==============================================================================
cat("Loading intermediate object and filtering doublets...\n")
seurat <- readRDS(paste0(path, '/seurat_ATAC.tmp.rds'))
ncells_orig <- ncol(seurat)

if (length(which(colnames(seurat@meta.data) == "scDblFinder.combined.p")) > 0) {
  seurat <- subset(seurat, subset = scDblFinder.combined.p >= 0.05)
  ncells_woDoublet <- ncol(seurat)
  dRemoved <- ncells_orig - ncells_woDoublet
  dRemovedPerc <- round(dRemoved / ncells_orig * 100, digits = 1)
} else {
  dRemoved <- 0
  dRemovedPerc <- 0
  ncells_woDoublet <- ncells_orig
}

cat(paste(dRemoved, paste0("(", dRemovedPerc, "%)"), "doublets removed\n", sep = " "))

# ==============================================================================
# Step 2: Calculate QC metrics
# ==============================================================================
cat("Calculating QC metrics...\n")
seurat <- NucleosomeSignal(seurat)
seurat <- TSSEnrichment(seurat)

# Count fragments per cell
frag_counts <- CountFragments(
  seurat@assays$ATAC@fragments[[1]]@path,
  colnames(seurat)
) %>% column_to_rownames('CB')

seurat@meta.data <- cbind(
  seurat@meta.data,
  frag_counts[colnames(seurat), ]
)

seurat$pct_reads_in_peaks <- seurat$nCount_peak_macs2 / seurat$frequency_count * 100
seurat$blacklist_fraction <- FractionCountsInRegion(
  object = seurat,
  assay = 'ATAC',
  regions = blacklist_hg38_unified
)

# ==============================================================================
# Step 3: Visualize QC metrics before filtering
# ==============================================================================
pdf(paste0(path, '/ATAC_densityScatter.pdf'), width = 8, height = 8)
DensityScatter(
  seurat, 
  x = 'nCount_peak_macs2', 
  y = 'TSS.enrichment', 
  log_x = TRUE, 
  quantiles = TRUE
)
dev.off()

p1 <- VlnPlot(seurat, features = "nCount_peak_macs2") + 
  geom_hline(yintercept = 1000, color = 'red') + 
  geom_hline(yintercept = 50000, color = 'red')
p2 <- VlnPlot(seurat, features = "TSS.enrichment") + 
  geom_hline(yintercept = 3, color = 'red')
p3 <- VlnPlot(seurat, features = "nucleosome_signal") + 
  geom_hline(yintercept = 4, color = 'red')
p4 <- VlnPlot(seurat, features = "pct_reads_in_peaks") + 
  geom_hline(yintercept = 15, color = 'red')
p5 <- VlnPlot(seurat, features = "blacklist_fraction") + 
  geom_hline(yintercept = 0.05, color = 'red')

pdf(paste0(path, '/ATAC_vlnPlotBefore.pdf'), width = 8, height = 15)
patchwork::wrap_plots(p1, p2, p3, p4, p5, ncol = 5) + NoLegend()
dev.off()

# ==============================================================================
# Step 4: Filter cells based on QC metrics
# ==============================================================================
cat("Filtering cells based on QC metrics...\n")
seurat <- subset(
  seurat,
  subset = nucleosome_signal < 2 & TSS.enrichment > 1 &
    nCount_peak_macs2 > 1000 & nCount_peak_macs2 < 50000 &
    blacklist_fraction < 0.05
)

ncells_final <- ncol(seurat)
cFiltered <- ncells_woDoublet - ncells_final
cFilteredPerc <- round(cFiltered / ncells_woDoublet * 100, digits = 1)

cat(paste(cFiltered, paste0("(", cFilteredPerc, "%)"), "cells filtered\n", sep = " "))

pdf(paste0(path, '/ATAC_vlnPlotAfter.pdf'), width = 8, height = 15)
VlnPlot(
  seurat, 
  features = c("nCount_peak_macs2", "TSS.enrichment", "nucleosome_signal",
               "pct_reads_in_peaks", "blacklist_fraction"), 
  ncol = 5
)
dev.off()

# ==============================================================================
# Step 5: LSI dimensionality reduction
# ==============================================================================
cat("Running LSI dimensionality reduction...\n")
DefaultAssay(seurat) <- "peak_macs2"
seurat <- seurat[rowSums(seurat) > 0, ]
seurat <- RunTFIDF(seurat)
seurat <- FindTopFeatures(seurat, min.cutoff = 5)
seurat <- RunSVD(
  seurat, 
  assay = 'peak_macs2', 
  reduction.name = 'peakMacs2.lsi', 
  verbose = FALSE, 
  reduction.key = 'PEAKMACS2LSI_'
)

pdf(paste0(path, '/ATAC_lsi_elbow.pdf'))
e <- ElbowPlot(seurat, ndims = 50, reduction = 'peakMacs2.lsi')
print(e)
dev.off()

pdf(paste0(path, '/ATAC_lsi_depthCor.pdf'))
DepthCor(seurat, reduction = 'peakMacs2.lsi')
dev.off()

# ==============================================================================
# Step 6: UMAP and clustering
# ==============================================================================
cat("Running UMAP and clustering...\n")
seurat <- RunUMAP(
  seurat, 
  reduction = "peakMacs2.lsi", 
  reduction.key = 'PEAKMACS2UMAP_',
  dims = 2:npc, 
  verbose = FALSE, 
  reduction.name = 'peakMacs2.umap'
)

seurat <- FindNeighbors(
  seurat, 
  reduction = "peakMacs2.lsi", 
  dims = 1:npc, 
  verbose = FALSE, 
  graph.name = c('peakMacs2.nn', 'peakMacs2.snn')
)

seurat <- FindClusters(
  seurat, 
  resolution = 0.5, 
  method = "igraph", 
  verbose = FALSE, 
  graph.name = 'peakMacs2.snn'
)
seurat$peakMacs2.clusters <- seurat$seurat_clusters

pdf(paste0(path, '/ATAC_peakMacs2_umap.pdf'))
d <- DimPlot(
  seurat, 
  group.by = "peakMacs2.clusters", 
  reduction = "peakMacs2.umap", 
  label = TRUE
)
print(d)
dev.off()

# ==============================================================================
# Step 7: Add metadata and save
# ==============================================================================
donorFull <- fs::path_file(path)
arr <- stringr::str_split(donorFull, "-")[[1]]
donor <- arr[1]

seurat@meta.data$donor <- donor
seurat@meta.data$modality <- 'multiomeATAC'
seurat <- RenameCells(seurat, new.names = paste0(donor, "_", colnames(seurat)))

data.table::fwrite(
  seurat@meta.data, 
  file = paste0(path, '/metadata_peakMacs2ATAC.csv'),
  row.names = TRUE,
  quote = FALSE
)

saveRDS(seurat, paste0(path, '/seurat_ATAC.rds'))

cat("ATAC sample processing complete!\n")

