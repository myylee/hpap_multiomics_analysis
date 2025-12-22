#!/usr/bin/env Rscript
# ==============================================================================
# Phase 1: Combine GEX and ATAC Data
#
# This script combines the integrated GEX and ATAC objects and adds metadata:
# 1. Loads unintegrated GEX object
# 2. Adds ATAC assays (Harmony and RLSI integrated)
# 3. Adds GEX integration metadata and reductions
# 4. Adds ATAC integration metadata and reductions
# 5. Calculates additional QC metrics
#
# This prepares the data for WNN integration.
#
# Usage: Rscript combine_gex_atac.R <input_path> [donor_metadata_path]
# ==============================================================================

library(Seurat)
library(stringr)
library(sctransform)
library(harmony)
library(tidyverse)
library(data.table)
library(Signac)
library(future)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript combine_gex_atac.R <input_path> [donor_metadata_path]")
}

inPath <- args[1]
donorMetadataPath <- ifelse(length(args) >= 2, args[2], NA)
n_workers <- ifelse(length(args) >= 3, as.numeric(args[3]), 4)

plan("multisession", workers = n_workers)
options(future.globals.maxSize = 200000 * 1024^2)  # 200GB

# ==============================================================================
# Step 1: Load GEX unintegrated object
# ==============================================================================
cat("Loading GEX unintegrated object...\n")
seurat_combined_v5 <- readRDS(paste0(inPath, '/seurat_sctNormGEX_unintegrated.rds'))

# ==============================================================================
# Step 2: Add ATAC assays
# ==============================================================================
cat("Adding ATAC assays...\n")

# Harmony integrated ATAC
atac_harmonyIntegrated <- readRDS(paste0(inPath, '/seurat_macs2Combined_wHarmony.rds'))
atac_harmonyIntegrated <- JoinLayers(atac_harmonyIntegrated)
seurat_combined_v5[["macs2_combined"]] <- atac_harmonyIntegrated[['macs2_combined']]

atac_harmonyIntegrated$macs2Combined_harmony_clusters <- factor(
  as.character(atac_harmonyIntegrated$macs2Combined_harmony_clusters),
  levels = as.character(sort(as.integer(levels(atac_harmonyIntegrated$macs2Combined_harmony_clusters))))
)

# RLSI integrated ATAC
atac_rlsiRefIntegrated <- readRDS(file = paste0(inPath, '/ATAC_macs2Combined_rlsiRefIntegrated.rds'))
atac_rlsiRefIntegrated$macs2Combined_rlsiRefIntegrated_clusters <- factor(
  as.character(atac_rlsiRefIntegrated$macs2Combined_rlsiRefIntegrated_clusters),
  levels = as.character(sort(as.integer(levels(atac_rlsiRefIntegrated$macs2Combined_rlsiRefIntegrated_clusters))))
)

# ==============================================================================
# Step 3: Add GEX integration metadata and reductions
# ==============================================================================
cat("Adding GEX integration metadata and reductions...\n")

# RPCA All
md_rpcaAll <- readRDS(paste0(inPath, '/metadata_sctNormGEX_rpcaAllIntegrated.rds'))
reduction_rpcaAll <- readRDS(paste0(inPath, '/reduction_sctNormGEX_rpcaAllIntegrated.rds'))
md_rpcaAll$combined_sct_rpcaAllIntegrated_clusters <- factor(
  as.character(md_rpcaAll$combined_sct_rpcaAllIntegrated_clusters),
  levels = as.character(sort(as.integer(levels(md_rpcaAll$combined_sct_rpcaAllIntegrated_clusters))))
)

# RPCA Ref
md_rpcaRef <- readRDS(paste0(inPath, '/metadata_sctNormGEX_rpcaRefIntegrated.rds'))
reduction_rpcaRef <- readRDS(paste0(inPath, '/reduction_sctNormGEX_rpcaRefIntegrated.rds'))
md_rpcaRef$combined_sct_rpcaRefIntegrated_clusters <- factor(
  as.character(md_rpcaRef$combined_sct_rpcaRefIntegrated_clusters),
  levels = as.character(sort(as.integer(levels(md_rpcaRef$combined_sct_rpcaRefIntegrated_clusters))))
)

# Harmony
md_harmony <- readRDS(paste0(inPath, '/metadata_sctNormGEX_harmonyIntegrated.rds'))
reduction_harmony <- readRDS(paste0(inPath, '/reduction_sctNormGEX_harmonyIntegrated.rds'))
md_harmony$combined_sct_harmony_clusters <- factor(
  as.character(md_harmony$combined_sct_harmony_clusters),
  levels = as.character(sort(as.integer(levels(md_harmony$combined_sct_harmony_clusters))))
)

# ==============================================================================
# Step 4: Combine all metadata
# ==============================================================================
cat("Combining metadata...\n")
seurat <- seurat_combined_v5

# Add donor metadata if available
if (!is.na(donorMetadataPath) && file.exists(donorMetadataPath)) {
  donorMetadata <- read.csv(donorMetadataPath)
  seurat@meta.data <- cbind(
    seurat@meta.data,
    seurat@meta.data %>%
      left_join(donorMetadata, by = c("donor" = 'Sample')) %>%
      select(Gender:disease_wAge)
  )
}

# Add GEX cluster assignments
seurat@meta.data <- cbind(
  seurat@meta.data,
  "combined_sct_rpcaAllIntegrated_clusters" = md_rpcaAll[rownames(seurat@meta.data), "combined_sct_rpcaAllIntegrated_clusters"],
  "combined_sct_rpcaRefIntegrated_clusters" = md_rpcaRef[rownames(seurat@meta.data), "combined_sct_rpcaRefIntegrated_clusters"],
  "combined_sct_harmony_clusters" = md_harmony[rownames(seurat@meta.data), "combined_sct_harmony_clusters"],
  "macs2Combined_harmony_clusters" = atac_harmonyIntegrated$macs2Combined_harmony_clusters,
  "macs2Combined_rlsiRefIntegrated_clusters" = atac_rlsiRefIntegrated$macs2Combined_rlsiRefIntegrated_clusters
)

# Add reductions
seurat@reductions <- c(
  seurat@reductions,
  reduction_rpcaAll[3:4],
  reduction_rpcaRef[3:4],
  reduction_harmony[3:4],
  atac_harmonyIntegrated@reductions[3:4],
  atac_rlsiRefIntegrated@reductions[1:2]
)

# ==============================================================================
# Step 5: Add ATAC metadata and calculate QC metrics
# ==============================================================================
cat("Adding ATAC metadata and calculating QC metrics...\n")
md_rlsiRef <- atac_rlsiRefIntegrated@meta.data
atac_col_idx <- c(2, 3, 5:18)
atac_colnames_picked <- colnames(md_rlsiRef)[atac_col_idx]

md_wATACInfo <- cbind(
  seurat@meta.data,
  md_rlsiRef[rownames(seurat@meta.data), atac_colnames_picked]
)

seurat@meta.data <- md_wATACInfo

# Calculate log-transformed QC metrics
seurat@meta.data$log10nFeature_RNA <- log10(seurat$nFeature_RNA)
seurat@meta.data$log10nCount_RNA <- log10(seurat$nCount_RNA)
seurat@meta.data$log10GenesPerUMI <- log10(seurat$nFeature_RNA) / log10(seurat$nCount_RNA)
seurat@meta.data$log10FreqCount <- log10(seurat$frequency_count)

# ==============================================================================
# Step 6: Save combined object
# ==============================================================================
cat("Saving combined object...\n")
saveRDS(seurat, paste0(inPath, '/seurat_gexAndAtac.rds'))

cat("Combined GEX and ATAC object saved!\n")
cat("Object:", paste0(inPath, '/seurat_gexAndAtac.rds'), "\n")
cat("Next step: Run integrate_wnn_phase1.R for WNN integration.\n")

