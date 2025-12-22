#!/usr/bin/env Rscript
# ==============================================================================
# Phase 1: Evaluate WNN Integration Results
#
# This script generates comprehensive cluster quality plots for each WNN
# integration result (wnn2-wnn6) to facilitate evaluation and selection of
# the best integration method.
#
# The evaluation plots include:
# 1. Marker gene expression (DotPlots and FeaturePlots)
# 2. Cluster quality metrics (VlnPlots)
# 3. UMAP visualizations with per-modality annotations
# 4. Comparison of WNN clusters with per-modality cluster assignments
#
# These plots are used for manual annotation and to select which WNN result
# provides the best separation and biological coherence.
#
# IMPORTANT: This script requires that per-modality clusters have been
# annotated (run annotate_per_modality_clusters.R first).
#
# Usage: Rscript evaluate_wnn_results.R <input_path> [n_workers]
# ==============================================================================

library(Seurat)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript evaluate_wnn_results.R <input_path> [n_workers]")
}

inPath <- args[1]

# ==============================================================================
# Step 1: Hard-coded marker gene lists
# ==============================================================================
# Marker genes from Golnaz marker list (markerGenes_Golnaz_wC11.txt)
# Organized by cell type
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
    "epsilon",
    "c11", "c11"
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
    "GHRL",
    "TPH2", "IKZF3"
  ),
  stringsAsFactors = FALSE
)

# Colors for cell types (matching original script)
colors <- data.frame(
  Type = c(
    "acinar", "alpha", "beta", "ductal", "stellates.mesenchymal",
    "endothelial", "immune", "pp.gamma", "delta", "epsilon",
    "c11", "unknown"
  ),
  col = c(
    "#e30800", "#f56505", "#dec400", "#A3AABE", "#006630",
    "#0223c7", "#5b02c7", "#00b0e6", "#c40080", "#02f00a",
    "#7d3301", "#000000"
  ),
  stringsAsFactors = FALSE
)

# Cycling alpha cell genes (from Elgamal et al. Diabetes 2023)
# Common cell cycle markers for alpha cells
# Note: If you have the original cyclingAlphaGenes file, replace this list
cyclingGenes <- data.frame(
  gene = c(
    "MKI67", "TOP2A", "PCNA", "CCNB1", "CCNB2",
    "CCNA2", "CCNE1", "CDK1", "BUB1", "CENPF",
    "CDK4", "CDK6", "CCND1", "CCND2"
  ),
  stringsAsFactors = FALSE
)

# Marker gene subsets (matching original script organization)
# Group 1: rows 1-14 and 31-32 (acinar, alpha, beta, ductal, c11) + ARX
markers_group1 <- c(markers$Marker[c(1:14, 31, 32)], "ARX")

# Group 2: rows 15-30 (stellates, endothelial, immune, pp.gamma, delta, epsilon)
markers_group2 <- markers$Marker[15:30]

# ==============================================================================
# Step 2: Load object with WNN integrations and per-modality annotations
# ==============================================================================
cat("Loading object with WNN integrations...\n")
seurat <- readRDS(paste0(inPath, '/seurat_gexAndAtac.rds'))

# Check if per-modality annotations exist
required_annotations <- c("gexHarmony", "gexRpcaRef", "gexRpcaAll", "atacHarmony", "atacRlsiRef")
missing_annotations <- setdiff(required_annotations, colnames(seurat@meta.data))
if (length(missing_annotations) > 0) {
  warning("Missing per-modality annotations: ", paste(missing_annotations, collapse = ", "), "\n")
  warning("Please run annotate_per_modality_clusters.R first.\n")
}

# ==============================================================================
# Step 3: Define cluster quality plotting function for WNN results
# ==============================================================================
# This function generates comprehensive plots for evaluating WNN cluster quality
cluster_quality_plots_wnn <- function(seurat, reduction, clusterName, fileName) {
  cat("Generating quality plots for", clusterName, "...\n")
  
  # Plot 1: DotPlot - Marker genes group 1 (acinar, alpha, beta, ductal, c11)
  d1 <- DotPlot(
    object = seurat,
    group.by = clusterName,
    features = markers_group1
  ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  # Plot 2: DotPlot - Marker genes group 2 (stellates, endothelial, immune, pp.gamma, delta, epsilon)
  d2 <- DotPlot(
    object = seurat,
    group.by = clusterName,
    features = markers_group2
  ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  # Plot 3: DotPlot - Cycling genes
  d3 <- DotPlot(
    object = seurat,
    group.by = clusterName,
    features = cyclingGenes$gene
  ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  # Plot 4: FeaturePlot - Marker genes group 1
  d4 <- FeaturePlot(
    seurat,
    features = markers_group1,
    combine = TRUE,
    ncol = 4,
    pt.size = 0.5,
    reduction = reduction,
    cols = c("grey", "red")
  )
  
  # Plot 5: FeaturePlot - Marker genes group 2
  d5 <- FeaturePlot(
    seurat,
    features = markers_group2,
    combine = TRUE,
    ncol = 4,
    pt.size = 0.5,
    reduction = reduction,
    cols = c("grey", "red")
  )
  
  # Plot 6: VlnPlot - QC metrics
  d6 <- VlnPlot(
    seurat,
    features = c("nFeature_RNA", "nCount_RNA", "log10GenesPerUMI"),
    group.by = clusterName,
    pt.size = 0,
    ncol = 2
  )
  
  # Plot 7a: DimPlot - Donor
  d7a <- DimPlot(
    seurat,
    group.by = "donor",
    reduction = reduction,
    label = FALSE
  ) + NoLegend()
  
  # Plot 7b: DimPlot - scSorter annotation (if available)
  if ("donorScSorter.sct" %in% colnames(seurat@meta.data)) {
    # Get colors for cell types present in the data
    types_present <- levels(as.factor(seurat$donorScSorter.sct))
    colors_new <- colors[colors$Type %in% types_present, ]
    color_vector <- setNames(colors_new$col, colors_new$Type)
    
    d7b <- DimPlot(
      seurat,
      group.by = "donorScSorter.sct",
      reduction = reduction,
      label = TRUE,
      cols = color_vector
    )
  } else {
    d7b <- NULL
  }
  
  # Plot 7c: DimPlot - WNN clusters
  d7c <- DimPlot(
    seurat,
    group.by = clusterName,
    reduction = reduction,
    label = TRUE
  ) + NoLegend()
  
  # Plot 7d: Bar plot - Disease composition (if available)
  if ("disease_wAge" %in% colnames(seurat@meta.data)) {
    d7d <- seurat@meta.data %>%
      as.data.frame() %>%
      group_by(!!sym(clusterName), disease_wAge) %>%
      tally(name = "count") %>%
      ggplot(aes_string(x = clusterName, y = "count", fill = "disease_wAge")) +
      geom_bar(stat = "identity", position = "fill") +
      geom_text(aes(label = count), position = position_fill(vjust = 0.5), size = 3) +
      labs(x = clusterName, y = "Proportion") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  } else {
    d7d <- NULL
  }
  
  # Plot 8: FeaturePlot - QC metrics
  if ("log10nFeature_RNA" %in% colnames(seurat@meta.data)) {
    d8 <- FeaturePlot(
      seurat,
      features = c("log10nFeature_RNA", "log10nCount_RNA", "log10GenesPerUMI"),
      combine = TRUE,
      ncol = 2,
      pt.size = 0.5,
      reduction = reduction
    ) & scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
  } else {
    d8 <- NULL
  }
  
  # Plot 9: Comparison with per-modality annotations
  comparison_annotations <- c(
    clusterName,
    "gexHarmony", "gexRpcaRef", "gexRpcaAll",
    "atacHarmony", "atacRlsiRef"
  )
  # Filter to only include annotations that exist
  comparison_annotations <- comparison_annotations[
    comparison_annotations %in% colnames(seurat@meta.data)
  ]
  
  if (length(comparison_annotations) > 1) {
    d9 <- DimPlot(
      seurat,
      group.by = comparison_annotations,
      reduction = reduction,
      label = TRUE,
      ncol = 2
    )
  } else {
    d9 <- NULL
  }
  
  # Combine plots into PDF
  pdf(fileName, width = 15, height = 15)
  
  # Layout first page: Overview plots
  plot_list <- list(d7a)
  if (!is.null(d7b)) plot_list <- c(plot_list, list(d7b))
  plot_list <- c(plot_list, list(d7c))
  if (!is.null(d7d)) plot_list <- c(plot_list, list(d7d))
  
  if (length(plot_list) > 0) {
    gridExtra::grid.arrange(grobs = plot_list, nrow = 2)
  }
  
  # Comparison with per-modality annotations
  if (!is.null(d9)) {
    plot(d9)
  }
  
  # Marker gene DotPlots
  plot(d1)
  plot(d2)
  
  # QC FeaturePlots
  if (!is.null(d8)) {
    plot(d8)
  }
  
  # Cycling genes
  plot(d3)
  
  # Marker gene FeaturePlots
  plot(d4)
  plot(d5)
  
  # QC VlnPlots
  plot(d6)
  
  dev.off()
  
  cat("Quality plots saved to:", fileName, "\n")
}

# ==============================================================================
# Step 4: Generate evaluation plots for each WNN result
# ==============================================================================
cat("Generating evaluation plots for WNN integrations...\n")

# WNN integrations to evaluate (wnn2-wnn6)
wnn_results <- paste0("clusters_wnn", 2:6)

for (wnn_cluster in wnn_results) {
  # Check if this WNN result exists
  if (!wnn_cluster %in% colnames(seurat@meta.data)) {
    warning("WNN cluster ", wnn_cluster, " not found. Skipping.\n")
    next
  }
  
  # Get the corresponding reduction name (e.g., "clusters_wnn5" -> "wnn5.umap")
  wnn_num <- gsub("clusters_wnn", "", wnn_cluster)
  reduction_name <- paste0("wnn", wnn_num, ".umap")
  
  # Check if reduction exists
  if (!reduction_name %in% names(seurat@reductions)) {
    warning("Reduction ", reduction_name, " not found. Skipping ", wnn_cluster, ".\n")
    next
  }
  
  cat("\nEvaluating", wnn_cluster, "...\n")
  fileName <- file.path(inPath, paste0("clusterQualityPlots_", wnn_cluster, ".pdf"))
  
  cluster_quality_plots_wnn(
    seurat = seurat,
    reduction = reduction_name,
    clusterName = wnn_cluster,
    fileName = fileName
  )
}

cat("\nWNN evaluation complete!\n")
cat("Evaluation plots saved to:", inPath, "\n")
cat("\nNext steps:\n")
cat("1. Review the cluster quality plots for each WNN result\n")
cat("2. Compare WNN results based on cluster separation and biological coherence\n")
cat("3. Select the best WNN result (wnn5 was selected in the published analysis)\n")
cat("4. Run annotate_wnn_filter_phase2.R to annotate and filter the selected WNN result\n")
