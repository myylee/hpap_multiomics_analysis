# HPAP Multiomics Analysis

Analysis code for the analysis of scRNA-seq, snATAC-seq and multiome data from the Human Pancreas Analysis Program (HPAP).

## Repository Structure

This repository contains organized analysis scripts for:

- **Multiome data analysis** (`multiome/`) - Comprehensive pipeline for processing and integrating multiome (scRNA-seq + snATAC-seq) data
- Additional analysis modules (to be added)

## Quick Start

See the [multiome README](multiome/README.md) for detailed documentation on the multiome analysis pipeline.

## Key Outputs

The final integrated multiome object:
- **File**: `multiome40DonorsRef-v2.rds`
- **Cell type annotation**: `celltype2` column
- **UMAP reduction**: `wnn5filteredRun2.umap`

## Software Versions

- Seurat 5.1.0
- Signac 1.13.0
- CellRanger Arc 2.0.2

## Citation

If you use this code, please cite the relevant publications for:
- Seurat (Hao et al., Cell 2021)
- Signac (Stuart et al., Cell 2019)
- Harmony (Korsunsky et al., Nat Methods 2019)
