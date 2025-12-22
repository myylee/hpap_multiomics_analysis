# Multiome Data Analysis Pipeline

This directory contains scripts for processing and analyzing multiome (scRNA-seq + snATAC-seq) data from the Human Pancreas Analysis Program (HPAP).

## Overview

The analysis pipeline consists of two main phases:

1. **Phase 1**: Initial integration of GEX and ATAC data separately, followed by WNN integration
2. **Phase 2**: Re-integration after quality filtering, followed by final WNN integration and annotation

## Workflow

```
CellRanger Arc Output
    ↓
01_preprocessing/          # Individual sample processing
    ├── process_gex_sample.R
    ├── process_gex_sample_annotation.R
    ├── process_atac_sample.R
    └── process_atac_sample_qc.R
    ↓
02_integration_phase1/     # Phase 1: Modality-specific integration
    ├── integrate_gex_phase1.R          # Harmony + RPCA (all/ref)
    ├── integrate_gex_rpca_all.R
    ├── integrate_gex_rpca_ref.R
    └── integrate_atac_phase1.R         # Harmony + RLSI (ref)
    ↓
03_wnn_phase1/            # Phase 1: WNN integration
    └── integrate_wnn_phase1.R         # Multiple GEX+ATAC combinations
    ↓
[Evaluation and Filtering]
    ↓
04_integration_phase2/     # Phase 2: Re-integration after filtering
    ├── integrate_gex_phase2.R         # RPCA reference-based
    └── integrate_atac_phase2.R        # Harmony
    ↓
05_wnn_phase2/             # Phase 2: Final WNN integration
    └── integrate_wnn_phase2.R         # Final WNN with subclustering
    ↓
06_quality_control/        # Quality checks
    └── check_sample_quality.R
    ↓
07_reference_mapping/      # Map new samples to reference
    ├── map_scrnaseq.R
    └── map_snatacseq.R
```

## Directory Structure

```
multiome/
├── 01_preprocessing/          # Individual sample processing
│   ├── process_gex_sample.R
│   ├── process_gex_sample_annotation.R
│   ├── process_atac_sample.R
│   └── process_atac_sample_qc.R
│
├── 02_integration_phase1/     # Phase 1 GEX/ATAC integration
│   ├── integrate_gex_phase1.R
│   ├── integrate_gex_rpca_all.R
│   ├── integrate_gex_rpca_ref.R
│   └── integrate_atac_phase1.R
│
├── 03_wnn_phase1/             # Phase 1 WNN integration
│   └── integrate_wnn_phase1.R
│
├── 04_integration_phase2/      # Phase 2 re-integration
│   ├── integrate_gex_phase2.R
│   └── integrate_atac_phase2.R
│
├── 05_wnn_phase2/              # Phase 2 final WNN
│   └── integrate_wnn_phase2.R
│
├── 06_quality_control/         # QC scripts
│   └── check_sample_quality.R
│
└── 07_reference_mapping/        # Reference mapping
    ├── map_scrnaseq.R
    └── map_snatacseq.R
```

## Key Output Files

The final integrated object is:
- **File**: `multiome40DonorsRef-v2.rds`
- **Cell type annotation column**: `celltype2`
- **UMAP reduction**: `wnn5filteredRun2.umap`

For ATAC visualization with fragments:
- **File**: `multiome40DonorsRef-v2_wFrags.rds`

## Software Versions

- **CellRanger Arc**: 2.0.2
- **Seurat**: 5.1.0
- **Signac**: 1.13.0
- **R**: 4.x

## Dependencies

All scripts require the following R packages:
- Seurat (v5.1.0)
- Signac (v1.13.0)
- sctransform
- harmony
- tidyverse
- data.table
- future (for parallel processing)
- SoupX (for GEX preprocessing)
- scDblFinder (for doublet removal)
- scSorter (for cell type annotation)

## Usage Examples

### Preprocessing a single sample

```bash
# Process GEX data
Rscript 01_preprocessing/process_gex_sample.R \
  /path/to/cellranger/output \
  /path/to/output \
  30

Rscript 01_preprocessing/process_gex_sample_annotation.R \
  /path/to/output \
  30

# Process ATAC data
Rscript 01_preprocessing/process_atac_sample.R \
  /path/to/cellranger/output \
  /path/to/output \
  /path/to/macs2 \
  /path/to/annotation.rds

Rscript 01_preprocessing/process_atac_sample_qc.R \
  /path/to/output \
  30
```

### Phase 1 Integration

```bash
# GEX integration
Rscript 02_integration_phase1/integrate_gex_phase1.R \
  /path/to/input \
  /path/to/metadata.csv \
  30 \
  4

# RPCA integrations
Rscript 02_integration_phase1/integrate_gex_rpca_all.R \
  /path/to/input \
  30 \
  4

Rscript 02_integration_phase1/integrate_gex_rpca_ref.R \
  /path/to/input \
  30 \
  4 \
  "5,6,9,16,21,27,29,31,32"
```

## Notes

- All paths in scripts use command-line arguments or can be modified at the top of each script
- Parallel processing is configured using the `future` package
- Memory requirements vary by step; adjust `future.globals.maxSize` as needed
- Reference samples (indices 5,6,9,16,21,27,29,31,32) are high-quality control samples used for reference-based integration

## Citation

If you use this pipeline, please cite:
- Seurat: Hao et al., Cell 2021
- Signac: Stuart et al., Cell 2019
- Harmony: Korsunsky et al., Nat Methods 2019

