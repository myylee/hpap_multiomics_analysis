# Multiome Data Analysis Pipeline

This directory contains scripts for processing and analyzing multiome (scRNA-seq + snATAC-seq) data from the Human Pancreas Analysis Program (HPAP).

## Overview

The analysis pipeline processes multiome (snRNA-seq + snATAC-seq) data through two main phases:

1. **Phase 1**: Initial integration of GEX and ATAC data separately, followed by WNN integration and evaluation
2. **Phase 2**: Re-integration after quality filtering, followed by final WNN integration and annotation

This two-phase approach allows for quality control and filtering between integration rounds, improving the final integrated object quality.

## Detailed Workflow

### Step 1: Preprocessing Individual Samples

**GEX (Gene Expression) Processing:**
- `process_gex_sample.R`: Loads 10X data, corrects for ambient RNA (SoupX), removes doublets, applies QC filters, runs SCTransform
- `process_gex_sample_annotation.R`: Runs UMAP, clustering, and basic visualization

**ATAC (Chromatin Accessibility) Processing:**
- `process_atac_sample.R`: Loads ATAC data, calls peaks using MACS2 (creates `peak_macs2` assay with sample-specific peaks), removes doublets
- `process_atac_sample_qc.R`: Calculates QC metrics (TSS enrichment, nucleosome signal), filters cells, runs LSI on sample-specific peaks
- `combine_atac_peaks.R`: Combines peaks from all samples into a unified peak set (run once after all samples processed)
- `process_atac_unified_peaks.R`: Quantifies each sample's fragments against the unified peak set, creates `macs2_combined` assay (essential for integration)

**Quality Control:**
- `check_sample_quality.R`: Aggregates QC metrics from all samples, creates summary tables and combined metadata

### Step 2: Phase 1 - Modality-Specific Integration

**GEX Integration:**
- `integrate_gex_phase1.R`: Creates combined Seurat object with layers, runs SCTransform, performs Harmony integration
- `integrate_gex_rpca_all.R`: RPCA integration using all samples
- `integrate_gex_rpca_ref.R`: RPCA integration with reference samples (higher quality control samples)

**ATAC Integration:**
- `integrate_atac_harmony.R`: Harmony integration of ATAC data using unified peaks
- `integrate_atac_rlsi.R`: Reciprocal LSI (RLSI) integration with reference samples

### Step 3: Phase 1 - WNN Integration

- `combine_gex_atac.R`: Combines GEX and ATAC objects, adds integration metadata and reductions
- `integrate_wnn_phase1.R`: Performs WNN integration using different GEX+ATAC reduction combinations:
  - WNN2: Harmony GEX + RLSI ATAC
  - WNN3: RPCA All GEX + Harmony ATAC
  - WNN4: RPCA All GEX + RLSI ATAC
  - WNN5: RPCA Ref GEX + Harmony ATAC
  - WNN6: RPCA Ref GEX + RLSI ATAC

### Step 4: Evaluation and Filtering

Evaluate WNN results (typically WNN5 or WNN6) and filter out low-quality cells (e.g., mixed populations, low-quality clusters). This filtering step prepares the data for Phase 2 re-integration.

### Step 5: Phase 2 - Re-integration

After filtering, re-run integration on the cleaned dataset:
- `integrate_gex_phase2.R`: Re-runs RPCA integration (reference-based) on filtered cells
- `integrate_atac_phase2.R`: Re-runs Harmony integration on filtered cells

### Step 6: Phase 2 - Final WNN Integration

- `integrate_wnn_phase2.R`: Final WNN integration using Phase 2 reductions, creates final integrated object (`multiome40DonorsRef-v2.rds`)

### Step 7: Reference Mapping (Optional)

For mapping new samples to the integrated reference:
- `map_scrnaseq.R`: Maps scRNA-seq query data to multiome reference
- `map_snatacseq.R`: Maps snATAC-seq query data to multiome reference

## Workflow Diagram

```
CellRanger Arc Output
    ↓
01_preprocessing/          # Individual sample processing
    ├── process_gex_sample.R
    ├── process_gex_sample_annotation.R
    ├── process_atac_sample.R
    ├── process_atac_sample_qc.R
    ├── combine_atac_peaks.R          # Combine peaks across all samples
    └── process_atac_unified_peaks.R  # Quantify each sample against unified peaks
    ↓
06_quality_control/        # Quality checks (run anytime)
    └── check_sample_quality.R         # Aggregate QC metrics
    ↓
02_integration_phase1/     # Phase 1: Modality-specific integration
    ├── integrate_gex_phase1.R          # Main GEX integration (Harmony)
    ├── integrate_gex_rpca_all.R        # RPCA integration (all samples)
    ├── integrate_gex_rpca_ref.R        # RPCA integration (reference-based)
    ├── integrate_atac_harmony.R        # ATAC Harmony integration
    └── integrate_atac_rlsi.R           # ATAC RLSI integration (reference-based)
    ↓
03_wnn_phase1/            # Phase 1: WNN integration
    ├── combine_gex_atac.R              # Combine GEX and ATAC data
    └── integrate_wnn_phase1.R          # Multiple GEX+ATAC combinations
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
07_reference_mapping/      # Map new samples to reference (optional)
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
│   ├── process_atac_sample_qc.R
│   ├── combine_atac_peaks.R
│   └── process_atac_unified_peaks.R
│
├── 02_integration_phase1/     # Phase 1 GEX/ATAC integration
│   ├── integrate_gex_phase1.R
│   ├── integrate_gex_rpca_all.R
│   ├── integrate_gex_rpca_ref.R
│   ├── integrate_atac_harmony.R
│   └── integrate_atac_rlsi.R
│
├── 03_wnn_phase1/             # Phase 1 WNN integration
│   ├── combine_gex_atac.R
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

# After processing all samples, combine peaks and create unified peak set
Rscript 01_preprocessing/combine_atac_peaks.R \
  /path/to/input

# Then quantify each sample against unified peaks
Rscript 01_preprocessing/process_atac_unified_peaks.R \
  /path/to/sample/output \
  30 \
  /path/to/peak_macs2_combined_across_samples.csv
```

### Quality Control Summary

```bash
# After processing all samples, generate QC summary
Rscript 06_quality_control/check_sample_quality.R \
  /path/to/input
```

### Phase 1 Integration

```bash
# GEX integration - Main script (creates unintegrated + Harmony)
Rscript 02_integration_phase1/integrate_gex_phase1.R \
  /path/to/input \
  /path/to/metadata.csv \
  30 \
  4

# GEX RPCA integrations (run after integrate_gex_phase1.R)
Rscript 02_integration_phase1/integrate_gex_rpca_all.R \
  /path/to/input \
  30 \
  4

Rscript 02_integration_phase1/integrate_gex_rpca_ref.R \
  /path/to/input \
  30 \
  4 \
  "5,6,9,16,21,27,29,31,32"

# ATAC integration
Rscript 02_integration_phase1/integrate_atac_harmony.R \
  /path/to/input \
  /path/to/metadata.csv \
  30 \
  4

Rscript 02_integration_phase1/integrate_atac_rlsi.R \
  /path/to/input \
  /path/to/metadata.csv \
  30 \
  4 \
  "5,6,9,16,21,27,29,31,32"
```

### Phase 1 WNN Integration

```bash
# Combine GEX and ATAC data
Rscript 03_wnn_phase1/combine_gex_atac.R \
  /path/to/input \
  /path/to/donor_metadata.csv

# Perform WNN integration
Rscript 03_wnn_phase1/integrate_wnn_phase1.R \
  /path/to/input \
  4
```

### Phase 2 Re-integration (After Filtering)

```bash
# Filter cells based on Phase 1 WNN results (manual step - filter object in R)
# Then re-integrate:

# GEX re-integration
Rscript 04_integration_phase2/integrate_gex_phase2.R \
  /path/to/input \
  /path/to/filtered_object.rds \
  30 \
  4

# ATAC re-integration
Rscript 04_integration_phase2/integrate_atac_phase2.R \
  /path/to/input \
  /path/to/filtered_object.rds \
  4
```

### Phase 2 Final WNN Integration

```bash
# Final WNN integration
Rscript 05_wnn_phase2/integrate_wnn_phase2.R \
  /path/to/input \
  /path/to/filtered_object.rds \
  4
```

### Reference Mapping

```bash
# Map scRNA-seq to reference
Rscript 07_reference_mapping/map_scrnaseq.R \
  /path/to/reference.rds \
  /path/to/query.rds \
  /path/to/output

# Map snATAC-seq to reference
Rscript 07_reference_mapping/map_snatacseq.R \
  /path/to/reference.rds \
  /path/to/query.rds \
  /path/to/output
```

## Important Notes

### Execution Order
1. **Preprocessing**: Process all individual samples first (GEX and ATAC)
2. **Peak Combination**: Run `combine_atac_peaks.R` once after all ATAC samples are processed
3. **Unified Peaks**: Run `process_atac_unified_peaks.R` for each sample using the combined peaks
4. **Quality Control**: Run `check_sample_quality.R` anytime to generate QC summaries
5. **Phase 1 Integration**: Run GEX integrations, then ATAC integrations
6. **Phase 1 WNN**: Combine GEX+ATAC, then run WNN integration
7. **Evaluation**: Manually evaluate WNN results and filter low-quality cells
8. **Phase 2**: Re-integrate filtered data, then final WNN

### Key Concepts

**Unified Peak Set**: For ATAC integration to work, all samples must use the same set of peaks. The workflow:
- Calls peaks per sample (`peak_macs2` assay)
- Combines all peaks into a unified set
- Re-quantifies each sample against unified peaks (`macs2_combined` assay)

**Two-Phase Approach**: 
- Phase 1: Initial integration with all cells to identify low-quality populations
- Filtering: Remove low-quality cells based on Phase 1 results
- Phase 2: Re-integration on cleaned dataset for improved results

**WNN Integration**: Weighted Nearest Neighbor integration optimally combines GEX and ATAC modalities. Multiple combinations are tested (WNN2-6) to find the best approach.

### Technical Details
- All paths in scripts use command-line arguments
- Parallel processing is configured using the `future` package
- Memory requirements vary by step; adjust `future.globals.maxSize` as needed
- Reference samples (indices 5,6,9,16,21,27,29,31,32) are high-quality control samples used for reference-based integration
- The final object uses `celltype2` for cell type annotations and `wnn5filteredRun2.umap` for visualization

## Citation

If you use this pipeline, please cite:
- **Seurat**: Hao et al., Cell 2021
- **Signac**: Stuart et al., Cell 2019
- **Harmony**: Korsunsky et al., Nat Methods 2019
- **SoupX**: Young & Behjati, GigaScience 2020 (for GEX ambient RNA correction)
- **scDblFinder**: Germain et al., Genome Biology 2021 (for doublet detection; Amulet method for ATAC-seq doublet detection)

