# Minimal Required CellRanger Arc Input Files

This document describes the minimal set of files needed from CellRanger Arc output to run the preprocessing scripts in this directory.

## CellRanger Arc Output Structure

The scripts expect the CellRanger Arc output directory structure:

```
<sample_id>/
└── outs/
    ├── raw_feature_bc_matrix.h5          # Required for SoupX (GEX)
    ├── filtered_feature_bc_matrix.h5     # Required for both GEX and ATAC
    ├── atac_fragments.tsv.gz             # Required for ATAC processing
    └── [other files - not required]
```

## File Requirements by Script

### For GEX Processing (`process_gex_sample.R`)

**Required files:**
- `outs/raw_feature_bc_matrix.h5` - Used by SoupX to estimate ambient RNA contamination
  - Alternative: `outs/raw_feature_bc_matrix/` directory containing:
    - `matrix.mtx.gz`
    - `barcodes.tsv.gz`
    - `features.tsv.gz`
- `outs/filtered_feature_bc_matrix.h5` - Filtered count matrix
  - Alternative: `outs/filtered_feature_bc_matrix/` directory containing:
    - `matrix.mtx.gz`
    - `barcodes.tsv.gz`
    - `features.tsv.gz`

**Note:** SoupX's `load10X()` function automatically detects and uses both raw and filtered matrices. The raw matrix is essential for ambient RNA correction.

### For ATAC Processing (`process_atac_sample.R`)

**Required files:**
- `outs/filtered_feature_bc_matrix.h5` - Contains both RNA and ATAC peak counts
  - The script extracts `counts$Peaks` from this file to get ATAC peak counts
- `outs/atac_fragments.tsv.gz` - Fragment file for ATAC-seq data
  - Required for peak calling (MACS2)
  - Required for fragment counting and QC metrics
  - Required for downstream analysis

### Files NOT Required

The following files from CellRanger Arc output are not used by these scripts:
- `outs/raw_feature_bc_matrix/` directory (if `.h5` file exists)
- `outs/filtered_feature_bc_matrix/` directory (if `.h5` file exists)
- `outs/web_summary.html`
- `outs/molecule_info.h5`
- `outs/summary.csv`
- `outs/peak_annotation.tsv`
- `outs/singlecell.csv`
- Any files in `outs/analysis/` directory

## Minimal Reconstruction

To create a minimal input directory that will work with all preprocessing scripts:

```bash
# Create directory structure
mkdir -p <sample_id>/outs

# Copy required files (assuming CellRanger Arc has been run)
cp <cellranger_output>/outs/raw_feature_bc_matrix.h5 <sample_id>/outs/
cp <cellranger_output>/outs/filtered_feature_bc_matrix.h5 <sample_id>/outs/
cp <cellranger_output>/outs/atac_fragments.tsv.gz <sample_id>/outs/

# Optionally create index file for fragments (if not already present)
# This may be needed for some Signac functions
# tabix -p bed <sample_id>/outs/atac_fragments.tsv.gz
```

## Important Notes

1. **SoupX requirement**: The `raw_feature_bc_matrix.h5` file is essential for SoupX's ambient RNA correction. Without it, SoupX cannot estimate the contamination profile.

2. **Fragment file**: The `atac_fragments.tsv.gz` file must be the compressed TSV format. The scripts create fragment objects directly from this file path.

3. **H5 format preferred**: While alternative directory structures (matrix.mtx.gz, etc.) may work, the `.h5` format is more efficient and is what the scripts are designed to use.

4. **Sample naming**: The scripts extract donor/sample information from the directory name, expecting format: `<donor>-<source>-<sample>` (e.g., `HPAP079-islets-101146-101148`)

## Validation

To verify your minimal input structure is correct:

```bash
# Check that all required files exist
test -f <sample_id>/outs/raw_feature_bc_matrix.h5 && echo "✓ raw matrix exists" || echo "✗ missing raw matrix"
test -f <sample_id>/outs/filtered_feature_bc_matrix.h5 && echo "✓ filtered matrix exists" || echo "✗ missing filtered matrix"
test -f <sample_id>/outs/atac_fragments.tsv.gz && echo "✓ fragments file exists" || echo "✗ missing fragments file"

# Test that h5 files are readable (requires h5dump or similar)
h5dump -n <sample_id>/outs/filtered_feature_bc_matrix.h5 > /dev/null 2>&1 && echo "✓ h5 file is readable" || echo "✗ h5 file corrupted"
```

