# HPAP scRNA-seq Analysis Pipeline

This directory contains scripts for processing and analyzing single-cell RNA sequencing (scRNA-seq) data from the Human Pancreas Analysis Program (HPAP).

## Overview

This workflow processes individual scRNA-seq samples from CellRanger output through quality control, normalization, clustering, and cell type annotation. Unlike the multiome workflow, this focuses on per-sample processing and does not include integration across donors.

## Workflow

```
CellRanger Output
    ↓
01_preprocessing/
    ├── run_cellranger.sh              # Example CellRanger script
    ├── process_sample.R                # SoupX correction, doublet removal, QC, normalization
    └── cluster_and_annotate.R          # Clustering, UMAP, scSorter annotation
```

## Directory Structure

```
scRNA-seq/
├── 01_preprocessing/
│   ├── run_cellranger.sh
│   ├── process_sample.R
│   └── cluster_and_annotate.R
└── README.md
```

## Software Versions

- **R**: 4.2.2
- **CellRanger**: 7.1.0
- **Seurat**: 4.9.9.9041
- **SoupX**: For ambient RNA correction
- **scDblFinder**: For doublet detection
- **scSorter**: For cell type annotation

## Dependencies

All scripts require the following R packages:
- Seurat (v4.9.9.9041)
- SoupX (for ambient RNA correction)
- scDblFinder (for doublet removal)
- scSorter (for cell type annotation)
- tidyverse
- SingleCellExperiment

## Workflow Steps

### Step 1: CellRanger Processing

Run CellRanger count on raw FASTQ files to generate count matrices.

**Example script**: `run_cellranger.sh`

This script provides an example of how to run CellRanger. Adjust paths and parameters as needed for your system.

```bash
bash 01_preprocessing/run_cellranger.sh \
  HPAP079-islets-101146-101148 \
  /path/to/fastq/files \
  /path/to/output \
  4 \
  128
```

**Key parameters**:
- Reference: GRCh38 (with `--include-introns` by default)
- Adjust `localcores` and `localmem` based on your system resources

### Step 2: Process Sample

This script processes the CellRanger output:
1. Loads 10X data and corrects for ambient RNA using SoupX
2. Removes doublets using scDblFinder
3. Filters cells based on quality metrics:
   - `nFeature_RNA`: 200 < nFeature < 10,000
   - `percent.mt`: < 25%
   - `nCount_RNA`: 500 < nCount < 100,000
4. Normalizes data using LogNormalize
5. Runs PCA

```bash
Rscript 01_preprocessing/process_sample.R \
  /path/to/cellranger/output \
  /path/to/output \
  40
```

**Outputs**:
- `seurat.tmp.rds`: Intermediate Seurat object
- QC plots: `autoSoupX.pdf`, `vlnPlotBefore.pdf`, `vlnPlotAfter.pdf`
- Elbow plot: `rna_elbow.pdf`

### Step 3: Cluster and Annotate

This script performs clustering and cell type annotation:
1. Runs clustering and UMAP using RNA (LogNormalize) assay
2. Performs cell type annotation using scSorter with marker genes
   (Note: scSorter annotation is a first pass to get a general idea of cell types)
3. Visualizes annotations

```bash
Rscript 01_preprocessing/cluster_and_annotate.R \
  /path/to/output \
  15
```

**Marker genes**: The marker genes are hardcoded in the script (based on Golnaz marker list) and include markers for:
- Acinar, Alpha, Beta, Ductal cells
- Stellate/mesenchymal, Endothelial, Immune cells
- PP (gamma), Delta, Epsilon cells

**Outputs**:
- `seurat.rds`: Final Seurat object with annotations
- UMAP plot: `rna_umap.pdf` (colored by scSorter annotations)

## Key Output Files

The final processed object (`seurat.rds`) contains:
- **Assays**: `RNA` (LogNormalize)
- **Reductions**: `rna.pca`, `rna.umap`
- **Clusters**: `rna.clusters`
- **Annotations**: `donorScSorter.rna`
- **QC metrics**: `nFeature_RNA`, `nCount_RNA`, `percent.mt`, `scDblFinder.class`

## Important Notes

1. **Seurat Assay Version**: Scripts use `options(Seurat.object.assay.version = "v3")` for compatibility with Seurat v4.

2. **Normalization**: Data is normalized using LogNormalize, followed by variance stabilization, scaling, and PCA.

3. **scSorter Annotation**: Cell type annotation is performed using scSorter with hardcoded marker genes on the RNA assay. This is a first pass annotation to get a general idea of cell types.

4. **No Integration**: This workflow focuses on per-sample processing. Integration across donors is not included in this pipeline (see multiome workflow for integration approaches).

5. **Quality Filters**: The filtering thresholds were selected based on typical pancreas islet cell characteristics. Adjust as needed for your data:
   - `nFeature_RNA`: 200-10,000
   - `percent.mt`: < 25%
   - `nCount_RNA`: 500-100,000

## Citation

If you use this pipeline, please cite:
- **Seurat**: Hao et al., Cell 2021
- **SoupX**: Young & Behjati, GigaScience 2020
- **scDblFinder**: Germain et al., Genome Biology 2021
- **scSorter**: Guo et al., Genome Biology 2021

## Acknowledgments

These scripts were based on the scRNA-seq analysis workflow available at:
[https://hpap.pmacs.upenn.edu/explore/workflow/scRNAseq-analysis](https://hpap.pmacs.upenn.edu/explore/workflow/scRNAseq-analysis)

The original scripts were constructed by **Dr. Elisabetta Manduchi**.

