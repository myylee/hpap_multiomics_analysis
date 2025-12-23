#!/bin/bash
# ==============================================================================
# Example CellRanger Script for scRNA-seq Data
#
# This is an example script for running CellRanger count on scRNA-seq data.
# Adjust paths and parameters as needed for your specific data.
#
# CellRanger version: 7.1.0
# Reference: GRCh38 (refdata-gex-GRCh38-2020-A)
#
# Usage: bash run_cellranger.sh <sample_id> <fastq_paths> <output_dir> [localcores] [localmem]
#   sample_id: Sample identifier (e.g., HPAP124-islets-111817)
#   fastq_paths: Comma-separated paths to directories containing FASTQ files (e.g., /path/to/FGC2468_2,/path/to/FGC2468_1)
#   output_dir: Output directory for CellRanger results
#   localcores: Number of cores to use (default: 16)
#   localmem: Memory in GB to use (default: 64)
# ==============================================================================

set -e  # Exit on error

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: bash run_cellranger.sh <sample_id> <fastq_paths> <output_dir> [localcores] [localmem]"
    echo "  fastq_paths: Comma-separated paths to FASTQ directories"
    exit 1
fi

SAMPLE_ID=$1
FASTQ_PATHS=$2
OUTPUT_DIR=$3
LOCALCORES=${4:-16}
LOCALMEM=${5:-64}

# Set paths - ADJUST THESE TO YOUR SYSTEM
CELLRANGER_PATH="/path/to/cellranger-7.1.0/cellranger"
REFERENCE_PATH="/path/to/refdata-gex-GRCh38-2020-A"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Run CellRanger count
echo "Running CellRanger count for ${SAMPLE_ID}..."
echo "FASTQ paths: ${FASTQ_PATHS}"
echo "Output directory: ${OUTPUT_DIR}"

"${CELLRANGER_PATH}" count \
    --transcriptome="${REFERENCE_PATH}" \
    --id="${SAMPLE_ID}" \
    --fastqs="${FASTQ_PATHS}" \
    --chemistry=auto \
    --include-introns=true \
    --localcores="${LOCALCORES}" \
    --localmem="${LOCALMEM}"

echo "CellRanger count complete for ${SAMPLE_ID}"
echo "Results saved to: ${OUTPUT_DIR}/${SAMPLE_ID}/outs"

