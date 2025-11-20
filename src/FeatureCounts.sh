#!/bin/bash

# Script to run FeatureCounts on aligned BAM files
# Usage: ./FeatureCounts.sh

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ALIGN_DIR="${ROOT_DIR}/results/alignment"
OUTPUT_DIR="${ROOT_DIR}/results/counts"
GTF_FILE="${ROOT_DIR}/data/Saccharomyces_cerevisiae.R64-1-1.gtf"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run FeatureCounts
featureCounts -T 4 -a "$GTF_FILE" -o "${OUTPUT_DIR}/gene_counts.txt" \
    ${ALIGN_DIR}/*_sorted.bam

echo "FeatureCounts completed. Gene counts saved to ${OUTPUT_DIR}/gene_counts.txt"
