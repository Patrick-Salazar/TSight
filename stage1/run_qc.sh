#!/bin/bash
# run_qc.sh: Runs FastQC and MultiQC safely via Bash

set -e  # Exit on error

QC_DIR="qc/fastqc_reports"
SAMPLES_DIR="data/samples"

mkdir -p "$QC_DIR"  # Ensure QC directory exists

# Enable nullglob to handle unmatched globs safely
shopt -s nullglob

# Find FASTQ/FQ files
FILES=("$SAMPLES_DIR"/*.fastq "$SAMPLES_DIR"/*.fq)

if [ ${#FILES[@]} -eq 0 ]; then
  echo "No FASTQ/FQ files found in $SAMPLES_DIR. Aborting."
  exit 1
fi

# Clean old QC reports before running new ones
rm -rf "$QC_DIR"/*_fastqc* "$QC_DIR"/multiqc_report*.html "$QC_DIR"/multiqc_data*

# Run FastQC (parallel if desired)
for f in "${FILES[@]}"; do
  fastqc "$f" -o "$QC_DIR" &
done
wait

# Run MultiQC
multiqc "$QC_DIR" -o "$QC_DIR"
