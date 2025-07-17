#!/usr/bin/env bash
set -e

echo "=== TSight Stage 2 Pipeline Starting ==="

# run from this script’s directory
cd "$(dirname "$0")"

# 1. Activate environment
echo "Activating rnaseq-env-qc environment..."
source ~/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq-env-qc

# 2. Generate trimming plan
echo "→ Generating trimming plan..."
python trimming_plan_generator.py

# 3. Run trimming executor
echo "→ Running trimming executor..."
python trimmer.py

# 4. Run FastQC on trimmed FASTQ files
echo "→ Running FastQC on trimmed FASTQ files..."
mkdir -p ../qc/fastqc_reports/
fastqc ../data/trimmed/*.fastq -o ../qc/fastqc_reports/

# 5. Run MultiQC
echo "→ Running MultiQC..."
multiqc ../qc/fastqc_reports/ -o ../qc/fastqc_reports/

# 6. Run post-trim QC checkpoint
echo "→ Running post-trimming QC checkpoint..."
python post_trim_qc_check.py

echo "=== TSight Stage 2 Pipeline Complete ==="
