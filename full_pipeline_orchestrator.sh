#!/usr/bin/env bash
set -euo pipefail

echo "=== TSight Full Pipeline Starting ==="

# Stage 1: QC & Batching
echo ">>> Running Stage 1 Orchestrator..."
python stage1/run_pipeline.py

# Stage 2: Trimming & QC
echo ">>> Running Stage 2 Orchestrator..."
bash stage2/stage2_orchestrator.sh

# Stage 3: Alignment
echo ">>> Running Stage 3 Alignment..."
bash stage3/run_alignment.sh

# Stage 4: Quantification
 echo ">>> Running Stage 4 Quantification..."
 bash stage4/run_quantification.sh

# Stage 4: DEG
echo ">>> Running Stage 5: DESeq2"
bash stage5/run_deseq2.sh


echo "=== TSight Full Pipeline Complete ==="
