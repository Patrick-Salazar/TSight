#!/usr/bin/env bash
set -euo pipefail

# Resolve this script’s directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# 1) Check for Rscript
if ! command -v Rscript &> /dev/null; then
  echo "[ERROR] Rscript not found. Install R and DESeq2." >&2
  exit 1
fi

# 2) Move to project root so that R’s working dir is the root
cd "$PROJECT_ROOT"

echo "[INFO] Running DESeq2 from $PROJECT_ROOT"
Rscript "$SCRIPT_DIR/run_deseq2.R"

echo "[INFO] DESeq2 complete. Results in stage5/results/"
