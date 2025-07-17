#!/usr/bin/env bash
set -euo pipefail

# Resolve script directory
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Paths
BAM_DIR="$BASE_DIR/../stage3/bam"
GTF="$BASE_DIR/../stage3/config/genes.gtf"
OUT_DIR="$BASE_DIR/counts"
LOG_DIR="$BASE_DIR/logs"
THREADS=4

# 1) Ensure featureCounts is available
if ! command -v featureCounts &> /dev/null; then
  echo "[ERROR] featureCounts not found. Install with: conda install -c bioconda subread" >&2
  exit 1
fi

# 2) Validate GTF
if [[ ! -s "$GTF" ]]; then
  echo "[ERROR] GTF not found at $GTF" >&2
  exit 1
fi
if ! grep -qP "\texon\t" "$GTF"; then
  echo "[ERROR] GTF at $GTF contains no exon entries" >&2
  exit 1
fi

# 3) Prepare output directories
mkdir -p "$OUT_DIR" "$LOG_DIR"

# 4) Collect BAM files
mapfile -t BAM_FILES < <(ls "$BAM_DIR"/*_sorted.bam 2>/dev/null || true)
if (( ${#BAM_FILES[@]} == 0 )); then
  echo "[ERROR] No BAM files found in $BAM_DIR" >&2
  exit 1
fi

echo "[INFO] Running featureCounts on ${#BAM_FILES[@]} BAM file(s)"
echo "[INFO] Annotation: $GTF"

# 5) Run featureCounts with exon counting by gene_id
featureCounts \
  -T "$THREADS" \
  -t exon \
  -g gene_id \
  -a "$GTF" \
  -o "$OUT_DIR/gene_counts.txt" \
  "${BAM_FILES[@]}" \
  2> "$LOG_DIR/featureCounts.log"

echo "[INFO] featureCounts complete"
echo "[INFO] Counts → $OUT_DIR/gene_counts.txt"
echo "[INFO] Log    → $LOG_DIR/featureCounts.log"

# 6) Convert to CSV, excluding summary lines
DATA="$OUT_DIR/gene_counts.txt"
CSV="$OUT_DIR/gene_counts.csv"

echo "[INFO] Converting gene_counts.txt → gene_counts.csv (genes only)"

awk '
  BEGIN {
    FS = "\t"; OFS = ",";
    header_printed = 0;
  }
  # skip comment lines and blank lines
  /^#/    { next }
  /^$/    { next }
  # skip summary section
  $1 == "Status"        { next }
  $1 == "Assigned"      { next }
  $1 ~ /^Unassigned_/   { next }

  {
    if (header_printed == 0) {
      # print header: Geneid + cleaned sample names
      printf "%s", $1
      for (i = 7; i <= NF; i++) {
        name = $i
        sub(".*/", "", name)                 # strip path
        sub(/_sorted\.bam$/, "", name)       # strip suffix
        printf "%s%s", OFS, name
      }
      print ""
      header_printed = 1
    } else {
      # print gene counts
      printf "%s", $1
      for (i = 7; i <= NF; i++) {
        printf "%s%s", OFS, $i
      }
      print ""
    }
  }
' "$DATA" > "$CSV"

echo "[INFO] CSV → $CSV"
