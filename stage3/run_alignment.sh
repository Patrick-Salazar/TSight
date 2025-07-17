#!/usr/bin/env bash
set -euo pipefail

# Resolve script’s own directory
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Correct paths relative to stage3/
INPUT_DIR="$BASE_DIR/../data/trimmed"
OUTPUT_DIR="$BASE_DIR/bam"
LOG_DIR="$BASE_DIR/logs"
INDEX_BASE="$BASE_DIR/index/genome"
THREADS=4

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

declare -A PROCESSED
function process_bam {
  local s=$1
  samtools view -@ "$THREADS" -bS "$OUTPUT_DIR/$s.sam" \
    | samtools sort -@ "$THREADS" -o "$OUTPUT_DIR/${s}_sorted.bam"
  samtools index  "$OUTPUT_DIR/${s}_sorted.bam"
  rm "$OUTPUT_DIR/$s.sam"
  echo "[OK] $s → ${s}_sorted.bam"
}

# Paired-end pass
for fq1 in "$INPUT_DIR"/*_1_trimmed.fastq "$INPUT_DIR"/*_1_trimmed.fastq.gz; do
  [[ -e $fq1 ]] || continue
  s=$(basename "$fq1" _1_trimmed.fastq* )
  [[ -n ${PROCESSED[$s]:-} ]] && continue

  # detect mate
  if [[ -e "$INPUT_DIR/${s}_2_trimmed.fastq" ]]; then
    fq2="$INPUT_DIR/${s}_2_trimmed.fastq"
  elif [[ -e "$INPUT_DIR/${s}_2_trimmed.fastq.gz" ]]; then
    fq2="$INPUT_DIR/${s}_2_trimmed.fastq.gz"
  else
    fq2=""
  fi

  echo "[INFO] $s → paired? $([[ -n $fq2 ]] && echo yes || echo no)"
  if [[ -n $fq2 ]]; then
    hisat2 -x "$INDEX_BASE" -1 "$fq1" -2 "$fq2" --dta -p "$THREADS" \
      -S "$OUTPUT_DIR/$s.sam" 2> "$LOG_DIR/${s}.log"
  else
    hisat2 -x "$INDEX_BASE" -U "$fq1" --dta -p "$THREADS" \
      -S "$OUTPUT_DIR/$s.sam" 2> "$LOG_DIR/${s}.log"
  fi

  process_bam "$s"
  PROCESSED[$s]=1
done

# Single-end leftover
for fq in "$INPUT_DIR"/*_trimmed.fastq "$INPUT_DIR"/*_trimmed.fastq.gz; do
  [[ -e $fq ]] || continue
  b=$(basename "$fq")
  [[ $b =~ _[12]_trimmed\.fastq ]] && continue
  s=${b%%_trimmed.fastq*}
  [[ -n ${PROCESSED[$s]:-} ]] && continue

  echo "[INFO] $s → single"
  hisat2 -x "$INDEX_BASE" -U "$fq" --dta -p "$THREADS" \
    -S "$OUTPUT_DIR/$s.sam" 2> "$LOG_DIR/${s}.log"

  process_bam "$s"
  PROCESSED[$s]=1
done
