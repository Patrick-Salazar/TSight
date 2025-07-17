import json
import yaml
import sys
from pathlib import Path

# CONFIGURATION
MULTIQC_JSON_PATH   = Path("../qc/fastqc_reports/multiqc_data/multiqc_data.json")
QC_THRESHOLDS_PATH  = Path("config/qc_thresholds.yaml")
TRIMMED_DIR         = Path("../data/trimmed")
DEFER_DIR           = Path("../batches/defer")
REJECTED_SAMPLE_DIR = TRIMMED_DIR / "rejected"

# SANITY CHECKS
for p in (MULTIQC_JSON_PATH, QC_THRESHOLDS_PATH):
    if not p.is_file():
        sys.exit(f"ERROR: file not found: {p}")
if not TRIMMED_DIR.is_dir():
    sys.exit(f"ERROR: trimmed-data folder not found: {TRIMMED_DIR}")

# CREATE OUTPUT FOLDERS
DEFER_DIR.mkdir(parents=True, exist_ok=True)
REJECTED_SAMPLE_DIR.mkdir(parents=True, exist_ok=True)

# LOAD INPUTS
def load_multiqc_data(path):
    with open(path) as f:
        return json.load(f)

def load_qc_thresholds(path):
    with open(path) as f:
        return yaml.safe_load(f)

multiqc_data  = load_multiqc_data(MULTIQC_JSON_PATH)
qc_thresholds = load_qc_thresholds(QC_THRESHOLDS_PATH)
fastqc_stats  = multiqc_data["report_general_stats_data"]["fastqc"]

# QC CHECK FUNCTION
def check_sample_qc(metrics, thresholds):
    fails = []
    if metrics.get("percent_fails", 0) > thresholds["max_percent_fails"]:
        fails.append("High FastQC Fails")
    if metrics.get("percent_duplicates", 0) > thresholds["max_percent_duplicates"]:
        fails.append("High Duplication")
    if metrics.get("total_sequences", 0) < thresholds["min_read_count"]:
        fails.append("Low Read Count")
    return fails

# BUILD BATCH LISTS
batches = {}
for sample, metrics in fastqc_stats.items():
    for reason in check_sample_qc(metrics, qc_thresholds):
        batches.setdefault(reason, []).append(sample)

if not batches:
    print("All samples passed post-trim QC.")
    sys.exit(0)

# PROCESS EACH BATCH
for reason, samples in batches.items():
    tag = reason.replace(" ", "_")
    batch_file = DEFER_DIR / f"postqc_{tag}.txt"
    batch_file.write_text("\n".join(samples))

    print(f"\nBatch: {reason}")
    print("Samples:", ", ".join(samples))
    keep = input("Proceed with these? (yes/no): ").strip().lower()

    if keep != "yes":
        print("Moving rejected samples to", REJECTED_SAMPLE_DIR)
        for sample in samples:
            matches = list(TRIMMED_DIR.rglob(f"*{sample}*"))
            if not matches:
                print(f"  WARNING: no files found for sample '{sample}'")
            for f in matches:
                if f.is_file():
                    dest = REJECTED_SAMPLE_DIR / f.name
                    print(f"  {f.relative_to(TRIMMED_DIR)} → rejected/{f.name}")
                    f.rename(dest)
    else:
        print("Retaining this batch for alignment.")

print("Post-trim QC complete.")
