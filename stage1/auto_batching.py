# TSight Pipeline: Final Production-Ready Batching Script
# Features: Sample Tracking, Auto-Batching, Safe Reprocessing, Logging

import json
import csv
import pandas as pd
import os
import glob
import re
from pathlib import Path
from datetime import datetime

# Force working directory to pipeline root
os.chdir(os.path.dirname(os.path.dirname(__file__)))


# === SETTINGS ===
PIPELINE_VERSION = "v1.0"

# === Define Project Root ===
project_root = Path(".").resolve()

# === Path Adjustments ===
json_path = project_root / "qc/fastqc_reports/multiqc_data/multiqc_data.json"
sample_tracker_path = project_root / "tracking/sample_tracker.csv"
pipeline_log = project_root / "tracking/pipeline_run.log"

batches_dir = project_root / "batches"
good_dir = batches_dir / "good"
defer_dir = batches_dir / "defer"

# === Safety Check for JSON File ===
if not os.path.exists(json_path):
    print(f"Error: JSON file '{json_path}' not found.")
    print("Tip: Run FastQC + MultiQC first to generate this file.")
    print("This script requires completed QC reports to proceed.")
    exit(1)

# === Auto-Create Folders ===
good_dir.mkdir(parents=True, exist_ok=True)
defer_dir.mkdir(parents=True, exist_ok=True)

# === Load or Initialize Sample Tracker ===
if os.path.exists(sample_tracker_path):
    sample_tracker_df = pd.read_csv(sample_tracker_path)
else:
    sample_tracker_df = pd.DataFrame(columns=["Sample", "QC_Status", "Reason", "Batch_File", "Pipeline_Version"])

# === Load MultiQC JSON ===
with open(json_path) as f:
    data = json.load(f)
qc_stats = data["report_general_stats_data"]["fastqc"]

# === Identify New Samples ===
new_samples = [s for s in qc_stats if s not in sample_tracker_df["Sample"].values]
if not new_samples:
    print("No new samples to process. Skipping run.")
    exit(0)
print(f"Found {len(new_samples)} new samples to process.")

# === Process New Samples ===
batch_records = []
for sample in new_samples:
    metrics = qc_stats[sample]
    fails = metrics.get("percent_fails", 0)
    dupes = metrics.get("percent_duplicates", 0)
    reads = metrics.get("total_sequences", 0)

    reasons = []
    if fails > 10.0:
        reasons.append("High FastQC Fails")
    if dupes > 50.0:
        reasons.append("High Duplication")
    if reads < 500:
        reasons.append("Low Read Count")

    qc_status = "Defer" if reasons else "Good"
    batch_records.append({
        "Sample": sample,
        "QC_Status": qc_status,
        "Reason": "; ".join(reasons),
        "Batch_File": "",  # Will be assigned below
        "Pipeline_Version": PIPELINE_VERSION
    })

# === Append New Samples to Sample Tracker ===
new_df = pd.DataFrame(batch_records)
sample_tracker_df = pd.concat([sample_tracker_df, new_df], ignore_index=True)

# === Batch Good Samples ===
good_samples = new_df[new_df["QC_Status"] == "Good"]["Sample"].tolist()
batch_number = 1

existing_batches = sorted(good_dir.glob("g_batch*.txt"))
if existing_batches:
    last_batch = int(existing_batches[-1].stem.split("g_batch")[-1])
    batch_number = last_batch + 1

if good_samples:
    batch_file = good_dir / f"g_batch{batch_number:02d}.txt"
    with open(batch_file, "a") as f:
        for sample in good_samples:
            f.write(f"{sample}\n")
    sample_tracker_df.loc[sample_tracker_df["Sample"].isin(good_samples), "Batch_File"] = batch_file.name
    print(f"Added {len(good_samples)} Good samples to {batch_file.name}")

# === Batch Deferred Samples by Reason ===
defer_df = new_df[new_df["QC_Status"] == "Defer"]
if not defer_df.empty:
    for reason in defer_df["Reason"].unique():
        reason_samples = defer_df[defer_df["Reason"] == reason]["Sample"].tolist()
        safe_reason = re.sub(r'[^A-Za-z0-9_]+', '_', reason.strip())
        batch_file = defer_dir / f"d_{safe_reason}_batch{batch_number:02d}.txt"
        with open(batch_file, "a") as f:
            for sample in reason_samples:
                f.write(f"{sample}\n")
        sample_tracker_df.loc[sample_tracker_df["Sample"].isin(reason_samples), "Batch_File"] = batch_file.name
        print(f"Added {len(reason_samples)} Deferred samples to {batch_file.name}")

# === Save Updated Sample Tracker ===
sample_tracker_df.to_csv(sample_tracker_path, index=False)

# === Log Run ===
with open(pipeline_log, "a") as log:
    log.write(f"Run on {datetime.now().isoformat()} | {len(new_samples)} new samples processed.\n")

print("Pipeline run complete. Tracker and logs updated.")