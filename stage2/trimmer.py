import json
import subprocess
import sys
import yaml
from pathlib import Path

# === PATH SETUP ===
SCRIPT_DIR            = Path(__file__).parent.resolve()
BASE_DIR              = SCRIPT_DIR.parent.resolve()
TRIMMING_PLAN_PATH    = SCRIPT_DIR / "trimming_plan.json"
TRIMMING_CONFIG_PATH  = SCRIPT_DIR / "config" / "trimming_config.yaml"
SAMPLES_DIR           = BASE_DIR / "data" / "samples"
OUTPUT_DIR            = BASE_DIR / "data" / "trimmed"
FASTP_REPORT_DIR      = BASE_DIR / "qc" / "fastp"

# === SANITY CHECKS ===
for p in (TRIMMING_PLAN_PATH, TRIMMING_CONFIG_PATH):
    if not p.is_file():
        sys.exit(f"ERROR: required file not found: {p}")
if not SAMPLES_DIR.is_dir():
    sys.exit(f"ERROR: samples directory not found: {SAMPLES_DIR}")

# === LOAD CONFIGS ===
def load_json(path):
    with open(path, 'r') as f:
        return json.load(f)

def load_yaml(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

trimming_plan   = load_json(TRIMMING_PLAN_PATH)
trimming_config = load_yaml(TRIMMING_CONFIG_PATH)

# === PREPARE OUTPUTS ===
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FASTP_REPORT_DIR.mkdir(parents=True, exist_ok=True)

# === FASTP WRAPPER ===
def run_fastp(input_path, output_path, report_prefix, cfg):
    cmd = [
        "fastp",
        "-i",  str(input_path),
        "-o",  str(output_path),
        "-j",  f"{report_prefix}.json",
        "-h",  f"{report_prefix}.html",
        "--adapter_sequence", cfg["adapter_sequence"],
        "--cut_window_size",    str(cfg["cut_window_size"]),
        "--cut_mean_quality",   str(cfg["cut_mean_quality"]),
        "--qualified_quality_phred", str(cfg["qualified_quality_phred"]),
        "--length_required",    str(cfg["min_length"])
    ]
    if cfg.get("cut_front"):
        cmd.append("--cut_front")
    if cfg.get("cut_tail"):
        cmd.append("--cut_tail")
    subprocess.run(cmd, check=True)

# === EXECUTE TRIMMING ===
for sample, action in trimming_plan.items():
    matches = list(SAMPLES_DIR.glob(f"{sample}.*q"))
    if not matches:
        print(f"Sample not found, skipping: {sample}")
        continue

    if action == "skip":
        print(f"Skipping {sample} (QC clean).")
        continue

    cfg = trimming_config.get(action)
    if not cfg:
        print(f"Unknown action '{action}' for {sample}, skipping.")
        continue

    inp  = matches[0]
    outp = OUTPUT_DIR / f"{sample}_trimmed.fastq"
    rpt  = FASTP_REPORT_DIR / sample

    print(f"Trimming {sample} with '{action}' settings...")
    run_fastp(inp, outp, rpt, cfg)
    print(f"Completed: {sample}")

print("All samples processed.")
