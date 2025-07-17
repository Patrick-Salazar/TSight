import json
from pathlib import Path

# === CONFIGURATION ===
# All paths are relative to TSight/stage2/
MULTIQC_JSON_PATH = Path("../qc/fastqc_reports/multiqc_data/multiqc_data.json")
GOOD_BATCH_PATH   = Path("../batches/good/")
DEFER_BATCH_PATH  = Path("../batches/defer/")
OUTPUT_PLAN_PATH  = Path("trimming_plan.json")

# === LOAD MULTIQC DATA ===
def load_multiqc_data(json_path):
    if not json_path.is_file():
        raise FileNotFoundError(f"No such file: {json_path}")
    with open(json_path, 'r') as f:
        return json.load(f)

multiqc_data = load_multiqc_data(MULTIQC_JSON_PATH)

# === LOAD BATCH LISTS ===
def read_sample_list(txt_path):
    with open(txt_path) as f:
        return [line.strip() for line in f if line.strip()]

def load_samples_by_priority():
    samples = []
    # 1) good samples
    for txt in sorted(GOOD_BATCH_PATH.glob("g_batch*.txt")):
        samples.extend(read_sample_list(txt))
    # 2) deferred samples
    for txt in sorted(DEFER_BATCH_PATH.glob("postqc_*.txt")):
        samples.extend(read_sample_list(txt))
    return samples

priority_list = load_samples_by_priority()

# === PARSE FASTQC MODULE STATUSES ===
fastqc_modules = multiqc_data.get("report_saved_raw_data", {}) \
                              .get("multiqc_fastqc", {})

statuses = {
    sample: {
        'Per Base Sequence Quality': mods.get('per_base_sequence_quality', 'pass'),
        'Adapter Content':              mods.get('adapter_content', 'pass'),
        'Overrepresented Sequences':    mods.get('overrepresented_sequences', 'pass'),
    }
    for sample, mods in fastqc_modules.items()
}

# === DECIDE TRIMMING STRATEGY ===
def decide_trimming(mods):
    if all(v == 'pass' for v in mods.values()):
        return "skip"
    if any(mods[k] in ('fail','warn') for k in mods):
        return "normal"
    return "skip"

# === BUILD PLAN ===
trimming_plan = {}
for sample in priority_list:
    match = next((s for s in statuses if s.startswith(sample)), None)
    trimming_plan[sample] = decide_trimming(statuses[match]) if match else "skip"

# === WRITE OUT ===
with open(OUTPUT_PLAN_PATH, 'w') as f:
    json.dump(trimming_plan, f, indent=4)

print(f"Trimming plan saved to {OUTPUT_PLAN_PATH}")
