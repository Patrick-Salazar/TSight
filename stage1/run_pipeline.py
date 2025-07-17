import subprocess
import os

os.chdir(os.path.dirname(os.path.dirname(__file__)))

print("Running QC Stage...")
subprocess.run(
    ["conda", "run", "-n", "rnaseq-env-qc", "bash", "stage1/run_qc.sh"],
    check=True
)


print("Running Batching Stage...")
subprocess.run(
    ["conda", "run", "-n", "rnaseq-env-batching", "python", "auto_batching.py"],
    check=True,
    cwd=os.path.dirname(__file__)
)

print("Pipeline Completed.")
