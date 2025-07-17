
TSight Stage 1 Pipeline: RNA-Seq QC + Batching

###################################################
This pipeline performs:
1. Quality control (QC) with FastQC + MultiQC.
2. Auto-batching based on QC results.

###################################################
Folder Structure

/pipeline/
│
├── run_qc.sh             # Bash script for QC (FastQC + MultiQC)
├── auto_batching.py      # Python script for batching (requires pandas)
├── run_pipeline.py       # Main automation script (QC + Batching)
│
/data/samples/            # Input FASTQ/FQ samples
/qc/fastqc_reports/       # QC results
/batches/                 # Batch outputs
/tracking/                # Sample tracker & pipeline logs

###################################################
Environments (Separated to Avoid Dependency Conflicts):

QC Environment: 
conda create -n rnaseq-env-qc python=3.9 fastqc multiqc numpy=1.19 networkx=2.5 -c bioconda -c conda-forge

Batching Environment:
conda create -n rnaseq-env-batching python=3.9 pandas numpy=1.24 -c conda-forge

###################################################
Dependencies:

| Tool     | Environment         | Version / Notes                    |
| -------- | ------------------- | ---------------------------------- |
| FastQC   | rnaseq-env-qc       | Bioconda version                   |
| MultiQC  | rnaseq-env-qc       | Bioconda version                   |
| numpy    | rnaseq-env-qc       | 1.19.x (for MultiQC)               |
| networkx | rnaseq-env-qc       | 2.5                                |
| pandas   | rnaseq-env-batching | Latest compatible with numpy ≥1.24 |
| numpy    | rnaseq-env-batching | ≥1.24 (compatible with pandas)     |

Environment: rnaseq-env-qc
FastQC v0.12.1
multiqc, version 1.29
numpy 1.19.5
networkx 2.5.1

Environment: rnaseq-env-batching
pandas 2.3.0
numpy 1.24.4

###################################################
How to Run

1. python run_pipeline.py


###################################################

Notes:
The pipeline automatically removes old QC reports before each run.

No more redundant MultiQC reports or batch files.

Safe for cloud pipelines and ready for Firebase backend integration.