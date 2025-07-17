TSight Pipeline - Stage 2: Trimming & Secondary QC  
===================================================

Overview
--------
This module automates **Stage 2** of the TSight pipeline:
- QC-driven adapter trimming using **fastp** based on Stage 1 QC results.
- Post-trimming quality check via **FastQC** and **MultiQC**.
- Batch-aware post-trimming QC checkpoint for deferred sample review.

===================================================
Folder Structure
----------------
| Folder                          | Description                                       |
|---------------------------------|---------------------------------------------------|
| `data/trimmed/`                 | Trimmed FASTQ files.                              |
| `qc/fastp/`                     | fastp reports (JSON & HTML) for trimming runs.    |
| `qc/fastqc_reports/`            | FastQC + MultiQC post-trimming QC reports.        |
| `batches/good/`                 | Good sample batches from Stage 1 (input to Stage 2). |
| `batches/defer/`                | Deferred sample lists flagged after trimming QC.  |
| `pipeline/stage2/config/`       | YAML configuration files for trimming and QC thresholds. |
| `tracking/`                     | Sample tracker and pipeline logs.                 |

===================================================
Configuration Files
-------------------
Located in `pipeline/stage2/config/`:
| File                    | Purpose                                             |
|-------------------------|-----------------------------------------------------|
| `trimming_config.yaml`  | fastp trimming parameters per strategy (`normal`, `aggressive`). |
| `qc_thresholds.yaml`    | QC thresholds for post-trimming QC checkpoint (e.g., percent fails, duplication, read count). |

===================================================
Software & Environments
-----------------------
| Environment         | Software                                | Versions                             |
|---------------------|----------------------------------------|--------------------------------------|
| rnaseq-env-qc       | FastQC                                  | v0.12.1                              |
|                     | MultiQC                                 | v1.29                                |
|                     | numpy                                   | 1.19.5                               |
|                     | networkx                                | 2.5.1                                |
| rnaseq-env-batching | pandas                                  | 2.3.0                                |
|                     | numpy                                   | 1.24.4                               |

===================================================
Dependencies
------------
- Python 3.x
- fastp
- fastqc
- multiqc
- Python packages: `pyyaml`, `pandas`

===================================================
How to Run
----------
1. Update `trimming_config.yaml` and `qc_thresholds.yaml` as needed.
2. Execute the Stage 2 Orchestrator:
```bash
bash stage2_orchestrator.sh
