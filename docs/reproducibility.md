# Reproducibility Checklist

This repository provides all scripts and configuration necessary to re-create the hydrogen-peroxide stress analysis described in the README. Follow the steps below on a Linux machine with Docker and Nextflow installed.

## 1. Clone and Inspect

```bash
git clone https://github.com/<your-user>/yeast-oxidative-stress-analysis.git
cd yeast-oxidative-stress-analysis
```

The `metadata/samples.csv` file defines every sample, condition, and expected FASTQ path. The example entries correspond to SRR636633–SRR636638 (PRJNA184040).

## 2. Prepare the Conda Environment

Install [Mamba](https://github.com/mamba-org/mamba) or Conda, then create the shared environment:

```bash
mamba env create -f environment.yml
conda activate yeast-oxidative-stress-analysis
```

This adds Nextflow, RNA-seq tools, and R packages to your PATH so every downstream command works without Docker.

## 3. Download Raw FASTQs

Requires the [SRA Toolkit](https://github.com/ncbi/sra-tools) (`fasterq-dump`) and `pigz`, both already included in `environment.yml`. The helper script reads `metadata/samples.csv` and places all files in `data/`.

```bash
PARALLEL=3 THREADS=6 make download-data
```

This command automatically computes SHA256 checksums (`metadata/fastq_checksums.txt`) so you can verify file integrity across runs.

## 4. Prepare the sacCer3 Reference

Requires `curl`, `gunzip`, `hisat2-build`, and (optionally) `gffread` for GTF conversion (all provided by the Conda env).

```bash
BUILD_THREADS=12 make prepare-reference
```

Outputs:
- HISAT2 index: `data/genome/sacCer3.*.ht2`
- Annotation: `data/Saccharomyces_cerevisiae.R64-1-1.gtf`

## 5. Execute the Pipeline

With all inputs available:

```bash
nextflow run src/main.nf -profile conda \
    --data_dir ./data \
    --results_dir ./results \
    --samplesheet ./metadata/samples.csv
```

Nextflow orchestrates FastQC → fastp → HISAT2 → samtools → featureCounts → edgeR. Key outputs land under `results/` (see README for locations).

## 6. Capture Environment Details

For full reproducibility, record the versions of each tool:

```bash
nextflow -version
docker --version
fasterq-dump --version
hisat2 --version
featureCounts -v
Rscript -e 'sessionInfo()'
```

Store the resulting log alongside `metadata/fastq_checksums.txt` when publishing analyses.

## 7. Clean and Re-run

Use `make clean-results` to remove pipeline outputs without touching downloaded data, then re-run to confirm deterministic results.

These steps ensure anyone with the same repository can regenerate the figures and tables accompanying this study.
