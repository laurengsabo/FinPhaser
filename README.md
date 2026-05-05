# FinPhaser
<p align="center">
  <b>A reproducible framework for local ancestry inference and IBD detection in hybrid populations</b>
</p>
<p align="center">
  <img src="https://img.shields.io/badge/python-3.12-blue.svg" />
  <img src="https://img.shields.io/badge/R-4.3+-blue.svg" />
  <img src="https://img.shields.io/badge/nextflow-24.04-brightgreen.svg" />
  <img src="https://img.shields.io/badge/conda-environment-green.svg" />
  <img src="https://img.shields.io/badge/status-active-success.svg" />
  <img src="https://img.shields.io/badge/license-MIT-lightgrey.svg" />
</p>

---

## Overview

**FinPhaser** is a bioinformatics pipeline for estimating **local ancestry** and **Identity-By-Descent (IBD)** in hybrid populations, with a focus on African cichlid genomics.

The framework integrates:

- Variant filtering and phasing (PhaseParents_VCF.py)
- Hidden Markov Models for local ancestry inference (ancestry_hmm)
- Admixture analysis using SPORE
- IBD detection and ranked relatedness scoring

It is particularly suited for studying **complex sex determination systems**, including ZW/XY dynamics on the LG10 linkage group in *Aulonocara* 'Yellow Head'.

The entire pipeline runs in **a single command** via [Nextflow](https://www.nextflow.io/), requiring only two input files: a filtered VCF and a samples YAML.

---

## Quick Start

### Prerequisites

| Tool | Minimum version | Install |
|------|----------------|---------|
| Java | 11+ | `brew install openjdk` / system package manager |
| Nextflow | 23.04+ | See below |
| Conda or Mamba | any recent | [Miniforge](https://github.com/conda-forge/miniforge) recommended |

**Install Nextflow** (if not already installed):
```bash
# Option A — direct install (puts nextflow in your PATH)
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Option B — via conda (already included in environment.yml)
conda install -c conda-forge nextflow
```

### 1. Clone the repository

```bash
git clone https://github.com/your-username/FinPhaser.git
cd FinPhaser
```

### 2. Create the conda environment

```bash
conda env create -f environment.yml
conda activate FinPhaser
```

> **Faster alternative:** replace `conda` with `mamba` in the command above, or use the `mamba` profile (see [Execution Profiles](#execution-profiles)).

### 3. Configure your samples

Edit `config/samples.yml` to match your data.  Two things to set:

**a) Population assignments** — assign every VCF sample to a reference population (0–3) or mark it as admixed (-1):

```yaml
populations:
  - sample: "YH_P1"
    population: 0     # Paternal Haplotype 1 (reference)
  - sample: "YH_016"
    population: -1    # Admixed offspring — ancestry inferred by HMM
```

**b) SPORE sex metadata path:**

```yaml
spore:
  sex_metadata: "data/raw/Genomics_Sex.tsv"
  output_prefix: "brood1_final"
```

### 4. Run the pipeline

```bash
nextflow run main.nf -profile conda \
    --vcf         data/raw/YHPedigree1_FilteredSNVs.recode.vcf \
    --samples_yml config/samples.yml
```

That's it. Nextflow manages every step — no manual script invocations needed.

#### Quick test (bundled data)

```bash
nextflow run main.nf -profile test
```

---

## Execution Profiles

Select a profile with `-profile <name>`:

| Profile | Description |
|---------|-------------|
| `conda` | Local conda environment — recommended default |
| `mamba` | Same as conda but uses mamba for faster solves |
| `docker` | Docker container (see [Docker](#docker)) |
| `singularity` | Singularity container — best for HPC clusters |
| `test` | Runs bundled test dataset; completes in minutes |

### Docker

If you prefer a fully self-contained container:

```bash
# Build the image once
docker build -t finphaser:latest .

# Run the pipeline
nextflow run main.nf -profile docker \
    --vcf         data/raw/YHPedigree1_FilteredSNVs.recode.vcf \
    --samples_yml config/samples.yml
```

### HPC / Singularity

```bash
# Convert the Docker image to a Singularity image
singularity build finphaser.sif docker://finphaser:latest

# Run
nextflow run main.nf -profile singularity \
    --vcf         data/raw/YHPedigree1_FilteredSNVs.recode.vcf \
    --samples_yml config/samples.yml
```

---

## Pipeline Workflow

The pipeline runs all steps in sequence, with automatic dependency resolution:

```
YHPedigree1_FilteredSNVs.recode.vcf
samples.yml
        │
        ▼
┌─────────────────────────────────┐
│  Step 1 — PHASE_PARENTS         │  PhaseParents_VCF.py
│  Identify conserved SNVs        │  → YHPed1_conserved_phased.vcf
└────────────────┬────────────────┘
                 │
                 ▼
┌─────────────────────────────────┐
│  Steps 2+3 — PREPARE_HMM_INPUT  │  mod_vcf2ahmm.py --pop-yaml
│  Pop assignments from YAML      │  → ancestry_input.txt
│  VCF → ancestry_hmm format      │  → ahmm.ploidy
└────────────────┬────────────────┘
                 │
                 ▼
┌─────────────────────────────────┐
│  Step 4 — RUN_ANCESTRY_HMM      │  ancestry_hmm
│  Local ancestry posteriors      │  → *.posterior (per sample)
└────────────────┬────────────────┘
                 │
                 ▼
┌─────────────────────────────────┐
│  Step 5 — ANCESTRY_SUMMARY      │  SNV_count.py
│  Readable summary CSVs          │  → ancestry_summary_report/
└────────────────┬────────────────┘
                 │
                 ▼
┌─────────────────────────────────┐
│  Step 6 — RUN_SPORE             │  SPORE.R
│  IBD segment detection          │  → *.ibd
└────────────────┬────────────────┘
                 │
                 ▼
┌─────────────────────────────────┐
│  Step 7 — RANK_IBD              │  rank_ibd.py
│  Ranked relatedness scores      │  → inbreeding_rankings.txt/.tsv
└─────────────────────────────────┘
```

Steps 2 and 3 are combined into a single process: `mod_vcf2ahmm.py` now accepts `--pop-yaml` and builds the internal population info file automatically from `samples.yml`.

---

## Inputs

| File | Description |
|------|-------------|
| `data/raw/YHPedigree1_FilteredSNVs.recode.vcf` | Raw filtered SNV VCF — the only genomic input required |
| `config/samples.yml` | Population assignments + SPORE settings |
| `data/raw/Genomics_Sex.tsv` | Sex metadata referenced inside `samples.yml` |

### samples.yml structure

```yaml
spore:
  sex_metadata: "data/raw/Genomics_Sex.tsv"
  output_prefix: "brood1_final"

populations:
  - sample: "YH_P1"
    population: 0       # 0–3 = reference population index
  - sample: "YH_016"
    population: -1      # -1 = admixed; ancestry inferred by HMM
```

Population codes follow `ancestry_hmm` convention:
- `0` = Paternal Haplotype 1
- `1` = Paternal Haplotype 2
- `2` = Maternal Haplotype 1
- `3` = Maternal Haplotype 2
- `-1` = Admixed offspring

---

## Outputs

All outputs land in `results/` (or the directory specified by `--outdir`):

```
results/
├── ancestry/
│   ├── 01_phased/
│   │   └── YHPed1_conserved_phased.vcf
│   ├── 02_hmm_input/
│   │   ├── ancestry_input.txt
│   │   └── ahmm.ploidy
│   ├── 03_posteriors/
│   │   └── YH_016.posterior ... YH_042.posterior
│   └── ancestry_summary_report/
│       ├── calls_YH_016.csv
│       ├── ...
│       └── master_ancestry_summary.csv
├── ibd/
│   ├── 06_spore/
│   │   ├── brood1_final.vcf.gz-truffle.ibd
│   │   └── SPORE_output.log
│   └── 07_rankings/
│       ├── inbreeding_rankings.txt
│       └── inbreeding_rankings.tsv
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── pipeline_dag.html
```

| File | Description |
|------|-------------|
| `*.posterior` | Local ancestry probabilities per genomic site, per sample |
| `master_ancestry_summary.csv` | Aggregated ancestry calls across all samples |
| `*.ibd` | IBD segment calls from SPORE/TRUFFLE |
| `inbreeding_rankings.tsv` | Ranked relatedness scores (TSV for downstream plotting) |
| `execution_report.html` | Nextflow HTML report with runtime and resource usage |

---

## Repository Structure

```
.
├── main.nf                  # Pipeline entry point
├── nextflow.config          # Profiles, resource limits, reporting
├── environment.yml          # Pinned conda environment
├── Dockerfile               # Container definition (for -profile docker)
├── config/
│   ├── samples.yml          # User-facing config: pop assignments + SPORE settings
│   └── test_samples.yml     # Config for bundled test run
├── bin/                     # Pipeline helper scripts (called by Nextflow)
│   ├── generate_spore_settings.py
│   └── txt_to_tsv.py
├── workflows/               # Nextflow subworkflow definitions
│   ├── ancestry.nf          # Steps 1–5
│   └── ibd.nf               # Steps 6–7
├── src/
│   ├── ancestry/
│   │   ├── PhaseParents_VCF.py
│   │   ├── mod_vcf2ahmm.py  # Modified — now accepts --pop-yaml
│   │   └── SNV_count.py
│   ├── reporting/
│   │   └── rank_ibd.py
│   └── SPORE.R
├── data/
│   ├── raw/                 # Input files (VCF, sex metadata TSV)
│   └── test/                # Bundled test dataset
├── docs/                    # Architecture diagrams and notes
└── results/                 # Pipeline outputs (git-ignored)
```

---

## Troubleshooting

**Sample name formatting (SPORE)**
SPORE does not allow underscores in sample IDs. FinPhaser handles this automatically — the pipeline sanitises names before SPORE runs. You do not need to rename samples in your VCF or YAML.

**Unexpected LG10 results**
Verify the HMM pulse parameters in `nextflow.config` match your expected admixture history. The default (`-p 0..3 -2 0.25`) assumes ~2 generations of admixture. Adjust if your crossing design differs.

**Dependency or environment issues**
Ensure the conda environment is active (`conda activate FinPhaser`) if you are running scripts directly outside Nextflow. Within a Nextflow run the environment is activated automatically by the selected profile.

**Resuming a failed run**
Nextflow caches completed steps. If a run fails midway, fix the issue and resume without re-running successful steps:
```bash
nextflow run main.nf -profile conda -resume
```

**Nextflow not found**
If `nextflow` is not in your PATH after activating the conda environment, install it directly:
```bash
curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/
```

**Viewing logs**
Each step writes its stdout/stderr to `.nextflow/` and to `results/pipeline_info/`. Check `execution_report.html` for a summary of which processes ran, their duration, and exit codes.

---

## Contact

For questions or issues, please open a GitHub Issue or contact the project maintainer.
