# FinPhaser
<p align="center">
  <b>A reproducible framework for local ancestry inference and IBD detection in hybrid populations</b>
</p>
<p align="center">
  <img src="https://img.shields.io/badge/python-3.10-blue.svg" />
  <img src="https://img.shields.io/badge/R-4.3+-blue.svg" />
  <img src="https://img.shields.io/badge/nextflow-24.04-brightgreen.svg" />
  <img src="https://img.shields.io/badge/conda-environment-green.svg" />
  <img src="https://img.shields.io/badge/platform-Linux-orange.svg" />
  <img src="https://img.shields.io/badge/status-active-success.svg" />
  <img src="https://img.shields.io/badge/license-MIT-lightgrey.svg" />
</p>

---

<p align="center">
  <img src="docs/figures/finphaser.png" width="400" alt="FinPhaser Pipeline Overview"/>
</p>

---

## Overview

**FinPhaser** is a bioinformatics pipeline for estimating **local ancestry** and **Identity-By-Descent (IBD)** in hybrid populations, with a focus on African cichlid genomics.

The framework integrates:

- Variant filtering and phasing (`PhaseParents_VCF.py`)
- Hidden Markov Models for local ancestry inference (`ancestry_hmm`)
- IBD detection using SPORE + TRUFFLE
- Ranked relatedness scoring

It is particularly suited for studying **complex sex determination systems**, including ZW/XY dynamics on the LG10 linkage group in *Aulonocara* 'Yellow Head'. By backcrossing a male YH with *Mchenga conophoros*, parental haplotypes on LG10 can be isolated and tracked in F1 offspring. SPORE is used to assess IBD among offspring, quantifying inbreeding accumulated across many generations of lab-bred stock.

The entire pipeline runs in **a single command** via [Nextflow](https://www.nextflow.io/), requiring only **one input file**: `config/samples.yml`, which contains the VCF path and all pipeline settings.

---

## Platform Compatibility

> ⚠️ **Linux is required for full pipeline execution.**

| Step | Linux | macOS |
|------|-------|-------|
| PHASE_PARENTS | ✔ | ✔ |
| PREPARE_HMM_INPUT | ✔ | ✔ |
| RUN_ANCESTRY_HMM | ✔ | ✔ (slow — single-threaded, memory-intensive) |
| ANCESTRY_SUMMARY | ✔ | ✔ |
| COMPRESS_VCF | ✔ | ✔ |
| RUN_SPORE / TRUFFLE | ✔ | ✘ — TRUFFLE binary is Linux x86-64 only |
| RANK_IBD | ✔ | ✔ |

The TRUFFLE binary (`tools/truffle/truffle`) is compiled for Linux x86-64 and **cannot run on macOS ARM (Apple Silicon)** or macOS x86-64 without a Linux x86-64 environment. Run the full pipeline on a Linux server or HPC cluster.

**Recommended:** Linux x86-64, 16+ GB RAM, 8+ cores.

---

## Quick Start

### Prerequisites

| Tool | Minimum version | Platform | Install |
|------|----------------|----------|---------|
| Java | 11+ (17 recommended) | Linux / macOS | system package manager |
| Nextflow | 23.04+ | Linux / macOS | See below |
| Conda or Mamba | any recent | Linux / macOS | [Miniforge](https://github.com/conda-forge/miniforge) recommended |

**Install Nextflow:**
```bash
# Option A — direct install
curl -s https://get.nextflow.io | bash
mkdir -p ~/bin && mv nextflow ~/bin/
export PATH=$HOME/bin:$PATH

# Option B — via conda
conda install -c conda-forge nextflow
```

---

### Option A — Quick test (recommended first run)

No input files needed — everything is bundled. Expected runtime: **30–60 minutes on Linux**.

```bash
# 1. Clone the repository
git clone https://github.com/laurengsabo/FinPhaser.git
cd FinPhaser

# 2. Run the bundled test dataset
nextflow run main.nf -profile test
```

Results will appear in `test_results/`. See [Test Dataset](#test-dataset) for expected outputs.

---

### Option B — Full run on your own data

#### 1. Clone the repository

```bash
git clone https://github.com/laurengsabo/FinPhaser.git
cd FinPhaser
```

#### 2. Place your input VCF

```bash
mkdir -p data/raw
cp /path/to/YHPedigree1_FilteredSNVs.recode.vcf.gz data/raw/
```

#### 3. Conda environments

FinPhaser uses **process-specific Conda environments** — you do not need to manually create or activate them. Nextflow creates, activates, and caches them automatically per step:

- `envs/finphaser.yml` — Python + ancestry_hmm + bcftools
- `envs/spore.yml` — R 4.3.3 + SPORE dependencies

#### 4. Configure your samples

Edit `config/samples.yml`:

**a) Set your VCF path:**
```yaml
vcf: "data/raw/YHPedigree1_FilteredSNVs.recode.vcf.gz"
```

**b) Population assignments** — each diploid parent appears **twice** (once per haplotype):
```yaml
populations:
  - sample: "YH_006_f"
    population: 0
    sex: F
  - sample: "YH_006_f"
    population: 1
    sex: F
  - sample: "YH_011_m"
    population: 2
    sex: M
  - sample: "YH_011_m"
    population: 3
    sex: M
  - { sample: "YH_016", population: admixed, sex: F }
```

**c) TRUFFLE path** — use the absolute path to the binary on your system:
```yaml
spore:
  truffle_path: "/absolute/path/to/FinPhaser/tools/truffle/truffle"
```

> **Note:** Relative paths for `truffle_path` resolve to the Nextflow process working directory, not the repo root. Always use an absolute path. Run `realpath tools/truffle/truffle` from the repo root to get the correct path.

#### 5. Run the pipeline

```bash
nextflow run main.nf -profile conda
```

**Resume after a failure:**
```bash
nextflow run main.nf -profile conda -resume
```

---

## Execution Profiles

| Profile | Description |
|---------|-------------|
| `conda` | Per-process conda environments — recommended default |
| `mamba` | Same as conda but uses mamba for faster environment solves |
| `test` | Runs bundled 50k-SNP test dataset; completes in ~30–60 min on Linux |

## Reproducibility

FinPhaser uses **process-level environment isolation** via Nextflow:

- Each pipeline step runs in its own Conda environment
- SPORE's R dependencies are isolated from the modern Python/bioinformatics environment
- All environments are fully version-pinned in `envs/`
- The entire pipeline configuration is contained in a single `samples.yml` file

---

## Pipeline Workflow

One starting VCF (declared in `samples.yml`) feeds both the ancestry and IBD branches after phasing:

```
config/samples.yml  ←─ vcf: path + all settings
        │
        │  VCF path read at startup
        ▼
YHPedigree1_FilteredSNVs.recode.vcf.gz
        │
        ▼
┌──────────────────────────────────────┐
│  Step 1 — PHASE_PARENTS              │  PhaseParents_VCF.py
│  Identify conserved informative SNVs │  → YHPed1_conserved_phased.vcf
└────────────┬─────────────────────────┘
             │
     ┌───────┴────────┐
     │                │
     ▼                ▼
 AHMM branch      IBD branch
     │                │
     ▼                ▼
┌─────────────────┐  ┌──────────────────────────┐
│  Steps 2+3      │  │  Step 6a — COMPRESS_VCF   │
│  PREPARE_HMM    │  │  bgzip + CSI index        │
│  INPUT          │  │                           │
│                 │  │  Step 6b — WRITE_SEX_TSV  │
│  yaml_to_       │  │  yaml_to_sex_tsv.py       │
│  popinfo.py     │  │  → Genomics_Sex.tsv       │
│  → mod_         │  └──────────┬────────────────┘
│    popinfo.txt  │             │
│                 │             ▼
│  mod_vcf2ahmm   │  ┌──────────────────────────┐
│  .py (-s flag)  │  │  Step 6c — RUN_SPORE      │
│  → ancestry_    │  │  SPORE.R + TRUFFLE        │
│    input.txt    │  │  → *-truffle.ibd.iqr.gz   │
│  → ahmm.ploidy  │  └──────────┬────────────────┘
└──────┬──────────┘             │
       │                        ▼
       ▼             ┌──────────────────────────┐
┌──────────────────┐ │  Step 7 — RANK_IBD       │
│  Step 4          │ │  rank_ibd.py             │
│  RUN_ANCESTRY_   │ │  → inbreeding_rankings   │
│  HMM             │ │    .txt / .tsv           │
│  → *.posterior   │ └──────────────────────────┘
└──────┬───────────┘
       ▼
┌──────────────────┐
│  Step 5          │
│  ANCESTRY_       │
│  SUMMARY         │
│  SNV_count.py    │
│  → summary CSVs  │
└──────────────────┘
```

**Key design decisions:**
- The VCF path is declared in `samples.yml` under `vcf:` — the single source of truth. Change one line to switch input files.
- Steps 2 and 3 are combined: `bin/yaml_to_popinfo.py` converts `samples.yml` into the population info format that `mod_vcf2ahmm.py` reads via its `-s` flag. `mod_vcf2ahmm.py` is **not modified**.
- The chromosome → linkage group mapping lives in `samples.yml` under `linkage_groups:` and is extracted at runtime by `bin/yaml_to_linkage_map.py`. Edit there to support a different reference assembly.
- Sex metadata is generated at runtime by `bin/yaml_to_sex_tsv.py` from the `sex` field in `samples.yml`. Sample names are written **with underscores preserved** to match VCF column headers exactly.
- `SPORE-Settings.R` is fully auto-generated from `samples.yml` — no manual editing required.
- `NoTrioCalculation` is set to `TRUE` in `SPORE.R` — Mendelian trio errors require a known pedigree which is not available in this pipeline.

---

## Test Dataset

The bundled test dataset lives in `data/test/` and is used by `-profile test`.

| File | Description |
|------|-------------|
| `data/test/test_FilteredSNVs.recode.vcf.gz` | 50,000 SNPs from the first chromosome of the full dataset, bgzipped (~13 MB) |
| `config/test_samples.yml` | Config file pointing to the test VCF with appropriate settings |

**Expected test run outputs** (in `test_results/`):
- 24 `.posterior` files (one per admixed sample, YH_016–YH_042)
- `master_ancestry_summary.csv` with ancestry state calls
- `*-truffle.ibd.iqr.gz` IBD segment calls (compressed by SPORE)
- `inbreeding_rankings.txt` with ranked relatedness scores

**Expected runtime:** 30–60 minutes on a Linux server with 8 cores and 16 GB RAM (verified: 50 minutes on Georgia Tech ICE cluster).

---

## Inputs

Only **one file** is required to start the pipeline:

| File | Description |
|------|-------------|
| `config/samples.yml` | All configuration: VCF path, population assignments, sex, linkage groups, HMM params, SPORE settings |

The VCF itself is referenced inside `samples.yml` under the `vcf:` key and must be present at that path.

### samples.yml structure

```yaml
# ── Input VCF ─────────────────────────────────────────────────
vcf: "data/raw/YHPedigree1_FilteredSNVs_quarter.recode.vcf.gz"

# ── Chromosome → linkage group mapping ────────────────────────
linkage_groups:
  NC_036780.1: LG1
  NC_036789.1: LG10
  # ... (see config/samples.yml for full list)

# ── HMM parameters ────────────────────────────────────────────
hmm:
  recombination_rate: 1.0e-8
  min_distance_bp: 1000
  min_allele_freq_diff: 0.1
  use_genotypes: 0

# ── SPORE / TRUFFLE settings ───────────────────────────────────
spore:
  truffle_path: "/absolute/path/to/tools/truffle/truffle"
  truffle_maf: 0.05
  truffle_missing: 0.95
  APO: 24
  output_label: "YHPed1_brood1"
  max_memory: "4G"
  max_cores: 4
  plots: TRUE

# ── Sample assignments ─────────────────────────────────────────
populations:
  - sample: "YH_006_f"
    population: 0
    sex: F
  - sample: "YH_006_f"
    population: 1
    sex: F
  - sample: "YH_011_m"
    population: 2
    sex: M
  - sample: "YH_011_m"
    population: 3
    sex: M
  - { sample: "YH_016", population: admixed, sex: F }
```

**Population codes:**

| Code | Meaning |
|------|---------|
| `0` | Maternal haplotype 1 |
| `1` | Maternal haplotype 2 |
| `2` | Paternal haplotype 1 |
| `3` | Paternal haplotype 2 |
| `admixed` | Offspring — ancestry inferred by HMM |

> Parents appear **twice** — once per haplotype — because each diploid parent contributes two distinct haplotypes to the 4-population HMM model.

---

## Outputs

All outputs land in `results/` (or `test_results/` for `-profile test`):

```
results/
├── ancestry/
│   ├── 01_phased/
│   │   ├── YHPed1_conserved_phased.vcf
│   │   └── linkage_map.tsv
│   ├── 02_hmm_input/
│   │   ├── ancestry_input.txt
│   │   ├── ahmm.ploidy
│   │   └── mod_popinfo.txt
│   ├── 03_posteriors/
│   │   └── YH_016.posterior ... YH_042.posterior
│   └── ancestry_summary_report/
│       ├── calls_YH_016.posterior.csv
│       ├── ...
│       └── master_ancestry_summary.csv
├── ibd/
│   ├── 06_spore_input/
│   │   ├── YHPed1_conserved_phased.vcf.gz
│   │   ├── YHPed1_conserved_phased.vcf.gz.csi
│   │   └── Genomics_Sex.tsv
│   ├── 06_spore_output/
│   │   ├── YHPed1_conserved_phased.vcf.gz-truffle.ibd.iqr.gz
│   │   ├── YHPed1_conserved_phased.vcf.gz--IBD0_threshold.txt
│   │   ├── SPORE_output.log
│   │   └── SPORE-Settings.R
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
| `*-truffle.ibd.iqr.gz` | IBD segment calls from SPORE/TRUFFLE (compressed) |
| `inbreeding_rankings.tsv` | Ranked relatedness scores (TSV for downstream plotting) |
| `SPORE-Settings.R` | Auto-generated settings file — saved for reproducibility |
| `execution_report.html` | Nextflow HTML report — runtime, resource usage, per-step logs |

---

## Repository Structure

```
.
├── bin/                            # Helper scripts called by Nextflow processes
│   ├── generate_spore_settings.py  # Writes complete SPORE-Settings.R from samples.yml
│   ├── txt_to_tsv.py               # Converts rank_ibd .txt output to .tsv
│   ├── yaml_to_linkage_map.py      # Extracts linkage_groups → TSV for PhaseParents_VCF.py
│   ├── yaml_to_popinfo.py          # Converts samples.yml → mod_popinfo.txt for mod_vcf2ahmm.py
│   └── yaml_to_sex_tsv.py          # Converts samples.yml → Genomics_Sex.tsv for SPORE
├── config/
│   ├── samples.yml                 # Main config: VCF path + all pipeline settings
│   └── test_samples.yml            # Config for bundled test run (-profile test)
├── data/
│   ├── raw/
│   │   └── YHPedigree1_FilteredSNVs_quarter.recode.vcf.gz   # Quarter VCF (via Git LFS)
│   └── test/
│       └── test_FilteredSNVs.recode.vcf.gz                  # Bundled test VCF (13 MB)
├── docs/
│   └── figures/                    # Pipeline figures for README
├── envs/
│   ├── finphaser.yml               # Python + ancestry_hmm + bcftools
│   └── spore.yml                   # R 4.3.3 + SPORE dependencies
├── main.nf                         # Pipeline entry point
├── nextflow.config                 # Profiles, per-process env assignments, resource limits
├── src/
│   ├── ancestry/
│   │   ├── PhaseParents_VCF.py     # Step 1 — phase parental genotypes (.vcf and .vcf.gz)
│   │   ├── mod_vcf2ahmm.py         # Step 3 — convert phased VCF to ancestry_hmm format
│   │   └── SNV_count.py            # Step 5 — summarise posterior probability files
│   ├── reporting/
│   │   └── rank_ibd.py             # Step 7 — rank IBD segments by relatedness
│   └── spore/
│       ├── SPORE.R                 # Step 6 — IBD detection (NoTrioCalculation=TRUE)
│       ├── scripts/
│       │   └── Mendel1.R           # Sourced internally by SPORE.R
│       ├── LICENSE
│       └── README.MD
├── tools/
│   └── truffle/
│       ├── truffle                 # TRUFFLE binary — Linux x86-64 only
│       └── LICENSE.txt
├── workflows/
│   ├── ancestry.nf                 # Nextflow subworkflow: steps 1–5
│   └── ibd.nf                      # Nextflow subworkflow: steps 6–7
└── results/                        # Pipeline outputs (not committed — generated at runtime)
```

---

## Troubleshooting

**TRUFFLE cannot run on macOS**
The `tools/truffle/truffle` binary is compiled for Linux x86-64. It will not run on macOS ARM (Apple Silicon) or macOS Intel. Run the full pipeline on a Linux server or HPC cluster.

**TRUFFLE binary not found at runtime**
SPORE resolves `truffle_path` relative to the Nextflow process working directory, not the repo root. Always use an absolute path in `config/samples.yml` or `config/test_samples.yml`:
```bash
# Get the absolute path from the repo root
realpath tools/truffle/truffle
```
Then set that value in your config:
```yaml
spore:
  truffle_path: "/absolute/path/to/FinPhaser/tools/truffle/truffle"
```

**`import allel` error**
`PhaseParents_VCF.py` requires `scikit-allel`, pinned in `envs/finphaser.yml`. Ensure Nextflow is using `-profile conda` so the correct per-process environment is activated.

**TRUFFLE binary not executable**
```bash
chmod +x tools/truffle/truffle
```

**Sample name matching in SPORE**
Sample names in `Genomics_Sex.tsv` are written exactly as they appear in `samples.yml`, with underscores preserved. This is required for SPORE to match names against VCF column headers.

**Switching to a different VCF**
Edit the `vcf:` key in `config/samples.yml`:
```yaml
vcf: "data/raw/my_other_file.vcf.gz"
```
No other files need to be changed.

**ancestry_hmm is slow**
ancestry_hmm is single-threaded and scales with sites × admixed samples. On a Mac with limited RAM it can take 2+ hours even on a quarter VCF. On a Linux server with 16+ GB RAM, the same run completes in 20–40 minutes. Use `-profile test` for a fast end-to-end check.

**Conda environments fail to build**
If conda environments were built on a different OS (e.g. macOS) and transferred to Linux, they will fail with `Exec format error`. Delete the cached environments and let Nextflow rebuild them:
```bash
rm -rf work/conda/
nextflow run main.nf -profile test
```

**Resuming a failed run**
```bash
nextflow run main.nf -profile conda -resume
```
Nextflow caches completed steps — only the failed step and anything downstream will re-run.

**Viewing step-level logs**
Check `results/pipeline_info/execution_report.html` for per-process status, duration, and exit codes. Raw logs are in `.nextflow/` (generated at runtime, not committed).

---

## Contact

For questions or issues, please open a GitHub Issue or contact the project maintainer.