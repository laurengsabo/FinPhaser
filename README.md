# FinPhaser
<p align="center">
  <b>A reproducible framework for local ancestry inference and IBD detection in hybrid populations</b>
</p>
<p align="center">
  <img src="https://img.shields.io/badge/python-3.10-blue.svg" />
  <img src="https://img.shields.io/badge/R-multi--version-blue.svg" />
  <img src="https://img.shields.io/badge/nextflow-24.04-brightgreen.svg" />
  <img src="https://img.shields.io/badge/conda-environment-green.svg" />
  <img src="https://img.shields.io/badge/status-active-success.svg" />
  <img src="https://img.shields.io/badge/license-MIT-lightgrey.svg" />
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

The entire pipeline runs in **a single command** via [Nextflow](https://www.nextflow.io/), requiring only two input files: a filtered VCF and a samples YAML.

---

## Quick Start

### Prerequisites

| Tool | Minimum version | Install |
|------|----------------|---------|
| Java | 11+ (Java 17 recommended) | `brew install openjdk` / system package manager |
| Nextflow | 23.04+ | See below |
| Conda or Mamba | any recent | [Miniforge](https://github.com/conda-forge/miniforge) recommended |

**Install Nextflow** (if not already installed):
```bash
# Option A — direct install
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Option B — via conda (already included in environment.yml)
conda install -c conda-forge nextflow
```

---

### 1. Clone the repository

```bash
git clone https://github.com/laurengsabo/FinPhaser.git
cd FinPhaser
```

---

### 2. Place your input files

```bash
# Your starting VCF — the only genomic input needed
mkdir -p data/raw
cp /path/to/YHPedigree1_FilteredSNVs.recode.vcf data/raw/
```

---

### 3. Build conda environments

FinPhaser uses **process-specific Conda environments** for full reproducibility.

- `envs/finphaser.yml` — modern Python + bioinformatics tools
- `envs/spore.yml` — legacy R environment required by SPORE

You **do not need to manually create or activate environments**.

Nextflow will automatically:
- create environments as needed
- activate the correct environment per step
- cache them for reuse

---

### 4. Configure your samples

All pipeline configuration lives in `config/samples.yml`. The two things to set for your data:

**a) Population assignments** — each diploid parent appears **twice** (once per haplotype it contributes). Offspring are marked `admixed`.

```yaml
populations:
  # Female parent — contributes maternal haplotypes 0 and 1
  - sample: "YH_006_f"
    population: 0
    sex: F
  - sample: "YH_006_f"
    population: 1
    sex: F

  # Male parent — contributes paternal haplotypes 2 and 3
  - sample: "YH_011_m"
    population: 2
    sex: M
  - sample: "YH_011_m"
    population: 3
    sex: M

  # Admixed offspring — ancestry inferred by HMM
  - { sample: "YH_016", population: admixed, sex: F }
  - { sample: "YH_017", population: admixed, sex: M }
```

**b) TRUFFLE path** — point to the binary placed in `tools/truffle/`:

```yaml
spore:
  truffle_path: "tools/truffle/truffle"
```

Sex metadata (`Genomics_Sex.tsv`) is **generated automatically** at runtime from the `sex` field in `samples.yml` — you do not need to maintain a separate file.

---

### 5. Run the pipeline

```bash
nextflow run main.nf -profile conda
```

Nextflow will automatically handle all software dependencies using Conda.

**Resume after a failure** (completed steps are not re-run):
```bash
nextflow run main.nf -profile conda -resume
```

**Quick test on bundled data:**
```bash
nextflow run main.nf -profile test
```

---

## Execution Profiles

| Profile | Description |
|---------|-------------|
| `conda` | Local conda environment — recommended default |
| `mamba` | Same as conda but uses mamba for faster solves |
| `docker` | Docker container (build image first — see below) |
| `singularity` | Singularity container — best for HPC clusters |
| `test` | Runs bundled test dataset; completes in minutes |

## Reproducibility

FinPhaser uses **process-level environment isolation** via Nextflow:

- Each pipeline step runs in its own Conda environment
- Legacy tools (e.g., SPORE with R 3.6) are isolated from modern dependencies
- All environments are fully version-pinned in `envs/`

This design avoids dependency conflicts and ensures consistent results across systems.

### Docker

```bash
docker build -t finphaser:latest .
nextflow run main.nf -profile docker
```

### HPC / Singularity

```bash
singularity build finphaser.sif docker://finphaser:latest
nextflow run main.nf -profile singularity
```

---

## Pipeline Workflow

One starting VCF feeds both the ancestry and IBD branches after phasing:

```
YHPedigree1_FilteredSNVs.recode.vcf
config/samples.yml
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
│    input.txt    │  │  → *-truffle.ibd          │
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
- Steps 2 and 3 are combined: `bin/yaml_to_popinfo.py` converts `samples.yml` into the tab-separated population info file that `mod_vcf2ahmm.py` reads via its `-s` flag. `mod_vcf2ahmm.py` is **not modified**.
- Sex metadata is generated at runtime by `bin/yaml_to_sex_tsv.py` — no separate `Genomics_Sex.tsv` file needed as input.
- Underscores are stripped from sample names automatically before SPORE runs.
- SPORE-Settings.R is fully auto-generated from `samples.yml` — no manual editing required.

---

## Inputs

Only two files are required to start the pipeline:

| File | Description |
|------|-------------|
| `data/raw/YHPedigree1_FilteredSNVs.recode.vcf` | Raw filtered SNV VCF — the only genomic input |
| `config/samples.yml` | All configuration: population assignments, sex, HMM params, SPORE settings |

### samples.yml structure

```yaml
# HMM parameters — forwarded to mod_vcf2ahmm.py
hmm:
  recombination_rate: 1.0e-8
  min_distance_bp: 1000
  min_allele_freq_diff: 0.1
  use_genotypes: 0

# SPORE / TRUFFLE settings
spore:
  truffle_path: "tools/truffle/truffle"
  truffle_maf: 0.05
  truffle_missing: 0.95
  APO: 24
  output_label: "YHPed1_brood1"
  max_memory: "4G"
  max_cores: 4
  plots: TRUE
  # ... (see config/samples.yml for full parameter list)

# Sample assignments
# population → used by AHMM only
# sex        → used by SPORE only
populations:
  - sample: "YH_006_f"
    population: 0          # maternal haplotype 1
    sex: F
  - sample: "YH_006_f"
    population: 1          # maternal haplotype 2
    sex: F
  - sample: "YH_011_m"
    population: 2          # paternal haplotype 1
    sex: M
  - sample: "YH_011_m"
    population: 3          # paternal haplotype 2
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

> Parents appear **twice** in the list — once per haplotype — because each diploid parent contributes two distinct haplotypes to the 4-population HMM model.

---

## Outputs

All outputs land in `results/`:

```
results/
├── ancestry/
│   ├── 01_phased/
│   │   └── YHPed1_conserved_phased.vcf
│   ├── 02_hmm_input/
│   │   ├── ancestry_input.txt
│   │   ├── ahmm.ploidy
│   │   └── mod_popinfo.txt          # generated from samples.yml; saved for inspection
│   ├── 03_posteriors/
│   │   └── YH_016.posterior ... YH_042.posterior
│   └── ancestry_summary_report/
│       ├── calls_YH_016.posterior.csv
│       ├── ...
│       └── master_ancestry_summary.csv
├── ibd/
│   ├── 06_spore_input/
│   │   ├── YHPed1_conserved_phased.vcf.gz     # bgzipped phased VCF
│   │   ├── YHPed1_conserved_phased.vcf.gz.csi # CSI index
│   │   └── Genomics_Sex.tsv                   # auto-generated from samples.yml
│   ├── 06_spore_output/
│   │   ├── YHPed1_conserved_phased.vcf.gz-truffle.ibd
│   │   ├── SPORE_output.log
│   │   └── SPORE-Settings.R         # saved for full reproducibility
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
| `*-truffle.ibd` | IBD segment calls from SPORE/TRUFFLE |
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
│   ├── yaml_to_popinfo.py          # Converts samples.yml → mod_popinfo.txt for mod_vcf2ahmm.py
│   └── yaml_to_sex_tsv.py          # Converts samples.yml → Genomics_Sex.tsv for SPORE
├── config/
│   └── samples.yml                 # All user config: populations, sex, HMM + SPORE params
├── data/
│   ├── raw/
│   │   └── YHPedigree1_FilteredSNVs.recode.vcf   # Starting VCF (only genomic input)
│   └── test/                       # Bundled test dataset
├── docs/                           # Architecture diagrams and notes
├── envs/
│   ├── finphaser.yml               # finphaser conda init
│   └── spore.yml                   # spore conda init
├── main.nf                         # Pipeline entry point
├── nextflow.config                 # Profiles, resource limits, reporting
├── src/
│   ├── ancestry/
│   │   ├── PhaseParents_VCF.py     # Step 1 — phase parental genotypes
│   │   ├── mod_vcf2ahmm.py         # Step 3 — convert phased VCF to ancestry_hmm format
│   │   └── SNV_count.py            # Step 5 — summarise posterior probability files
│   ├── reporting/
│   │   └── rank_ibd.py             # Step 7 — rank IBD segments by relatedness
│   └── spore/
│       ├── SPORE.R                 # Step 6 — IBD detection (SPORE + TRUFFLE)
│       ├── scripts/
│       │   └── Mendel1.R           # Sourced internally by SPORE.R
│       ├── LICENSE
│       └── README.MD
├── tools/
│   └── truffle/
│       ├── truffle                 # TRUFFLE binary (not available via conda)
│       └── LICENSE.txt
├── workflows/
│   ├── ancestry.nf                 # Nextflow subworkflow: steps 1–5
│   └── ibd.nf                      # Nextflow subworkflow: steps 6–7
└── results/                        # Pipeline outputs (git-ignored)
```

---

## Troubleshooting

**`import allel` error**
`PhaseParents_VCF.py` requires `scikit-allel`. Ensure the conda environment was built from the current `environment.yml` (`scikit-allel=1.3.7` is pinned there) and is active.

**TRUFFLE binary not found**
Ensure `tools/truffle/truffle` exists and is executable:
```bash
chmod +x tools/truffle/truffle
```
The path in `config/samples.yml` under `spore.truffle_path` must match exactly.

**Sample name formatting (SPORE)**
SPORE cannot handle underscores in sample IDs. FinPhaser strips them automatically when generating `Genomics_Sex.tsv` — no manual renaming needed in your VCF or `samples.yml`.

**`data/raw/Genomics_Sex.tsv` not needed**
If you have this file from a previous manual run, you can safely remove it. Sex metadata is now read from `samples.yml` and the TSV is generated automatically at runtime.

**Unexpected LG10 results**
Verify the HMM pulse parameters in `config/samples.yml` match your crossing design. The default (`-p 0..3 -2 0.25`) assumes ~2 generations of admixture.

**Resuming a failed run**
```bash
nextflow run main.nf -profile conda -resume
```

**Viewing step-level logs**
Check `results/pipeline_info/execution_report.html` for per-process status, duration, and exit codes. Raw logs are in `.nextflow/`.

---

## Contact

For questions or issues, please open a GitHub Issue or contact the project maintainer.