# FinPhaser

<p align="center">
  <b>A reproducible framework for local ancestry inference and IBD detection in hybrid populations</b>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/python-3.9+-blue.svg" />
  <img src="https://img.shields.io/badge/R-4.0+-blue.svg" />
  <img src="https://img.shields.io/badge/conda-environment-green.svg" />
  <img src="https://img.shields.io/badge/status-active-success.svg" />
  <img src="https://img.shields.io/badge/license-MIT-lightgrey.svg" />
</p>

---

## Overview

**FinPhaser** is a bioinformatics pipeline for estimating **local ancestry** and **Identity-By-Descent (IBD)** in hybrid populations, with a focus on African cichlid genomics.

The framework integrates:
- Variant filtering and phasing
- Hidden Markov Models for ancestry inference
- Admixture analysis using SPORE
- IBD detection and ranking

It is particularly suited for studying **complex sex determination systems**, including ZW/XY dynamics on the LG10 linkage group in *Aulonocara* 'Yellow Head'.

---

## Quick Start

```bash id="qs1x2a"
# Clone repository
git clone https://github.com/your-username/FinPhaser.git
cd FinPhaser

# Create environment
conda env create -f environment.yml
conda activate FinPhaser

# Run core pipeline (example)
python src/ancestry/PhaseParents_VCF.py data/raw/YHPedigree1_FilteredSNVs.recode.vcf
ancestry_hmm -i ancestry_input.txt -s ahmm.ploidy -a 4 0.25 0.25 0.25 0.25
Rscript src/SPORE.R config/SPORE-Settings.R
```

---

## Features
- Reproducible pipeline with unified environment
- Local ancestry inference via HMM
- IBD detection and ranking
- Custom SNV filtering and phasing
- Modular design for easy extension

---

## Repository Structure

```
.
├── config/             # Parameter files (e.g., SPORE-Settings.R)
├── data/
│   ├── raw/            # Input VCF files
│   └── processed/      # Phased and filtered outputs
├── docs/               # Architecture diagrams and notes
├── results/            # Final outputs and logs
├── src/
│   ├── ancestry/       # Phasing + HMM preprocessing
│   ├── reporting/      # IBD ranking + summaries
│   └── SPORE.R         # Admixture analysis
└── environment.yml     # Conda environment
```

---

## Pipeline Workflow
**Phasing**
- Input: FilteredSNVs.vcf
- Output: conserved_phased.vcf
*Identifies conserved informative SNVs*

**Ancestry Inference (HMM)**
- Input: ancestry_input.txt
- Output: .posterior
*Computes local ancestry probabilities*

**Admixture Analysis (SPORE)**
- Input: Genomics_Sex.tsv
- Output: truffle.ibd
*Detects IBD segments*

**IBD Ranking**
- Input: truffle.ibd
- Output: inbreeding_rankings.txt
*Produces ranked relatedness scores*

---

## Example Output
.posterior — Local ancestry probabilities per site

truffle.ibd — IBD segment calls

inbreeding_rankings.txt — Ranked relatedness metrics

---

## Troubleshooting
*Sample Name Formatting*
SPORE requires sample IDs without underscores.

*Unexpected LG10 Results*
Verify HMM pulse parameters (-p) match expected admixture history.

*Dependency Issues*
Ensure the Conda environment is active before running any scripts.

---

## Contributing
Contributions are welcome. Please open an issue to discuss proposed changes or submit a pull request.

---

## License
This project is licensed under the MIT License.

---

## Contact

Lauren Sabo

Bioinformatics MS Student, McGrath Lab

Email: lsabo8@gatech.edu

GitHub: https://github.com/laurengsabo