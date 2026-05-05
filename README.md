FinPhaser
A Reproducible Pipeline for Genomic Inversion Analysis in Cichlids

Overview
FinPhaser is an automated bioinformatics pipeline designed to analyze the LG10 inversion in Aulonocara (Yellow Head) cichlids. By processing backcrosses with Mchenga conoforos, the pipeline identifies ancestral mosaicism through high-confidence SNP discovery, statistical phasing, and haplotype inference.

Repository Structure
Following the model of high-standard repositories like scikit-learn, FinPhaser is organized as follows:

bin/: Contains core analysis scripts, including the scikit-allel phasing script and the SPORE R-script.

config/: Nextflow profiles and software configurations.

data/: Bundled test dataset (subsampled BAMs for LG10) to exercise pipeline features.

docs/: Documentation, including VCF_generation.pdf and the project prospectus.  

results/: Automated output directory for VCFs, phased haplotypes, and LOD plots.

Pipeline Workflow
The pipeline executes the following stages in a single command:

Variant Calling: Implements GATK HaplotypeCaller and GenotypeGVCFs logic to produce a multisample VCF.  

Filtration: Standardizes SNP quality by filtering for biallelic sites and depth (DP) thresholds.  

Statistical Phasing: Uses scikit-allel to resolve long-range haplotypes across the chromosome.

Haplotype Inference: Runs Ancestry_HMM followed by SPORE for ancestral tracking and data visualization.

Getting Started
1. Environment Setup

Bash
conda env create -f environment.yml
conda activate finphaser_env
2. Running the Pipeline
To run the end-to-end analysis on your test data:

Bash
nextflow run main.nf --input 'data/*.bam' --ref 'data/reference.fasta'
Validation & Quality Control
Execution: The workflow is managed via Nextflow for full parallelization and dependency resolution.

Automated Checks: Integrated validation scripts confirm output correctness by checking row counts and checksums against known-result comparisons.