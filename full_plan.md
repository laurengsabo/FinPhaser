Ancestry HMM Pipeline: From VCF to Posterior Probabilities
This pipeline processes filtered SNVs to estimate local ancestry using a Hidden Markov Model (HMM). It involves phasing parental genotypes, formatting genomic data for ancestry analysis, and running ancestry_hmm.

Pipeline Overview
Parental Phasing: Generate a phased VCF with masked genotypes.

Input Preparation: Generate sample metadata and ancestry-specific input formats.

HMM Formatting: Convert VCF and metadata into ancestry_hmm compatible files.

Ancestry Estimation: Run the HMM to produce posterior probabilities for each individual.

1. Parental Phasing and Masking
Input: YHPedigree1_FilteredSNVs.recode.vcf

Script: PhaseParents_VCF.py (Modified from Dr. McGrath’s PhaseParents_Scikit.py)

The original script was adapted to work directly with VCF formats. This step identifies conserved informative sites and produces a phased VCF file.

Output: YHPed1_conserved_phased.vcf

2. Ancestry Metadata Generation
Script: ancestry_input.py

This script processes the sample metadata to create the required input structure for downstream conversion.

Output: YHPed1_conv_input.txt

3. VCF to HMM Input Conversion
Input: YHPed1_conv_input.txt (in phased_vcf2ahmm) OR mod_popinfo.txt (in AC_mod_vcfahmm), YHPed1_conserved_phased.vcf

Script: mod_vcf2ahmm.py (Modified)

Note: The vcf2ahmm.py script was specifically modified to allow the inclusion of phased genotypes, ensuring the HMM utilizes the linkage information established in Step 1.

Outputs: ancestry_input.txt (Genotype likelihoods/counts) and ahmm.ploidy (Ploidy information for each individual)

4. Running Ancestry_HMM
The final estimation was performed using the following command:

Bash
ancestry_hmm -i ancestry_input.txt -s ahmm.ploidy \
 -a 4 0.25 0.25 0.25 0.25 \
 -p 0 -2 0.25 \
 -p 1 -2 0.25 \
 -p 2 -2 0.25 \
 -p 3 -2 0.25

Parameter Justification
-a 2 0.25 0.25 0.25 0.25: Specifies 4 ancestral populations (haplotypes) with an initial starting admixture proportion of 25% for each haplotype. This assumes an even contribution from both parental sources and their haplotypes.

-p 0 -2 0.25 (Pulse 0): Represents an ancient admixture event (~2 generations ago) (Paternal Haplotype 1).
-p 1 -2 0.25 (Pulse 1): Represents an ancient admixture event (~2 generations ago) (Paternal Haplotype 2).
-p 2 -2 0.25 (Pulse 2): Represents an ancient admixture event (~2 generations ago) (Maternal Haplotype 1).
-p 3 -2 0.25 (Pulse 3): Represents an ancient admixture event (~2 generations ago) (Maternal Haplotype 2).

Results
The pipeline generates .posterior files for samples YH_016 through YH_042. These files contain the probabilities of each ancestry state at every site across the genome.

5. Creating an AHMM Summary Report
Input: dir path to AHMM posterior files

Script: SNV_count.py

Note: This script allows for the AHMM output to become readable, can be used for future visualization (might need tweaks)

Outputs: ancestry_summary_report/ + calls_[sample].csv + master_ancestry_summary.csv

6. Running SPORE
Input: modify the SPORE-Settings.R file with the pathnames and desired parameters + Genomics_Sex.tsv

Bash
Rscript ../SPORE.R ../SPORE-Settings.R > SPORE_output.log 2>&1

Note: SPORE doesn't allow underscores in sample names, so the pipeline will need to check for this after the 

Outputs: brood1_final.vcf.gz-truffle.ibd + more

7. IBD Rankings Report
Input: brood1_final.vcf.gz-truffle.ibd

Script: rank_ibd.py

Note: The output is txt, so in order to use it for future visualization, the type will need to be ._sv

Outputs: inbreeding_rankings.txt

8. Visualizations