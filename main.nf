nextflow.enable.dsl=2

// Define input parameters
params.input_bam = "data/*.bam"
params.reference = "data/reference.fasta"
params.outdir = "results"

process VARIANT_CALLING {
    tag "GATK on $bam"
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    path bam
    path ref

    output:
    path "raw_variants.vcf.gz"

    script:
    """
    # Logic derived from VCF_generation.pdf
    gatk HaplotypeCaller -R $ref -I $bam -O raw_variants.vcf.gz
    bcftools view -v snps -m2 -M2 raw_variants.vcf.gz -Oz -o filtered_biallelic.vcf.gz
    """
}

process PHASE_VARIANTS {
    tag "Scikit-Allel Phasing"
    publishDir "${params.outdir}/phased", mode: 'copy'

    input:
    path vcf

    output:
    path "phased_variants.vcf.gz"

    script:
    """
    # Runs your custom scikit-allel phasing script
    python bin/phase_script.py --input $vcf --output phased_variants.vcf.gz
    """
}

process ANCESTRY_ANALYSIS {
    tag "Ancestry_HMM & SPORE"
    publishDir "${params.outdir}/final", mode: 'copy'

    input:
    path phased_vcf

    output:
    path "lod_plots.png"
    path "spore_results.csv"

    script:
    """
    # Ancestry_HMM followed by SPORE R-script
    Ancestry_HMM -i $phased_vcf -o ancestry_output.txt
    Rscript bin/spore_analysis.R ancestry_output.txt
    """
}

workflow {
    vcf_ch = VARIANT_CALLING(params.input_bam, params.reference)
    phased_ch = PHASE_VARIANTS(vcf_ch)
    ANCESTRY_ANALYSIS(phased_ch)
}