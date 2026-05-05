#!/usr/bin/env nextflow
/*
 * ============================================================
 *  FinPhaser — Main Pipeline Entry Point
 * ============================================================
 *  Local ancestry inference and IBD detection in hybrid populations.
 *
 *  Two inputs required:
 *    1. A raw filtered VCF  (e.g. YHPedigree1_FilteredSNVs.recode.vcf)
 *    2. A samples YAML      (e.g. config/samples.yml)
 *
 *  Usage:
 *    nextflow run main.nf -profile conda
 *    nextflow run main.nf -profile conda \
 *        --vcf data/raw/YHPedigree1_FilteredSNVs.recode.vcf \
 *        --samples_yml config/samples.yml
 * ============================================================
 */

nextflow.enable.dsl = 2

// ── Help message ──────────────────────────────────────────────────────────────
if (params.help) {
    log.info """
    ┌─────────────────────────────────────────────────────────┐
    │           F I N P H A S E R  v1.0.0                    │
    │  Local Ancestry Inference & IBD Detection Pipeline      │
    └─────────────────────────────────────────────────────────┘

    USAGE:
        nextflow run main.nf -profile conda [options]

    REQUIRED INPUTS (set in nextflow.config or pass via CLI):
        --vcf           Path to raw filtered VCF file
        --samples_yml   Path to sample/population config YAML

    OPTIONS:
        --outdir        Output directory           [default: results/]
        --help          Show this help message

    PROFILES:
        conda           Local conda environment    [recommended]
        mamba           Local mamba (faster solves)
        docker          Docker container
        singularity     Singularity container (HPC)
        test            Run bundled test dataset

    EXAMPLES:
        # Standard run
        nextflow run main.nf -profile conda

        # Custom inputs
        nextflow run main.nf -profile conda \\
            --vcf data/raw/MyFile.vcf \\
            --samples_yml config/samples.yml

        # HPC with Singularity
        nextflow run main.nf -profile singularity

        # Quick test
        nextflow run main.nf -profile test
    """.stripIndent()
    exit 0
}

// ── Import subworkflows ───────────────────────────────────────────────────────
include { ANCESTRY_WORKFLOW } from './workflows/ancestry'
include { IBD_WORKFLOW       } from './workflows/ibd'

// ── Startup log ───────────────────────────────────────────────────────────────
log.info """
┌─────────────────────────────────────────────────────────┐
│           F I N P H A S E R  v1.0.0                    │
│  Local Ancestry Inference & IBD Detection Pipeline      │
└─────────────────────────────────────────────────────────┘
vcf          : ${params.vcf}
samples_yml  : ${params.samples_yml}
outdir       : ${params.outdir}
profile      : ${workflow.profile}
──────────────────────────────────────────────────────────
""".stripIndent()

// ── Main workflow ─────────────────────────────────────────────────────────────
workflow {

    // Input channels — both files must exist before the pipeline starts
    ch_vcf         = Channel.fromPath(params.vcf,         checkIfExists: true)
    ch_samples_yml = Channel.fromPath(params.samples_yml, checkIfExists: true)

    // Steps 1–5: Phasing → HMM input prep → ancestry_hmm → summary report
    ANCESTRY_WORKFLOW(ch_vcf, ch_samples_yml)

    // Steps 6–7: SPORE admixture → IBD ranking
    IBD_WORKFLOW(
        ANCESTRY_WORKFLOW.out.phased_vcf,
        ch_samples_yml
    )
}

// ── Completion handler ────────────────────────────────────────────────────────
workflow.onComplete {
    def status = workflow.success ? "SUCCESS ✔" : "FAILED ✘"
    log.info """
    ──────────────────────────────────────────────────────────
    Pipeline complete!
    Status   : ${status}
    Duration : ${workflow.duration}
    Outputs  : ${params.outdir}/
    Report   : ${params.tracedir}/execution_report.html
    ──────────────────────────────────────────────────────────
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline failed: ${workflow.errorReport}"
}
