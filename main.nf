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


// ── Import subworkflows ───────────────────────────────────────────────────────
include { ANCESTRY_WORKFLOW } from './workflows/ancestry'
include { IBD_WORKFLOW       } from './workflows/ibd'

// ── Main workflow ─────────────────────────────────────────────────────────────
workflow {

    // ── Startup log ───────────────────────────────────────────────────────────────
        log.info """
        ┌─────────────────────────────────────────────────────────┐
        │           F I N P H A S E R  v1.0.0                     │
        │  Local Ancestry Inference & IBD Detection Pipeline      │
        └─────────────────────────────────────────────────────────┘
        vcf          : ${params.vcf}
        samples_yml  : ${params.samples_yml}
        outdir       : ${params.outdir}
        profile      : ${workflow.profile}
        ──────────────────────────────────────────────────────────
        """.stripIndent()

    // ── Help message ──────────────────────────────────────────
    if (params.help) {
        log.info """
        ┌─────────────────────────────────────────────────────────┐
        │           F I N P H A S E R  v1.0.0                     │
        │  Local Ancestry Inference & IBD Detection Pipeline      │
        └─────────────────────────────────────────────────────────┘

        USAGE:
            nextflow run main.nf -profile conda [options]

        REQUIRED INPUTS:
            --vcf           Path to raw filtered VCF file
            --samples_yml   Path to sample/population config YAML

        OPTIONS:
            --outdir        Output directory           [default: results/]
            --help          Show this help message
        """.stripIndent()

        return
    }

    // Input channels — both files must exist before the pipeline starts
    ch_vcf         = Channel.fromPath(params.vcf,         checkIfExists: true)
    ch_samples_yml = Channel.fromPath(params.samples_yml, checkIfExists: true)

    // Phasing → HMM input prep → ancestry_hmm → summary report
    ANCESTRY_WORKFLOW(ch_vcf, ch_samples_yml)

    // SPORE admixture → IBD ranking
    IBD_WORKFLOW(
        ANCESTRY_WORKFLOW.out.phased_vcf,
        ch_samples_yml
    )
}
