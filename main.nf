#!/usr/bin/env nextflow
/*
 * ============================================================
 *  FinPhaser — Main Pipeline Entry Point
 * ============================================================
 *  Local ancestry inference and IBD detection in hybrid populations.
 *
 *  ONE input required:
 *    config/samples.yml — contains the vcf: path and all settings.
 *
 *  Usage:
 *    nextflow run main.nf -profile conda
 *    nextflow run main.nf -profile conda \
 *        --samples_yml config/samples.yml
 * ============================================================
 */

nextflow.enable.dsl = 2

// ── Import subworkflows ───────────────────────────────────────────────────────
include { ANCESTRY_WORKFLOW } from './workflows/ancestry'
include { IBD_WORKFLOW       } from './workflows/ibd'

// ── Main workflow ─────────────────────────────────────────────────────────────
workflow {

    // ── Read VCF path from samples.yml ───────────────────────────────────────
    // The vcf: key in samples.yml is the single source of truth for the input
    // file. This means the full run is reproducible from samples.yml alone —
    // no CLI flags or config edits needed when switching VCF files.
    def yml_text  = new File(params.samples_yml).text
    def vcf_match = yml_text =~ /(?m)^vcf:\s*["']?([^"'\s#\r\n]+)["']?/
    if (!vcf_match) {
        error "ERROR: 'vcf:' key not found in ${params.samples_yml}.\n" +
              "Add a line like:  vcf: \"data/raw/your_file.vcf\""
    }
    def vcf_path = vcf_match[0][1]

    // ── Startup log ──────────────────────────────────────────────────────────
    log.info """
    ┌─────────────────────────────────────────────────────────┐
    │           F I N P H A S E R  v1.0.0                    │
    │  Local Ancestry Inference & IBD Detection Pipeline      │
    └─────────────────────────────────────────────────────────┘
    vcf          : ${vcf_path}  (from samples.yml)
    samples_yml  : ${params.samples_yml}
    outdir       : ${params.outdir}
    profile      : ${workflow.profile}
    ──────────────────────────────────────────────────────────
    """.stripIndent()

    // ── Help message ─────────────────────────────────────────────────────────
    if (params.help) {
        log.info """
        USAGE:
            nextflow run main.nf -profile conda [options]

        REQUIRED INPUT:
            --samples_yml   Path to sample config YAML  [default: config/samples.yml]
                            The VCF path is read from the vcf: key inside this file.

        OPTIONS:
            --outdir        Output directory  [default: results/]
            --help          Show this message
        """.stripIndent()
        return
    }

    // ── Input channels ────────────────────────────────────────────────────────
    ch_samples_yml = Channel.fromPath(params.samples_yml, checkIfExists: true)
    ch_vcf         = Channel.fromPath(vcf_path,           checkIfExists: true)

    // ── Phasing → HMM input prep → ancestry_hmm → summary ───────────────────
    ANCESTRY_WORKFLOW(ch_vcf, ch_samples_yml)

    // ── SPORE admixture → IBD ranking ─────────────────────────────────────────
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
