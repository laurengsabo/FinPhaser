/*
 * ============================================================
 *  workflows/ancestry.nf
 * ============================================================
 *  Ancestry inference subworkflow.
 *
 *  Step 1 — PHASE_PARENTS
 *    Script : src/ancestry/PhaseParents_VCF.py
 *    In     : raw filtered VCF
 *    Out    : YHPed1_conserved_phased.vcf
 *
 *  Step 2+3 — PREPARE_HMM_INPUT  (combined)
 *    Script : src/ancestry/mod_vcf2ahmm.py --pop-yaml
 *    In     : phased VCF + samples.yml (pop assignments read directly)
 *    Out    : ancestry_input.txt, ahmm.ploidy
 *
 *  Step 4 — RUN_ANCESTRY_HMM
 *    Tool   : ancestry_hmm (CLI)
 *    In     : ancestry_input.txt, ahmm.ploidy
 *    Out    : *.posterior files
 *
 *  Step 5 — ANCESTRY_SUMMARY
 *    Script : src/ancestry/SNV_count.py
 *    In     : directory of *.posterior files
 *    Out    : ancestry_summary_report/ directory
 * ============================================================
 */

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS: PHASE_PARENTS
// ─────────────────────────────────────────────────────────────────────────────
process PHASE_PARENTS {
    tag "PhaseParents"

    // Publish the phased VCF so users can inspect it
    publishDir "${params.outdir}/ancestry/01_phased", mode: 'copy'

    input:
    path vcf            // raw filtered VCF

    output:
    path "YHPed1_conserved_phased.vcf", emit: phased_vcf

    script:
    /*
     * PhaseParents_VCF.py (modified from Dr. McGrath's PhaseParents_Scikit.py)
     * Identifies conserved informative SNVs and masks parental genotypes
     * to produce a phased VCF.
     */
    """
    python ${projectDir}/src/ancestry/PhaseParents_VCF.py ${vcf}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS: PREPARE_HMM_INPUT  (steps 2 + 3 combined)
// ─────────────────────────────────────────────────────────────────────────────
process PREPARE_HMM_INPUT {
    tag "PrepareHMMInput"

    publishDir "${params.outdir}/ancestry/02_hmm_input", mode: 'copy'

    input:
    path phased_vcf     // conserved phased VCF from step 1
    path samples_yml    // user-supplied YAML with pop assignments

    output:
    path "ancestry_input.txt", emit: ancestry_input
    path "ahmm.ploidy",        emit: ploidy_file

    script:
    /*
     * mod_vcf2ahmm.py --pop-yaml reads population assignments (0/1/2/3/-1)
     * directly from samples.yml, building the internal popinfo structure and
     * then producing both ancestry_hmm input files in one pass.
     *
     * This combines the original ancestry_input.py (step 2) and
     * mod_vcf2ahmm.py (step 3) into a single invocation.
     */
    """
    python ${projectDir}/src/ancestry/mod_vcf2ahmm.py \\
        --vcf        ${phased_vcf}     \\
        --pop-yaml   ${samples_yml}    \\
        --out-input  ancestry_input.txt \\
        --out-ploidy ahmm.ploidy
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS: RUN_ANCESTRY_HMM
// ─────────────────────────────────────────────────────────────────────────────
process RUN_ANCESTRY_HMM {
    tag "AncestryHMM"

    publishDir "${params.outdir}/ancestry/03_posteriors", mode: 'copy'

    // Give this process more CPUs and memory — it is the most compute-heavy step
    cpus   4
    memory '16 GB'
    time   '6 h'

    input:
    path ancestry_input
    path ploidy_file

    output:
    // Capture every .posterior file produced (one per admixed sample)
    path "*.posterior", emit: posterior_files

    script:
    /*
     * ancestry_hmm parameters:
     *   -a 4 0.25 0.25 0.25 0.25  — 4 ancestral populations, equal priors
     *   -p 0..3 -2 0.25           — each pulse ~2 generations ago, equal weight
     *
     * These defaults come from the samples.yml [hmm] block; override them
     * there rather than editing this file.
     */
    """
    ancestry_hmm \\
        -i ${ancestry_input} \\
        -s ${ploidy_file}    \\
        -a 4 0.25 0.25 0.25 0.25 \\
        -p 0 -2 0.25 \\
        -p 1 -2 0.25 \\
        -p 2 -2 0.25 \\
        -p 3 -2 0.25
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS: ANCESTRY_SUMMARY
// ─────────────────────────────────────────────────────────────────────────────
process ANCESTRY_SUMMARY {
    tag "AncestrySummary"

    publishDir "${params.outdir}/ancestry", mode: 'copy'

    input:
    // .collect() gathers all individual .posterior paths into one list,
    // staging them all into the same working directory before the script runs.
    path posterior_files

    output:
    path "ancestry_summary_report/", emit: summary_dir

    script:
    /*
     * SNV_count.py reads every .posterior file in the current working directory
     * and writes per-sample calls_<sample>.csv files plus a master summary CSV.
     */
    """
    mkdir -p ancestry_summary_report
    python ${projectDir}/src/ancestry/SNV_count.py \\
        --input-dir  . \\
        --output-dir ancestry_summary_report
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  SUBWORKFLOW: ANCESTRY_WORKFLOW
// ─────────────────────────────────────────────────────────────────────────────
workflow ANCESTRY_WORKFLOW {

    take:
    ch_vcf          // Channel<Path> — raw VCF
    ch_samples_yml  // Channel<Path> — user YAML

    main:

    // Step 1
    PHASE_PARENTS(ch_vcf)

    // Steps 2 + 3
    PREPARE_HMM_INPUT(
        PHASE_PARENTS.out.phased_vcf,
        ch_samples_yml
    )

    // Step 4
    RUN_ANCESTRY_HMM(
        PREPARE_HMM_INPUT.out.ancestry_input,
        PREPARE_HMM_INPUT.out.ploidy_file
    )

    // Step 5 — collect all .posterior files before passing to summary
    ANCESTRY_SUMMARY(
        RUN_ANCESTRY_HMM.out.posterior_files.collect()
    )

    emit:
    phased_vcf      = PHASE_PARENTS.out.phased_vcf
    posterior_files = RUN_ANCESTRY_HMM.out.posterior_files
    summary_dir     = ANCESTRY_SUMMARY.out.summary_dir
}
