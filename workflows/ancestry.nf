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
 *    Helper : bin/yaml_to_popinfo.py
 *      Reads samples.yml → writes mod_popinfo.txt in the exact
 *      tab-separated format that mod_vcf2ahmm.py expects
 *      (sample<TAB>pop, one row per haplotype, "admixed" literal).
 *    Script : src/ancestry/mod_vcf2ahmm.py
 *      Reads the generated mod_popinfo.txt and the phased VCF.
 *      All tunable parameters (recombination rate, min distance,
 *      allele freq diff threshold) are pulled from samples.yml
 *      and forwarded as CLI flags.
 *    In     : phased VCF + samples.yml
 *    Out    : ancestry_input.txt, ahmm.ploidy
 *
 *  Step 4 — RUN_ANCESTRY_HMM
 *    Tool   : ancestry_hmm (CLI)
 *    In     : ancestry_input.txt, ahmm.ploidy
 *    Out    : *.posterior files (one per admixed sample)
 *
 *  Step 5 — ANCESTRY_SUMMARY
 *    Script : src/ancestry/SNV_count.py
 *    In     : directory of *.posterior files
 *    Out    : ancestry_summary_report/
 * ============================================================
 */

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS: PHASE_PARENTS
// ─────────────────────────────────────────────────────────────────────────────
process PHASE_PARENTS {
    tag "PhaseParents"
    publishDir "${params.outdir}/ancestry/01_phased", mode: 'copy'

    input:
    path vcf

    output:
    path "YHPed1_conserved_phased.vcf", emit: phased_vcf

    script:
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
    path phased_vcf
    path samples_yml

    output:
    path "ancestry_input.txt", emit: ancestry_input
    path "ahmm.ploidy",        emit: ploidy_file
    path "mod_popinfo.txt",    emit: popinfo       // saved for reproducibility / inspection

    script:
    /*
     * Step 2 (combined into one process):
     *   bin/yaml_to_popinfo.py converts the populations block of
     *   samples.yml into the tab-separated mod_popinfo.txt that
     *   mod_vcf2ahmm.py expects:
     *       YH_006_f<TAB>0
     *       YH_006_f<TAB>1
     *       YH_011_m<TAB>2
     *       YH_011_m<TAB>3
     *       YH_016<TAB>admixed
     *       ...
     *
     * Step 3:
     *   mod_vcf2ahmm.py uses its original CLI flags:
     *     -v  phased VCF
     *     -s  mod_popinfo.txt (generated above)
     *     -g  use_genotypes (0 = read counts, 1 = genotypes)
     *     -r  recombination rate in Morgans/bp
     *     -m  min distance between SNPs in bp
     *     --min_diff  min allele frequency difference between haplotypes
     */
    """
    # Convert samples.yml → mod_popinfo.txt
    python ${projectDir}/bin/yaml_to_popinfo.py \
        --samples-yml ${samples_yml} \
        --out         mod_popinfo.txt

    # Convert phased VCF + popinfo → ancestry_hmm input files
    python ${projectDir}/src/ancestry/mod_vcf2ahmm.py \
        -v       ${phased_vcf}      \
        -s       mod_popinfo.txt    \
        -o_txt   ancestry_input.txt \
        -o_ploidy ahmm.ploidy       \
        -g       ${params.hmm_use_genotypes}     \
        -r       ${params.hmm_recombination_rate} \
        -m       ${params.hmm_min_distance_bp}    \
        --min_diff ${params.hmm_min_allele_freq_diff}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS: RUN_ANCESTRY_HMM
// ─────────────────────────────────────────────────────────────────────────────
process RUN_ANCESTRY_HMM {
    tag "AncestryHMM"
    publishDir "${params.outdir}/ancestry/03_posteriors", mode: 'copy'

    cpus   4
    memory '16 GB'
    time   '6 h'

    input:
    path ancestry_input
    path ploidy_file

    output:
    path "*.posterior", emit: posterior_files

    script:
    /*
     * 4-population model:
     *   -a 4 0.25 0.25 0.25 0.25
     *       Four ancestral haplotypes; equal prior on each.
     *   -p 0..3  -2  0.25
     *       Each haplotype pulse placed ~2 generations back with
     *       equal initial proportion. The negative sign for
     *       generations is ancestry_hmm convention for discrete
     *       admixture pulses (not continuous migration).
     */
    """
    ancestry_hmm \
        -i ${ancestry_input} \
        -s ${ploidy_file}    \
        -a 4 0.25 0.25 0.25 0.25 \
        -p 0 -2 0.25 \
        -p 1 -2 0.25 \
        -p 2 -2 0.25 \
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
    path posterior_files    // .collect() gathers all into one staging dir

    output:
    path "ancestry_summary_report/", emit: summary_dir

    script:
    """
    mkdir -p ancestry_summary_report
    python ${projectDir}/src/ancestry/SNV_count.py \
        --input-dir  . \
        --output-dir ancestry_summary_report
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  SUBWORKFLOW
// ─────────────────────────────────────────────────────────────────────────────
workflow ANCESTRY_WORKFLOW {

    take:
    ch_vcf
    ch_samples_yml

    main:
    PHASE_PARENTS(ch_vcf)

    PREPARE_HMM_INPUT(
        PHASE_PARENTS.out.phased_vcf,
        ch_samples_yml
    )

    RUN_ANCESTRY_HMM(
        PREPARE_HMM_INPUT.out.ancestry_input,
        PREPARE_HMM_INPUT.out.ploidy_file
    )

    ANCESTRY_SUMMARY(
        RUN_ANCESTRY_HMM.out.posterior_files.collect()
    )

    emit:
    phased_vcf      = PHASE_PARENTS.out.phased_vcf
    posterior_files = RUN_ANCESTRY_HMM.out.posterior_files
    summary_dir     = ANCESTRY_SUMMARY.out.summary_dir
}
