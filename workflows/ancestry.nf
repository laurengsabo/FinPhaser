/*
 * ============================================================
 *  workflows/ancestry.nf
 * ============================================================
 *  Ancestry inference subworkflow.
 *
 *  Step 1 — PHASE_PARENTS
 *    Helper : bin/yaml_to_linkage_map.py
 *      Extracts the linkage_groups block from samples.yml and
 *      writes linkage_map.tsv — passed to PhaseParents_VCF.py
 *      via --linkage-map so the mapping is versioned in the repo
 *      and not hardcoded in the script.
 *    Script : src/ancestry/PhaseParents_VCF.py
 *    CLI    : python PhaseParents_VCF.py <VCFFile>
 *                 --out <outname> --linkage-map linkage_map.tsv
 *    In     : raw filtered VCF + samples.yml
 *    Out    : YHPed1_conserved_phased.vcf
 *             linkage_map.tsv  (saved for reproducibility)
 *
 *  Step 2+3 — PREPARE_HMM_INPUT  (combined into one process)
 *    Helper : bin/yaml_to_popinfo.py
 *      Reads samples.yml → writes mod_popinfo.txt in the exact
 *      format mod_vcf2ahmm.py expects (-s flag input).
 *    Script : src/ancestry/mod_vcf2ahmm.py  (original, unmodified)
 *      CLI  : python mod_vcf2ahmm.py -v <vcf> -s <popinfo>
 *                 -o_txt ancestry_input.txt -o_ploidy ahmm.ploidy
 *                 -g <int> -r <float> -m <int> --min_diff <float>
 *    In     : phased VCF + samples.yml
 *    Out    : ancestry_input.txt, ahmm.ploidy, mod_popinfo.txt
 *
 *  Step 4 — RUN_ANCESTRY_HMM
 *    Tool   : ancestry_hmm
 *    In     : ancestry_input.txt, ahmm.ploidy
 *    Out    : *.posterior (one file per admixed sample)
 *
 *  Step 5 — ANCESTRY_SUMMARY
 *    Script : src/ancestry/SNV_count.py
 *    CLI    : python SNV_count.py <directory>
 *             (always writes to ./ancestry_summary_report/ in cwd)
 *    In     : directory containing *.posterior files
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
    path samples_yml

    output:
    path "YHPed1_conserved_phased.vcf", emit: phased_vcf
    path "linkage_map.tsv",             emit: linkage_map  // saved for reproducibility

    script:
    /*
     * yaml_to_linkage_map.py reads the linkage_groups: block from
     * samples.yml and writes a two-column TSV:
     *     NC_036780.1<TAB>LG1
     *     NC_036781.1<TAB>LG2
     *     ...
     *
     * That TSV is then passed to PhaseParents_VCF.py via --linkage-map,
     * replacing the hardcoded dict in the original script. If the
     * linkage_groups block is ever absent from samples.yml, the script
     * exits with a clear error rather than silently using old values.
     */
    """
    # Extract linkage group mapping from samples.yml
    python ${projectDir}/bin/yaml_to_linkage_map.py \
        --samples-yml ${samples_yml} \
        --out         linkage_map.tsv

    # Phase and filter the VCF using the extracted mapping
    python ${projectDir}/src/ancestry/PhaseParents_VCF.py \
        ${vcf} \
        --out          YHPed1_conserved_phased.vcf \
        --linkage-map  linkage_map.tsv
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
    path "mod_popinfo.txt",    emit: popinfo       // saved for inspection / reproducibility

    script:
    /*
     * Step 2 (yaml_to_popinfo.py):
     *   Reads the populations block of samples.yml and writes mod_popinfo.txt
     *   in the exact tab-separated format mod_vcf2ahmm.py expects:
     *       YH_006_f<TAB>0
     *       YH_006_f<TAB>1
     *       YH_011_m<TAB>2
     *       YH_011_m<TAB>3
     *       YH_016<TAB>admixed
     *       ...
     *
     * Step 3 (mod_vcf2ahmm.py — original script, not modified):
     *   Uses its real CLI flags exactly as documented in the script header.
     *   HMM parameters come from nextflow.config (params.hmm_*).
     */
    """
    # Step 2: convert samples.yml → mod_popinfo.txt
    python ${projectDir}/bin/yaml_to_popinfo.py \
        --samples-yml ${samples_yml} \
        --out         mod_popinfo.txt

    # Step 3: convert phased VCF + popinfo → ancestry_hmm input files
    python ${projectDir}/src/ancestry/mod_vcf2ahmm.py \
        -v        ${phased_vcf}                    \
        -s        mod_popinfo.txt                  \
        -o_txt    ancestry_input.txt               \
        -o_ploidy ahmm.ploidy                      \
        -g        ${params.hmm_use_genotypes}      \
        -r        ${params.hmm_recombination_rate} \
        -m        ${params.hmm_min_distance_bp}    \
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
     * 4-population model matching the parental haplotype structure:
     *   pop 0 = maternal haplotype 1  (YH_006_f, allele 1)
     *   pop 1 = maternal haplotype 2  (YH_006_f, allele 2)
     *   pop 2 = paternal haplotype 1  (YH_011_m, allele 1)
     *   pop 3 = paternal haplotype 2  (YH_011_m, allele 2)
     *
     * -a 4 0.25 0.25 0.25 0.25   equal prior on all 4 haplotypes
     * -p 0..3  -2  0.25          each pulse ~2 generations back
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
    path posterior_files    // .collect() stages all *.posterior into cwd

    output:
    path "ancestry_summary_report/", emit: summary_dir

    script:
    /*
     * SNV_count.py takes one positional arg: the directory with *.posterior
     * files. It always writes to ./ancestry_summary_report/ in cwd.
     * We pass "." because Nextflow staged all files here.
     */
    """
    python ${projectDir}/src/ancestry/SNV_count.py .
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

    // Step 1 — phase parents (linkage map extracted from samples.yml)
    PHASE_PARENTS(ch_vcf, ch_samples_yml)

    // Steps 2 + 3 — build popinfo then run mod_vcf2ahmm.py
    PREPARE_HMM_INPUT(
        PHASE_PARENTS.out.phased_vcf,
        ch_samples_yml
    )

    // Step 4 — run ancestry_hmm
    RUN_ANCESTRY_HMM(
        PREPARE_HMM_INPUT.out.ancestry_input,
        PREPARE_HMM_INPUT.out.ploidy_file
    )

    // Step 5 — summarise posteriors (collect all files first)
    ANCESTRY_SUMMARY(
        RUN_ANCESTRY_HMM.out.posterior_files.collect()
    )

    emit:
    phased_vcf      = PHASE_PARENTS.out.phased_vcf
    posterior_files = RUN_ANCESTRY_HMM.out.posterior_files
    summary_dir     = ANCESTRY_SUMMARY.out.summary_dir
}
