/*
 * ============================================================
 *  workflows/ibd.nf
 * ============================================================
 *  IBD detection subworkflow.
 *
 *  The phased VCF produced by PhaseParents_VCF.py is used as
 *  the single VCF input for both the AHMM branch and this
 *  branch. SPORE requires the VCF to be bgzipped and
 *  tabix-indexed, so COMPRESS_VCF handles that before
 *  generating the SPORE settings.
 *
 *  Step 6a — COMPRESS_VCF
 *    Tools  : bgzip, tabix (from htslib, already in environment)
 *    In     : YHPed1_conserved_phased.vcf
 *    Out    : YHPed1_conserved_phased.vcf.gz + .tbi index
 *
 *  Step 6b — WRITE_SEX_TSV
 *    Helper : bin/yaml_to_sex_tsv.py
 *    In     : samples.yml (sex field per sample)
 *    Out    : Genomics_Sex.tsv
 *             Columns: indv, GenomicsSex
 *             Underscores stripped from names (SPORE requirement):
 *               YH_006_f → YH006f,  YH_016 → YH016
 *
 *  Step 6c — RUN_SPORE
 *    Helper : bin/generate_spore_settings.py  (no yaml import)
 *    Script : src/SPORE.R
 *    In     : compressed VCF + Genomics_Sex.tsv + samples.yml
 *    Out    : *.ibd (truffle IBD segment calls), SPORE log
 *
 *  Step 7 — RANK_IBD
 *    Script : src/reporting/rank_ibd.py
 *    In     : *.ibd
 *    Out    : inbreeding_rankings.txt
 *             inbreeding_rankings.tsv  (TSV copy for plotting)
 * ============================================================
 */

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS: COMPRESS_VCF
// ─────────────────────────────────────────────────────────────────────────────
process COMPRESS_VCF {
    tag "CompressVCF"
    // Publish both the .gz and its index so they are available for inspection
    publishDir "${params.outdir}/ibd/06_spore_input", mode: 'copy'

    input:
    path phased_vcf     // plain .vcf from PhaseParents_VCF.py

    output:
    path "*.vcf.gz",    emit: vcf_gz
    path "*.vcf.gz.tbi", emit: vcf_tbi

    script:
    /*
     * bgzip compresses the VCF in a way that tabix can index.
     * tabix -p vcf builds the .tbi index that SPORE/bcftools need.
     * Both htslib tools are already installed via environment.yml.
     */
    """
    bgzip -c ${phased_vcf} > ${phased_vcf}.gz
    tabix -p vcf ${phased_vcf}.gz
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS: WRITE_SEX_TSV
// ─────────────────────────────────────────────────────────────────────────────
process WRITE_SEX_TSV {
    tag "WriteSexTSV"
    publishDir "${params.outdir}/ibd/06_spore_input", mode: 'copy'

    input:
    path samples_yml

    output:
    path "Genomics_Sex.tsv", emit: sex_tsv

    script:
    /*
     * yaml_to_sex_tsv.py reads the sex field from every sample entry
     * in samples.yml and writes Genomics_Sex.tsv with:
     *   - header: indv<TAB>GenomicsSex
     *   - underscores stripped from sample names (SPORE requirement)
     *   - one unique row per sample (parents deduplicated)
     */
    """
    python ${projectDir}/bin/yaml_to_sex_tsv.py \
        --samples-yml ${samples_yml} \
        --out         Genomics_Sex.tsv
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS: RUN_SPORE
// ─────────────────────────────────────────────────────────────────────────────
process RUN_SPORE {
    tag "SPORE"
    publishDir "${params.outdir}/ibd/06_spore_output", mode: 'copy'

    input:
    path vcf_gz          // bgzipped + tabix-indexed phased VCF
    path vcf_tbi         // .tbi index (must be staged alongside .gz)
    path sex_tsv         // Genomics_Sex.tsv from WRITE_SEX_TSV
    path samples_yml     // for reading spore.output_prefix

    output:
    path "*.ibd",            emit: ibd_file
    path "SPORE_output.log", emit: spore_log
    path "SPORE-Settings.R", emit: spore_settings    // saved for reproducibility

    script:
    /*
     * generate_spore_settings.py (no yaml import — args via CLI only):
     *   --vcf      path to the .vcf.gz file
     *   --sex-tsv  path to Genomics_Sex.tsv
     *   --prefix   output file prefix (from nextflow.config params)
     *   --out      SPORE-Settings.R
     *
     * SPORE.R then reads SPORE-Settings.R to find its inputs.
     */
    """
    python ${projectDir}/bin/generate_spore_settings.py \
        --vcf     ${vcf_gz}           \
        --sex-tsv ${sex_tsv}          \
        --prefix  ${params.spore_output_prefix} \
        --out     SPORE-Settings.R

    Rscript ${projectDir}/src/SPORE.R SPORE-Settings.R > SPORE_output.log 2>&1
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS: RANK_IBD
// ─────────────────────────────────────────────────────────────────────────────
process RANK_IBD {
    tag "RankIBD"
    publishDir "${params.outdir}/ibd/07_rankings", mode: 'copy'

    input:
    path ibd_file

    output:
    path "inbreeding_rankings.txt", emit: rankings_txt
    path "inbreeding_rankings.tsv", emit: rankings_tsv

    script:
    /*
     * rank_ibd.py produces the plain-text rankings file.
     * txt_to_tsv.py converts it to .tsv for downstream R/Python plotting.
     */
    """
    python ${projectDir}/src/reporting/rank_ibd.py \
        --input  ${ibd_file} \
        --output inbreeding_rankings.txt

    python ${projectDir}/bin/txt_to_tsv.py \
        --input  inbreeding_rankings.txt \
        --output inbreeding_rankings.tsv
    """
}

// ─────────────────────────────────────────────────────────────────────────────
//  SUBWORKFLOW
// ─────────────────────────────────────────────────────────────────────────────
workflow IBD_WORKFLOW {

    take:
    ch_phased_vcf   // Channel<Path> — plain .vcf from ANCESTRY_WORKFLOW
    ch_samples_yml  // Channel<Path> — samples.yml

    main:

    // Step 6a — compress phased VCF for SPORE
    COMPRESS_VCF(ch_phased_vcf)

    // Step 6b — write Genomics_Sex.tsv from samples.yml
    WRITE_SEX_TSV(ch_samples_yml)

    // Step 6c — run SPORE
    RUN_SPORE(
        COMPRESS_VCF.out.vcf_gz,
        COMPRESS_VCF.out.vcf_tbi,
        WRITE_SEX_TSV.out.sex_tsv,
        ch_samples_yml
    )

    // Step 7 — rank IBD segments
    RANK_IBD(RUN_SPORE.out.ibd_file)

    emit:
    ibd_file     = RUN_SPORE.out.ibd_file
    rankings_txt = RANK_IBD.out.rankings_txt
    rankings_tsv = RANK_IBD.out.rankings_tsv
}
