/*
 * ============================================================
 *  workflows/ibd.nf
 * ============================================================
 *  IBD detection subworkflow.
 *
 *  Step 6 — RUN_SPORE
 *    Script : src/SPORE.R  (with auto-generated SPORE-Settings.R)
 *    In     : phased VCF, samples.yml (provides sex metadata path)
 *    Out    : *.ibd  (truffle IBD segment calls)
 *
 *    Two naming-format issues are handled automatically:
 *      a) The VCF sample names use underscores (e.g. YH_006_f,
 *         YH_016) but SPORE cannot handle underscores.
 *         bin/generate_spore_settings.py documents this and the
 *         SPORE.R script is expected to sanitise as needed.
 *      b) Genomics_Sex.tsv uses a different format (YH006f, YH016)
 *         from the VCF/popinfo names.  generate_spore_settings.py
 *         notes this discrepancy in the generated settings file.
 *
 *  Step 7 — RANK_IBD
 *    Script : src/reporting/rank_ibd.py
 *    In     : *.ibd file from SPORE
 *    Out    : inbreeding_rankings.txt
 *             inbreeding_rankings.tsv  (for downstream plotting)
 * ============================================================
 */

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS: RUN_SPORE
// ─────────────────────────────────────────────────────────────────────────────
process RUN_SPORE {
    tag "SPORE"
    publishDir "${params.outdir}/ibd/06_spore", mode: 'copy'

    input:
    path phased_vcf
    path samples_yml

    output:
    path "*.ibd",            emit: ibd_file
    path "SPORE_output.log", emit: spore_log
    path "SPORE-Settings.R", emit: spore_settings   // saved for reproducibility

    script:
    /*
     * generate_spore_settings.py reads samples.yml and writes a
     * SPORE-Settings.R pointing to the correct VCF and sex metadata.
     *
     * Important naming notes captured in the generated settings file:
     *   - VCF sample names contain underscores (SPORE limitation)
     *   - Genomics_Sex.tsv uses a different naming convention
     *     (no underscores, no separator: YH006f vs YH_006_f)
     * Both issues are documented in SPORE-Settings.R for transparency.
     */
    """
    python ${projectDir}/bin/generate_spore_settings.py \
        --vcf         ${phased_vcf} \
        --samples-yml ${samples_yml} \
        --out         SPORE-Settings.R

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
     * rank_ibd.py produces a ranked relatedness table.
     * txt_to_tsv.py writes a .tsv copy for use in R/Python
     * visualisation scripts without manual reformatting.
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
    ch_phased_vcf
    ch_samples_yml

    main:
    RUN_SPORE(ch_phased_vcf, ch_samples_yml)
    RANK_IBD(RUN_SPORE.out.ibd_file)

    emit:
    ibd_file     = RUN_SPORE.out.ibd_file
    rankings_txt = RANK_IBD.out.rankings_txt
    rankings_tsv = RANK_IBD.out.rankings_tsv
}
