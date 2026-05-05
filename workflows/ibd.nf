/*
 * ============================================================
 *  workflows/ibd.nf
 * ============================================================
 *  IBD detection subworkflow.
 *
 *  The phased VCF from PhaseParents_VCF.py is the single VCF
 *  source for this branch. SPORE requires it to be bgzipped
 *  and tabix-indexed with a CSI index (not TBI), matching the
 *  working manual run which produced .csi files.
 *
 *  Step 6a — COMPRESS_VCF
 *    Tools  : bgzip, bcftools index (CSI)
 *    In     : YHPed1_conserved_phased.vcf
 *    Out    : YHPed1_conserved_phased.vcf.gz  +  .vcf.gz.csi
 *
 *  Step 6b — WRITE_SEX_TSV
 *    Helper : bin/yaml_to_sex_tsv.py
 *    In     : samples.yml
 *    Out    : Genomics_Sex.tsv  (indv / GenomicsSex columns,
 *             underscores stripped from names for SPORE)
 *
 *  Step 6c — RUN_SPORE
 *    Helper : bin/generate_spore_settings.py
 *      Writes a complete SPORE-Settings.R with every variable
 *      SPORE.R reads — truffle path, filters, APO, memory, etc.
 *      `folder` = Nextflow working directory (where VCF is staged)
 *      `vcf`    = just the filename (SPORE builds the full path)
 *    Script : src/spore/SPORE.R
 *    Dep    : src/spore/scripts/Mendel1.R  (sourced by SPORE.R)
 *             tools/truffle               (TRUFFLE binary)
 *    In     : compressed VCF + CSI index + Genomics_Sex.tsv
 *    Out    : *.ibd, *.ibd.iqr.gz, threshold files, plots, log
 *
 *  Step 7 — RANK_IBD
 *    Script : src/reporting/rank_ibd.py
 *    In     : *.ibd
 *    Out    : inbreeding_rankings.txt + .tsv
 * ============================================================
 */

nextflow.enable.dsl = 2

// ─────────────────────────────────────────────────────────────────────────────
//  PROCESS: COMPRESS_VCF
// ─────────────────────────────────────────────────────────────────────────────
process COMPRESS_VCF {
    tag "CompressVCF"
    publishDir "${params.outdir}/ibd/06_spore_input", mode: 'copy'

    input:
    path phased_vcf

    output:
    path "*.vcf.gz",     emit: vcf_gz
    path "*.vcf.gz.csi", emit: vcf_csi

    script:
    /*
     * bgzip compresses the plain VCF, keeping its original filename.
     * bcftools index --csi builds the CSI index SPORE/TRUFFLE expects.
     * The VCF name is never changed — SPORE's own output files will
     * be prefixed with whatever name this compressed VCF has.
     */
    """
    bgzip -c ${phased_vcf} > ${phased_vcf}.gz
    bcftools index --csi ${phased_vcf}.gz
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
     *   - underscores stripped from names (SPORE requirement):
     *       YH_006_f → YH006f,  YH_016 → YH016
     *   - parents deduplicated (one row per unique sample)
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

    // SPORE is memory-intensive — allocate generously
    cpus   { params.spore_max_cores }
    memory { params.spore_max_memory }
    time   '8 h'

    input:
    path vcf_gz          // bgzipped VCF
    path vcf_csi         // CSI index — must be staged alongside .gz
    path sex_tsv         // Genomics_Sex.tsv from WRITE_SEX_TSV
    path samples_yml     // for reading the full spore: block

    output:
    path "*-truffle.ibd",    emit: ibd_file          // main IBD output
    path "*",                emit: all_outputs        // all SPORE outputs
    path "SPORE_output.log", emit: spore_log
    path "SPORE-Settings.R", emit: spore_settings     // saved for reproducibility

    script:
    /*
     * generate_spore_settings.py writes a complete SPORE-Settings.R:
     *   --vcf-filename  just the filename (SPORE builds folder+vcf)
     *   --workdir       the Nextflow process working dir where VCF is staged
     *   --sex-tsv       absolute path to Genomics_Sex.tsv
     *   --samples-yml   the user config
     *
     * SPORE.R then sources scripts/Mendel1.R internally.
     * We pass the project-relative paths via ${projectDir}.
     */
    """
    # Generate complete SPORE-Settings.R from samples.yml
    python ${projectDir}/bin/generate_spore_settings.py \
        --vcf-filename ${vcf_gz}                           \
        --workdir      \$(pwd)/                            \
        --sex-tsv      \$(pwd)/${sex_tsv}                  \
        --samples-yml  ${samples_yml}                      \
        --out          SPORE-Settings.R

    # Run SPORE — it sources scripts/Mendel1.R relative to SPORE.R's location
    Rscript ${projectDir}/src/spore/SPORE.R SPORE-Settings.R > SPORE_output.log 2>&1
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
     * txt_to_tsv.py writes a .tsv copy for downstream R/Python plotting.
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

    // Step 6a — compress and index phased VCF for SPORE
    COMPRESS_VCF(ch_phased_vcf)

    // Step 6b — generate Genomics_Sex.tsv from samples.yml
    WRITE_SEX_TSV(ch_samples_yml)

    // Step 6c — run SPORE + TRUFFLE
    RUN_SPORE(
        COMPRESS_VCF.out.vcf_gz,
        COMPRESS_VCF.out.vcf_csi,
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
