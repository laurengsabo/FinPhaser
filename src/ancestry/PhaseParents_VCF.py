import allel
import argparse
import numpy as np

# Default chromosome → linkage group mapping for Aulonocara 'Yellow Head'.
# This is used as a fallback only if --linkage-map is not supplied.
# To change the mapping, edit the linkage_groups: block in config/samples.yml
# and the pipeline will pass a generated TSV via --linkage-map automatically.
DEFAULT_LINKAGE_GROUPS = {
    'NC_036780.1': 'LG1',  'NC_036781.1': 'LG2',  'NC_036782.1': 'LG3',
    'NC_036783.1': 'LG4',  'NC_036784.1': 'LG5',  'NC_036785.1': 'LG6',
    'NC_036786.1': 'LG7',  'NC_036787.1': 'LG8',  'NC_036788.1': 'LG9',
    'NC_036789.1': 'LG10', 'NC_036790.1': 'LG11', 'NC_036791.1': 'LG12',
    'NC_036792.1': 'LG13', 'NC_036793.1': 'LG14', 'NC_036794.1': 'LG15',
    'NC_036796.1': 'LG17', 'NC_036797.1': 'LG18', 'NC_036798.1': 'LG19',
    'NC_036799.1': 'LG20', 'NC_036800.1': 'LG22', 'NC_036801.1': 'LG23',
}


def load_linkage_map(tsv_path: str) -> dict:
    """
    Load a chromosome → linkage group mapping from a two-column TSV file.

    File format (no header, tab-separated):
        NC_036780.1<TAB>LG1
        NC_036781.1<TAB>LG2
        ...

    Generated at runtime by bin/yaml_to_linkage_map.py from samples.yml.
    Any line starting with '#' is treated as a comment and skipped.
    """
    mapping = {}
    with open(tsv_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) != 2:
                raise ValueError(
                    f"Bad line in linkage map (expected 2 tab-separated fields): {line!r}"
                )
            chrom, lg = parts
            mapping[chrom.strip()] = lg.strip()
    print(f"[PhaseParents_VCF] Loaded {len(mapping)} chromosome mappings from {tsv_path}")
    return mapping


def phase_and_filter_vcf(input_path, output_path, linkage_groups):
    # 1. READ VCF FOR MASKING LOGIC
    print(f"[PhaseParents_VCF] Reading {input_path} for masking logic...")
    vcf_dict = allel.read_vcf(input_path, fields=['variants/CHROM', 'calldata/GT'])

    genotypes = allel.GenotypeArray(vcf_dict['calldata/GT'])
    paternal_gt = genotypes[:, 0]
    maternal_gt = genotypes[:, 1]

    # Masks: keep sites where one parent is het and the other is hom.
    # These are the informative sites for phasing by transmission.
    pat_mask = paternal_gt.is_het() & maternal_gt.is_hom()
    mat_mask = maternal_gt.is_het() & paternal_gt.is_hom()
    combined_mask = pat_mask | mat_mask

    valid_indices = np.where(combined_mask)[0]
    print(f"[PhaseParents_VCF] {len(valid_indices)} conserved informative SNPs pass filter")

    # 2. PHASE ONLY THE FILTERED GENOTYPES
    print(f"[PhaseParents_VCF] Phasing {len(valid_indices)} conserved SNPs...")
    masked_genotypes = genotypes[combined_mask]
    phased_genotypes = allel.phase_by_transmission(masked_genotypes, window_size=1000)

    # Map original VCF row index → index in the phased array
    idx_map = {original_idx: i for i, original_idx in enumerate(valid_indices)}

    # 3. WRITE THE NEW VCF
    print(f"[PhaseParents_VCF] Writing phased and filtered data to {output_path}...")
    written = 0
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        data_row_count = 0

        for line in infile:
            # Preserve all header lines unchanged
            if line.startswith('#'):
                outfile.write(line)
                continue

            if data_row_count in idx_map:
                phased_idx = idx_map[data_row_count]
                fields = line.strip().split('\t')

                # Rename chromosome using the linkage group mapping.
                # If the CHROM value is not in the map, it is kept as-is.
                original_chrom = fields[0]
                fields[0] = linkage_groups.get(original_chrom, original_chrom)

                # Replace unphased GT (0/1) with phased GT (0|1) for each sample.
                # All other FORMAT fields (AD, DP, GQ, etc.) are preserved.
                for sample_col in range(9, len(fields)):
                    sample_data = fields[sample_col].split(':')
                    alleles = phased_genotypes[phased_idx, sample_col - 9]
                    sample_data[0] = f"{alleles[0]}|{alleles[1]}"
                    fields[sample_col] = ':'.join(sample_data)

                outfile.write('\t'.join(fields) + '\n')
                written += 1

            data_row_count += 1

    print(f"[PhaseParents_VCF] Done — wrote {written} phased SNVs to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Phase parental genotypes and filter conserved informative SNVs."
    )
    parser.add_argument(
        'VCFFile', type=str,
        help='Input VCF file (filtered SNVs)'
    )
    parser.add_argument(
        '--out', type=str, default='conserved_phased.vcf',
        help='Output phased VCF filename (default: conserved_phased.vcf)'
    )
    parser.add_argument(
        '--linkage-map', type=str, default=None,
        help=(
            'Path to a two-column tab-separated file mapping chromosome IDs '
            'to linkage group names (e.g. NC_036780.1<TAB>LG1). '
            'Generated automatically from samples.yml by the pipeline. '
            'If omitted, the built-in Aulonocara YH default mapping is used.'
        )
    )
    args = parser.parse_args()

    # Load linkage group mapping — from file if provided, else use default
    if args.linkage_map:
        linkage_groups = load_linkage_map(args.linkage_map)
    else:
        print(
            "[PhaseParents_VCF] No --linkage-map provided — "
            "using built-in default mapping."
        )
        linkage_groups = DEFAULT_LINKAGE_GROUPS

    phase_and_filter_vcf(args.VCFFile, args.out, linkage_groups)
