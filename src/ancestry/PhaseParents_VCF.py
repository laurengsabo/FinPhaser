import allel
import argparse
import numpy as np

# Mapping for Chromosome renaming if needed
linkageGroups = {
    'NC_036780.1':'LG1', 'NC_036781.1':'LG2', 'NC_036782.1':'LG3', 
    'NC_036783.1':'LG4', 'NC_036784.1':'LG5', 'NC_036785.1':'LG6', 
    'NC_036786.1':'LG7', 'NC_036787.1':'LG8', 'NC_036788.1':'LG9', 
    'NC_036789.1':'LG10', 'NC_036790.1':'LG11', 'NC_036791.1':'LG12', 
    'NC_036792.1':'LG13', 'NC_036793.1':'LG14', 'NC_036794.1':'LG15', 
    'NC_036796.1':'LG17', 'NC_036797.1':'LG18', 'NC_036798.1':'LG19', 
    'NC_036799.1':'LG20', 'NC_036800.1':'LG22', 'NC_036801.1':'LG23'
}

def phase_and_filter_vcf(input_path, output_path):
    # 1. READ VCF FOR MASKING LOGIC
    print(f"Reading {input_path} for masking logic...")
    vcf_dict = allel.read_vcf(input_path, fields=['variants/CHROM', 'calldata/GT'])
    
    genotypes = allel.GenotypeArray(vcf_dict['calldata/GT'])
    paternal_gt = genotypes[:, 0]
    maternal_gt = genotypes[:, 1]

    # Calculate masks: (Pat Het & Mat Hom) OR (Mat Het & Pat Hom)
    pat_mask = paternal_gt.is_het() & maternal_gt.is_hom()
    mat_mask = maternal_gt.is_het() & paternal_gt.is_hom()
    combined_mask = pat_mask | mat_mask
    
    # Get indices of rows that pass the mask
    valid_indices = np.where(combined_mask)[0]
    
    # 2. PHASE ONLY THE FILTERED GENOTYPES
    print(f"Phasing {len(valid_indices)} conserved SNPs...")
    masked_genotypes = genotypes[combined_mask]
    phased_genotypes = allel.phase_by_transmission(masked_genotypes, window_size=1000)
    
    # Create a map for quick lookup: {Original_Row_Index: Phased_Row_Index}
    idx_map = {original_idx: i for i, original_idx in enumerate(valid_indices)}

    # 3. WRITE THE NEW VCF
    print(f"Writing phased and filtered data to {output_path}...")
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        data_row_count = 0
        
        for line in infile:
            # Conserve all header and informational rows starting with #
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            # Check if this data row is one we want to keep
            if data_row_count in idx_map:
                phased_idx = idx_map[data_row_count]
                fields = line.strip().split('\t')
                
                # Optional: Rename Chromosome using your linkageGroups map
                original_chrom = fields[0]
                fields[0] = linkageGroups.get(original_chrom, original_chrom)
                
                # Update Sample Columns (Genotypes start at index 9)
                # We replace the unphased GT (0/1) with phased GT (0|1)
                for sample_col in range(9, len(fields)):
                    sample_data = fields[sample_col].split(':')
                    
                    # Get the phased alleles for this sample
                    alleles = phased_genotypes[phased_idx, sample_col - 9]
                    # Update only the GT portion, keeping AD, DP, etc. intact
                    sample_data[0] = f"{alleles[0]}|{alleles[1]}"
                    
                    fields[sample_col] = ':'.join(sample_data)
                
                outfile.write('\t'.join(fields) + '\n')
            
            data_row_count += 1

    print("Done!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('VCFFile', type=str, help='Input VCF file')
    parser.add_argument('--out', type=str, default='conserved_phased.vcf', help='Output VCF name')
    args = parser.parse_args()
    
    phase_and_filter_vcf(args.VCFFile, args.out)