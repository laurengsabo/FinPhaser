import pandas as pd

def rank_inbreeding_potential(input_file, output_file):
    # Read the space-delimited IBD file
    # TRUFFLE output usually has multiple spaces as delimiters
    try:
        df = pd.read_csv(input_file, sep=r'\s+')
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Calculate a simple Kinship Coefficient (Phi)
    # Phi = (0.5 * IBD2) + (0.25 * IBD1)
    # High IBD2 is the strongest indicator of potential inbreeding in future crosses
    df['Relatedness_Score'] = (df['IBD2'] * 0.5) + (df['IBD1'] * 0.25)

    # Sort by Relatedness_Score descending (highest sharing at the top)
    ranked_df = df.sort_values(by='Relatedness_Score', ascending=False)

    # Select and rename columns for the final report
    report = ranked_df[[
        'ID1', 'ID2', 'IBD0', 'IBD1', 'IBD2', 'Relatedness_Score'
    ]].copy()

    # Add a column to flag pairs that are likely related based on your plot threshold
    report['Likely_Related'] = report['IBD0'].apply(lambda x: 'YES' if x < 0.95 else 'NO')

    # Save to a text file
    with open(output_file, 'w') as f:
        f.write("RANKED PAIRS BY INBREEDING POTENTIAL\n")
        f.write("====================================\n")
        f.write("Relatedness Score calculated as: (0.5 * IBD2) + (0.25 * IBD1)\n\n")
        f.write(report.to_string(index=False))

    print(f"Ranking complete. Results saved to {output_file}")

# Usage
rank_inbreeding_potential('spore_brood1_fullsex/brood1_final.vcf.gz-truffle.ibd', 'inbreeding_rankings.txt')