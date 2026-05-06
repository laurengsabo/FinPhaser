import argparse
import pandas as pd


def rank_inbreeding_potential(input_file, output_file):
    """
    Read a TRUFFLE .ibd file and produce a ranked relatedness report.

    Relatedness Score = (0.5 * IBD2) + (0.25 * IBD1)
    High IBD2 indicates close inbreeding; high IBD1 indicates
    more distant shared ancestry.
    """
    try:
        df = pd.read_csv(input_file, sep=r'\s+')
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Calculate kinship coefficient proxy
    df['Relatedness_Score'] = (df['IBD2'] * 0.5) + (df['IBD1'] * 0.25)

    # Sort descending — most related pairs first
    ranked_df = df.sort_values(by='Relatedness_Score', ascending=False)

    report = ranked_df[[
        'ID1', 'ID2', 'IBD0', 'IBD1', 'IBD2', 'Relatedness_Score'
    ]].copy()

    # Flag pairs with IBD0 < 0.95 as likely related
    report['Likely_Related'] = report['IBD0'].apply(
        lambda x: 'YES' if x < 0.95 else 'NO'
    )

    with open(output_file, 'w') as f:
        f.write("RANKED PAIRS BY INBREEDING POTENTIAL\n")
        f.write("====================================\n")
        f.write("Relatedness Score calculated as: (0.5 * IBD2) + (0.25 * IBD1)\n\n")
        f.write(report.to_string(index=False))

    print(f"Ranking complete. Results saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Rank IBD pairs by inbreeding potential from TRUFFLE output."
    )
    parser.add_argument(
        "--input",  required=True,
        help="Path to TRUFFLE .ibd file"
    )
    parser.add_argument(
        "--output", required=True,
        help="Path to output rankings text file"
    )
    args = parser.parse_args()
    rank_inbreeding_potential(args.input, args.output)