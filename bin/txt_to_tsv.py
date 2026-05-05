#!/usr/bin/env python3
"""
bin/txt_to_tsv.py
==================
Converts the whitespace-delimited inbreeding_rankings.txt produced by
rank_ibd.py into a proper tab-separated .tsv file.

The .tsv copy is more convenient for downstream visualisation in R or
Python (pandas.read_csv with sep='\\t').

Usage (called by Nextflow — not normally run directly):
    python bin/txt_to_tsv.py \\
        --input  inbreeding_rankings.txt \\
        --output inbreeding_rankings.tsv
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Convert rank_ibd .txt output to tab-separated .tsv"
    )
    p.add_argument("--input",  required=True, help="Input .txt file")
    p.add_argument("--output", required=True, help="Output .tsv file")
    return p.parse_args()


def main():
    args = parse_args()
    in_path  = Path(args.input)
    out_path = Path(args.output)

    if not in_path.exists():
        sys.exit(f"ERROR: Input file not found: {in_path}")

    lines = in_path.read_text().splitlines()
    if not lines:
        sys.exit(f"ERROR: Input file is empty: {in_path}")

    out_lines = []
    for line in lines:
        if not line.strip():
            out_lines.append("")
            continue
        if "\t" in line:
            out_lines.append(line)
        else:
            # Collapse any run of whitespace into a single tab
            out_lines.append("\t".join(line.split()))

    out_path.write_text("\n".join(out_lines) + "\n")
    print(
        f"[txt_to_tsv] Converted {len(out_lines)} lines: "
        f"{in_path.name} → {out_path.name}"
    )


if __name__ == "__main__":
    main()
