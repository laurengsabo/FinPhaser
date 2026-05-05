#!/usr/bin/env python3
"""
bin/yaml_to_linkage_map.py
===========================
Extracts the linkage_groups block from samples.yml and writes a
two-column tab-separated file that PhaseParents_VCF.py reads via
its --linkage-map flag.

Output format (no header, tab-separated):
    NC_036780.1<TAB>LG1
    NC_036781.1<TAB>LG2
    ...

If linkage_groups is absent from samples.yml, this script exits with
a clear message — the Nextflow process only calls it when the block
is present, so PhaseParents_VCF.py falls back to its built-in default
when the block is missing.

Uses only stdlib — no pyyaml dependency.

Usage (called automatically by Nextflow — not normally run directly):
    python bin/yaml_to_linkage_map.py \\
        --samples-yml config/samples.yml \\
        --out         linkage_map.tsv
"""

import argparse
import sys
from pathlib import Path


def parse_linkage_groups(yml_path: Path) -> dict:
    """
    Parse the linkage_groups: block from samples.yml.

    Expected format in samples.yml:
        linkage_groups:
          NC_036780.1: LG1
          NC_036781.1: LG2
          ...

    Returns a dict {chrom_id: lg_name} in file order.
    Returns an empty dict if the block is absent.
    """
    text = yml_path.read_text()
    groups = {}
    in_block = False

    for raw_line in text.splitlines():
        line = raw_line.strip()

        if line == "linkage_groups:":
            in_block = True
            continue

        # Stop at next top-level key (zero-indent line ending with colon)
        if in_block and line and not line.startswith("#"):
            indent = len(raw_line) - len(raw_line.lstrip())
            if indent == 0 and line.endswith(":"):
                break

        if not in_block or not line or line.startswith("#"):
            continue

        if ":" in line:
            k, v = line.split(":", 1)
            chrom = k.strip().strip('"').strip("'")
            lg    = v.strip().strip('"').strip("'")
            if chrom and lg:
                groups[chrom] = lg

    return groups


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract linkage_groups from samples.yml → TSV for PhaseParents_VCF.py"
    )
    p.add_argument("--samples-yml", required=True, help="Path to samples.yml")
    p.add_argument("--out", default="linkage_map.tsv", help="Output TSV file")
    return p.parse_args()


def main():
    args = parse_args()
    yml_path = Path(args.samples_yml)

    if not yml_path.exists():
        sys.exit(f"ERROR: samples.yml not found: {yml_path}")

    groups = parse_linkage_groups(yml_path)

    if not groups:
        sys.exit(
            "ERROR: No linkage_groups block found in samples.yml.\n"
            "Add a linkage_groups: section or remove --linkage-map from "
            "the PhaseParents_VCF.py call to use the built-in default."
        )

    out_path = Path(args.out)
    with open(out_path, "w") as fh:
        fh.write("# Chromosome ID → Linkage Group mapping\n")
        fh.write("# Generated from samples.yml linkage_groups block\n")
        for chrom, lg in groups.items():
            fh.write(f"{chrom}\t{lg}\n")

    print(f"[yaml_to_linkage_map] Wrote {len(groups)} mappings to {out_path}")


if __name__ == "__main__":
    main()
