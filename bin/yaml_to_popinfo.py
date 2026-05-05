#!/usr/bin/env python3
"""
bin/yaml_to_popinfo.py
=======================
Converts the `populations` block of samples.yml into the tab-separated
mod_popinfo.txt format that mod_vcf2ahmm.py expects.

This is the bridge that allows the Nextflow pipeline to use a single
samples.yml as its only config input, rather than requiring the user to
manually maintain a separate popinfo file.

Output format (matches the original mod_popinfo.txt exactly):
    YH_006_f<TAB>0
    YH_006_f<TAB>1
    YH_011_m<TAB>2
    YH_011_m<TAB>3
    YH_016<TAB>admixed
    YH_017<TAB>admixed
    ...

Key rules reflected here:
  - A parent sample appears ONCE per haplotype it contributes, so a
    diploid parent appears TWICE (population 0 and 1, or 2 and 3).
  - Admixed samples appear once with the literal string "admixed"
    (not -1 — that's what mod_vcf2ahmm.py checks for).
  - Output row order follows the order in samples.yml.

Usage (called automatically by Nextflow — not normally run directly):
    python bin/yaml_to_popinfo.py \\
        --samples-yml config/samples.yml \\
        --out         mod_popinfo.txt
"""

import argparse
import sys
from pathlib import Path

import yaml


def parse_args():
    p = argparse.ArgumentParser(
        description="Convert samples.yml populations block → mod_popinfo.txt"
    )
    p.add_argument(
        "--samples-yml", required=True,
        help="Path to samples.yml"
    )
    p.add_argument(
        "--out", default="mod_popinfo.txt",
        help="Output popinfo file (default: mod_popinfo.txt)"
    )
    return p.parse_args()


def main():
    args = parse_args()

    yml_path = Path(args.samples_yml)
    if not yml_path.exists():
        sys.exit(f"ERROR: samples.yml not found: {yml_path}")

    with open(yml_path) as fh:
        cfg = yaml.safe_load(fh)

    populations = cfg.get("populations", [])
    if not populations:
        sys.exit(
            "ERROR: No entries found under 'populations' in samples.yml. "
            "Check that the key exists and has at least one entry."
        )

    # Validate: every entry must have both 'sample' and 'population' keys
    for i, entry in enumerate(populations):
        if "sample" not in entry or "population" not in entry:
            sys.exit(
                f"ERROR: populations entry #{i+1} is missing 'sample' or "
                f"'population' key: {entry}"
            )

    # Count haplotype entries per parent (numeric population values)
    # for informational output — helps catch misconfigured YAMLs
    ref_entries = [e for e in populations if str(e["population"]) != "admixed"]
    adm_entries = [e for e in populations if str(e["population"]) == "admixed"]

    # Write popinfo — one line per entry, preserving YAML order
    out_path = Path(args.out)
    with open(out_path, "w") as fh:
        for entry in populations:
            sample = entry["sample"]
            pop    = entry["population"]   # 0, 1, 2, 3, or "admixed"
            fh.write(f"{sample}\t{pop}\n")

    print(
        f"[yaml_to_popinfo] Wrote {len(populations)} entries to {out_path}  "
        f"({len(ref_entries)} reference haplotype rows, "
        f"{len(adm_entries)} admixed rows)"
    )


if __name__ == "__main__":
    main()