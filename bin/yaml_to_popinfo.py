#!/usr/bin/env python3
"""
bin/yaml_to_popinfo.py
=======================
Converts the populations block of samples.yml into the tab-separated
mod_popinfo.txt format that mod_vcf2ahmm.py expects via its -s flag.

This replaces ancestry_input.py, which assigned populations based only
on name suffixes (_f → 0, _m → 1) and only supported 2 populations.
yaml_to_popinfo.py reads the explicit 4-haplotype assignments from
samples.yml and writes them in the same format mod_vcf2ahmm.py requires.

Output format (identical to the original mod_popinfo.txt):
    YH_006_f<TAB>0
    YH_006_f<TAB>1
    YH_011_m<TAB>2
    YH_011_m<TAB>3
    YH_016<TAB>admixed
    YH_017<TAB>admixed
    ...

Key rules:
  - A diploid parent appears TWICE (one row per haplotype it contributes).
  - Admixed offspring appear ONCE with the literal string "admixed".
    mod_vcf2ahmm.py checks for the string "admixed", not -1.
  - Row order follows the order in samples.yml.
  - The sex field in samples.yml is intentionally ignored here —
    it is used separately by yaml_to_sex_tsv.py for SPORE.

Uses only stdlib (no pyyaml) so it works even before the conda
environment is fully resolved.

Usage (called automatically by Nextflow — not normally run directly):
    python bin/yaml_to_popinfo.py \\
        --samples-yml config/samples.yml \\
        --out         mod_popinfo.txt
"""

import argparse
import sys
from pathlib import Path


# ── Minimal YAML parser (stdlib only) ────────────────────────────────────────
# samples.yml is simple enough that we can parse it without pyyaml.
# We only need the populations list, which is a flat sequence of mappings.

def parse_populations(yml_path: Path) -> list:
    """
    Parse only the `populations:` block from samples.yml.

    Returns a list of dicts, each with at least 'sample' and 'population'.
    Supports both block style and flow style (inline {}) entries.

    This parser is intentionally minimal — it handles the specific structure
    of FinPhaser's samples.yml and is not a general YAML parser.
    """
    text = yml_path.read_text()
    populations = []
    in_populations = False

    for raw_line in text.splitlines():
        line = raw_line.strip()

        # Detect the populations block
        if line == "populations:":
            in_populations = True
            continue

        # Stop at the next top-level key (unindented, ends with colon)
        if in_populations and line and not line.startswith("-") and not line.startswith("#"):
            if ":" in line and not line.startswith(" "):
                in_populations = False
                continue

        if not in_populations:
            continue

        # Skip blank lines and comments
        if not line or line.startswith("#"):
            continue

        # ── Flow style: - { sample: "X", population: 0, sex: F } ─────────────
        if line.startswith("- {") and "}" in line:
            inner = line[line.index("{") + 1: line.index("}")]
            entry = {}
            for part in inner.split(","):
                if ":" in part:
                    k, v = part.split(":", 1)
                    entry[k.strip()] = v.strip().strip('"').strip("'")
            if "sample" in entry and "population" in entry:
                populations.append(entry)
            continue

        # ── Block style: starts a new list item ───────────────────────────────
        if line.startswith("- sample:"):
            val = line.split(":", 1)[1].strip().strip('"').strip("'")
            populations.append({"sample": val})
            continue

        # ── Block style: continuation key under a list item ──────────────────
        if populations and ":" in line and not line.startswith("-"):
            k, v = line.split(":", 1)
            populations[-1][k.strip()] = v.strip().strip('"').strip("'")

    return populations


# ── Main ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Convert samples.yml populations block → mod_popinfo.txt"
    )
    p.add_argument("--samples-yml", required=True, help="Path to samples.yml")
    p.add_argument("--out", default="mod_popinfo.txt", help="Output popinfo file")
    return p.parse_args()


def main():
    args = parse_args()
    yml_path = Path(args.samples_yml)

    if not yml_path.exists():
        sys.exit(f"ERROR: samples.yml not found: {yml_path}")

    populations = parse_populations(yml_path)

    if not populations:
        sys.exit(
            "ERROR: No entries found under 'populations:' in samples.yml.\n"
            "Check indentation and that the key exists."
        )

    # Validate every entry has the required fields
    for i, entry in enumerate(populations):
        if "sample" not in entry:
            sys.exit(f"ERROR: populations entry #{i+1} is missing 'sample': {entry}")
        if "population" not in entry:
            sys.exit(f"ERROR: populations entry #{i+1} is missing 'population': {entry}")

    # Tally for informational output
    ref_rows = [e for e in populations if str(e["population"]) != "admixed"]
    adm_rows = [e for e in populations if str(e["population"]) == "admixed"]

    out_path = Path(args.out)
    with open(out_path, "w") as fh:
        for entry in populations:
            fh.write(f"{entry['sample']}\t{entry['population']}\n")

    print(
        f"[yaml_to_popinfo] Wrote {len(populations)} rows to {out_path}  "
        f"({len(ref_rows)} reference haplotype rows, {len(adm_rows)} admixed rows)"
    )


if __name__ == "__main__":
    main()
