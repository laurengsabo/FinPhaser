#!/usr/bin/env python3
"""
bin/yaml_to_sex_tsv.py
=======================
Generates the Genomics_Sex.tsv file that SPORE needs from the sex
fields embedded in samples.yml.

Output format (matches original Genomics_Sex.tsv exactly):
    indv<TAB>GenomicsSex
    YH_011_m<TAB>M
    YH_006_f<TAB>F
    YH_016<TAB>F
    ...

IMPORTANT: Sample names are written exactly as they appear in samples.yml,
with underscores preserved. SPORE matches these names against the VCF
sample names, which also contain underscores — stripping them causes
mismatches and breaks the sex-based analysis.

Each unique sample is written ONCE (parents appear twice in the
populations block for haplotype tracking, but only once in the TSV).

Uses only stdlib — no pyyaml dependency.

Usage (called automatically by Nextflow — not normally run directly):
    python bin/yaml_to_sex_tsv.py \\
        --samples-yml config/samples.yml \\
        --out         Genomics_Sex.tsv
"""

import argparse
import sys
from pathlib import Path


def parse_populations(yml_path: Path) -> list:
    """
    Parse the populations block from samples.yml.
    Returns a list of dicts with at least 'sample', 'population', and 'sex'.
    """
    text = yml_path.read_text()
    populations = []
    in_populations = False

    for raw_line in text.splitlines():
        line = raw_line.strip()

        if line == "populations:":
            in_populations = True
            continue

        if in_populations and line and not line.startswith("-") and not line.startswith("#"):
            if ":" in line and not line.startswith(" "):
                in_populations = False
                continue

        if not in_populations:
            continue
        if not line or line.startswith("#"):
            continue

        # Flow style: - { sample: "X", population: admixed, sex: F }
        if line.startswith("- {") and "}" in line:
            inner = line[line.index("{") + 1: line.index("}")]
            entry = {}
            for part in inner.split(","):
                if ":" in part:
                    k, v = part.split(":", 1)
                    entry[k.strip()] = v.strip().strip('"').strip("'")
            if "sample" in entry and "sex" in entry:
                populations.append(entry)
            continue

        # Block style — start new item
        if line.startswith("- sample:"):
            val = line.split(":", 1)[1].strip().strip('"').strip("'")
            populations.append({"sample": val})
            continue

        # Block style — continuation
        if populations and ":" in line and not line.startswith("-"):
            k, v = line.split(":", 1)
            # Strip inline comments
            v_clean = v.split("#")[0].strip().strip('"').strip("'")
            populations[-1][k.strip()] = v_clean

    return populations


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate Genomics_Sex.tsv for SPORE from samples.yml"
    )
    p.add_argument("--samples-yml", required=True, help="Path to samples.yml")
    p.add_argument("--out", default="Genomics_Sex.tsv", help="Output TSV file")
    return p.parse_args()


def main():
    args = parse_args()
    yml_path = Path(args.samples_yml)

    if not yml_path.exists():
        sys.exit(f"ERROR: samples.yml not found: {yml_path}")

    populations = parse_populations(yml_path)

    if not populations:
        sys.exit("ERROR: No entries found under 'populations:' in samples.yml.")

    # Deduplicate by sample name — parents appear twice in populations
    # (once per haplotype) but should appear only once in the sex TSV.
    seen = {}
    for entry in populations:
        sample = entry.get("sample", "")
        sex    = entry.get("sex", "")
        if not sample:
            continue
        if sample not in seen:
            if not sex:
                print(
                    f"[yaml_to_sex_tsv] WARNING: no sex specified for {sample} — "
                    f"omitting from TSV",
                    file=sys.stderr,
                )
                continue
            seen[sample] = sex
        elif seen[sample] != sex:
            sys.exit(
                f"ERROR: conflicting sex values for sample '{sample}': "
                f"'{seen[sample]}' vs '{sex}'. Check samples.yml."
            )

    # Write TSV — preserve sample names exactly as in samples.yml
    # (underscores kept — SPORE matches against VCF names which also have underscores)
    out_path = Path(args.out)
    with open(out_path, "w") as fh:
        fh.write("indv\tGenomicsSex\n")
        for sample, sex in seen.items():
            fh.write(f"{sample}\t{sex}\n")

    print(
        f"[yaml_to_sex_tsv] Wrote {len(seen)} samples to {out_path} "
        f"(sample names preserved with underscores)"
    )


if __name__ == "__main__":
    main()