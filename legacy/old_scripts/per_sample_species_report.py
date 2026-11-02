#!/usr/bin/env python3
"""
per_sample_species_report.py

Generates one markdown table per sample from final_report.tsv (species-level,
NOT the genus-collapsed version) — species name + % of that sample's
ITS-region depth, sorted highest first. This is species-level detail (the
genus_rollup.py / pivot_genus_matrix.py scripts collapse tie-breaking noise
within a genus; this script deliberately does NOT, so you see every
distinct BLAST-called species label as-is, for a full-detail per-sample view).

The "(from title, taxdb not installed)" suffix on names (see
summarize_blast_hits.py) is stripped for display and explained once in a
footnote instead, to keep the tables readable.

Usage:
    python3 per_sample_species_report.py --report final_report.tsv --out per_sample_species.md
    python3 per_sample_species_report.py --report final_report.tsv --out per_sample_species.md --min-pct 0.1
"""
import argparse
import csv
import re
from collections import defaultdict

_SUFFIX_RE = re.compile(r"\s*\(from title, taxdb not installed\)\s*$")


def clean_name(species: str) -> tuple[str, bool]:
    """Returns (display_name, was_inferred_from_title)."""
    stripped = _SUFFIX_RE.sub("", species)
    return stripped, (stripped != species)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", required=True, help="Path to final_report.tsv (species-level)")
    ap.add_argument("--out", required=True, help="Path to write the markdown report")
    ap.add_argument("--min-pct", type=float, default=0.0, help="Omit species below this %% of sample depth (default: show all)")
    args = ap.parse_args()

    rows_by_sample: dict[str, list[dict]] = defaultdict(list)
    any_inferred = False

    with open(args.report, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pct = float(row["pct_of_sample_its_depth"])
            if pct < args.min_pct:
                continue
            name, inferred = clean_name(row["species"])
            any_inferred = any_inferred or inferred
            rows_by_sample[row["sample"]].append(
                {
                    "species": name,
                    "pct": pct,
                    "depth": float(row["summed_depth"]),
                    "n_contigs": row["n_contigs"],
                    "pident": row["best_pident"],
                    "bitscore": row["best_bitscore"],
                }
            )

    lines = ["# Per-sample species composition\n"]
    if any_inferred:
        lines.append(
            "_Species names were resolved from BLAST record titles "
            "(local `taxdb` not installed on the machine that ran BLAST) "
            "rather than from NCBI's taxonomy database directly. See "
            "`summarize_blast_hits.py` for how names were extracted._\n"
        )

    for sample in sorted(rows_by_sample):
        hits = sorted(rows_by_sample[sample], key=lambda r: r["pct"], reverse=True)
        lines.append(f"\n## {sample}\n")
        lines.append("| Species | % of sample ITS depth | Depth | # contigs | %identity | bitscore |")
        lines.append("|---|---|---|---|---|---|")
        for h in hits:
            lines.append(
                f"| {h['species']} | {h['pct']:.2f}% | {h['depth']:.2f} | "
                f"{h['n_contigs']} | {h['pident']}% | {h['bitscore']} |"
            )

    with open(args.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"Wrote {len(rows_by_sample)} sample table(s) -> {args.out}")


if __name__ == "__main__":
    main()
