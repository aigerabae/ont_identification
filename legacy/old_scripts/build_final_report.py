#!/usr/bin/env python3
"""
build_final_report.py

Joins BLAST species calls (blast_results/summary.tsv, from
summarize_blast_hits.py) with per-contig read depth
(pipeline_out/<sample>/qc/depth_its_mean.tsv, from the Snakemake pipeline)
into one master table: per sample, per species, total supporting depth and
% of that sample's total ITS-region depth — the same manual cross-reference
done by hand for the pilot sample (ONplants-1_S40), now automated across
all samples.

Two things this script handles that a naive join would get wrong:

1. Contig ID suffix mismatch: extract_its_regions.py appends
   "_ITS1_ITS2_combined" to the header of any contig whose ITS1+ITS2 were
   concatenated (see that script), but the depth file's contig IDs come
   from the original assembly headers, which never carry that suffix. This
   script strips it before joining, so those contigs aren't silently
   dropped as "no depth found".

2. Multiple contigs -> one species: as with the pilot's Hippophae rhamnoides
   (detected via 3 separate contigs/sequence variants), depths for every
   contig assigned to the same top-hit species in a sample are SUMMED
   before computing % of total, rather than reporting each contig as if it
   were an independent species.

Quality filtering: hits below --min-bitscore or above --max-evalue are
excluded from the report entirely (default: bitscore >= 100, evalue <=
1e-15) — adjust these based on your own confirmed hits if needed. Excluded
low-confidence hits are written to a separate file so nothing is silently
discarded without a record.

Usage:
    python3 build_final_report.py \
        --blast-summary blast_results/summary.tsv \
        --pipeline-out pipeline_out \
        --out final_report.tsv
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

_COMBINED_SUFFIX_RE = re.compile(r"_ITS1_ITS2_combined$")


def strip_combined_suffix(contig: str) -> str:
    return _COMBINED_SUFFIX_RE.sub("", contig)


def load_blast_summary(path: str, min_bitscore: float, max_evalue: float):
    """Returns {sample: {contig(depth-key): top_hit_row}}, plus a list of
    excluded (below-threshold) rows for the record."""
    top_hit_per_contig: dict[str, dict[str, dict]] = defaultdict(dict)
    excluded_rows: list[dict] = []

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["rank"] != "1":
                continue  # only the best hit per contig feeds the report
            try:
                bitscore = float(row["bitscore"])
                evalue = float(row["evalue"])
            except ValueError:
                continue

            sample = row["sample"]
            depth_key = strip_combined_suffix(row["contig"])

            if bitscore < min_bitscore or evalue > max_evalue:
                excluded_rows.append(row)
                continue

            top_hit_per_contig[sample][depth_key] = row

    return top_hit_per_contig, excluded_rows


def load_depth_file(path: str) -> dict[str, float]:
    """depth_its_mean.tsv has no header: columns are contig<TAB>mean_depth."""
    depths: dict[str, float] = {}
    if not os.path.exists(path):
        return depths
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2:
                continue
            contig, depth_str = parts
            try:
                depths[contig] = float(depth_str)
            except ValueError:
                continue
    return depths


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--blast-summary", required=True, help="Path to blast_results/summary.tsv")
    ap.add_argument("--pipeline-out", required=True, help="Path to the Snakemake pipeline_out directory")
    ap.add_argument("--out", required=True, help="Path to write the master report TSV")
    ap.add_argument("--min-bitscore", type=float, default=100.0)
    ap.add_argument("--max-evalue", type=float, default=1e-15)
    args = ap.parse_args()

    top_hits, excluded = load_blast_summary(args.blast_summary, args.min_bitscore, args.max_evalue)
    if not top_hits:
        sys.exit(f"No hits passed the quality filter (bitscore>={args.min_bitscore}, evalue<={args.max_evalue}) in {args.blast_summary}")

    out_rows = []
    samples_missing_depth = []

    for sample, contig_hits in sorted(top_hits.items()):
        depth_path = os.path.join(args.pipeline_out, sample, "qc", "depth_its_mean.tsv")
        depths = load_depth_file(depth_path)
        if not depths:
            samples_missing_depth.append(sample)

        # group contigs by species for this sample
        species_depth: dict[str, float] = defaultdict(float)
        species_contigs: dict[str, list[str]] = defaultdict(list)
        species_best_hit: dict[str, dict] = {}

        for depth_key, hit in contig_hits.items():
            species = hit["scientific_name"]
            depth = depths.get(depth_key, 0.0)
            species_depth[species] += depth
            species_contigs[species].append(depth_key)
            # keep the highest-bitscore hit as the representative record
            if species not in species_best_hit or float(hit["bitscore"]) > float(species_best_hit[species]["bitscore"]):
                species_best_hit[species] = hit

        total_depth = sum(species_depth.values())

        for species, depth in sorted(species_depth.items(), key=lambda kv: kv[1], reverse=True):
            pct = (depth / total_depth * 100) if total_depth > 0 else 0.0
            best = species_best_hit[species]
            out_rows.append(
                {
                    "sample": sample,
                    "species": species,
                    "n_contigs": len(species_contigs[species]),
                    "contigs": ";".join(species_contigs[species]),
                    "summed_depth": round(depth, 2),
                    "pct_of_sample_its_depth": round(pct, 3),
                    "best_pident": best["pident"],
                    "best_bitscore": best["bitscore"],
                    "best_evalue": best["evalue"],
                    "best_accession": best["accession"],
                }
            )

    out_cols = [
        "sample", "species", "n_contigs", "contigs", "summed_depth",
        "pct_of_sample_its_depth", "best_pident", "best_bitscore",
        "best_evalue", "best_accession",
    ]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    excluded_path = args.out + ".excluded_low_confidence.tsv"
    if excluded:
        with open(excluded_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=excluded[0].keys(), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(excluded)

    print(f"{len(top_hits)} sample(s) with passing hits.")
    print(f"Wrote {len(out_rows)} sample-species rows -> {args.out}")
    if excluded:
        print(f"{len(excluded)} low-confidence hit(s) excluded (below bitscore={args.min_bitscore} or above evalue={args.max_evalue}) -> {excluded_path}")
    if samples_missing_depth:
        print(f"WARNING: no depth_its_mean.tsv found/empty for: {', '.join(samples_missing_depth)}"
              " (species calls included with depth=0; check these samples' qc/ output)")


if __name__ == "__main__":
    main()
