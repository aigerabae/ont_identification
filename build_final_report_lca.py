#!/usr/bin/env python3
"""
build_final_report_lca.py

Same purpose as build_final_report.py (join species/taxon calls with
per-contig depth to get relative abundance per sample), but consumes
lca_assign.py's LCA-resolved output instead of summarize_blast_hits.py's
top-hit-only output. Where a contig's top hits were genuinely tied across
species (or genera), this reports the taxon at the RANK actually supported
by the evidence (e.g. "Centaureinae (subtribe)") rather than an arbitrary
species name — see lca_assign.py.

Usage:
    python3 build_final_report_lca.py \
        --lca-summary blast_results/lca_summary.tsv \
        --pipeline-out pipeline_out \
        --out final_report_lca.tsv
"""
import argparse
import csv
import os
import re
from collections import defaultdict

_COMBINED_SUFFIX_RE = re.compile(r"_ITS1_ITS2_combined$")

# Coarse -> fine numeric ordering of taxonomic ranks. Ranks not listed
# (including "no rank" and "clade", which can sit almost anywhere in the
# tree and are frequently the product of a very broad, low-information
# tie) default to -1, i.e. excluded by --min-informative-rank's default.
RANK_LEVELS = {
    "domain": 0, "superkingdom": 0,
    "kingdom": 1,
    "phylum": 2, "subphylum": 2,
    "class": 3, "subclass": 3,
    "order": 4, "suborder": 4,
    "family": 5, "subfamily": 5,
    "tribe": 6, "subtribe": 6,
    "genus": 7, "subgenus": 7,
    "species": 8, "subspecies": 8,
}


def strip_combined_suffix(contig: str) -> str:
    return _COMBINED_SUFFIX_RE.sub("", contig)


def load_depth_file(path: str) -> dict[str, float]:
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
    ap.add_argument("--lca-summary", required=True, help="Path to lca_assign.py's output")
    ap.add_argument("--pipeline-out", required=True, help="Path to the Snakemake pipeline_out directory")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-informative-rank", default="family",
                     help="Calls resolved BROADER than this rank (e.g. kingdom, phylum, 'no rank', 'clade') "
                          "are excluded from the composition report and written to a separate audit file "
                          "instead, since they carry little to no real taxonomic information. Default: family.")
    args = ap.parse_args()

    if args.min_informative_rank not in RANK_LEVELS:
        valid = ", ".join(sorted(RANK_LEVELS, key=lambda r: RANK_LEVELS[r]))
        raise SystemExit(f"--min-informative-rank must be one of: {valid}")
    min_level = RANK_LEVELS[args.min_informative_rank]

    calls_by_sample: dict[str, list[dict]] = defaultdict(list)
    excluded_uninformative: list[dict] = []
    with open(args.lca_summary, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            row_level = RANK_LEVELS.get(row["assigned_rank"], -1)
            if row_level < min_level:
                excluded_uninformative.append(row)
                continue
            calls_by_sample[row["sample"]].append(row)

    out_rows = []
    samples_missing_depth = []

    for sample, calls in sorted(calls_by_sample.items()):
        depth_path = os.path.join(args.pipeline_out, sample, "qc", "depth_its_mean.tsv")
        depths = load_depth_file(depth_path)
        if not depths:
            samples_missing_depth.append(sample)

        taxon_depth: dict[str, float] = defaultdict(float)
        taxon_contigs: dict[str, list[str]] = defaultdict(list)
        taxon_rank: dict[str, str] = {}
        taxon_tie_flag: dict[str, bool] = {}

        for call in calls:
            # taxon label carries its rank so e.g. genus-level and species-level
            # calls of a similarly-named taxon never silently merge
            taxon_label = f"{call['assigned_name']} ({call['assigned_rank']})"
            depth_key = strip_combined_suffix(call["contig"])
            depth = depths.get(depth_key, 0.0)

            taxon_depth[taxon_label] += depth
            taxon_contigs[taxon_label].append(depth_key)
            taxon_rank[taxon_label] = call["assigned_rank"]
            taxon_tie_flag[taxon_label] = taxon_tie_flag.get(taxon_label, False) or (int(call["tie_group_size"]) > 1)

        total_depth = sum(taxon_depth.values())
        for taxon, depth in sorted(taxon_depth.items(), key=lambda kv: kv[1], reverse=True):
            pct = (depth / total_depth * 100) if total_depth > 0 else 0.0
            out_rows.append({
                "sample": sample,
                "taxon": taxon,
                "rank": taxon_rank[taxon],
                "was_tie_resolved": "yes" if taxon_tie_flag[taxon] else "no",
                "n_contigs": len(taxon_contigs[taxon]),
                "contigs": ";".join(taxon_contigs[taxon]),
                "summed_depth": round(depth, 2),
                "pct_of_sample_its_depth": round(pct, 3),
            })

    out_cols = [
        "sample", "taxon", "rank", "was_tie_resolved", "n_contigs",
        "contigs", "summed_depth", "pct_of_sample_its_depth",
    ]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    n_tie_resolved = sum(1 for r in out_rows if r["was_tie_resolved"] == "yes")
    print(f"Wrote {len(out_rows)} sample-taxon rows -> {args.out}")
    print(f"{n_tie_resolved} row(s) were tie-resolved (rank backed off above species due to a genuine BLAST tie) -- these are flagged, not hidden.")
    if samples_missing_depth:
        print(f"WARNING: no depth data for: {', '.join(samples_missing_depth)}")

    if excluded_uninformative:
        excluded_path = args.out + ".excluded_uninformative.tsv"
        with open(excluded_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=excluded_uninformative[0].keys(), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(excluded_uninformative)
        print(f"{len(excluded_uninformative)} contig(s) resolved BROADER than '{args.min_informative_rank}' "
              f"(e.g. kingdom/phylum/domain/clade/no-rank) were excluded from the composition report "
              f"and written to {excluded_path} for audit -- they were NOT counted toward any sample's "
              f"depth total or percentages.")


if __name__ == "__main__":
    main()
