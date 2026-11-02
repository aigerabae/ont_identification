#!/usr/bin/env python3
"""
pivot_genus_matrix.py

Turns genus_summary.tsv (long format: one row per sample-genus) into the
actual final deliverable: one row per SAMPLE, one column per genus, values
= % of that sample's ITS-region read depth. This is the per-sample
biodiversity + proportions table that was the original goal of the whole
pipeline.

Genera absent from a sample show 0.0, not a blank — so every row's numeric
columns are directly comparable and (roughly) sum to 100% per sample,
modulo whatever was excluded by the quality filter in build_final_report.py
(check that sample's row in final_report.tsv.excluded_low_confidence.tsv if
a row's total is noticeably below 100%).

Usage:
    python3 pivot_genus_matrix.py --genus-summary genus_summary.tsv --out final_matrix.tsv
"""
import argparse
import csv
from collections import defaultdict


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genus-summary", required=True, help="Path to genus_summary.tsv")
    ap.add_argument("--out", required=True, help="Path to write the pivoted matrix TSV")
    args = ap.parse_args()

    pct_by_sample_genus: dict[str, dict[str, float]] = defaultdict(dict)
    all_genera: set = set()
    all_samples: list = []
    seen_samples: set = set()

    with open(args.genus_summary, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            genus = row["genus"]
            pct = float(row["pct_of_sample_its_depth"])
            pct_by_sample_genus[sample][genus] = pct
            all_genera.add(genus)
            if sample not in seen_samples:
                seen_samples.add(sample)
                all_samples.append(sample)

    # order genera by overall prevalence (most samples first), then alphabetically
    genus_prevalence = defaultdict(int)
    for sample_data in pct_by_sample_genus.values():
        for genus in sample_data:
            genus_prevalence[genus] += 1
    ordered_genera = sorted(all_genera, key=lambda g: (-genus_prevalence[g], g))

    out_cols = ["sample"] + ordered_genera + ["row_total_pct"]
    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(out_cols)
        for sample in all_samples:
            row_data = pct_by_sample_genus[sample]
            values = [round(row_data.get(g, 0.0), 3) for g in ordered_genera]
            row_total = round(sum(values), 3)
            writer.writerow([sample] + values + [row_total])

    print(f"Wrote {len(all_samples)} sample(s) x {len(ordered_genera)} genera -> {args.out}")
    print("row_total_pct should be ~100 for each sample; a lower total means some of")
    print("that sample's ITS depth came from hits excluded by the quality filter")
    print("(check final_report.tsv.excluded_low_confidence.tsv for that sample).")


if __name__ == "__main__":
    main()
