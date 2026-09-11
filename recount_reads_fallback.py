#!/usr/bin/env python3
"""
recount_reads_fallback.py

Computes per-taxon read counts for the 8 no-ITSx fallback samples, using
reads mapped against each sample's top-N-longest contigs
(no_itsx_fallback/{sample}_mapped_fallback.sorted.bam) and the taxon calls
in blast_results_no_itsx/lca_summary_no_itsx.tsv.

IMPORTANT -- these numbers are NOT directly comparable to Table 2's ITS-anchored
read counts. A read counted here maps to a *whole contig* that may be dominated
by conserved 18S/28S flanking sequence (never passed through ITS extraction),
so the taxon it's counted toward reflects that contig's overall BLAST/LCA call,
not a validated ITS-region match. Every output row is explicitly flagged via
the `confidence_note` column so this distinction survives into any downstream
table or report and can't be silently merged with Table 2's numbers.

Usage:
    python3 recount_reads_fallback.py \
        --lca-summary blast_results_no_itsx/lca_summary_no_itsx.tsv \
        --fallback-dir no_itsx_fallback \
        --out read_counts_report_fallback.tsv
"""
import argparse
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict

CONFIDENCE_NOTE = (
    "LOWER CONFIDENCE -- contig-level, no ITS extraction (see report Section 4.5/5.2)"
)


def get_idxstats_counts(bam_path):
    """Returns {contig_name: mapped_read_count} via plain samtools idxstats,
    no MAPQ filtering, consistent with the main read-counting method."""
    if not os.path.exists(bam_path) or os.path.getsize(bam_path) == 0:
        return {}
    result = subprocess.run(
        ["samtools", "idxstats", bam_path],
        capture_output=True, text=True, check=True,
    )
    counts = {}
    for line in result.stdout.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) != 4 or parts[0] == "*":
            continue
        contig, _length, mapped, _unmapped = parts
        counts[contig] = int(mapped)
    return counts


def get_total_and_mapped_from_flagstat(flagstat_path):
    """Returns (total_reads, mapped_reads), same parsing convention used
    for the main ITS-anchored read-counting method."""
    if not os.path.exists(flagstat_path):
        return (0, 0)
    total, mapped = 0, 0
    with open(flagstat_path) as fh:
        for line in fh:
            if " in total " in line or line.strip().endswith("in total"):
                total = int(line.split("+")[0].strip())
            elif re.search(r"\bmapped\b", line) and "primary mapped" not in line and "mate mapped" not in line:
                m = re.match(r"(\d+)\s*\+", line)
                if m and mapped == 0:
                    mapped = int(m.group(1))
    return (total, mapped)


def load_lca_summary(path):
    """Group LCA calls by sample. Schema: sample, contig, assigned_name,
    assigned_rank, tie_group_size, best_bitscore, best_pident."""
    by_sample = defaultdict(list)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            by_sample[row["sample"]].append(row)
    return by_sample


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lca-summary", required=True)
    ap.add_argument("--fallback-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    lca_by_sample = load_lca_summary(args.lca_summary)

    rows_out = []
    samples_missing_bam = []

    for sample, calls in lca_by_sample.items():
        bam_path = os.path.join(args.fallback_dir, f"{sample}_mapped_fallback.sorted.bam")
        flagstat_path = os.path.join(args.fallback_dir, f"{sample}_flagstat_fallback.txt")
        if not os.path.exists(bam_path):
            samples_missing_bam.append(sample)
            continue

        contig_counts = get_idxstats_counts(bam_path)
        total_reads, mapped_reads = get_total_and_mapped_from_flagstat(flagstat_path)

        taxon_reads = defaultdict(int)
        taxon_level = {}
        unmatched_contigs = []
        for call in calls:
            contig_key = call["contig"]
            n = contig_counts.get(contig_key, None)
            if n is None:
                # contig from lca_summary not found in idxstats output --
                # flag rather than silently treat as 0, since a naming
                # mismatch here would silently under-count (exactly the
                # failure mode ruled out for the main samples earlier).
                unmatched_contigs.append(contig_key)
                n = 0
            taxon_label = call["assigned_name"]
            taxon_reads[taxon_label] += n
            taxon_level[taxon_label] = call["assigned_rank"]

        if unmatched_contigs:
            print(f"WARNING [{sample}]: {len(unmatched_contigs)} contig(s) in lca_summary "
                  f"not found in BAM reference -- check contig naming: {unmatched_contigs}",
                  file=sys.stderr)

        unclassified = max(total_reads - mapped_reads, 0)

        for taxon, n in sorted(taxon_reads.items(), key=lambda kv: -kv[1]):
            pct = (n / total_reads * 100) if total_reads > 0 else 0.0
            rows_out.append({
                "sample": sample,
                "taxon": taxon,
                "read_count": n,
                "pct_of_total_reads": round(pct, 2),
                "identification_level": taxon_level[taxon],
                "confidence_note": CONFIDENCE_NOTE,
            })

        pct_unclass = (unclassified / total_reads * 100) if total_reads > 0 else 0.0
        rows_out.append({
            "sample": sample,
            "taxon": "Unclassified",
            "read_count": unclassified,
            "pct_of_total_reads": round(pct_unclass, 2),
            "identification_level": "-",
            "confidence_note": CONFIDENCE_NOTE,
        })

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["sample", "taxon", "read_count", "pct_of_total_reads",
                           "identification_level", "confidence_note"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Wrote {len(rows_out)} rows to {args.out}")
    if samples_missing_bam:
        print(f"WARNING: no fallback BAM found for: {', '.join(samples_missing_bam)}", file=sys.stderr)
        print("Run map_fallback_contigs.sh first.", file=sys.stderr)


if __name__ == "__main__":
    main()
