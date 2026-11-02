#!/usr/bin/env python3
"""
Recompute per-taxon read counts using pre-deduplication reads mapped against
each sample's existing ITS-positive contigs (mapped_its_prededup.sorted.bam),
instead of the deduplicated reads used in the original read_counts_report.tsv.

Rationale: deduplication is appropriate before assembly (removes PCR duplicates
that distort coverage-based graph resolution) but distorts relative abundance
when reused for quantification, since collision probability during dedup rises
with local read depth -- high-abundance taxa lose disproportionately more reads
than rare ones. See conversation notes / report Section 4.4 for full rationale.

Usage:
    python3 recount_reads_prededup.py \
        --lca-summary blast_results/lca_summary.tsv \
        --pipeline-out pipeline_out \
        --build-script scripts/build_final_report_lca.py \
        --mapq 10 \
        --out read_counts_report_prededup.tsv
"""
import argparse
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict

# Matches the exact, already-validated logic in scripts/extract_read_counts_for_report.py
# -- confirmed by reading that script directly, not assumed from the depth rule.
_COMBINED_SUFFIX_RE = re.compile(r"_ITS1_ITS2_combined$")


def strip_combined_suffix(contig: str) -> str:
    return _COMBINED_SUFFIX_RE.sub("", contig)


def get_idxstats_counts(bam_path):
    """Returns {contig_name: mapped_read_count} via plain `samtools idxstats`,
    with NO MAPQ filtering -- matching extract_read_counts_for_report.py exactly.
    (An earlier version of this script incorrectly added a MAPQ>=10 prefilter,
    borrowed from the unrelated depth_its_mean rule -- that was a bug.)"""
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
    """Returns (total_reads, mapped_reads), matching
    extract_read_counts_for_report.py's own flagstat parsing exactly."""
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
    """Group LCA calls by sample. Assumes the same schema used by lca_assign.py's
    output: sample, contig, assigned_name, assigned_rank (confirmed from
    lca_summary_no_itsx.tsv's header in this conversation)."""
    by_sample = defaultdict(list)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            by_sample[row["sample"]].append(row)
    return by_sample


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lca-summary", required=True)
    ap.add_argument("--pipeline-out", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    lca_by_sample = load_lca_summary(args.lca_summary)

    rows_out = []
    samples_missing_bam = []

    for sample, calls in lca_by_sample.items():
        bam_path = os.path.join(args.pipeline_out, sample, "qc", "mapped_its_prededup.sorted.bam")
        flagstat_path = os.path.join(args.pipeline_out, sample, "qc", "flagstat_its_prededup.txt")
        if not os.path.exists(bam_path):
            samples_missing_bam.append(sample)
            continue

        contig_counts = get_idxstats_counts(bam_path)
        total_reads, mapped_reads = get_total_and_mapped_from_flagstat(flagstat_path)

        taxon_reads = defaultdict(int)
        taxon_level = {}
        for call in calls:
            contig_key = strip_combined_suffix(call["contig"])
            n = contig_counts.get(contig_key, 0)
            taxon_label = call["assigned_name"]
            taxon_reads[taxon_label] += n
            taxon_level[taxon_label] = call["assigned_rank"]

        # Unclassified = total - mapped, from flagstat directly, matching
        # extract_read_counts_for_report.py exactly (not total - sum(taxon counts),
        # though the two should agree closely if the contig join is complete).
        unclassified = max(total_reads - mapped_reads, 0)

        for taxon, n in sorted(taxon_reads.items(), key=lambda kv: -kv[1]):
            pct = (n / total_reads * 100) if total_reads > 0 else 0.0
            rows_out.append({
                "sample": sample,
                "taxon": taxon,
                "read_count": n,
                "pct_of_total_reads": round(pct, 2),
                "identification_level": taxon_level[taxon],
            })

        pct_unclass = (unclassified / total_reads * 100) if total_reads > 0 else 0.0
        rows_out.append({
            "sample": sample,
            "taxon": "Unclassified",
            "read_count": unclassified,
            "pct_of_total_reads": round(pct_unclass, 2),
            "identification_level": "-",
        })

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["sample", "taxon", "read_count", "pct_of_total_reads", "identification_level"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Wrote {len(rows_out)} rows to {args.out}")
    if samples_missing_bam:
        print(f"WARNING: no pre-dedup BAM found for: {', '.join(samples_missing_bam)}", file=sys.stderr)
        print("Run the map_its_only_prededup Snakemake rule for these samples first.", file=sys.stderr)


if __name__ == "__main__":
    main()
