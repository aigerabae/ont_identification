#!/usr/bin/env python3
"""
extract_read_counts_for_report.py

Produces the exact table format required: per sample, per taxon, LITERAL
read counts (not mean depth) and the count of reads that could not be
classified -- using data your pipeline already has on disk.

Read counts per taxon: uses `samtools idxstats` on each sample's
mapped_its.sorted.bam, which reports the number of mapped reads PER
REFERENCE CONTIG directly (not a depth average). Reads for every contig
assigned to the same taxon (from lca_assign.py's output) are summed.

Unclassified reads: (reads that entered ITS-region mapping) minus (reads
that mapped to any ITS-flagged contig at all) = reads with no taxonomic
signal at this locus. Read from each sample's flagstat_its.txt.

Usage:
    python3 extract_read_counts_for_report.py \
        --pipeline-out pipeline_out \
        --lca-summary blast_results/lca_summary.tsv \
        --out read_counts_report.tsv
"""
import argparse
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict

_COMBINED_SUFFIX_RE = re.compile(r"_ITS1_ITS2_combined$")


def strip_combined_suffix(contig: str) -> str:
    return _COMBINED_SUFFIX_RE.sub("", contig)


def get_idxstats_counts(bam_path: str) -> dict[str, int]:
    """Returns {contig_name: mapped_read_count} via `samtools idxstats`."""
    if not os.path.exists(bam_path):
        print(f"WARNING: BAM not found (skipping, treated as 0 reads): {bam_path}", file=sys.stderr)
        return {}
    if os.path.getsize(bam_path) == 0:
        print(f"NOTE: empty placeholder BAM (sample had no ITS-flagged contigs): {bam_path}", file=sys.stderr)
        return {}
    try:
        result = subprocess.run(
            ["samtools", "idxstats", bam_path],
            capture_output=True, text=True, check=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"WARNING: samtools idxstats failed for {bam_path}: {e}", file=sys.stderr)
        return {}

    counts = {}
    for line in result.stdout.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) != 4 or parts[0] == "*":
            continue
        contig, _length, mapped, _unmapped = parts
        counts[contig] = int(mapped)
    return counts


def get_total_and_mapped_from_flagstat(flagstat_path: str) -> tuple[int, int]:
    """Returns (total_reads, mapped_reads) parsed from samtools flagstat output."""
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pipeline-out", required=True)
    ap.add_argument("--lca-summary", required=True, help="Path to lca_assign.py's per-contig output")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # load contig -> taxon assignment per sample from lca_summary.tsv
    contig_taxon: dict[str, dict[str, tuple[str, str]]] = defaultdict(dict)  # sample -> contig -> (name, rank)
    with open(args.lca_summary, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            depth_key = strip_combined_suffix(row["contig"])
            contig_taxon[row["sample"]][depth_key] = (row["assigned_name"], row["assigned_rank"])

    out_rows = []
    n_bam_missing = 0
    n_bam_empty = 0
    n_bam_ok = 0
    for sample, contigs in sorted(contig_taxon.items()):
        bam_path = os.path.join(args.pipeline_out, sample, "qc", "mapped_its.sorted.bam")
        flagstat_path = os.path.join(args.pipeline_out, sample, "qc", "flagstat_its.txt")

        if not os.path.exists(bam_path):
            n_bam_missing += 1
        elif os.path.getsize(bam_path) == 0:
            n_bam_empty += 1
        else:
            n_bam_ok += 1

        idx_counts = get_idxstats_counts(bam_path)
        total_reads, mapped_reads = get_total_and_mapped_from_flagstat(flagstat_path)
        unclassified = max(total_reads - mapped_reads, 0)

        taxon_reads: dict[str, int] = defaultdict(int)
        taxon_rank: dict[str, str] = {}
        for contig, (name, rank) in contigs.items():
            reads = idx_counts.get(contig, 0)
            taxon_reads[name] += reads
            taxon_rank[name] = rank

        total_classified = sum(taxon_reads.values())
        grand_total = total_classified + unclassified

        for taxon, reads in sorted(taxon_reads.items(), key=lambda kv: kv[1], reverse=True):
            pct = (reads / grand_total * 100) if grand_total > 0 else 0.0
            out_rows.append({
                "sample": sample, "taxon": taxon, "read_count": reads,
                "pct_of_total_reads": round(pct, 2), "identification_level": taxon_rank[taxon],
            })

        # unclassified row
        pct_unclassified = (unclassified / grand_total * 100) if grand_total > 0 else 0.0
        out_rows.append({
            "sample": sample, "taxon": "Unclassified", "read_count": unclassified,
            "pct_of_total_reads": round(pct_unclassified, 2), "identification_level": "-",
        })

    out_cols = ["sample", "taxon", "read_count", "pct_of_total_reads", "identification_level"]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Wrote {len(out_rows)} rows -> {args.out}")
    print(f"\nBAM file status across {len(contig_taxon)} sample(s) in lca_summary.tsv:")
    print(f"  {n_bam_ok} sample(s): real BAM found, read counts extracted normally")
    print(f"  {n_bam_empty} sample(s): empty placeholder BAM (expected -- no ITS-flagged contigs for these)")
    print(f"  {n_bam_missing} sample(s): BAM file MISSING ENTIRELY -- if this number is high, something is")
    print(f"                wrong with --pipeline-out path or these samples never completed map_its_only.")
    print("NOTE: 'Unclassified' = reads that did not map to any ITS-flagged contig for that sample")
    print("      (computed as flagstat total minus flagstat mapped). This does NOT include reads")
    print("      that mapped to a contig later excluded as taxonomically uninformative")
    print("      (kingdom/phylum-level LCA calls) -- check final_report_lca.tsv.excluded_uninformative.tsv")
    print("      separately if you want those folded into 'unclassified' as well.")


if __name__ == "__main__":
    main()
