#!/usr/bin/env python3
"""
plot_read_counts.py

Scans every sample's fastp.json (from the Snakemake pipeline) and produces
two plots:
  1. A sorted bar chart of total reads per sample -- makes low-yield
     outliers (like the four samples that crashed SPAdes) immediately
     visible, especially against everything else in the batch.
  2. A histogram of the read-count distribution across all samples -- shows
     whether your batch is roughly uniform or has distinct clusters (e.g.
     a low-yield cluster vs. a normal-yield cluster).

Also prints a plain-text sorted table, and flags samples below
--crash-threshold (default 5000, chosen from the observed pattern where
assemblies below ~1000 reads crashed SPAdes outright) directly in the
console output and with a red bar/line in both plots.

Usage:
    python3 plot_read_counts.py --pipeline-out pipeline_out --out read_counts.png
    python3 plot_read_counts.py --pipeline-out pipeline_out --out read_counts.png --crash-threshold 1000
"""
import argparse
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")  # no display needed, just save to file
import matplotlib.pyplot as plt


def load_read_counts(pipeline_out: str) -> dict[str, dict]:
    """Returns {sample: {"raw": int, "after_fastp": int, "duplication_pct": float}}."""
    counts: dict[str, dict] = {}
    pattern = os.path.join(pipeline_out, "*", "qc_reads", "fastp.json")
    for path in sorted(glob.glob(pattern)):
        sample = path.split(os.sep)[-3]
        try:
            with open(path) as fh:
                data = json.load(fh)
            summary = data["summary"]
            dup_pct = data.get("duplication", {}).get("rate", None)
            counts[sample] = {
                "raw": summary["before_filtering"]["total_reads"],
                "after_fastp": summary["after_filtering"]["total_reads"],
                "duplication_pct": (dup_pct * 100) if dup_pct is not None else None,
            }
        except (json.JSONDecodeError, KeyError) as e:
            print(f"WARNING: could not read expected fields from {path}: {e}", file=sys.stderr)
    return counts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pipeline-out", required=True)
    ap.add_argument("--out", required=True, help="Output image path, e.g. read_counts.png")
    ap.add_argument("--crash-threshold", type=int, default=5000,
                     help="Samples at or below this read count are flagged/colored red (default 5000)")
    args = ap.parse_args()

    counts = load_read_counts(args.pipeline_out)
    if not counts:
        sys.exit(f"No fastp.json files found under {args.pipeline_out}/*/qc_reads/")

    # sort ascending by RAW read count -- lowest-yield samples first
    items = sorted(counts.items(), key=lambda kv: kv[1]["raw"])
    samples = [s for s, _ in items]
    values = [d["raw"] for _, d in items]  # plots still use raw counts
    colors = ["#d62728" if v <= args.crash_threshold else "#2ca02c" for v in values]

    print(f"{'sample':<20}{'raw_reads':>12}{'after_fastp':>14}{'survival%':>11}{'dup%':>8}")
    print("-" * 65)
    low_yield = []
    for s, d in items:
        raw, after = d["raw"], d["after_fastp"]
        survival = (after / raw * 100) if raw else 0
        dup = f"{d['duplication_pct']:.1f}" if d["duplication_pct"] is not None else "n/a"
        flag = "  <-- LOW YIELD" if raw <= args.crash_threshold else ""
        if flag:
            low_yield.append(s)
        print(f"{s:<20}{raw:>12,}{after:>14,}{survival:>10.1f}%{dup:>8}{flag}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(max(10, len(samples) * 0.4), 5))

    # Plot 1: sorted bar chart, log scale (read counts span orders of magnitude)
    ax1.bar(range(len(samples)), values, color=colors)
    ax1.set_xticks(range(len(samples)))
    ax1.set_xticklabels(samples, rotation=90, fontsize=7)
    ax1.set_yscale("log")
    ax1.set_ylabel("Raw reads, before fastp (log scale)")
    ax1.set_title("Raw reads per sample (sorted)")
    ax1.axhline(args.crash_threshold, color="gray", linestyle="--", linewidth=1,
                label=f"threshold ({args.crash_threshold:,})")
    ax1.legend(fontsize=8)

    # Plot 2: distribution histogram
    ax2.hist(values, bins=15, color="#1f77b4", edgecolor="black")
    ax2.axvline(args.crash_threshold, color="gray", linestyle="--", linewidth=1)
    ax2.set_xlabel("Total reads")
    ax2.set_ylabel("Number of samples")
    ax2.set_title("Distribution of read counts across samples")

    plt.tight_layout()
    plt.savefig(args.out, dpi=150)
    print(f"\nSaved plot -> {args.out}")

    if low_yield:
        print(f"\n{len(low_yield)} sample(s) at or below the {args.crash_threshold:,}-read threshold:")
        for s in low_yield:
            print(f"  - {s}")


if __name__ == "__main__":
    main()
