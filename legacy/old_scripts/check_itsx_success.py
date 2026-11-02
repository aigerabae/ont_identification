#!/usr/bin/env python3
"""
check_itsx_success.py

Scans every sample's itsx/its.summary.txt and reports, per sample:
  - did ITSx run to completion at all (summary file exists and ends with
    "ITSx run finished")
  - how many contigs it searched (input sequences)
  - how many were actually flagged as containing ITS sequence
  - detection rate (ITS-flagged / input)

This distinguishes three different situations that otherwise look similar
from the outside: (1) ITSx never ran or crashed for a sample, (2) ITSx ran
fine but found zero ITS content (e.g. because the assembly was too
fragmented -- see plot_read_counts.py / check_assembly_yield.sh for why),
and (3) ITSx ran and found ITS content normally.

Usage:
    python3 check_itsx_success.py --pipeline-out pipeline_out
    python3 check_itsx_success.py --pipeline-out pipeline_out --sort-by rate
"""
import argparse
import glob
import os
import re
import sys

INPUT_RE = re.compile(r"Number of sequences in input file:\s*(\d+)")
DETECTED_RE = re.compile(r"Sequences detected as ITS by ITSx:\s*(\d+)")
FINISHED_RE = re.compile(r"ITSx run finished")


def check_sample(summary_path: str) -> dict:
    if not os.path.exists(summary_path):
        return {"ran": False, "finished": False, "n_input": None, "n_detected": None}

    with open(summary_path) as fh:
        text = fh.read()

    finished = bool(FINISHED_RE.search(text))
    input_match = INPUT_RE.search(text)
    detected_match = DETECTED_RE.search(text)

    return {
        "ran": True,
        "finished": finished,
        "n_input": int(input_match.group(1)) if input_match else None,
        "n_detected": int(detected_match.group(1)) if detected_match else None,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pipeline-out", required=True)
    ap.add_argument("--sort-by", choices=["sample", "rate", "n_detected"], default="rate",
                     help="Sort order for the table (default: rate, lowest first)")
    args = ap.parse_args()

    sample_dirs = sorted(glob.glob(os.path.join(args.pipeline_out, "*")))
    if not sample_dirs:
        sys.exit(f"No sample directories found under {args.pipeline_out}/")

    rows = []
    for sdir in sample_dirs:
        if not os.path.isdir(sdir):
            continue
        sample = os.path.basename(sdir)
        summary_path = os.path.join(sdir, "itsx", "its.summary.txt")
        result = check_sample(summary_path)
        rate = (result["n_detected"] / result["n_input"] * 100
                if result["n_input"] else None)
        rows.append({"sample": sample, "rate": rate, **result})

    def sort_key(r):
        if args.sort_by == "sample":
            return r["sample"]
        if args.sort_by == "n_detected":
            return r["n_detected"] if r["n_detected"] is not None else -1
        # rate: put never-ran/crashed samples first (most concerning), then by rate ascending
        if not r["ran"] or not r["finished"] or r["rate"] is None:
            return -1
        return r["rate"]

    rows.sort(key=sort_key)

    print(f"{'sample':<20}{'ITSx ran?':<12}{'input_seqs':>12}{'its_detected':>14}{'rate':>8}")
    print("-" * 68)
    n_never_ran, n_zero_detected, n_normal = 0, 0, 0
    for r in rows:
        if not r["ran"]:
            print(f"{r['sample']:<20}{'NO (missing)':<12}{'-':>12}{'-':>14}{'-':>8}")
            n_never_ran += 1
        elif not r["finished"]:
            print(f"{r['sample']:<20}{'CRASHED?':<12}{'-':>12}{'-':>14}{'-':>8}")
            n_never_ran += 1
        else:
            rate_str = f"{r['rate']:.1f}%" if r["rate"] is not None else "n/a"
            flag = "  <-- ZERO ITS FOUND" if r["n_detected"] == 0 else ""
            if r["n_detected"] == 0:
                n_zero_detected += 1
            else:
                n_normal += 1
            print(f"{r['sample']:<20}{'yes':<12}{r['n_input']:>12}{r['n_detected']:>14}{rate_str:>8}{flag}")

    print(f"\n{n_never_ran} sample(s) where ITSx never ran or didn't finish")
    print(f"{n_zero_detected} sample(s) where ITSx ran fine but found zero ITS-flagged contigs")
    print(f"{n_normal} sample(s) with normal ITS detection")


if __name__ == "__main__":
    main()
