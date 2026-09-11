#!/usr/bin/env python3
"""
diagnose_failed_samples.py

For the 4 samples that failed assembly (18S-1_S37, 18S-3_S41, 18S-9_S38,
18S-31_S33), extracts the most-duplicated unique sequences from each
sample's post-fastp, pre-dedup reads and writes them as a small FASTA per
sample -- small enough to paste into NCBI web BLAST by hand, or run
through blastn locally, as a quick check of what's actually in these
samples before investing in a full reanalysis.

Uses merged.fastq.gz preferentially (full-length, highest quality single
reads); falls back to unmerged R1 if a sample's merged file is too sparse
to be informative on its own.

Usage:
    python3 diagnose_failed_samples.py \
        --pipeline-out pipeline_out \
        --samples 18S-1_S37 18S-3_S41 18S-9_S38 18S-31_S33 \
        --top-n 20 \
        --out-dir failed_samples_diagnostic
"""
import argparse
import gzip
import os
from collections import Counter


def read_fastq_seqs(path):
    """Yields sequences from a (possibly gzipped) FASTQ file."""
    if not os.path.exists(path):
        return
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # sequence line
                yield line.strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipeline-out", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--top-n", type=int, default=20)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    for sample in args.samples:
        qc_dir = os.path.join(args.pipeline_out, sample, "qc_reads")
        merged_path = os.path.join(qc_dir, "merged.fastq.gz")
        r1_path = os.path.join(qc_dir, "unmerged_R1.fastq.gz")

        counts = Counter(read_fastq_seqs(merged_path))
        source = "merged"
        total_seqs = sum(counts.values())

        # If merged reads are too sparse to be informative, fall back to R1
        if total_seqs < 5:
            counts = Counter(read_fastq_seqs(r1_path))
            source = "unmerged_R1"
            total_seqs = sum(counts.values())

        if total_seqs == 0:
            print(f"WARNING [{sample}]: no reads found in {qc_dir} -- check the path.")
            continue

        top = counts.most_common(args.top_n)
        unique_seqs = len(counts)

        out_path = os.path.join(args.out_dir, f"{sample}_top_duplicated.fasta")
        with open(out_path, "w") as f:
            for i, (seq, count) in enumerate(top, 1):
                f.write(f">{sample}_{source}_rank{i}_count{count}\n{seq}\n")

        print(f"{sample}: {total_seqs} reads ({source}), {unique_seqs} unique sequences, "
              f"top sequence duplicated {top[0][1]}x -- wrote {out_path}")


if __name__ == "__main__":
    main()
