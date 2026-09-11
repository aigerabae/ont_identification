#!/usr/bin/env python3
"""
export_unique_reads_fasta.py

Dereplicates ALL post-fastp reads (not just the top-N most duplicated) for
each failed-assembly sample and writes one FASTA file per sample, ready to
paste directly into the NCBI BLAST web form (blast.ncbi.nlm.nih.gov -> 
blastn -> core_nt).

Each sequence's FASTA header encodes its sample and duplication count, so
the companion script (parse_ncbi_text_results.py) can recover exact
read-weighted percentages from the downloaded results later --
no separate lookup table needed.

Sequences are written most-duplicated-first, in case a sample's unique
sequence count is large enough that you only want to submit the top
portion of it manually.

Usage:
    python3 export_unique_reads_fasta.py \
        --pipeline-out pipeline_out \
        --samples 18S-1_S37 18S-3_S41 18S-9_S38 18S-31_S33 \
        --out-dir failed_samples_diagnostic
"""
import argparse
import gzip
import os
from collections import Counter


def read_fastq_seqs(path):
    if not os.path.exists(path):
        return
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                yield line.strip()


def dereplicate_all(pipeline_out, sample):
    qc_dir = os.path.join(pipeline_out, sample, "qc_reads")
    merged_path = os.path.join(qc_dir, "merged.fastq.gz")
    r1_path = os.path.join(qc_dir, "unmerged_R1.fastq.gz")

    counts = Counter(read_fastq_seqs(merged_path))
    source = "merged"
    total = sum(counts.values())

    if total < 5:
        counts = Counter(read_fastq_seqs(r1_path))
        source = "unmerged_R1"
        total = sum(counts.values())

    return counts, source, total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipeline-out", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    for sample in args.samples:
        counts, source, total = dereplicate_all(args.pipeline_out, sample)
        if total == 0:
            print(f"WARNING [{sample}]: no reads found, skipping.")
            continue

        out_path = os.path.join(args.out_dir, f"{sample}_all_unique.fasta")
        with open(out_path, "w") as f:
            for i, (seq, count) in enumerate(counts.most_common(), 1):
                # ID format: {sample}__u{rank}__count{N} -- the parser
                # extracts sample name and count directly from this string.
                f.write(f">{sample}__u{i}__count{count}\n{seq}\n")

        print(f"{sample}: {total} reads ({source}), {len(counts)} unique sequences "
              f"-> {out_path}")
        print(f"  Total reads represented in this file: {sum(counts.values())} "
              f"(should equal {total})")

    print("\nPaste each *_all_unique.fasta file's contents into "
          "https://blast.ncbi.nlm.nih.gov/Blast.cgi (blastn, core_nt database), "
          "one sample per submission. After each job finishes, use "
          "'Download All' -> 'Text' to save the full results as a single "
          ".txt file per sample -- that's what parse_ncbi_text_results.py reads.")


if __name__ == "__main__":
    main()
