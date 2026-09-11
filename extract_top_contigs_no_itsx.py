#!/usr/bin/env python3
"""
extract_top_contigs_no_itsx.py

For samples where ITSx found zero ITS-flagged contigs (assembly worked,
but contigs were too short/fragmented for ITSx's HMM boundaries to
anchor), this extracts the longest available contigs directly from the
chimera-filtered assembly and prepares them for BLAST -- skipping ITSx
entirely. Resolution will be coarser (these often won't isolate a clean
ITS-only region the way extract_its_regions.py does), but paired with
LCA-based assignment, ambiguous/conserved contigs will honestly back off
to genus/family/order rather than produce a false species claim.

Usage:
    python3 extract_top_contigs_no_itsx.py \
        --pipeline-out pipeline_out \
        --samples 18S-13_S39 18S-17_S40 18S-23_S45 18S-25_S46 18S-32_S44 18S-34_S43 18S-37_S42 18S-4_S47 \
        --out-dir no_itsx_fallback \
        --top-n 10 --min-length 100
"""
import argparse
import os
import sys


def read_fasta(path: str) -> list[tuple[str, str]]:
    records = []
    header = None
    chunks: list[str] = []
    if not os.path.exists(path):
        return records
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            records.append((header, "".join(chunks)))
    return records


def write_fasta(path: str, records: list[tuple[str, str]], wrap: int = 80) -> None:
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), wrap):
                fh.write(seq[i : i + wrap] + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pipeline-out", required=True)
    ap.add_argument("--samples", required=True, nargs="+", help="Sample names to process")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--top-n", type=int, default=10, help="Max contigs per sample (default 10; fewer if the sample has less)")
    ap.add_argument("--min-length", type=int, default=100, help="Skip contigs shorter than this (default 100bp)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    for sample in args.samples:
        fasta_path = os.path.join(args.pipeline_out, sample, "chimera", "combined_contigs_nonchimeric.fasta")
        records = read_fasta(fasta_path)
        if not records:
            print(f"[warn] {sample}: no contigs found at {fasta_path} -- skipping", file=sys.stderr)
            continue

        filtered = [(h, s) for h, s in records if len(s) >= args.min_length]
        filtered.sort(key=lambda r: len(r[1]), reverse=True)
        top = filtered[: args.top_n]

        if not top:
            print(f"[warn] {sample}: 0 contigs pass min-length={args.min_length} "
                  f"(had {len(records)} total, longest was {max((len(s) for _, s in records), default=0)}bp) -- skipping", file=sys.stderr)
            continue

        out_path = os.path.join(args.out_dir, f"{sample}_top_contigs.fasta")
        write_fasta(out_path, top)
        lengths = [len(s) for _, s in top]
        print(f"[ok] {sample}: {len(top)} contig(s) written (lengths {min(lengths)}-{max(lengths)}bp) -> {out_path}")


if __name__ == "__main__":
    main()
