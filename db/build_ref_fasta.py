#!/usr/bin/env python3
"""
build_ref_fasta.py

Turns ITSx output on downloaded GenBank records into the final reference set:
  * one ITS reference per accession, with the SAME priority as
    extract_its_regions.py:  full (ITS1-5.8S-ITS2) > ITS1+ITS2 concatenated
    > ITS1 only > ITS2 only
  * exact-duplicate sequences removed ONLY within the same taxid (identical
    sequences shared between different taxa are kept on purpose: they are
    what makes the LCA step report genus instead of a false species)
  * short / N-rich records dropped
  * an id -> taxid map in makeblastdb -taxid_map format
  * a manifest (id, accession, taxid, region, length, source group)

Usage:
  build_ref_fasta.py --itsx-prefix build/itsx_out --accession-map build/accession_taxid_group.tsv \
                     --out-fasta build/refs.fasta --out-taxid-map build/refs.taxid.map \
                     --out-manifest build/refs_manifest.tsv [--min-len 100]
"""
import argparse
import re
import sys
from collections import Counter
from pathlib import Path

ACC_RE = re.compile(r"([A-Z]{1,6}_?\d+\.\d+)")


def read_fasta(path):
    seqs = {}
    p = Path(path)
    if not p.exists():
        return seqs
    hdr, buf = None, []
    for line in p.open():
        line = line.rstrip("\n")
        if line.startswith(">"):
            if hdr is not None:
                seqs[hdr] = "".join(buf)
            m = ACC_RE.match(line[1:].strip())
            hdr = m.group(1) if m else None
            buf = []
        elif hdr is not None:
            buf.append(line.strip())
    if hdr is not None:
        seqs[hdr] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--itsx-prefix", required=True)
    ap.add_argument("--accession-map", required=True,
                    help="TSV: accession.version <TAB> taxid <TAB> source_group")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-taxid-map", required=True)
    ap.add_argument("--out-manifest", required=True)
    ap.add_argument("--min-len", type=int, default=100)
    ap.add_argument("--max-n-frac", type=float, default=0.05)
    a = ap.parse_args()

    acc_info = {}
    for line in open(a.accession_map):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3 and f[0] not in acc_info:      # first group wins
            acc_info[f[0]] = (f[1], f[2])

    pre = a.itsx_prefix
    full = read_fasta(f"{pre}.full.fasta")
    its1 = read_fasta(f"{pre}.ITS1.fasta")
    its2 = read_fasta(f"{pre}.ITS2.fasta")
    candidates = set(full) | set(its1) | set(its2)

    seen = set()
    stats = Counter()
    with open(a.out_fasta, "w") as fa, open(a.out_taxid_map, "w") as tm, \
            open(a.out_manifest, "w") as mf:
        mf.write("ref_id\taccession\ttaxid\tregion\tlength\tsource_group\n")
        for acc in sorted(candidates):
            if acc not in acc_info:
                stats["no_taxid_dropped"] += 1
                continue
            taxid, group = acc_info[acc]
            if acc in full:
                region, seq = "full", full[acc]
            elif acc in its1 and acc in its2:
                region, seq = "its12", its1[acc] + its2[acc]
            elif acc in its1:
                region, seq = "its1", its1[acc]
            else:
                region, seq = "its2", its2[acc]
            seq = seq.upper().replace("U", "T")
            if len(seq) < a.min_len:
                stats["too_short_dropped"] += 1
                continue
            if seq.count("N") / len(seq) > a.max_n_frac:
                stats["n_rich_dropped"] += 1
                continue
            key = (taxid, seq)
            if key in seen:
                stats["dup_same_taxid_dropped"] += 1
                continue
            seen.add(key)
            rid = f"{acc}_{region}"
            fa.write(f">{rid}\n")
            for i in range(0, len(seq), 80):
                fa.write(seq[i:i + 80] + "\n")
            tm.write(f"{rid}\t{taxid}\n")
            mf.write(f"{rid}\t{acc}\t{taxid}\t{region}\t{len(seq)}\t{group}\n")
            stats[f"kept_{region}"] += 1
            stats[f"kept_group:{group}"] += 1

    stats["records_downloaded"] = len(acc_info)
    stats["records_with_ITS_detected"] = len(candidates)
    print("== build_ref_fasta summary ==", file=sys.stderr)
    for k in sorted(stats):
        print(f"  {k}\t{stats[k]}", file=sys.stderr)


if __name__ == "__main__":
    main()
