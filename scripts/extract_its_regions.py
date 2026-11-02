#!/usr/bin/env python3
"""
extract_its_regions.py -- turn ITSx output into (a) one BLAST query per ITS-flagged
contig and (b) the subset of full contigs to map reads onto.

Per flagged contig the best available region is used:
  full ITS1-5.8S-ITS2  >  ITS1+ITS2 concatenated  >  ITS1 only  >  ITS2 only
Query ids are the contig ids unchanged, so BLAST calls, read counts and
mapping all join on the same name.

Outputs (always created, possibly empty):
  <outdir>/its_regions.fasta   one record per flagged contig
  <outdir>/its_contigs.fasta   those contigs, full length (mapping reference)
  <outdir>/extraction_log.tsv  contig, source, length
"""
import argparse, os, re, sys

ACC = re.compile(r"^(\S+?)(?:\||$)")


def read_fasta(path):
    d, hdr, buf = {}, None, []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return d
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if hdr is not None:
                d[hdr] = "".join(buf)
            hdr = ACC.match(line[1:].strip().split()[0]).group(1)
            buf = []
        elif hdr is not None:
            buf.append(line.strip())
    if hdr is not None:
        d[hdr] = "".join(buf)
    return d


def write_rec(fh, h, s, wrap=80):
    fh.write(f">{h}\n")
    for i in range(0, len(s), wrap):
        fh.write(s[i:i + wrap] + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True, help="ITSx output prefix")
    ap.add_argument("--contigs", required=True, help="contigs given to ITSx")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    full = read_fasta(f"{a.prefix}.full.fasta")
    its1 = read_fasta(f"{a.prefix}.ITS1.fasta")
    its2 = read_fasta(f"{a.prefix}.ITS2.fasta")
    contigs = read_fasta(a.contigs)

    flagged = []
    pos = f"{a.prefix}.positions.txt"
    if os.path.exists(pos):
        for line in open(pos):
            if line.strip():
                flagged.append(line.split("\t")[0].strip())

    log = []
    with open(os.path.join(a.outdir, "its_regions.fasta"), "w") as fr, \
            open(os.path.join(a.outdir, "its_contigs.fasta"), "w") as fc:
        for c in dict.fromkeys(flagged):          # unique, order kept
            if c in full:
                seq, src = full[c], "full"
            elif c in its1 and c in its2:
                seq, src = its1[c] + its2[c], "ITS1+ITS2"
            elif c in its1:
                seq, src = its1[c], "ITS1_only"
            elif c in its2:
                seq, src = its2[c], "ITS2_only"
            else:
                log.append((c, "not_extracted", 0))
                continue
            if c not in contigs:
                log.append((c, "contig_missing", 0))
                continue
            write_rec(fr, c, seq)
            write_rec(fc, c, contigs[c])
            log.append((c, src, len(seq)))
    with open(os.path.join(a.outdir, "extraction_log.tsv"), "w") as fh:
        fh.write("contig\tsource\tlength\n")
        for r in log:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"{sum(1 for r in log if r[2] > 0)}/{len(set(flagged))} flagged contigs extracted", file=sys.stderr)


if __name__ == "__main__":
    main()
