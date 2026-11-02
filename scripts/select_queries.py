#!/usr/bin/env python3
"""
select_queries.py -- decide, per sample, what represents the assembly-based
evidence. Three modes (written to mode.txt, carried into the report):

  its      ITSx found ITS regions      -> query = ITS regions, map reads to those contigs
  contigs  assembly OK but no ITS      -> query = top-N longest contigs (lower confidence)
  none     no usable assembly          -> nothing here; every read goes to read-level
                                          classification (residual_reads.py)

Always creates queries.fasta, mapping_ref.fasta, mode.txt (possibly empty).
"""
import argparse, os, sys


def read_fasta(path):
    recs, h, b = [], None, []
    if os.path.exists(path) and os.path.getsize(path) > 0:
        for line in open(path):
            line = line.rstrip("\n")
            if line.startswith(">"):
                if h is not None:
                    recs.append((h, "".join(b)))
                h, b = line[1:].split()[0], []
            elif h is not None:
                b.append(line.strip())
        if h is not None:
            recs.append((h, "".join(b)))
    return recs


def write(path, recs):
    with open(path, "w") as fh:
        for h, s in recs:
            fh.write(f">{h}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")


ap = argparse.ArgumentParser()
ap.add_argument("--its-regions", required=True)
ap.add_argument("--its-contigs", required=True)
ap.add_argument("--contigs", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--top-n", type=int, default=10)
ap.add_argument("--min-len", type=int, default=100)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

regions = read_fasta(a.its_regions)
if regions:
    mode = "its"
    q, ref = regions, read_fasta(a.its_contigs)
else:
    c = [r for r in read_fasta(a.contigs) if len(r[1]) >= a.min_len]
    c.sort(key=lambda r: len(r[1]), reverse=True)
    c = c[: a.top_n]
    if c:
        mode, q, ref = "contigs", c, c
    else:
        mode, q, ref = "none", [], []
write(os.path.join(a.outdir, "queries.fasta"), q)
write(os.path.join(a.outdir, "mapping_ref.fasta"), ref)
open(os.path.join(a.outdir, "mode.txt"), "w").write(mode + "\n")
print(f"mode={mode} n_queries={len(q)}", file=sys.stderr)
