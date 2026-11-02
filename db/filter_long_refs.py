#!/usr/bin/env python3
"""filter_long_refs.py -- drop reference records whose length is impossible for the ITS region they claim.

Why: a few GenBank records keep 18S/28S flanks after ITSx trimming (e.g. AF250821.1_its12, 2108 bp).
Plant 18S/28S reads and contigs then match those flanks at 78-90 % identity and are reported as fungi.
In the 2026-10-28 build, 3 such records attracted ~75 000 reads in the 32 project samples.

Rule (record id suffix comes from build_ref_fasta.py: _full, _its12, _its1, _its2):
  * single-spacer records (_its1, _its2)      longer than --max-spacer bp (default 450)  -> dropped
  * ITS1-5.8S-ITS2 records (_full, _its12)    longer than --max-its bp    (default 1000) -> dropped
  * records whose lineage contains an exempt taxid are always kept (default: Metazoa 33208 and
    gymnosperms 1437180 -- animal decoy ITS and Ginkgo/conifer ITS are genuinely long).

Output order and format are identical to the input (records are only removed), so the result is
deterministic: the same input + taxdump gives the same md5.

Usage:
  python3 db/filter_long_refs.py --in-fasta refs.unfiltered.fasta --in-map refs.unfiltered.taxid.map \
      --taxdump-dir db/taxdump --out-fasta refs.fasta --out-map refs.taxid.map \
      --out-dropped dropped_records.tsv [--max-spacer 450] [--max-its 1000] [--exempt-taxids 33208,1437180]
"""
import argparse
import collections
import os
import sys

# NCBI taxids: 1437180 = Acrogymnospermae (gymnosperms), 3398 = Magnoliopsida (flowering plants).
# Not 58024: that is Spermatophyta (all seed plants) and would exempt every flowering plant too.
GROUPS = {"4751": "Fungi", "1437180": "Gymnosperms", "3398": "Angiosperms",
          "33090": "Other plants", "33208": "Metazoa"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-fasta", required=True)
    ap.add_argument("--in-map", required=True)
    ap.add_argument("--taxdump-dir", required=True, help="directory with NCBI nodes.dmp and names.dmp")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-map", required=True)
    ap.add_argument("--out-dropped", required=True)
    ap.add_argument("--max-spacer", type=int, default=450)
    ap.add_argument("--max-its", type=int, default=1000)
    ap.add_argument("--exempt-taxids", default="33208,1437180")
    a = ap.parse_args()

    nodes = os.path.join(a.taxdump_dir, "nodes.dmp")
    names = os.path.join(a.taxdump_dir, "names.dmp")
    for f in (a.in_fasta, a.in_map, nodes, names):
        if not os.path.isfile(f):
            sys.exit(f"ERROR: file not found: {f}")
    exempt = {t.strip() for t in a.exempt_taxids.split(",") if t.strip()}

    tax = {}
    for line in open(a.in_map):
        p = line.split()
        if len(p) >= 2:
            tax[p[0]] = p[1]
    parent, rank, name = {}, {}, {}
    for line in open(nodes):
        p = [x.strip() for x in line.split("|")]
        parent[p[0]], rank[p[0]] = p[1], p[2]
    for line in open(names):
        p = [x.strip() for x in line.split("|")]
        if len(p) > 3 and p[3] == "scientific name":
            name[p[0]] = p[1]

    def lineage(t):
        out = []
        while t in parent and t != "1":
            out.append(t)
            t = parent[t]
        return out

    def group(lin):
        return next((GROUPS[x] for x in lin if x in GROUPS), "Other")

    def genus(t, lin):
        return next((name.get(x, x) for x in lin if rank.get(x) == "genus"), name.get(t, t))

    def reason(acc, n, lin):
        if exempt.intersection(lin):
            return None
        region = acc.rsplit("_", 1)[-1]
        if region in ("its1", "its2"):
            return f"single spacer >{a.max_spacer} bp" if n > a.max_spacer else None
        return f"ITS1-5.8S-ITS2 >{a.max_its} bp" if n > a.max_its else None

    recs, cur = [], None
    for line in open(a.in_fasta):
        if line.startswith(">"):
            cur = [line[1:].split()[0], line, []]
            recs.append(cur)
        elif cur is not None:
            cur[2].append(line)

    kept, dropped = 0, collections.Counter()
    with open(a.out_fasta, "w") as fa, open(a.out_map, "w") as fm, open(a.out_dropped, "w") as fd:
        fd.write("accession\tlength_bp\ttaxid\tgroup\tgenus\treason\n")
        for acc, hdr, seq in recs:
            n = sum(len(s.strip()) for s in seq)
            t = tax.get(acc, "0")
            lin = lineage(t)
            r = reason(acc, n, lin)
            if r:
                fd.write(f"{acc}\t{n}\t{t}\t{group(lin)}\t{genus(t, lin)}\t{r}\n")
                dropped[group(lin)] += 1
                continue
            fa.write(hdr)
            fa.writelines(seq)
            fm.write(f"{acc}\t{t}\n")
            kept += 1
    print(f"length filter: kept {kept}, dropped {sum(dropped.values())} {dict(dropped)} "
          f"(max_spacer={a.max_spacer}, max_its={a.max_its}, exempt={','.join(sorted(exempt))})")


if __name__ == "__main__":
    main()
