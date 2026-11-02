#!/usr/bin/env python3
"""
concordance.py -- custom-DB call vs core_nt call, per ITS contig.

Needs the core_nt control to have been run with --select all (config:
core_nt.validate_all_its_contigs: true). For every query with both calls:
  concordant              same taxon
  compatible_custom_finer one is an ancestor of the other, custom is more specific
  compatible_core_finer   ... core_nt is more specific
  discordant              different lineages
  no_core_nt_call         core_nt gave no confident call
Discordant rows are the ones to inspect by hand.
"""
import argparse, csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from taxdump import Taxonomy  # noqa: E402


def load(path):
    d = {}
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                d[r["query"]] = r
    return d


ap = argparse.ArgumentParser()
ap.add_argument("--taxdump-dir", required=True)
ap.add_argument("--pairs", nargs="+", required=True, help="sample=custom_lca.tsv=control_lca.tsv")
ap.add_argument("--out", required=True)
a = ap.parse_args()
tax = Taxonomy(a.taxdump_dir)

counts = {}
with open(a.out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["sample", "query", "custom_taxon", "custom_rank", "core_nt_taxon", "core_nt_rank", "status"])
    for item in a.pairs:
        sample, cpath, kpath = item.split("=")
        cust, ctrl = load(cpath), load(kpath)
        for q, c in cust.items():
            k = ctrl.get(q)
            if k is None:
                st = "no_core_nt_call"
            else:
                ct = int(c["assigned_taxid"]) if c["assigned_taxid"] else None
                kt = int(k["assigned_taxid"]) if k["assigned_taxid"] else None
                if ct is None or kt is None:
                    st = "discordant"
                elif ct == kt:
                    st = "concordant"
                elif ct in tax.lineage(kt):
                    st = "compatible_core_finer"
                elif kt in tax.lineage(ct):
                    st = "compatible_custom_finer"
                else:
                    st = "discordant"
            counts[st] = counts.get(st, 0) + 1
            w.writerow([sample, q, c["assigned_name"], c["assigned_rank"],
                        k["assigned_name"] if k else "", k["assigned_rank"] if k else "", st])
print("concordance:", counts, file=sys.stderr)
