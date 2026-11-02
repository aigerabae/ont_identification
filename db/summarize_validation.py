#!/usr/bin/env python3
"""
summarize_validation.py -- break the leave-one-out result down by source group, and list the
'wrong' and 'too_broad' cases so they can be inspected.

Usage: db/summarize_validation.py --loo db/build/validation/validation_loo.tsv \
           --manifest db/build/refs_manifest.tsv [--out db/build/validation/by_group.tsv]
"""
import argparse, csv
from collections import Counter, defaultdict

V = ["correct_genus", "correct_family", "too_broad", "wrong"]
ap = argparse.ArgumentParser()
ap.add_argument("--loo", required=True); ap.add_argument("--manifest", required=True)
ap.add_argument("--out", default="by_group.tsv")
a = ap.parse_args()

grp = {r["ref_id"]: r["source_group"] for r in csv.DictReader(open(a.manifest), delimiter="\t")}
rows = list(csv.DictReader(open(a.loo), delimiter="\t"))
by = defaultdict(Counter)
for r in rows:
    by[grp.get(r["reference"], "?")][r["verdict"]] += 1

with open(a.out, "w") as fo:
    hdr = ["group", "n"] + V + ["pct_genus_or_family"]
    fo.write("\t".join(hdr) + "\n")
    print(f"{'group':38s}{'n':>6s}" + "".join(f"{v:>16s}" for v in V) + f"{'genus+family %':>16s}")
    for g, c in sorted(by.items(), key=lambda kv: -sum(kv[1].values())):
        n = sum(c.values()); ok = 100 * (c["correct_genus"] + c["correct_family"]) / n
        fo.write("\t".join([g, str(n)] + [str(c[v]) for v in V] + [f"{ok:.1f}"]) + "\n")
        print(f"{g:38s}{n:6d}" + "".join(f"{c[v]:16d}" for v in V) + f"{ok:16.1f}")

print("\n--- most frequent 'wrong' calls (true species -> called name) ---")
w = Counter((r["true_species"], r["called_name"]) for r in rows if r["verdict"] == "wrong")
for (t, c), n in w.most_common(15): print(f"{n:4d}  {t}  ->  {c}")
print("\n--- 'too_broad': what rank/name was called ---")
tb = Counter((r["called_rank"] or "none", r["called_name"] or "-") for r in rows if r["verdict"] == "too_broad")
for (rk, nm), n in tb.most_common(12): print(f"{n:4d}  {rk:12s} {nm}")
