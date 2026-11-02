#!/usr/bin/env python3
"""Count reads after fastp. total = merged + R1 + R2 (every read counted once,
exactly what bowtie2/flagstat later sees when -1/-2/-U are all supplied)."""
import argparse, gzip, json

ap = argparse.ArgumentParser()
ap.add_argument("--merged", required=True)
ap.add_argument("--r1", required=True)
ap.add_argument("--r2", required=True)
ap.add_argument("--fastp-json", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()


def n(path):
    with gzip.open(path, "rt") as f:
        return sum(1 for _ in f) // 4


raw = 0
try:
    j = json.load(open(a.fastp_json))
    raw = int(j["summary"]["before_filtering"]["total_reads"])
except Exception:
    pass
m, r1, r2 = n(a.merged), n(a.r1), n(a.r2)
json.dump({"raw_reads": raw, "merged_reads": m, "unmerged_R1": r1, "unmerged_R2": r2,
           "post_fastp_reads": m + r1 + r2}, open(a.out, "w"), indent=1)
