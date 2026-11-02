#!/usr/bin/env python3
"""
residual_reads.py -- reads NOT explained by the assembly-based evidence.

Residual reads = reads that did not map to the ITS contigs / fallback contigs
(or ALL reads when there is no usable assembly). From them:
  * reads that map to the human reference (human_screen BAM) are counted and
    removed (they are reported as Homo sapiens, not sent to BLAST)
  * the rest are dereplicated (exact sequence) and written as BLAST queries
    with their read counts, most abundant first, capped at --max-unique

Read identity key: QNAME (+ '/1' or '/2' for mates), matching the BAMs.
"""
import argparse, gzip, json, os, subprocess, sys
from collections import Counter


def bam_keys(bam, include_flags="", exclude_flags=""):
    """Set of read keys from a BAM (empty set if the BAM is empty/missing)."""
    keys = set()
    if not bam or not os.path.exists(bam) or os.path.getsize(bam) == 0:
        return keys
    cmd = ["samtools", "view"]
    if include_flags:
        cmd += ["-f", include_flags]
    if exclude_flags:
        cmd += ["-F", exclude_flags]
    cmd.append(bam)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.split("\t", 2)
        flag = int(f[1])
        suf = "/1" if flag & 64 else ("/2" if flag & 128 else "")
        keys.add(f[0] + suf)
    p.wait()
    if p.returncode != 0:
        sys.exit(f"samtools view failed on {bam}")
    return keys


def fastq(path, suffix):
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            seq = fh.readline().strip()
            fh.readline(); fh.readline()
            name = h[1:].split()[0]
            if name.endswith("/1") or name.endswith("/2"):
                name = name[:-2]
            yield name + suffix, seq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--quant-bam", default="")
    ap.add_argument("--human-bam", default="")
    ap.add_argument("--max-unique", type=int, default=50000)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-counts", required=True)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    have_quant = bool(a.quant_bam) and os.path.exists(a.quant_bam) and os.path.getsize(a.quant_bam) > 0
    unmapped = bam_keys(a.quant_bam, include_flags="4", exclude_flags="2304") if have_quant else None
    human = bam_keys(a.human_bam, exclude_flags="2308")

    counter, residual, n_human = Counter(), 0, 0
    for path, suf in ((a.merged, ""), (a.r1, "/1"), (a.r2, "/2")):
        for key, seq in fastq(path, suf):
            if unmapped is not None and key not in unmapped:
                continue
            residual += 1
            if key in human:
                n_human += 1
                continue
            counter[seq] += 1

    ranked = sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))
    keep, rest = ranked[: a.max_unique], ranked[a.max_unique:]
    with open(a.out_fasta, "w") as ff, open(a.out_counts, "w") as fc:
        for i, (seq, n) in enumerate(keep, 1):
            ff.write(f">rd{i}\n{seq}\n")
            fc.write(f"rd{i}\t{n}\n")
    summary = {
        "assembly_mapping_used": have_quant,
        "residual_reads_total": residual,
        "human_screen_reads": n_human,
        "unique_sequences": len(counter),
        "queried_unique": len(keep),
        "queried_reads": sum(n for _, n in keep),
        "unqueried_reads": sum(n for _, n in rest),
    }
    json.dump(summary, open(a.out_json, "w"), indent=1)
    print(json.dumps(summary), file=sys.stderr)


if __name__ == "__main__":
    main()
