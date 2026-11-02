#!/usr/bin/env python3
"""
diagnose_residual.py -- why are residual reads unclassified?

For one sample folder (results/<sample>) it takes the residual reads (reads that did not map to the ITS
contigs and are not human), and says, weighted by read count, what happens to them:
  * called by the custom ITS database (confident LCA call)
  * no BLAST hit at all                       -> not ITS: 18S/28S flank, other loci, host, junk
  * hit, but covers <80 % of the read         -> partial overlap (read spans an ITS boundary / contig end)
  * hit covers the read but identity <90 %    -> divergent relative missing from the database
  * hit passes, but LCA is above family       -> conserved (5.8S) region shared by many taxa
and how many of them align (bowtie2 local mode) to the assembled ITS contigs / to ALL assembled contigs,
i.e. could be rescued by a more tolerant mapping. The 5 most abundant unclassified sequences are printed
so they can be checked by hand (e.g. web BLAST).

Usage: diagnose_residual.py --sample-dir results/ONplants-1_S40 [--threads 4] [--top 5]
"""
import argparse, csv, os, subprocess, sys, tempfile
from collections import defaultdict, Counter


def fasta(path):
    out, h, b = {}, None, []
    if not os.path.exists(path): return out
    for l in open(path):
        l = l.rstrip("\n")
        if l.startswith(">"):
            if h: out[h] = "".join(b)
            h, b = l[1:].split()[0], []
        elif h: b.append(l.strip())
    if h: out[h] = "".join(b)
    return out


def local_map(ref, reads, threads):
    """set of read ids that align (local mode) to ref"""
    if not os.path.exists(ref) or os.path.getsize(ref) == 0 or os.path.getsize(reads) == 0:
        return set()
    d = tempfile.mkdtemp()
    subprocess.run(["bowtie2-build", "--quiet", ref, f"{d}/i"], check=True)
    p = subprocess.run(["bowtie2", "--very-sensitive-local", "-f", "-U", reads, "-x", f"{d}/i",
                        "-p", str(threads), "--no-unal", "--no-hd", "--no-sq", "-S", "-"],
                       capture_output=True, text=True, check=True)
    return {l.split("\t", 1)[0] for l in p.stdout.splitlines() if l}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-dir", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--min-qcov", type=float, default=80)
    ap.add_argument("--min-pident", type=float, default=90)
    a = ap.parse_args()
    S = a.sample_dir
    reads = fasta(f"{S}/residual/residual_reads.fasta")
    cnt = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(f"{S}/residual/residual_counts.tsv") if l.strip()}
    called = {r["query"] for r in csv.DictReader(open(f"{S}/classification/custom_residual.lca.tsv"), delimiter="\t")
              if r.get("assigned_taxid")} if os.path.exists(f"{S}/classification/custom_residual.lca.tsv") else set()
    best = {}
    bf = f"{S}/classification/custom_residual.blast.tsv"
    if os.path.exists(bf):
        for l in open(bf):
            f = l.rstrip("\n").split("\t")
            if len(f) < 8: continue
            q, pid, qc, bs = f[0], float(f[3]), float(f[5]), float(f[7])
            if q not in best or bs > best[q][2]: best[q] = (pid, qc, bs, f[8] if len(f) > 8 else "")

    cat = defaultdict(lambda: [0, 0]); cat_of = {}
    for q in reads:
        n = cnt.get(q, 1)
        if q in called: c = "1 called by the custom database"
        elif q not in best: c = "2 no BLAST hit (not ITS: flanks / other loci / junk)"
        else:
            pid, qc, bs, _ = best[q]
            if qc < a.min_qcov: c = f"3 hit covers <{a.min_qcov:.0f}% of the read (partial overlap)"
            elif pid < a.min_pident: c = f"4 hit identity <{a.min_pident:.0f}% (divergent relative)"
            else: c = "5 good hit, but LCA too broad / low score (conserved region)"
        cat_of[q] = c; cat[c][0] += 1; cat[c][1] += n
    total = sum(v[1] for v in cat.values()) or 1
    print(f"Residual reads for {S}: {sum(v[0] for v in cat.values())} unique sequences, {total} reads")
    print(f"{'category':68s}{'unique':>9s}{'reads':>10s}{'% resid':>9s}")
    for c in sorted(cat): print(f"{c:68s}{cat[c][0]:9d}{cat[c][1]:10d}{100*cat[c][1]/total:8.1f}%")

    unc = [q for q in reads if q not in called]
    tmp = tempfile.mkdtemp(); un = f"{tmp}/unclassified.fa"
    with open(un, "w") as fh:
        for q in unc: fh.write(f">{q}\n{reads[q]}\n")
    for label, ref in (("ITS contigs", f"{S}/queries/mapping_ref.fasta"),
                       ("ALL assembled contigs", f"{S}/chimera/nonchimeric.fasta")):
        m = local_map(ref, un, a.threads)
        r = sum(cnt.get(q, 1) for q in m); u = sum(cnt.get(q, 1) for q in unc) or 1
        print(f"unclassified residual reads that align locally to {label}: {len(m)} unique / {r} reads ({100*r/u:.1f}% of unclassified)")
    print(f"\nlengths of unclassified sequences (median): {sorted(len(reads[q]) for q in unc)[len(unc)//2] if unc else 0}")
    print(f"\nTop {a.top} most abundant unclassified sequences (paste into web BLAST):")
    for q in sorted(unc, key=lambda x: -cnt.get(x, 1))[:a.top]:
        b = best.get(q); info = f"best hit {b[0]:.1f}% id, {b[1]:.0f}% cov: {b[3][:70]}" if b else "no hit"
        print(f"  {cnt.get(q,1):7d} reads | {cat_of[q][:1]} | {info}\n          {reads[q]}")


if __name__ == "__main__":
    main()
