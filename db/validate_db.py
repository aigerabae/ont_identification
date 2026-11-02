#!/usr/bin/env python3
"""
validate_db.py -- leave-one-out (LOO) validation of the custom ITS database.

Every reference sequence is used as a query against the database; hits to the
SAME SPECIES are removed afterwards (so the species is "unknown" to the
classifier, as a real novel/under-represented species would be). The LCA call
that remains is compared with the truth:

  correct_genus     the call is the true genus (or a species of it)
  correct_family    the call is the true family (genus not resolvable)
  too_broad         call above the true family / no call
  wrong             the call lies in a different family

and, as a separate line, the plain (non-LOO) self-hit run: does each sequence
return its own species at 100%?

Usage:
  db/validate_db.py --fasta db/build/refs.fasta --taxid-map db/build/refs.taxid.map \
      --blast-db db/build/db/its_custom --taxdump-dir db/taxdump --out db/build/validation \
      [--max-queries 3000] [--threads 8]
Output: validation_loo.tsv (per query), validation_summary.txt
"""
import argparse, os, random, subprocess, sys, tempfile, csv
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
from taxdump import Taxonomy  # noqa: E402


def read_fasta(p):
    out, h, b = {}, None, []
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"):
            if h: out[h] = "".join(b)
            h, b = l[1:].split()[0], []
        elif h:
            b.append(l.strip())
    if h: out[h] = "".join(b)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True); ap.add_argument("--taxid-map", required=True)
    ap.add_argument("--blast-db", required=True); ap.add_argument("--taxdump-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-queries", type=int, default=3000)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--tie-fraction", type=float, default=0.02)
    ap.add_argument("--rank-caps", default="species:97,genus:90,family:80")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    tax = Taxonomy(a.taxdump_dir)
    seqs = read_fasta(a.fasta)
    tid = {}
    for l in open(a.taxid_map):
        i, t = l.split()[:2]; tid[i] = int(t)
    ids = [i for i in seqs if i in tid and tax.known(tid[i])]
    random.Random(a.seed).shuffle(ids)
    ids = ids[:a.max_queries]

    q = os.path.join(a.out, "loo_queries.fasta")
    with open(q, "w") as fh:
        for i in ids: fh.write(f">{i}\n{seqs[i]}\n")
    raw = os.path.join(a.out, "loo.blast.tsv")
    subprocess.run(["blastn", "-task", "blastn", "-db", a.blast_db, "-query", q, "-dust", "no",
                    "-evalue", "1e-5", "-max_target_seqs", "500", "-num_threads", str(a.threads),
                    "-outfmt", "6 qseqid sacc staxids pident length qcovs evalue bitscore stitle",
                    "-out", raw], check=True)

    # remove same-species hits
    qsp = {i: tax.to_species(tid[i]) for i in ids}
    filt, selfhit = os.path.join(a.out, "loo.filtered.blast.tsv"), Counter()
    with open(raw) as fi, open(filt, "w") as fo:
        for line in fi:
            f = line.rstrip("\n").split("\t")
            hit_t = [int(x) for x in f[2].replace(";", " ").split() if x.isdigit()]
            hit_sp = {tax.to_species(t) for t in hit_t if tax.known(t)}
            if f[0] in qsp and qsp[f[0]] in hit_sp:
                if float(f[3]) >= 99.9: selfhit[f[0]] += 1
                continue
            fo.write(line)
    lca = os.path.join(a.out, "loo.lca.tsv")
    subprocess.run([sys.executable, os.path.join(HERE, "..", "scripts", "lca_assign.py"), "--blast", filt,
                    "--sample", "LOO", "--taxdump-dir", a.taxdump_dir, "--out", lca,
                    "--tie-fraction", str(a.tie_fraction), "--min-bitscore", "100", "--max-evalue", "1e-15",
                    "--rank-caps", a.rank_caps], check=True)
    calls = {r["query"]: r for r in csv.DictReader(open(lca), delimiter="\t")}

    def anc(t, rank):
        for x in reversed(tax.lineage(t)):
            if tax.rank(x) == rank: return x
        return None

    res, cnt = [], Counter()
    for i in ids:
        truth = tid[i]; tg, tf = anc(truth, "genus"), anc(truth, "family")
        c = calls.get(i)
        if not c or not c.get("assigned_taxid"):
            verdict = "too_broad"
        else:
            ct = int(c["assigned_taxid"]); lin = set(tax.lineage(ct))
            if tg and tg in lin and tax.rank(ct) in ("genus", "species", "subspecies"):
                verdict = "correct_genus"
            elif tf and tf in lin and tax.rank(ct) in ("family", "tribe", "subfamily", "genus", "species"):
                verdict = "correct_family"
            elif ct in set(tax.lineage(truth)):
                verdict = "too_broad"
            else:
                verdict = "wrong"
        cnt[verdict] += 1
        res.append((i, tax.name(truth), c["assigned_rank"] if c else "", c["assigned_name"] if c else "", verdict))
    with open(os.path.join(a.out, "validation_loo.tsv"), "w") as fh:
        fh.write("reference\ttrue_species\tcalled_rank\tcalled_name\tverdict\n")
        for r in res: fh.write("\t".join(map(str, r)) + "\n")
    n = max(len(ids), 1)
    with open(os.path.join(a.out, "validation_summary.txt"), "w") as fh:
        fh.write(f"LOO queries: {len(ids)}\n")
        for k in ("correct_genus", "correct_family", "too_broad", "wrong"):
            fh.write(f"{k}\t{cnt[k]}\t{100*cnt[k]/n:.1f}%\n")
        fh.write(f"queries whose own species was found at >=99.9% before filtering: {len(selfhit)} / {len(ids)}\n")
    print(open(os.path.join(a.out, "validation_summary.txt")).read())


if __name__ == "__main__":
    main()
