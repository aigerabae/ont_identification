#!/usr/bin/env python3
"""
select_rdna_queries.py -- choose what goes to the SECONDARY 18S/28S (nuclear rDNA) database.

The custom ITS database stays the primary classifier. The rDNA database is asked only about evidence
that ITS cannot explain:

  --tier contigs   samples in mode "contigs" (assembly OK, but ITSx found no ITS): the same fallback
                   contigs that select_queries.py chose (queries.fasta). Other modes -> empty output.
  --tier residual  residual reads (not on a contig, not human) WITHOUT an informative custom-ITS call,
                   most abundant first, capped at --max-unique dereplicated sequences. Reads beyond the
                   cap stay "Unclassified" (their number is written to the summary JSON).

Always writes --out-fasta and --out-json (possibly empty), so the workflow never breaks.
"""
import argparse, csv, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from taxdump import rank_level  # noqa: E402

RANK_BY_NAME = {"family": 5, "order": 4, "genus": 7, "species": 8, "class": 3, "tribe": 6}


def read_fasta(path):
    recs, h, b = [], None, []
    if path and os.path.exists(path) and os.path.getsize(path) > 0:
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tier", choices=["contigs", "residual"], required=True)
    ap.add_argument("--enabled", default="True", help="False -> write empty outputs (rDNA pass switched off)")
    ap.add_argument("--mode-file", default="")
    ap.add_argument("--queries", default="", help="contigs tier: queries.fasta from select_queries.py")
    ap.add_argument("--fasta", default="", help="residual tier: residual_reads.fasta")
    ap.add_argument("--counts", default="", help="residual tier: residual_counts.tsv")
    ap.add_argument("--lca", default="", help="residual tier: custom_residual.lca.tsv (ITS calls)")
    ap.add_argument("--min-informative-rank", default="family")
    ap.add_argument("--max-unique", type=int, default=20000)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args()

    out, summary = [], {"tier": a.tier, "enabled": a.enabled == "True"}
    if a.enabled != "True":
        summary["note"] = "rDNA pass disabled"
    elif a.tier == "contigs":
        mode = open(a.mode_file).read().strip() if a.mode_file and os.path.exists(a.mode_file) else "none"
        summary["mode"] = mode
        if mode == "contigs":
            out = read_fasta(a.queries)
        summary["queried_contigs"] = len(out)
    else:
        min_level = RANK_BY_NAME.get(a.min_informative_rank, 5)
        called = set()
        if a.lca and os.path.exists(a.lca) and os.path.getsize(a.lca) > 0:
            for r in csv.DictReader(open(a.lca), delimiter="\t"):
                if rank_level(r.get("assigned_rank", "")) >= min_level:
                    called.add(r["query"])
        counts = {}
        if a.counts and os.path.exists(a.counts):
            for line in open(a.counts):
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2 and p[1].isdigit():
                    counts[p[0]] = int(p[1])
        cand = [(q, s) for q, s in read_fasta(a.fasta) if q not in called]
        cand.sort(key=lambda r: (-counts.get(r[0], 0), r[0]))     # most abundant first, deterministic
        out, rest = cand[: a.max_unique], cand[a.max_unique:]
        summary.update({
            "residual_unique_without_its_call": len(cand),
            "residual_reads_without_its_call": sum(counts.get(q, 0) for q, _ in cand),
            "queried_unique": len(out),
            "queried_reads": sum(counts.get(q, 0) for q, _ in out),
            "unqueried_unique": len(rest),
            "unqueried_reads": sum(counts.get(q, 0) for q, _ in rest),
        })

    with open(a.out_fasta, "w") as fh:
        for h, s in out:
            fh.write(f">{h}\n{s}\n")
    json.dump(summary, open(a.out_json, "w"), indent=1)
    print(json.dumps(summary), file=sys.stderr)


if __name__ == "__main__":
    main()
