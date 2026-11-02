#!/usr/bin/env python3
"""
core_nt_control.py -- the CONTROL search against the general NCBI core_nt database.

core_nt is not the primary classifier in this pipeline (the custom ITS
database is). It is used to look at sequences the custom database could not
explain: contamination (human, bacteria, ...) and unexpected taxa.

Query selection:
  --select unassigned : only queries without an informative custom-DB call
  --select all        : every query (used for the custom-vs-core_nt concordance check)
Highest read-count queries first, at most --max-queries.

Modes:  remote (NCBI web service, polite chunking), local (a local core_nt copy), skip.
Writes the same 9-column BLAST table as the primary search, plus a status JSON
recording how much of the sample's reads the control actually covered.
"""
import argparse, csv, json, os, subprocess, sys, tempfile, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from taxdump import rank_level  # noqa: E402

FMT = "6 qseqid sacc staxids pident length qcovs evalue bitscore stitle"


def read_fasta(path):
    d, h, b = {}, None, []
    if os.path.exists(path) and os.path.getsize(path) > 0:
        for line in open(path):
            line = line.rstrip("\n")
            if line.startswith(">"):
                if h is not None:
                    d[h] = "".join(b)
                h, b = line[1:].split()[0], []
            elif h is not None:
                b.append(line.strip())
        if h is not None:
            d[h] = "".join(b)
    return d


def read_counts(path):
    c = {}
    if path and os.path.exists(path):
        for line in open(path):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1].isdigit():
                c[p[0]] = int(p[1])
    return c


def informative_custom(path, min_level):
    ok = set()
    if path and os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if rank_level(r["assigned_rank"]) >= min_level:
                    ok.add(r["query"])
    return ok


def run_blast(chunk_fa, out_tmp, mode, db, threads, entrez, evalue):
    cmd = ["blastn", "-task", "blastn", "-db", db, "-query", chunk_fa, "-outfmt", FMT,
           "-max_target_seqs", "50", "-evalue", str(evalue), "-out", out_tmp]
    if mode == "remote":
        cmd.append("-remote")
        if entrez:
            cmd += ["-entrez_query", entrez]
    else:
        cmd += ["-num_threads", str(threads)]
    return subprocess.run(cmd, capture_output=True, text=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--queries", required=True)
    ap.add_argument("--counts", default="")
    ap.add_argument("--custom-lca", default="")
    ap.add_argument("--select", choices=["unassigned", "all"], default="unassigned")
    ap.add_argument("--min-informative-rank-level", type=int, default=5, help="5 = family")
    ap.add_argument("--mode", choices=["remote", "local", "skip"], default="skip")
    ap.add_argument("--db", default="core_nt")
    ap.add_argument("--max-queries", type=int, default=200)
    ap.add_argument("--chunk-size", type=int, default=25)
    ap.add_argument("--sleep", type=float, default=15)
    ap.add_argument("--entrez-query", default="")
    ap.add_argument("--evalue", type=float, default=1e-10)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", required=True)
    ap.add_argument("--status", required=True)
    a = ap.parse_args()

    queries = read_fasta(a.queries)
    counts = read_counts(a.counts)
    done = informative_custom(a.custom_lca, a.min_informative_rank_level)
    ids = [q for q in queries if a.select == "all" or q not in done]
    ids.sort(key=lambda q: (-counts.get(q, 0), q))
    ids = ids[: a.max_queries]
    total_reads = sum(counts.get(q, 0) for q in queries)
    sel_reads = sum(counts.get(q, 0) for q in ids)

    status = {"mode": a.mode, "db": a.db, "select": a.select, "n_queries_total": len(queries),
              "n_selected": len(ids), "n_done": 0, "n_failed_chunks": 0,
              "reads_in_query_set": total_reads, "reads_in_selected_queries": sel_reads,
              "errors": []}
    open(a.out, "w").close()
    if a.mode == "skip" or not ids:
        status["note"] = "control skipped" if a.mode == "skip" else "nothing to control"
        json.dump(status, open(a.status, "w"), indent=1)
        return

    size = a.chunk_size if a.mode == "remote" else max(len(ids), 1)
    for i in range(0, len(ids), size):
        chunk = ids[i:i + size]
        with tempfile.TemporaryDirectory() as td:
            fa, tmp = os.path.join(td, "q.fa"), os.path.join(td, "o.tsv")
            with open(fa, "w") as fh:
                for q in chunk:
                    fh.write(f">{q}\n{queries[q]}\n")
            res = run_blast(fa, tmp, a.mode, a.db, a.threads, a.entrez_query, a.evalue)
            if res.returncode != 0 and a.mode == "remote":
                time.sleep(60)                       # one polite retry
                res = run_blast(fa, tmp, a.mode, a.db, a.threads, a.entrez_query, a.evalue)
            if res.returncode == 0:
                if os.path.exists(tmp):
                    with open(a.out, "a") as fo, open(tmp) as fi:
                        fo.write(fi.read())
                status["n_done"] += len(chunk)
            else:
                status["n_failed_chunks"] += 1
                status["errors"].append(res.stderr.strip()[-300:])
        if a.mode == "remote" and i + size < len(ids):
            time.sleep(a.sleep)
    json.dump(status, open(a.status, "w"), indent=1)
    print(json.dumps({k: v for k, v in status.items() if k != "errors"}), file=sys.stderr)


if __name__ == "__main__":
    main()
