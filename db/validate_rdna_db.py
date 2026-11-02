#!/usr/bin/env python3
"""
validate_rdna_db.py -- leave-one-out (LOO) calibration of the secondary 18S/28S (nuclear rDNA) database.

The rDNA database is used only for contigs and reads that carry no ITS, so the queries here imitate those:
  * fragments, not whole references: by default 1000 bp (contig-like) and 150 bp (read-like);
  * a reference of >= --cistron-len bp is assumed to be a whole 18S-ITS-28S cistron, and its fragments are
    taken from the 5' end (18S) and the 3' end (28S), so the ITS in the middle cannot help the call;
    shorter references give one random window (seeded, reproducible).
Queries are drawn from plants only (Viridiplantae): a random sample plus up to --per-target references from
each --target-taxa group (so thinly covered families such as Elaeagnaceae are always tested).

BLAST runs once (same output format as the pipeline). Hits are then removed at two LOO levels:
  species  hits to the query's own species are removed (an unsampled species of a sampled genus)
  genus    hits to the query's whole genus are removed (a genus absent from GenBank rDNA, e.g. Leuzea 28S)
and lca_assign.py (the pipeline's own LCA) is run for every --grid rank-cap scheme, plus once without caps.

Verdicts:
  correct_genus        call is the true genus (or a species of it)
  correct_tribe        call is the true tribe or subtribe (an ancestor of the truth below the family)
  correct_family       call is the true family or subfamily
  too_broad            no call, or a correct call above the true family
  wrong_within_family  right family, but a WRONG genus / tribe / subfamily (an over-confident call;
                       validate_db.py counts this as correct_family -- here it is an error)
  wrong_family         call lies in a different family

Outputs in --out:
  calibration_grid.tsv      verdict % per (rank-cap scheme x LOO level x fragment size)
  identity_bins.tsv         uncapped calls: verdict % per best-identity bin (where do errors start?)
  target_groups.tsv         uncapped + chosen scheme, verdict % per target group
  calibration_queries.tsv   per query, per level/scheme: called rank/name, best identity, verdict
  calibration_summary.txt   human-readable summary

Usage (from the repo root):
  python3 db/validate_rdna_db.py --fasta db/build_rdna/refs.fasta --taxid-map db/build_rdna/refs.taxid.map \
      --blast-db db/build_rdna/db/rdna_custom --taxdump-dir db/taxdump --out db/build_rdna/validation \
      [--max-queries 1500] [--per-target 60] [--frag-sizes 1000,150] [--threads 8] [--max-target-seqs 500]
"""
import argparse, csv, os, random, subprocess, sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
from taxdump import Taxonomy  # noqa: E402

VIRIDIPLANTAE = 33090
DEFAULT_GRID = ";".join([
    "species:101,genus:101,family:97",
    "species:101,genus:101,family:95",
    "species:101,genus:101,tribe:97,family:95",
    "species:101,genus:101,tribe:95,family:95",
    "species:101,genus:101,family:93",
    "species:101,genus:101,tribe:95,family:93",
])
DEFAULT_TARGETS = ("Elaeagnaceae,Cardueae,Anthemideae,Cichorieae,Crassulaceae,Lamiaceae,Rosaceae,"
                   "Apiaceae,Fabaceae,Ranunculaceae,Cannabaceae,Zingiberaceae,Poaceae,Brassicaceae")
VERDICTS = ("correct_genus", "correct_tribe", "correct_family", "too_broad", "wrong_within_family", "wrong_family")
BINS = [(99.5, 101, ">=99.5"), (99, 99.5, "99-99.5"), (98, 99, "98-99"), (97, 98, "97-98"),
        (95, 97, "95-97"), (93, 95, "93-95"), (90, 93, "90-93"), (0, 90, "<90")]


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


def make_fragments(acc, seq, sizes, cistron_len, rng):
    """Return [(query_id, fragment)] -- 18S/28S ends of whole cistrons, else one random window per size."""
    out, n = [], len(seq)
    for L in sizes:
        if n <= L:
            out.append((f"{acc}__f{L}_whole", seq))
        elif n >= cistron_len:
            out.append((f"{acc}__f{L}_5end", seq[:L]))
            out.append((f"{acc}__f{L}_3end", seq[-L:]))
        else:
            s = rng.randrange(0, n - L + 1)
            out.append((f"{acc}__f{L}_at{s}", seq[s:s + L]))
    return out


def frag_size(qid):
    return int(qid.split("__f", 1)[1].split("_", 1)[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True); ap.add_argument("--taxid-map", required=True)
    ap.add_argument("--blast-db", required=True); ap.add_argument("--taxdump-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-queries", type=int, default=1500, help="random plant references (before targets)")
    ap.add_argument("--per-target", type=int, default=60, help="extra references per target group")
    ap.add_argument("--target-taxa", default=DEFAULT_TARGETS)
    ap.add_argument("--frag-sizes", default="1000,150")
    ap.add_argument("--cistron-len", type=int, default=4000,
                    help="references at least this long are treated as whole 18S-ITS-28S cistrons")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--task", default="blastn", help="blastn task, as in the pipeline (default blastn)")
    ap.add_argument("--max-target-seqs", type=int, default=500)
    ap.add_argument("--tie-fraction", type=float, default=0.02)
    ap.add_argument("--min-bitscore", default="100"); ap.add_argument("--max-evalue", default="1e-15")
    ap.add_argument("--grid", default=DEFAULT_GRID, help="rank-cap schemes separated by ';'")
    ap.add_argument("--report-scheme", default="species:101,genus:101,tribe:95,family:95",
                    help="scheme shown per target group")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--reuse-blast", action="store_true", help="skip BLAST if loo.blast.tsv exists")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)
    sizes = [int(x) for x in a.frag_sizes.split(",") if x.strip()]
    schemes = [s.strip() for s in a.grid.split(";") if s.strip()]
    if a.report_scheme not in schemes: schemes.append(a.report_scheme)

    tax = Taxonomy(a.taxdump_dir)
    seqs = read_fasta(a.fasta)
    tid = {}
    for l in open(a.taxid_map):
        p = l.split()
        if len(p) >= 2 and p[1].isdigit(): tid[p[0]] = int(p[1])

    def anc(t, rank):
        for x in reversed(tax.lineage(t)):
            if tax.rank(x) == rank: return x
        return None

    # ---- pick plant references: random sample + up to --per-target from each target group
    plants = sorted(i for i in seqs if i in tid and tax.known(tid[i]) and VIRIDIPLANTAE in set(tax.lineage(tid[i])))
    name2tid = {}
    targets = [t.strip() for t in a.target_taxa.split(",") if t.strip()]
    group_of = {}
    for i in plants:
        lin_names = {tax.name(x) for x in tax.lineage(tid[i])}
        for tg in targets:
            if tg in lin_names: group_of[i] = tg; break
    chosen = set(rng.sample(plants, min(a.max_queries, len(plants))))
    for tg in targets:
        members = sorted(i for i in plants if group_of.get(i) == tg)
        rng.shuffle(members)
        chosen.update(members[:a.per_target])
    chosen = sorted(chosen)
    # a reference without a genus cannot be scored at genus level; keep it (family scoring still works)

    q = os.path.join(a.out, "loo_queries.fasta")
    qinfo = {}
    with open(q, "w") as fh:
        for acc in chosen:
            for qid, frag in make_fragments(acc, seqs[acc], sizes, a.cistron_len, rng):
                fh.write(f">{qid}\n{frag}\n")
                qinfo[qid] = acc
    print(f"queries: {len(chosen)} plant references -> {len(qinfo)} fragments "
          f"({', '.join(f'{s} bp' for s in sizes)}); target groups: "
          + ", ".join(f"{tg} {sum(1 for i in chosen if group_of.get(i) == tg)}" for tg in targets), flush=True)

    raw = os.path.join(a.out, "loo.blast.tsv")
    if not (a.reuse_blast and os.path.exists(raw)):
        print("running BLAST (one run; this is the slow step)...", flush=True)
        subprocess.run(["blastn", "-task", a.task, "-db", a.blast_db, "-query", q, "-dust", "no",
                        "-evalue", "1e-5", "-max_target_seqs", str(a.max_target_seqs),
                        "-num_threads", str(a.threads),
                        "-outfmt", "6 qseqid sacc staxids pident length qcovs evalue bitscore stitle",
                        "-out", raw], check=True)

    # ---- LOO filtering at two levels
    qsp = {qid: tax.to_species(tid[acc]) for qid, acc in qinfo.items()}
    qgen = {qid: anc(tid[acc], "genus") for qid, acc in qinfo.items()}
    filt = {lvl: os.path.join(a.out, f"loo.{lvl}.blast.tsv") for lvl in ("species", "genus")}
    sp_cache, gen_cache = {}, {}
    with open(raw) as fi, open(filt["species"], "w") as fs, open(filt["genus"], "w") as fg:
        for line in fi:
            f = line.rstrip("\n").split("\t")
            if f[0] not in qinfo: continue
            hit_t = [int(x) for x in f[2].replace(";", " ").split() if x.isdigit()]
            hit_sp, hit_gen = set(), set()
            for t in hit_t:
                if not tax.known(t): continue
                if t not in sp_cache: sp_cache[t] = tax.to_species(t); gen_cache[t] = anc(t, "genus")
                hit_sp.add(sp_cache[t]); hit_gen.add(gen_cache[t])
            if qsp[f[0]] not in hit_sp: fs.write(line)
            if qgen[f[0]] is None or qgen[f[0]] not in hit_gen: fg.write(line)

    # ---- LCA for each scheme (+ uncapped) and level
    def run_lca(level, scheme):
        tag = (scheme or "uncapped").replace(":", "").replace(",", "_")
        out = os.path.join(a.out, f"loo.{level}.{tag}.lca.tsv")
        subprocess.run([sys.executable, os.path.join(HERE, "..", "scripts", "lca_assign.py"),
                        "--blast", filt[level], "--sample", "LOO", "--taxdump-dir", a.taxdump_dir, "--out", out,
                        "--tie-fraction", str(a.tie_fraction), "--min-bitscore", a.min_bitscore,
                        "--max-evalue", a.max_evalue, "--rank-caps", scheme], check=True)
        return {r["query"]: r for r in csv.DictReader(open(out), delimiter="\t")}

    def verdict(qid, c):
        truth = tid[qinfo[qid]]; tg, tf = anc(truth, "genus"), anc(truth, "family")
        if not c or not c.get("assigned_taxid"):
            return "too_broad"
        ct = int(c["assigned_taxid"]); lin = set(tax.lineage(ct)) | {ct}
        truth_lin = set(tax.lineage(truth)) | {truth}
        if tg and tg in lin and tax.rank(ct) in ("genus", "species", "subspecies", "varietas"):
            return "correct_genus"
        if ct in truth_lin:                       # an ancestor of the truth: family/subfamily/tribe or broader
            if not (tf and tf in lin):
                return "too_broad"
            return "correct_tribe" if tax.rank(ct) in ("tribe", "subtribe") else "correct_family"
        if tf and tf in lin:                      # right family, but a wrong genus / tribe / subfamily
            return "wrong_within_family"
        return "wrong_family"

    grid_rows, per_query = [], defaultdict(dict)
    uncapped = {}
    for level in ("species", "genus"):
        print(f"LCA, LOO level = {level}: uncapped + {len(schemes)} schemes", flush=True)
        uncapped[level] = run_lca(level, "")
        for qid in qinfo:
            c = uncapped[level].get(qid)
            per_query[qid][(level, "uncapped")] = (c, verdict(qid, c))
        for scheme in schemes:
            calls = run_lca(level, scheme)
            cnt = defaultdict(Counter)
            for qid in qinfo:
                v = verdict(qid, calls.get(qid))
                cnt[frag_size(qid)][v] += 1
                per_query[qid][(level, scheme)] = (calls.get(qid), v)
            for L in sizes:
                n = max(sum(cnt[L].values()), 1)
                grid_rows.append([scheme, level, L, n] + [f"{100 * cnt[L][k] / n:.1f}" for k in VERDICTS])

    with open(os.path.join(a.out, "calibration_grid.tsv"), "w") as fh:
        fh.write("rank_caps\tloo_level\tfragment_bp\tn\t" + "\t".join(f"{k}_pct" for k in VERDICTS) + "\n")
        for r in grid_rows: fh.write("\t".join(map(str, r)) + "\n")

    # ---- identity bins (uncapped): where do wrong calls start?
    bin_rows = []
    for level in ("species", "genus"):
        for L in sizes:
            cnt = defaultdict(Counter)
            for qid in qinfo:
                if frag_size(qid) != L: continue
                c, v = per_query[qid][(level, "uncapped")]
                if not c or not c.get("best_pident"): cnt["no hit"][v] += 1; continue
                pid = float(c["best_pident"])
                lab = next(lab for lo, hi, lab in BINS if lo <= pid < hi)
                cnt[lab][v] += 1
            for lab in [b[2] for b in BINS] + ["no hit"]:
                n = sum(cnt[lab].values())
                if n: bin_rows.append([level, L, lab, n] + [f"{100 * cnt[lab][k] / n:.1f}" for k in VERDICTS])
    with open(os.path.join(a.out, "identity_bins.tsv"), "w") as fh:
        fh.write("loo_level\tfragment_bp\tbest_identity\tn\t" + "\t".join(f"{k}_pct" for k in VERDICTS) + "\n")
        for r in bin_rows: fh.write("\t".join(map(str, r)) + "\n")

    # ---- per target group, report scheme
    tg_rows = []
    for level in ("species", "genus"):
        for tg in targets + ["(all plants)"]:
            for L in sizes:
                cnt = Counter()
                for qid, acc in qinfo.items():
                    if frag_size(qid) != L: continue
                    if tg != "(all plants)" and group_of.get(acc) != tg: continue
                    cnt[per_query[qid][(level, a.report_scheme)][1]] += 1
                n = sum(cnt.values())
                if n: tg_rows.append([level, tg, L, n] + [f"{100 * cnt[k] / n:.1f}" for k in VERDICTS])
    with open(os.path.join(a.out, "target_groups.tsv"), "w") as fh:
        fh.write(f"loo_level\tgroup\tfragment_bp\tn\t" + "\t".join(f"{k}_pct" for k in VERDICTS)
                 + f"\t(scheme: {a.report_scheme})\n")
        for r in tg_rows: fh.write("\t".join(map(str, r)) + "\n")

    with open(os.path.join(a.out, "calibration_queries.tsv"), "w") as fh:
        keys = [(lvl, s) for lvl in ("species", "genus") for s in ["uncapped"] + schemes]
        fh.write("query\treference\ttrue_species\tgroup\t" +
                 "\t".join(f"{lvl}|{s}|{x}" for lvl, s in keys for x in ("rank", "name", "pident", "verdict")) + "\n")
        for qid, acc in sorted(qinfo.items()):
            row = [qid, acc, tax.name(tid[acc]), group_of.get(acc, "")]
            for k in keys:
                c, v = per_query[qid][k]
                row += [c.get("assigned_rank", "") if c else "", c.get("assigned_name", "") if c else "",
                        c.get("best_pident", "") if c else "", v]
            fh.write("\t".join(map(str, row)) + "\n")

    # ---- summary
    with open(os.path.join(a.out, "calibration_summary.txt"), "w") as fh:
        fh.write(f"rDNA LOO calibration: {len(chosen)} plant references, {len(qinfo)} fragments, "
                 f"sizes {sizes}, blastn -task {a.task}, max_target_seqs {a.max_target_seqs}, "
                 f"tie {a.tie_fraction}\n\n")
        fh.write("Verdict % per rank-cap scheme (species = own species hidden; genus = own genus hidden)\n")
        fh.write(f"{'rank_caps':<36}{'level':<9}{'bp':>6}{'n':>7}" + "".join(f"{k:>21}" for k in VERDICTS) + "\n")
        for r in grid_rows:
            fh.write(f"{r[0]:<36}{r[1]:<9}{r[2]:>6}{r[3]:>7}" + "".join(f"{x:>21}" for x in r[4:]) + "\n")
        fh.write("\nUncapped calls by best identity (verdict %): see identity_bins.tsv\n")
        fh.write(f"Per target group with {a.report_scheme}: see target_groups.tsv\n")
    print(open(os.path.join(a.out, "calibration_summary.txt")).read())


if __name__ == "__main__":
    main()
