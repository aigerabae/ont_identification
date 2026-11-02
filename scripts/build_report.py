#!/usr/bin/env python3
"""
build_report.py -- per-sample accounting of EVERY read into taxon rows.

Each post-fastp read ends in exactly one row:

  A. reads mapped to an assembled contig (ITS contig, or fallback contig):
       taxon from the CUSTOM ITS DB if that call is informative (>= min rank),
       else from the core_nt CONTROL if informative,
       else from the secondary 18S/28S (rDNA) DB if informative (fallback contigs only),
       else row "Mapped to contig, no confident taxon".
  B. reads not mapped to any contig ("residual"):
       reads mapping to the human reference screen -> Homo sapiens
       other residual reads, dereplicated and classified the same way
       (custom ITS DB first, core_nt control second, rDNA DB third), else "Unclassified".
  The rDNA DB only ever fills gaps: a read that already has an informative ITS (or control) call keeps it.
  C. anything left (should be 0) -> "Unclassified" and a QC warning.

Evidence column records which line of evidence gave the taxon, so results of
different confidence are never silently merged.
"""
import argparse, csv, json, os, sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from taxdump import Taxonomy, rank_level  # noqa: E402

RANK_BY_NAME = {"family": 5, "order": 4, "genus": 7, "species": 8, "class": 3, "tribe": 6}


def load_lca(path):
    d = {}
    if path and os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                d[r["query"]] = r
    return d


def load_counts(path):
    d = {}
    if path and os.path.exists(path):
        for line in open(path):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1].isdigit():
                d[p[0]] = int(p[1])
    return d


def load_json(path):
    try:
        return json.load(open(path))
    except Exception:
        return {}


def divergence_class(pident, thr):
    hi, mid, lo = thr
    if pident >= hi:
        return "identical (>=%g%%)" % hi
    if pident >= mid:
        return "low divergence (%g-%g%%)" % (mid, hi)
    if pident >= lo:
        return "moderate divergence (%g-%g%%)" % (lo, mid)
    return "high divergence (<%g%%)" % lo


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--taxdump-dir", required=True)
    ap.add_argument("--read-counts", required=True)
    ap.add_argument("--mode-file", required=True)
    ap.add_argument("--contig-counts", required=True)
    ap.add_argument("--lca-custom-contigs", required=True)
    ap.add_argument("--lca-control-contigs", required=True)
    ap.add_argument("--control-status-contigs", default="")
    ap.add_argument("--residual-json", required=True)
    ap.add_argument("--residual-counts", required=True)
    ap.add_argument("--lca-custom-residual", required=True)
    ap.add_argument("--lca-control-residual", required=True)
    ap.add_argument("--control-status-residual", default="")
    ap.add_argument("--lca-rdna-contigs", default="", help="LCA of fallback contigs vs the 18S/28S DB")
    ap.add_argument("--lca-rdna-residual", default="", help="LCA of residual reads vs the 18S/28S DB")
    ap.add_argument("--rdna-residual-json", default="", help="select_rdna_queries.py summary (residual tier)")
    ap.add_argument("--assembly-status", nargs="*", default=[])
    ap.add_argument("--min-informative-rank", default="family")
    ap.add_argument("--divergence-thresholds", default="99.5,97,93")
    ap.add_argument("--out-rows", required=True)
    ap.add_argument("--out-qc", required=True)
    a = ap.parse_args()

    tax = Taxonomy(a.taxdump_dir)
    min_level = RANK_BY_NAME.get(a.min_informative_rank, 5)
    div_thr = tuple(float(x) for x in a.divergence_thresholds.split(","))

    rc = load_json(a.read_counts)
    total = int(rc.get("post_fastp_reads", 0))
    mode = open(a.mode_file).read().strip() or "none"
    contig_counts = load_counts(a.contig_counts)
    res_counts = load_counts(a.residual_counts)
    rj = load_json(a.residual_json)

    custom = {"contigs": load_lca(a.lca_custom_contigs), "residual": load_lca(a.lca_custom_residual)}
    control = {"contigs": load_lca(a.lca_control_contigs), "residual": load_lca(a.lca_control_residual)}
    ev_custom = {"contigs": "custom DB / " + ("ITS region" if mode == "its" else "contig, no ITS"),
                 "residual": "custom DB / read-level"}
    ev_ctrl = {"contigs": "core_nt control / contig", "residual": "core_nt control / read-level"}
    rdna = {"contigs": load_lca(a.lca_rdna_contigs), "residual": load_lca(a.lca_rdna_residual)}
    ev_rdna = {"contigs": "rDNA DB (18S/28S) / contig, no ITS", "residual": "rDNA DB (18S/28S) / read-level"}

    acc = defaultdict(int)       # (taxon, taxid, rank, group, evidence, div) -> reads
    stats = defaultdict(int)

    def add(taxon, taxid, rank, group, evidence, div, n):
        acc[(taxon, taxid, rank, group, evidence, div)] += n

    def informative(r):
        return r is not None and rank_level(r["assigned_rank"]) >= min_level

    def resolve(tier, q, n):
        c, k = custom[tier].get(q), control[tier].get(q)
        if informative(c):
            tid = int(c["assigned_taxid"]) if c["assigned_taxid"] else None
            add(c["assigned_name"], tid or "", c["assigned_rank"],
                tax.group(tid) if tid else "Unknown", ev_custom[tier],
                divergence_class(float(c["best_pident"]), div_thr), n)
            stats["custom_assigned_reads_" + tier] += n
            return True
        if informative(k):
            tid = int(k["assigned_taxid"]) if k["assigned_taxid"] else None
            add(k["assigned_name"], tid or "", k["assigned_rank"],
                tax.group(tid) if tid else "Unknown", ev_ctrl[tier], "-", n)
            stats["control_assigned_reads_" + tier] += n
            return True
        d = rdna[tier].get(q)
        if informative(d):
            tid = int(d["assigned_taxid"]) if d["assigned_taxid"] else None
            add(d["assigned_name"], tid or "", d["assigned_rank"],
                tax.group(tid) if tid else "Unknown", ev_rdna[tier],
                divergence_class(float(d["best_pident"]), div_thr), n)
            stats["rdna_assigned_reads_" + tier] += n
            return True
        return False

    # A. contig-mapped reads
    mapped_total = 0
    for q, n in contig_counts.items():
        mapped_total += n
        if n == 0:
            continue
        if not resolve("contigs", q, n):
            add("Mapped to contig, no confident taxon", "", "-", "-", "assembly", "-", n)
            stats["contig_no_taxon_reads"] += n

    # B. residual reads
    residual_total = int(rj.get("residual_reads_total", 0))
    human = int(rj.get("human_screen_reads", 0))
    if human:
        add("Homo sapiens", 9606, "species", "Metazoa", "human reference screen", "-", human)
    res_unclassified = int(rj.get("unqueried_reads", 0))
    for q, n in res_counts.items():
        if not resolve("residual", q, n):
            res_unclassified += n
    if res_unclassified:
        add("Unclassified", "", "-", "-", "no match", "-", res_unclassified)

    # C. reconciliation
    accounted = mapped_total + residual_total
    leftover = total - accounted
    if leftover > 0:
        add("Unclassified", "", "-", "-", "no match", "-", leftover)
    qc_warn = []
    if leftover != 0:
        qc_warn.append(f"read accounting mismatch: total={total}, mapped={mapped_total}, "
                       f"residual={residual_total} (difference {leftover})")

    rows = []
    for (taxon, taxid, rank, group, ev, div), n in acc.items():
        rows.append({"sample": a.sample, "taxon": taxon, "taxid": taxid, "rank": rank,
                     "group": group, "reads": n,
                     "pct_of_total_reads": round(n / total * 100, 4) if total else 0.0,
                     "evidence": ev, "divergence": div})
    rows.sort(key=lambda r: (-r["reads"]))
    cols = ["sample", "taxon", "taxid", "rank", "group", "reads", "pct_of_total_reads", "evidence", "divergence"]
    with open(a.out_rows, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    plant = sum(r["reads"] for r in rows if r["group"] == "Viridiplantae")
    other = sum(r["reads"] for r in rows if r["group"] not in ("Viridiplantae", "-", "Unknown"))
    cs_c, cs_r = load_json(a.control_status_contigs), load_json(a.control_status_residual)
    asm = {}
    for p in a.assembly_status:
        asm[os.path.basename(p).split(".")[0]] = open(p).read().strip() if os.path.exists(p) else "NA"
    qc = {"sample": a.sample, "mode": mode, "raw_reads": rc.get("raw_reads"),
          "post_fastp_reads": total, "mapped_to_contigs": mapped_total,
          "residual_reads": residual_total, "human_screen_reads": human,
          "plant_reads": plant, "plant_pct": round(plant / total * 100, 3) if total else 0,
          "non_plant_assigned_reads": other,
          "unclassified_reads": sum(r["reads"] for r in rows if r["taxon"] == "Unclassified"),
          "assembly_status": asm,
          "core_nt_control_contigs": cs_c, "core_nt_control_residual": cs_r,
          "residual_unqueried_reads": rj.get("unqueried_reads", 0),
          "rdna_assigned_reads": stats.get("rdna_assigned_reads_contigs", 0) + stats.get("rdna_assigned_reads_residual", 0),
          "rdna_plant_reads": sum(r["reads"] for r in rows
                                  if r["group"] == "Viridiplantae" and r["evidence"].startswith("rDNA DB")),
          "rdna_residual_unqueried_reads": load_json(a.rdna_residual_json).get("unqueried_reads", 0)
                                           if a.rdna_residual_json else 0,
          "warnings": qc_warn, **{k: v for k, v in stats.items()}}
    json.dump(qc, open(a.out_qc, "w"), indent=1)
    for w_ in qc_warn:
        print(f"WARNING [{a.sample}] {w_}", file=sys.stderr)


if __name__ == "__main__":
    main()
