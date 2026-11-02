#!/usr/bin/env python3
"""
merge_reports.py -- combine per-sample rows into the deliverable tables.

  final_table.tsv           exactly the requested 5 columns
                            (Sample | Taxon | Rank | Reads | % of reads),
                            reads of the same taxon from different evidence lines summed
  final_table_detailed.tsv  one row per (taxon, evidence line): group, taxid,
                            divergence class, % of plant reads
  qc_summary.tsv            one row per sample: mode, assembly status, read
                            accounting, human fraction, control coverage
--lang ru|en selects column headers and rank names of final_table.tsv.
"""
import argparse, csv, json, re, sys
from collections import defaultdict

RANK_RU = {"species": "вид", "genus": "род", "family": "семейство", "order": "порядок",
           "class": "класс", "phylum": "тип/отдел", "kingdom": "царство", "tribe": "триба",
           "subtribe": "подтриба", "subfamily": "подсемейство", "subclass": "подкласс",
           "suborder": "подпорядок", "subgenus": "подрод", "-": "-"}
HEAD = {"ru": ["Образец", "Таксон", "Таксономический ранг", "Количество прочтений", "Доля прочтений, %"],
        "en": ["Sample", "Taxon", "Rank", "Reads", "% of reads"]}


def natkey(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", nargs="+", required=True)
    ap.add_argument("--qc", nargs="+", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--lang", choices=["ru", "en"], default="ru")
    a = ap.parse_args()

    rows = []
    for p in a.rows:
        with open(p) as fh:
            rows.extend(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["reads"] = int(r["reads"])
    qcs = [json.load(open(p)) for p in a.qc]
    total = {q["sample"]: q["post_fastp_reads"] for q in qcs}
    plant_total = defaultdict(int)
    for r in rows:
        if r["group"] == "Viridiplantae":
            plant_total[r["sample"]] += r["reads"]

    SPECIAL = {"Mapped to contig, no confident taxon": 1, "Unclassified": 2}
    order = lambda r: (natkey(r["sample"]), SPECIAL.get(r["taxon"], 0), -r["reads"], r["taxon"])

    # detailed
    det = sorted(rows, key=order)
    with open(f"{a.out_dir}/final_table_detailed.tsv", "w", newline="") as fh:
        cols = ["sample", "taxon", "taxid", "rank", "group", "reads", "pct_of_total_reads",
                "pct_of_plant_reads", "evidence", "divergence"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in det:
            pp = ""
            if r["group"] == "Viridiplantae" and plant_total[r["sample"]]:
                pp = round(r["reads"] / plant_total[r["sample"]] * 100, 3)
            w.writerow({**r, "pct_of_plant_reads": pp})

    # collapsed, requested format
    agg = defaultdict(int)
    for r in rows:
        agg[(r["sample"], r["taxon"], r["rank"])] += r["reads"]
    flat = [{"sample": s, "taxon": t, "rank": k, "reads": n} for (s, t, k), n in agg.items()]
    flat.sort(key=order)
    with open(f"{a.out_dir}/final_table.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(HEAD[a.lang])
        for r in flat:
            tot = total.get(r["sample"], 0)
            rank = RANK_RU.get(r["rank"], r["rank"]) if a.lang == "ru" else r["rank"]
            w.writerow([r["sample"], r["taxon"], rank, r["reads"],
                        f"{r['reads'] / tot * 100:.2f}" if tot else "0.00"])

    # qc
    with open(f"{a.out_dir}/qc_summary.tsv", "w", newline="") as fh:
        cols = ["sample", "mode", "raw_reads", "post_fastp_reads", "mapped_to_contigs",
                "residual_reads", "human_screen_reads", "plant_reads", "plant_pct",
                "non_plant_assigned_reads", "unclassified_reads", "assembly_long",
                "assembly_short", "rdna_assigned_reads", "rdna_plant_reads", "warnings"]
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for q in sorted(qcs, key=lambda q: natkey(q["sample"])):
            asm = q.get("assembly_status", {})
            w.writerow([q["sample"], q["mode"], q.get("raw_reads"), q["post_fastp_reads"],
                        q["mapped_to_contigs"], q["residual_reads"], q["human_screen_reads"],
                        q["plant_reads"], q["plant_pct"], q["non_plant_assigned_reads"],
                        q["unclassified_reads"], asm.get("longinsert", "NA"),
                        asm.get("shortinsert", "NA"), q.get("rdna_assigned_reads", 0),
                        q.get("rdna_plant_reads", 0), "; ".join(q.get("warnings", []))])
    print(f"wrote final_table.tsv ({len(flat)} rows) for {len(qcs)} samples", file=sys.stderr)


if __name__ == "__main__":
    main()
