#!/usr/bin/env python3
"""
lca_assign.py

Replaces top-hit-only species assignment with LCA (lowest common ancestor)
resolution: for each contig, instead of trusting whichever hit happened to
rank #1, this groups together every hit within --tie-fraction of the top
bit score and reports the taxonomic rank actually supported by ALL of
them. A clean single winner still resolves to species; a genuine tie
(like Serratula coronata vs. four Rhaponticoides species, all at
identical bitscore/identity/e-value) correctly backs off to whatever
rank they share (e.g. subtribe), instead of arbitrarily naming one.

Supports TWO BLAST output schemas, auto-detected per file:

  NEW (9 columns, from the current run_blast_all_samples.sh):
    qseqid sacc staxids pident length qcovs evalue bitscore stitle
    -- taxids come directly from BLAST's own -outfmt staxids field, no
    separate accession->taxid lookup needed.

  OLD (10 columns, from samples BLASTed before staxids was added):
    qseqid sacc sscinames scomnames pident length qcovs evalue bitscore stitle
    -- requires --taxid-lookup (from fetch_staxids.py) to resolve taxids,
    since this schema has no staxids column.

This means you do NOT need to re-run BLAST for samples already completed
under the old format -- just keep passing --taxid-lookup for those, and
any new samples run under the updated format are handled automatically
without it.

Requires a local NCBI taxdump directory (nodes.dmp + names.dmp) --
download once from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

Usage:
    # all-new-format samples, no lookup file needed:
    python3 lca_assign.py \
        --blast-dir blast_results \
        --taxdump-dir taxdump \
        --out blast_results/lca_summary.tsv

    # mixed old+new format samples in the same directory:
    python3 lca_assign.py \
        --blast-dir blast_results \
        --taxid-lookup blast_results/accession_taxid_lookup.tsv \
        --taxdump-dir taxdump \
        --out blast_results/lca_summary.tsv
"""
import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from taxdump import Taxonomy

NEW_COLUMNS = ["qseqid", "sacc", "staxids", "pident", "length", "qcovs", "evalue", "bitscore", "stitle"]
OLD_COLUMNS = ["qseqid", "sacc", "sscinames", "scomnames", "pident", "length", "qcovs", "evalue", "bitscore", "stitle"]

NUMERIC_FIELDS = ["pident", "qcovs", "evalue", "bitscore"]


def load_accession_taxid_lookup(path: str) -> dict[str, int]:
    lookup: dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2:
                continue
            acc, taxid_str = parts
            try:
                lookup[acc] = int(taxid_str)
            except ValueError:
                continue
    return lookup


def _parse_row(fields: list[str], columns: list[str]) -> dict | None:
    if len(fields) != len(columns):
        if len(fields) > len(columns):
            fields = fields[: len(columns) - 1] + ["\t".join(fields[len(columns) - 1 :])]
        else:
            return None
    row = dict(zip(columns, fields))
    try:
        for f in NUMERIC_FIELDS:
            row[f] = float(row[f])
    except ValueError:
        return None
    return row


def parse_sample_tsv(path: str) -> tuple[list[dict], str]:
    """Returns (rows, schema) where schema is 'new' or 'old'. Detected from
    the first line that successfully parses under either schema."""
    schema = None
    columns = None
    rows = []
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for fields in reader:
            if schema is None:
                # try new schema first (9 cols), then old (10 cols)
                row = _parse_row(fields, NEW_COLUMNS)
                if row is not None:
                    schema, columns = "new", NEW_COLUMNS
                else:
                    row = _parse_row(fields, OLD_COLUMNS)
                    if row is not None:
                        schema, columns = "old", OLD_COLUMNS
                if row is None:
                    continue  # neither schema matched this line; skip and keep trying
            else:
                row = _parse_row(fields, columns)
                if row is None:
                    continue
            rows.append(row)
    return rows, (schema or "unknown")


def taxids_from_row(row: dict, schema: str, acc_to_taxid: dict[str, int] | None) -> list[int]:
    if schema == "new":
        raw = row.get("staxids", "").strip()
        if not raw or raw.upper() in ("N/A", "NA", "-"):
            return []
        out = []
        for part in raw.split(";"):
            part = part.strip()
            if part.isdigit():
                out.append(int(part))
        return out
    else:  # old schema: look up by accession
        if not acc_to_taxid:
            return []
        taxid = acc_to_taxid.get(row["sacc"])
        return [taxid] if taxid is not None else []


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--blast-dir", required=True)
    ap.add_argument("--taxid-lookup", default=None,
                     help="Accession->taxid lookup (from fetch_staxids.py). Only required if any "
                          "*_blast.tsv files are in the OLD (pre-staxids) 10-column format.")
    ap.add_argument("--taxdump-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tie-fraction", type=float, default=0.02,
                     help="Hits within this fraction of the top bitscore are treated as tied (default 0.02 = within 2%%)")
    ap.add_argument("--min-bitscore", type=float, default=100.0)
    ap.add_argument("--max-evalue", type=float, default=1e-15)
    args = ap.parse_args()

    tax = Taxonomy(args.taxdump_dir)
    acc_to_taxid = load_accession_taxid_lookup(args.taxid_lookup) if args.taxid_lookup else None

    tsv_files = sorted(glob.glob(os.path.join(args.blast_dir, "*_blast.tsv")))
    if not tsv_files:
        sys.exit(f"No *_blast.tsv files found in {args.blast_dir}")

    out_rows = []
    n_no_taxid = 0
    n_old_schema_no_lookup = 0
    unmapped_accessions: set[str] = set()
    schema_counts: dict[str, int] = defaultdict(int)

    for tsv_path in tsv_files:
        sample = os.path.basename(tsv_path).removesuffix("_blast.tsv")
        rows, schema = parse_sample_tsv(tsv_path)
        schema_counts[schema] += 1

        if schema == "old" and acc_to_taxid is None:
            print(f"WARNING: {sample} is in the OLD BLAST format (no staxids column) but no "
                  f"--taxid-lookup was given -- skipping this sample. Run fetch_staxids.py for it, "
                  f"or pass --taxid-lookup.", file=sys.stderr)
            n_old_schema_no_lookup += 1
            continue

        by_contig: dict[str, list[dict]] = defaultdict(list)
        for r in rows:
            by_contig[r["qseqid"]].append(r)

        for contig, hits in by_contig.items():
            hits.sort(key=lambda r: r["bitscore"], reverse=True)
            top = hits[0]
            if top["bitscore"] < args.min_bitscore or top["evalue"] > args.max_evalue:
                continue  # whole contig fails the confidence floor

            tie_threshold = top["bitscore"] * (1 - args.tie_fraction)
            tied_hits = [h for h in hits if h["bitscore"] >= tie_threshold]

            taxids: list[int] = []
            for h in tied_hits:
                hit_taxids = taxids_from_row(h, schema, acc_to_taxid)
                if not hit_taxids:
                    n_no_taxid += 1
                    unmapped_accessions.add(h["sacc"])
                    continue
                taxids.extend(hit_taxids)

            if not taxids:
                out_rows.append({
                    "sample": sample, "contig": contig,
                    "assigned_rank": "unresolved (no taxid)", "assigned_name": top["stitle"][:80],
                    "assigned_taxid": "", "tie_group_size": len(tied_hits),
                    "member_species": "", "best_pident": top["pident"],
                    "best_bitscore": top["bitscore"], "best_evalue": top["evalue"],
                    "best_accession": top["sacc"],
                })
                continue

            lca_taxid = tax.lca(taxids)
            member_names = sorted({tax.name(t) for t in taxids})

            out_rows.append({
                "sample": sample, "contig": contig,
                "assigned_rank": tax.rank(lca_taxid),
                "assigned_name": tax.name(lca_taxid),
                "assigned_taxid": lca_taxid if lca_taxid is not None else "",
                "tie_group_size": len(tied_hits),
                "member_species": ";".join(member_names),
                "best_pident": top["pident"], "best_bitscore": top["bitscore"],
                "best_evalue": top["evalue"], "best_accession": top["sacc"],
            })

    out_cols = [
        "sample", "contig", "assigned_rank", "assigned_name", "assigned_taxid",
        "tie_group_size", "member_species", "best_pident", "best_bitscore",
        "best_evalue", "best_accession",
    ]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    n_species = sum(1 for r in out_rows if r["assigned_rank"] == "species")
    n_backed_off = sum(1 for r in out_rows if r["tie_group_size"] > 1 and r["assigned_rank"] != "species")
    n_ties_total = sum(1 for r in out_rows if r["tie_group_size"] > 1)

    print(f"Wrote {len(out_rows)} LCA-resolved contig calls -> {args.out}")
    print(f"  Input schema breakdown: {dict(schema_counts)}")
    print(f"  {n_species} resolved to species level")
    print(f"  {n_ties_total} contigs had a tied/near-tied top hit group (tie-fraction={args.tie_fraction})")
    print(f"  {n_backed_off} of those ties were correctly backed off ABOVE species level (this is LCA doing its job)")
    if n_old_schema_no_lookup:
        print(f"  WARNING: {n_old_schema_no_lookup} sample(s) skipped entirely -- old-format BLAST file with no --taxid-lookup given")
    if n_no_taxid:
        print(f"  WARNING: {n_no_taxid} hit(s) across {len(unmapped_accessions)} accession(s) had no resolvable taxid")


if __name__ == "__main__":
    main()
