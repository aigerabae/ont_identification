#!/usr/bin/env python3
"""
summarize_blast_hits.py

Parses the tabular BLAST output produced by run_blast_all_samples.sh
(one <sample>_blast.tsv per sample, columns:
 qseqid sacc sscinames scomnames pident length qcovs evalue bitscore stitle)
and produces one master table: for every contig (qseqid) in every sample,
the top N hits ranked by bit score, deduplicated by scientific name (so a
contig that matches 10 near-identical accessions of the same species shows
up as one row, not ten) — mirroring how you were reading the NCBI web
"Descriptions" tab by eye.

Usage:
    python3 summarize_blast_hits.py --blast-dir blast_results --out blast_results/summary.tsv
    python3 summarize_blast_hits.py --blast-dir blast_results --out blast_results/summary.tsv --top-n 3
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict

# Matches a leading binomial (+ optional subsp./var.) at the start of a
# GenBank record title, e.g.:
#   "Hippophae rhamnoides subsp. mongolica voucher ..." -> "Hippophae rhamnoides subsp. mongolica"
#   "PREDICTED: Rhaponticum carthamoides from Russia ..." -> "Rhaponticum carthamoides"
#   "Leuzea carthamoides voucher ALT202210 ..."           -> "Leuzea carthamoides"
_ORGANISM_RE = re.compile(
    r"^(?:PREDICTED:\s*)?"
    r"([A-Z][a-z]+ [a-z\-]+(?: subsp\. [a-z\-]+| var\. [a-z\-]+)?)"
)


def organism_from_title(title: str) -> str:
    """Best-effort organism name extracted from a GenBank stitle, used as a
    fallback when sscinames is blank (i.e. local taxdb isn't installed)."""
    m = _ORGANISM_RE.match(title.strip())
    return m.group(1) if m else title.strip()[:60]


_MISSING_VALUES = {"", "n/a", "na", "-"}


def is_missing(value: str) -> bool:
    return value.strip().lower() in _MISSING_VALUES

COLUMNS = [
    "qseqid", "sacc", "sscinames", "scomnames",
    "pident", "length", "qcovs", "evalue", "bitscore", "stitle",
]


def parse_sample_tsv(path: str) -> list[dict]:
    rows = []
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for fields in reader:
            if len(fields) != len(COLUMNS):
                # tolerate stitle containing extra tabs by rejoining the tail
                if len(fields) > len(COLUMNS):
                    fields = fields[: len(COLUMNS) - 1] + ["\t".join(fields[len(COLUMNS) - 1 :])]
                else:
                    continue
            row = dict(zip(COLUMNS, fields))
            try:
                row["pident"] = float(row["pident"])
                row["qcovs"] = float(row["qcovs"])
                row["evalue"] = float(row["evalue"])
                row["bitscore"] = float(row["bitscore"])
            except ValueError:
                continue
            rows.append(row)
    return rows


def top_hits_per_contig(rows: list[dict], top_n: int) -> dict[str, list[dict]]:
    by_contig: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_contig[r["qseqid"]].append(r)

    result: dict[str, list[dict]] = {}
    for contig, hits in by_contig.items():
        hits.sort(key=lambda r: r["bitscore"], reverse=True)
        seen_species: set[str] = set()
        deduped: list[dict] = []
        for h in hits:
            # sscinames can be a ';'-separated list when one HSP matches
            # several identical accessions from different species records.
            # If it's blank OR the literal string "N/A" (what blastn prints
            # when local taxdb isn't installed), fall back to parsing the
            # organism name out of the GenBank title instead — and crucially,
            # make that fallback name the dedup key too, so hits to genuinely
            # different species don't all collapse under a shared "N/A" key.
            raw_sscinames = h["sscinames"].strip()
            if is_missing(raw_sscinames):
                species_key = organism_from_title(h["stitle"])
                h["sscinames"] = species_key + " (from title, taxdb not installed)"
            else:
                species_key = raw_sscinames
            if species_key in seen_species:
                continue
            seen_species.add(species_key)
            deduped.append(h)
            if len(deduped) >= top_n:
                break
        result[contig] = deduped
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--blast-dir", required=True, help="Directory containing <sample>_blast.tsv files")
    ap.add_argument("--out", required=True, help="Path to write the master summary TSV")
    ap.add_argument("--top-n", type=int, default=5, help="Top N distinct-species hits to keep per contig (default 5)")
    args = ap.parse_args()

    tsv_files = sorted(glob.glob(os.path.join(args.blast_dir, "*_blast.tsv")))
    if not tsv_files:
        sys.exit(f"No *_blast.tsv files found in {args.blast_dir}")

    out_rows = []
    n_contigs_total = 0
    n_contigs_no_hits = 0

    for tsv_path in tsv_files:
        sample = os.path.basename(tsv_path).removesuffix("_blast.tsv")
        rows = parse_sample_tsv(tsv_path)
        top_hits = top_hits_per_contig(rows, args.top_n)

        for contig, hits in top_hits.items():
            n_contigs_total += 1
            if not hits:
                n_contigs_no_hits += 1
                continue
            for rank, h in enumerate(hits, start=1):
                out_rows.append(
                    {
                        "sample": sample,
                        "contig": contig,
                        "rank": rank,
                        "scientific_name": h["sscinames"],
                        "common_name": h["scomnames"],
                        "pident": h["pident"],
                        "qcovs": h["qcovs"],
                        "evalue": h["evalue"],
                        "bitscore": h["bitscore"],
                        "accession": h["sacc"],
                        "title": h["stitle"],
                    }
                )

    out_cols = [
        "sample", "contig", "rank", "scientific_name", "common_name",
        "pident", "qcovs", "evalue", "bitscore", "accession", "title",
    ]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Parsed {len(tsv_files)} sample BLAST result file(s).")
    print(f"{n_contigs_total} contig(s) had BLAST hits; {n_contigs_no_hits} had none.")
    print(f"Wrote {len(out_rows)} ranked hit rows -> {args.out}")
    print("\nTip: filter rank==1 for a quick top-call-per-contig view, e.g.:")
    print(f"  awk -F'\\t' 'NR==1 || $3==1' {args.out} | column -t -s$'\\t' | less -S")


if __name__ == "__main__":
    main()
