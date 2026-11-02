#!/usr/bin/env python3
"""
genus_rollup.py

Collapses final_report.tsv's species-level rows into genus-level totals.
This matters because within a genus with many close congeners/subspecies
in GenBank and a divergent/degraded query sequence (exactly your Hippophae
and Rhaponticum situation), BLAST's single top hit can flip between several
near-tied reference records across different contigs/samples — producing
what LOOKS like many different species when it's really one biological
signal reported under different labels.

"Uncultured ..." / "environmental sample" hits are bucketed separately
(as "Uncultured/environmental (unresolved)") rather than by genus, since
there's no genus to roll up to and they're a distinct interpretive category
(background/contamination signal, not a named ingredient).

Usage:
    python3 genus_rollup.py --report final_report.tsv --out genus_summary.tsv
"""
import argparse
import csv
import re
from collections import defaultdict

_GENUS_RE = re.compile(r"^([A-Z][a-z]+)\b")
_UNCULTURED_RE = re.compile(r"^uncultured", re.IGNORECASE)

# Known taxonomic synonyms that would otherwise show up as separate genera.
# Leuzea carthamoides / Rhaponticum carthamoides are the SAME species under
# two generic names still both in circulation in GenBank records — this is
# documented synonymy, not a guess, so it's safe to auto-merge.
#
# Deliberately NOT included: Serratula, Stemmacantha, Acroptilon. These are
# genuinely distinct genera in the same tribe (Cardueae) as Rhaponticum —
# closely related enough that a short/divergent contig can tie-break onto
# one of them by BLAST score, but merging them automatically would assert a
# species identity the data doesn't actually establish. If these show up
# with meaningful abundance, check the underlying contig(s) by hand rather
# than trusting an automatic merge.
GENUS_SYNONYMS = {
    "Leuzea": "Rhaponticum",
}


def genus_of(species: str) -> str:
    name = re.sub(r"\s*\(from title.*?\)\s*$", "", species).strip()
    if _UNCULTURED_RE.match(name):
        return "Uncultured/environmental (unresolved)"
    m = _GENUS_RE.match(name)
    if not m:
        return name
    genus = m.group(1)
    return GENUS_SYNONYMS.get(genus, genus)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", required=True, help="Path to final_report.tsv")
    ap.add_argument("--out", required=True, help="Path to write the genus-level summary TSV")
    args = ap.parse_args()

    # per (sample, genus): summed depth, set of species labels seen, n contigs
    per_sample_genus_depth: dict[tuple[str, str], float] = defaultdict(float)
    per_sample_genus_species: dict[tuple[str, str], set] = defaultdict(set)
    per_sample_genus_contigs: dict[tuple[str, str], int] = defaultdict(int)
    sample_total_depth: dict[str, float] = defaultdict(float)

    with open(args.report, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            species = row["species"]
            genus = genus_of(species)
            depth = float(row["summed_depth"])
            n_contigs = int(row["n_contigs"])

            key = (sample, genus)
            per_sample_genus_depth[key] += depth
            per_sample_genus_species[key].add(species)
            per_sample_genus_contigs[key] += n_contigs
            sample_total_depth[sample] += depth

    out_rows = []
    for (sample, genus), depth in sorted(per_sample_genus_depth.items()):
        total = sample_total_depth[sample]
        pct = (depth / total * 100) if total > 0 else 0.0
        species_labels = per_sample_genus_species[(sample, genus)]
        out_rows.append(
            {
                "sample": sample,
                "genus": genus,
                "n_species_labels_collapsed": len(species_labels),
                "n_contigs_total": per_sample_genus_contigs[(sample, genus)],
                "summed_depth": round(depth, 2),
                "pct_of_sample_its_depth": round(pct, 3),
                "species_labels_seen": ";".join(sorted(species_labels)),
            }
        )

    out_cols = [
        "sample", "genus", "n_species_labels_collapsed", "n_contigs_total",
        "summed_depth", "pct_of_sample_its_depth", "species_labels_seen",
    ]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    # cross-sample prevalence summary, printed directly (most useful quick view)
    genus_sample_count: dict[str, int] = defaultdict(int)
    genus_avg_pct: dict[str, list] = defaultdict(list)
    for row in out_rows:
        genus_sample_count[row["genus"]] += 1
        genus_avg_pct[row["genus"]].append(row["pct_of_sample_its_depth"])

    print(f"Wrote {len(out_rows)} sample-genus rows -> {args.out}\n")
    print(f"{'Genus':<45} {'#samples':>9} {'avg % of ITS depth':>20}")
    print("-" * 76)
    for genus, count in sorted(genus_sample_count.items(), key=lambda kv: kv[1], reverse=True):
        avg_pct = sum(genus_avg_pct[genus]) / len(genus_avg_pct[genus])
        print(f"{genus:<45} {count:>9} {avg_pct:>19.3f}%")

    watch_genera = {"Serratula", "Stemmacantha", "Acroptilon"} & set(genus_sample_count)
    if watch_genera:
        print(
            f"\nNote: {', '.join(sorted(watch_genera))} appeared above. These are "
            "close Cardueae-tribe relatives of Rhaponticum, not confirmed synonyms — "
            "NOT auto-merged. If their abundance is non-trivial, check the underlying "
            "contig(s) (species_labels_seen / n_contigs_total columns) by hand before "
            "deciding whether they're a real separate taxon or another Rhaponticum "
            "tie-break."
        )


if __name__ == "__main__":
    main()
