#!/usr/bin/env python3
"""
human_fraction_kraken.py

From Kraken2 .kreport files, report how many reads (fragments) are human
vs. other-classified vs. unclassified, per sample.

Uses the Kraken2 report (NOT the Bracken output: Bracken re-distributes and
drops reads, so it cannot give an honest human fraction).

.kreport columns: percent, clade_reads, direct_reads, rank, taxid, name
  - "unclassified" row  -> reads Kraken2 could not classify
  - "root" row          -> all classified reads
  - "Homo" / "Homo sapiens" -> human (matched by NAME, because SILVA taxids
    are not NCBI taxids)

Usage:
    python3 human_fraction_kraken.py kraken2_fallback_results/kraken2_out/{18S-1_S37,18S-3_S41,18S-9_S38,18S-31_S33}.kreport
"""
import sys
from pathlib import Path


def parse(path):
    unclassified = classified = 0
    homo = homo_sapiens = 0
    found_homo = False
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        clade = int(f[1])
        name = f[5].strip()
        if name == "unclassified":
            unclassified = clade
        elif name == "root":
            classified = clade
        elif name == "Homo":
            homo = clade
            found_homo = True
        elif name == "Homo sapiens":
            homo_sapiens = clade
            found_homo = True
    human = homo if homo else homo_sapiens
    return unclassified, classified, human, found_homo


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    hdr = ("sample", "total", "human", "human_%total", "other_classified",
           "other_%total", "unclassified", "unclass_%total", "human_%classified")
    print("\t".join(hdr))
    for p in sys.argv[1:]:
        sample = Path(p).name.replace(".kreport", "")
        unc, cls, human, found = parse(p)
        total = unc + cls
        other = cls - human
        pct = lambda x, d: f"{100 * x / d:.2f}" if d else "NA"
        print("\t".join(map(str, (
            sample, total, human, pct(human, total), other, pct(other, total),
            unc, pct(unc, total), pct(human, cls)))))
        if not found:
            print(f"  note: no 'Homo' row in {p} (SILVA SSU is mostly rRNA of "
                  "all life; human may be absent or labelled differently)",
                  file=sys.stderr)


if __name__ == "__main__":
    main()
