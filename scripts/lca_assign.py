#!/usr/bin/env python3
"""
lca_assign.py -- lowest-common-ancestor taxon assignment from BLAST tabular output.

For every query: take all hits within --tie-fraction of the top bit score, find
the LCA of their taxids. A clean single winner stays at species; a genuine tie
between several species backs off to the shared genus (or higher). Then:

  * sub-specific nodes (subspecies/variety/strain) are reported as their species
  * optional identity caps (--rank-caps species:97,genus:90,family:80): if the
    best identity is below the cap for the assigned rank, the call is moved up
    one rank at a time. This is the guard against the "taxon is missing from
    the reference database, so it is forced onto its nearest relative" failure.

Input columns (outfmt 6):
  qseqid sacc staxids pident length qcovs evalue bitscore stitle
"""
from __future__ import annotations
import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from taxdump import Taxonomy, rank_level  # noqa: E402

COLS = ["qseqid", "sacc", "staxids", "pident", "length", "qcovs", "evalue", "bitscore", "stitle"]
OUT_COLS = ["sample", "query", "assigned_rank", "assigned_name", "assigned_taxid",
            "tie_group_size", "member_taxa", "best_pident", "best_qcov", "best_bitscore",
            "best_evalue", "best_accession", "capped_from"]


def parse_caps(s: str | None) -> dict[str, float]:
    caps = {}
    if s:
        for part in s.split(","):
            if ":" in part:
                k, v = part.split(":")
                caps[k.strip()] = float(v)
    return caps


def parse_line(line):
    f = line.rstrip("\n").split("\t")
    if len(f) < len(COLS):
        return None
    if len(f) > len(COLS):
        f = f[:len(COLS) - 1] + ["\t".join(f[len(COLS) - 1:])]
    r = dict(zip(COLS, f))
    try:
        for k in ("pident", "qcovs", "evalue", "bitscore"):
            r[k] = float(r[k])
    except ValueError:
        return None
    return r


def read_hits(path: str):
    """All hits in memory, keyed by query (insertion order = first appearance). Used as a fallback."""
    hits = defaultdict(list)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return hits
    with open(path) as fh:
        for line in fh:
            r = parse_line(line)
            if r is not None:
                hits[r["qseqid"]].append(r)
    return hits


class NotGrouped(Exception):
    pass


def iter_hit_groups(path: str):
    """Yield (query, hits) one query at a time. BLAST writes all hits of a query consecutively, so only
    one query's hits are held in memory (a 20-million-line file needs ~MBs instead of ~16 GB).
    Raises NotGrouped if a query reappears after its block ended (file not grouped by query)."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return
    done, cur, block = set(), None, []
    with open(path) as fh:
        for line in fh:
            r = parse_line(line)
            if r is None:
                continue
            q = r["qseqid"]
            if q != cur:
                if cur is not None:
                    yield cur, block
                    done.add(cur)
                if q in done:
                    raise NotGrouped(q)
                cur, block = q, []
            block.append(r)
    if cur is not None:
        yield cur, block


def taxids_of(row) -> list[int]:
    out = []
    for p in row["staxids"].replace(",", ";").split(";"):
        p = p.strip()
        if p.isdigit() and int(p) > 0:
            out.append(int(p))
    return out


def apply_caps(tax: Taxonomy, node: int, pident: float, caps: dict[str, float]):
    """Move the call up while the identity is too low for the claimed rank.
    Returns (node, original_rank_or_empty).

    Without a "tribe" cap (the ITS settings): a capped genus goes to subfamily/family (level 5), and the
    family cap applies to tribe and family level alike -- the original behaviour.
    With a "tribe" cap (used by the 18S/28S settings): a capped genus stops at the tribe (level 6; a
    subtribe is standardised to its tribe afterwards), and a tribe-level call needs >= the tribe cap,
    otherwise it moves to subfamily/family, where the family cap applies."""
    original = ""
    changed = True
    has_tribe = "tribe" in caps
    while changed and node is not None:
        changed = False
        rk = tax.rank(node)
        lv = rank_level(rk)
        if rk == "species" and "species" in caps and pident < caps["species"]:
            nxt = tax.climb_to_level(node, 7)
        elif lv == 7 and "genus" in caps and pident < caps["genus"]:
            nxt = tax.climb_to_level(node, 6 if has_tribe else 5)
        elif has_tribe and lv == 6 and pident < caps["tribe"]:
            nxt = tax.climb_to_level(node, 5)
        elif lv in ((5,) if has_tribe else (5, 6)) and "family" in caps and pident < caps["family"]:
            nxt = tax.climb_to_level(node, 4)
        else:
            break
        if nxt != node:
            original = original or rk
            node, changed = nxt, True
    return node, original


def assign_query(hits, tax, tie_fraction, min_bitscore, max_evalue,
                 min_pident, min_qcov, caps):
    """Return a dict (one output row without sample/query) or None."""
    ok = [h for h in hits if h["pident"] >= min_pident and h["qcovs"] >= min_qcov]
    if not ok:
        return None
    ok.sort(key=lambda h: h["bitscore"], reverse=True)
    top = ok[0]
    if top["bitscore"] < min_bitscore or top["evalue"] > max_evalue:
        return None
    tied = [h for h in ok if h["bitscore"] >= top["bitscore"] * (1 - tie_fraction)]
    taxids = []
    for h in tied:
        taxids.extend(taxids_of(h))
    taxids = [t for t in taxids if tax.known(t)]
    if not taxids:
        return {"assigned_rank": "unresolved", "assigned_name": top["stitle"][:80],
                "assigned_taxid": "", "tie_group_size": len(tied), "member_taxa": "",
                "best_pident": top["pident"], "best_qcov": top["qcovs"],
                "best_bitscore": top["bitscore"], "best_evalue": top["evalue"],
                "best_accession": top["sacc"], "capped_from": ""}
    node = tax.lca(taxids)
    node = tax.to_species(node) if node is not None else node
    node = tax.standardize(node)          # "unclassified X", clades, subtribes -> real ranks
    node, capped_from = apply_caps(tax, node, top["pident"], caps)
    node = tax.standardize(node)
    members = sorted({tax.name(tax.to_species(t)) for t in taxids})
    return {"assigned_rank": tax.rank(node), "assigned_name": tax.name(node),
            "assigned_taxid": node if node is not None else "",
            "tie_group_size": len(tied), "member_taxa": ";".join(members[:20]),
            "best_pident": top["pident"], "best_qcov": top["qcovs"],
            "best_bitscore": top["bitscore"], "best_evalue": top["evalue"],
            "best_accession": top["sacc"], "capped_from": capped_from}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--blast", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--taxdump-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tie-fraction", type=float, default=0.02)
    ap.add_argument("--min-bitscore", type=float, default=100)
    ap.add_argument("--max-evalue", type=float, default=1e-15)
    ap.add_argument("--min-pident", type=float, default=0)
    ap.add_argument("--min-qcov", type=float, default=0)
    ap.add_argument("--rank-caps", default="",
                    help="e.g. species:97,genus:90,family:80 (empty = off); optional tribe:NN makes a capped "
                         "genus stop at its tribe (used for the 18S/28S database)")
    a = ap.parse_args()

    tax = Taxonomy(a.taxdump_dir)
    caps = parse_caps(a.rank_caps)

    def run(groups):
        n_q = n_called = 0
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            for q, hs in groups:
                n_q += 1
                row = assign_query(hs, tax, a.tie_fraction, a.min_bitscore, a.max_evalue,
                                   a.min_pident, a.min_qcov, caps)
                if row is None:
                    continue
                row.update({"sample": a.sample, "query": q})
                w.writerow(row)
                n_called += 1
        return n_q, n_called

    try:                                   # streaming: one query in memory at a time
        n_q, n_called = run(iter_hit_groups(a.blast))
    except NotGrouped as e:                # input not grouped by query: fall back to loading it all
        print(f"[{a.sample}] BLAST file not grouped by query (e.g. {e}); using in-memory mode",
              file=sys.stderr)
        n_q, n_called = run(read_hits(a.blast).items())
    print(f"[{a.sample}] {n_q} queries with hits -> {n_called} confident LCA calls -> {a.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
