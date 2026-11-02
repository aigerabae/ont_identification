#!/usr/bin/env python3
"""
taxdump.py -- minimal NCBI taxonomy (nodes.dmp / names.dmp / merged.dmp) with
the helpers the classification steps need: LCA, rank normalisation,
climbing to a coarser rank, and a coarse "group" label (plants / fungi / ...).

Self-contained on purpose (no taxonkit dependency) so behaviour is inspectable.
"""
from __future__ import annotations
import os

# coarse -> fine; ranks not listed are "unranked" (clade, no rank, ...)
RANK_LEVELS = {
    "domain": 0, "superkingdom": 0, "kingdom": 1,
    "phylum": 2, "subphylum": 2, "class": 3, "subclass": 3,
    "order": 4, "suborder": 4, "family": 5, "subfamily": 5,
    "tribe": 6, "subtribe": 6, "genus": 7, "subgenus": 7,
    "species": 8,
}
# Ranks a result may be REPORTED at. NCBI also has nodes such as "unclassified Fusarium",
# "Salvia incertae sedis", "apioid superclade" (clade / no rank) or subtribes/sections; a call that
# lands on one of them is moved up to the nearest standard-rank ancestor.
STANDARD_RANKS = {"superkingdom", "kingdom", "phylum", "class", "order", "family",
                  "subfamily", "tribe", "genus", "species"}
SUBSPECIFIC = {"subspecies", "varietas", "forma", "strain", "isolate",
               "subvariety", "forma specialis", "serotype", "no rank"}

GROUP_TAXIDS = [(33090, "Viridiplantae"), (4751, "Fungi"), (33208, "Metazoa"),
                (2, "Bacteria"), (2157, "Archaea"), (10239, "Viruses")]


def rank_level(rank: str) -> int:
    """-1 for unranked / unknown (always treated as uninformative)."""
    return RANK_LEVELS.get(rank, -1)


class Taxonomy:
    def __init__(self, taxdump_dir: str):
        self.parent_of: dict[int, int] = {}
        self.rank_of: dict[int, str] = {}
        self.name_of: dict[int, str] = {}
        self.merged: dict[int, int] = {}
        self._load_nodes(os.path.join(taxdump_dir, "nodes.dmp"))
        self._load_names(os.path.join(taxdump_dir, "names.dmp"))
        merged = os.path.join(taxdump_dir, "merged.dmp")
        if os.path.exists(merged):
            with open(merged) as fh:
                for line in fh:
                    f = [x.strip() for x in line.split("|")]
                    if len(f) >= 2 and f[0].isdigit() and f[1].isdigit():
                        self.merged[int(f[0])] = int(f[1])

    def _load_nodes(self, path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"nodes.dmp not found at {path}")
        with open(path) as fh:
            for line in fh:
                f = [x.strip() for x in line.split("|")]
                t = int(f[0])
                self.parent_of[t] = int(f[1])
                self.rank_of[t] = f[2]

    def _load_names(self, path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"names.dmp not found at {path}")
        with open(path) as fh:
            for line in fh:
                f = [x.strip() for x in line.split("|")]
                if len(f) >= 4 and f[3] == "scientific name":
                    self.name_of[int(f[0])] = f[1]

    # -- basic ---------------------------------------------------------
    def resolve(self, taxid: int) -> int:
        """Follow merged.dmp (deprecated taxids) to the current taxid."""
        seen = set()
        while taxid in self.merged and taxid not in seen:
            seen.add(taxid)
            taxid = self.merged[taxid]
        return taxid

    def known(self, taxid: int) -> bool:
        return self.resolve(taxid) in self.parent_of

    def lineage(self, taxid: int) -> list[int]:
        """Leaf-first list up to the root."""
        taxid = self.resolve(taxid)
        path, seen, cur = [], set(), taxid
        while cur in self.parent_of and cur not in seen:
            seen.add(cur)
            path.append(cur)
            par = self.parent_of[cur]
            if par == cur:
                break
            cur = par
        return path

    def lca(self, taxids: list[int]) -> int | None:
        valid = [self.resolve(t) for t in taxids if self.known(t)]
        if not valid:
            return None
        lineages = [self.lineage(t) for t in valid]
        others = [set(p) for p in lineages[1:]]
        for node in lineages[0]:
            if all(node in s for s in others):
                return node
        return None

    def name(self, taxid):
        if taxid is None:
            return "Unresolved"
        return self.name_of.get(self.resolve(taxid), f"taxid:{taxid}")

    def rank(self, taxid):
        if taxid is None:
            return "unresolved"
        return self.rank_of.get(self.resolve(taxid), "unknown rank")

    # -- normalisation ------------------------------------------------------
    def to_species(self, taxid: int) -> int:
        """Sub-specific nodes (subspecies, varieties, strains, 'no rank'
        below a species) are reported at the species they belong to, because
        sub-specific labels in GenBank are not reliable identifications."""
        taxid = self.resolve(taxid)
        if self.rank(taxid) not in SUBSPECIFIC:
            return taxid
        for anc in self.lineage(taxid)[1:]:
            r = self.rank(anc)
            if r == "species":
                return anc
            if rank_level(r) >= 0:        # reached genus or coarser: stop
                break
        return taxid

    def standardize(self, taxid):
        """Nearest ancestor-or-self with a standard reporting rank (see STANDARD_RANKS)."""
        if taxid is None:
            return None
        for node in self.lineage(taxid):
            if self.rank(node) in STANDARD_RANKS:
                return node
        return taxid

    def climb_to_level(self, taxid: int, max_level: int) -> int:
        """First node on the path to the root whose (known) rank level is
        <= max_level (0=domain ... 7=genus, 8=species)."""
        for node in self.lineage(taxid):
            lv = rank_level(self.rank(node))
            if 0 <= lv <= max_level:
                return node
        return taxid

    def group(self, taxid: int) -> str:
        lin = set(self.lineage(taxid))
        for tid, label in GROUP_TAXIDS:
            if tid in lin:
                return label
        return "Other" if lin else "Unknown"
