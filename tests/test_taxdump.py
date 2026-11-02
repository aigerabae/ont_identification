#!/usr/bin/env python3
"""Unit check: odd NCBI nodes are reported at a standard rank. Run: python3 tests/test_taxdump.py"""
import os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from taxdump import Taxonomy

nodes = [(1,1,"no rank","root"),(10,1,"superkingdom","Eukaryota"),(20,10,"kingdom","Fungi"),
         (30,20,"family","Nectriaceae"),(31,30,"genus","Fusarium"),(32,31,"no rank","unclassified Fusarium"),
         (33,31,"species","Fusarium oxysporum"),(40,10,"family","Apiaceae"),(41,40,"subfamily","Apioideae"),
         (42,41,"clade","apioid superclade"),(43,42,"genus","Arracacia"),(50,40,"genus","Salvia"),
         (51,50,"no rank","Salvia incertae sedis"),(52,40,"tribe","Artemisieae"),(53,52,"subtribe","Artemisiinae")]
d = tempfile.mkdtemp()
with open(f"{d}/nodes.dmp","w") as fn, open(f"{d}/names.dmp","w") as fm:
    for t,p,r,n in nodes:
        fn.write(f"{t}\t|\t{p}\t|\t{r}\t|\t\t|\n"); fm.write(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n")
tax = Taxonomy(d)
exp = {32: 31, 42: 41, 51: 50, 53: 52, 33: 33, 31: 31, 10: 10}
for k, v in exp.items():
    got = tax.standardize(tax.to_species(k))
    assert got == v, f"{tax.name(k)} -> {tax.name(got)} (expected {tax.name(v)})"
print("taxdump standardize: OK")
