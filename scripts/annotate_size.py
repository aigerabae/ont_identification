#!/usr/bin/env python3
"""Append ';size=N' (from SPAdes cov_ value x scale) to each header, as vsearch
--uchime_denovo needs abundance annotations."""
import re, sys

inp, out, scale = sys.argv[1], sys.argv[2], float(sys.argv[3])
rx = re.compile(r"cov_([\d.]+)")
with open(inp) as fi, open(out, "w") as fo:
    for line in fi:
        if line.startswith(">"):
            m = rx.search(line)
            size = max(int(float(m.group(1)) * scale), 1) if m else 1
            line = line.rstrip("\n") + f";size={size}\n"
        fo.write(line)
