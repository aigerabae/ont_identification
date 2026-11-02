#!/bin/bash
# Download and unpack the NCBI taxonomy (nodes.dmp, names.dmp, merged.dmp).
# Usage: db/get_taxdump.sh [outdir=db/taxdump]
set -euo pipefail
OUT="${1:-db/taxdump}"
mkdir -p "$OUT"
URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
curl -fL --retry 3 -o "$OUT/taxdump.tar.gz" "$URL"
curl -fL --retry 3 -o "$OUT/taxdump.tar.gz.md5" "$URL.md5"
( cd "$OUT" && md5sum -c taxdump.tar.gz.md5 && tar -xzf taxdump.tar.gz nodes.dmp names.dmp merged.dmp )
echo "taxdump downloaded $(date -u +%F)" > "$OUT/TAXDUMP_DATE.txt"
echo "OK: $OUT/{nodes,names,merged}.dmp"
