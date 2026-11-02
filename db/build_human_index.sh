#!/bin/bash
# Build the bowtie2 index of the human reference used by the human-contamination screen.
# Usage: db/build_human_index.sh [outdir=db/human] [threads=16]
# For tests, set HUMAN_FASTA=/path/to/small.fasta to skip the download.
set -euo pipefail
OUT="${1:-db/human}"; THREADS="${2:-16}"
mkdir -p "$OUT"
FA="${HUMAN_FASTA:-}"
if [ -z "$FA" ]; then
    FA="$OUT/GRCh38.fna.gz"
    URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    # resumable: -C - continues a partial file; retried until gzip says the file is complete
    for attempt in $(seq 1 30); do
        curl -fL -C - --retry 5 --retry-all-errors --retry-delay 5 -o "$FA" "$URL" || true
        if gzip -t "$FA" 2>/dev/null; then break; fi
        echo "[download incomplete, attempt $attempt/30: resuming in 10 s]"; sleep 10
    done
    gzip -t "$FA" || { echo "ERROR: download still incomplete; run the script again to resume" >&2; exit 1; }
fi
bowtie2-build --threads "$THREADS" "$FA" "$OUT/GRCh38"
echo "human reference: $FA  built $(date -u +%F)" > "$OUT/HUMAN_INDEX.txt"
echo "OK: index prefix $OUT/GRCh38"
