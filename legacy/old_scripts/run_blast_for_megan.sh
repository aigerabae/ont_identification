#!/bin/bash
# run_blast_for_megan.sh
#
# Same submission pattern as run_blast_all_samples.sh (one combined query
# per sample, resumable, rate-limited against NCBI's remote server) but
# outputs the classical 12-column BLAST tabular format MEGAN's blast2lca
# requires -- NOT the custom-column format your other scripts use.
#
# This is a SEPARATE set of output files (different suffix) so it doesn't
# overwrite or interfere with your existing blast_results/*_blast.tsv used
# by summarize_blast_hits.py / lca_assign.py.
#
# Usage:
#   ./run_blast_for_megan.sh [pipeline_out_dir] [output_dir]

set -euo pipefail

PIPELINE_OUT="${1:-pipeline_out}"
OUTDIR="${2:-blast_results_megan}"
DB="core_nt"
SLEEP_BETWEEN_JOBS=30

mkdir -p "$OUTDIR"

shopt -s nullglob
sample_dirs=("$PIPELINE_OUT"/*/its_extracted)

if [[ ${#sample_dirs[@]} -eq 0 ]]; then
    echo "No */its_extracted directories found under $PIPELINE_OUT — check the path." >&2
    exit 1
fi

for sample_dir in "${sample_dirs[@]}"; do
    sample=$(basename "$(dirname "$sample_dir")")
    combined_query="$OUTDIR/${sample}_query.fasta"
    out_file="$OUTDIR/${sample}_blast_megan.tsv"

    if [[ -s "$out_file" ]]; then
        echo "[skip] $sample already has MEGAN-format results: $out_file"
        continue
    fi

    its_files=("$sample_dir"/*_its.fasta)
    if [[ ${#its_files[@]} -eq 0 ]]; then
        echo "[skip] $sample: no *_its.fasta files in $sample_dir"
        continue
    fi

    cat "${its_files[@]}" > "$combined_query"
    n_seqs=$(grep -c ">" "$combined_query" || true)
    echo "[submit] $sample: $n_seqs sequence(s) -> $DB (MEGAN-format, task=blastn, remote)"

    if blastn -remote -db "$DB" \
        -task blastn \
        -query "$combined_query" \
        -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -max_target_seqs 10 \
        -out "$out_file"
    then
        n_hits=$(wc -l < "$out_file")
        echo "[done]   $sample: $n_hits hit rows -> $out_file"
    else
        echo "[FAILED] $sample — re-run this script to retry (completed samples are skipped)." >&2
        rm -f "$out_file"
    fi

    sleep "$SLEEP_BETWEEN_JOBS"
done

echo ""
echo "All samples processed. MEGAN-format BLAST results in: $OUTDIR/"
echo "Next: run blast2lca on each ${OUTDIR}/*_blast_megan.tsv file."
