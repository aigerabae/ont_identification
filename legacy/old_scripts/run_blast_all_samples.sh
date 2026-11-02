#!/bin/bash
# run_blast_all_samples.sh
#
# Submits ONE combined BLASTN job per sample (all its_extracted/*_its.fasta
# contigs for that sample concatenated into one query) against NCBI's
# remote core_nt database — same database your manual web BLAST runs used.
#
# Why combined per-sample rather than per-contig: NCBI asks scripts not to
# hammer their servers with many small, back-to-back requests. Submitting
# one multi-sequence query per sample (18 jobs total, not ~100+) is both
# faster for you and respectful of their usage guidelines
# (https://www.ncbi.nlm.nih.gov/books/NBK25497/).
#
# -task blastn (not the default megablast) is used deliberately: megablast
# is tuned for near-identical/same-species matches and would likely miss
# the divergent ~86-87%-identity hits you got manually on the web (which
# defaults to blastn for "somewhat similar sequences").
#
# Usage:
#   ./run_blast_all_samples.sh [pipeline_out_dir] [output_dir]
#
# Resumable: already-completed samples (output file exists) are skipped,
# so a crashed/interrupted run can just be re-launched.

set -euo pipefail

PIPELINE_OUT="${1:-pipeline_out}"
OUTDIR="${2:-blast_results}"
DB="core_nt"
SLEEP_BETWEEN_JOBS=30   # seconds; be polite to NCBI's shared servers

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
    out_file="$OUTDIR/${sample}_blast.tsv"

    if [[ -s "$out_file" ]]; then
        echo "[skip] $sample already has results: $out_file"
        continue
    fi

    its_files=("$sample_dir"/*_its.fasta)
    if [[ ${#its_files[@]} -eq 0 ]]; then
        echo "[skip] $sample: no *_its.fasta files in $sample_dir"
        continue
    fi

    cat "${its_files[@]}" > "$combined_query"
    n_seqs=$(grep -c ">" "$combined_query" || true)
    echo "[submit] $sample: $n_seqs sequence(s) -> $DB (task=blastn, remote)"

    # Columns: contig id | subject accession | subject taxonomy ID(s) |
    #          %identity | alignment length | %query coverage per subject |
    #          e-value | bit score | subject title
    #
    # staxids gives us the NCBI taxonomy ID directly from BLAST -- no local
    # taxdb needed for this (only needed to resolve an ID to a readable
    # name, which lca_assign.py does separately via a local taxdump).
    # This removes the need for the separate fetch_staxids.py lookup step
    # for any sample processed with this updated format.
    if blastn -remote -db "$DB" \
        -task blastn \
        -query "$combined_query" \
        -outfmt "6 qseqid sacc staxids pident length qcovs evalue bitscore stitle" \
        -max_target_seqs 10 \
        -out "$out_file"
    then
        n_hits=$(wc -l < "$out_file")
        echo "[done]   $sample: $n_hits hit rows -> $out_file"
    else
        echo "[FAILED] $sample — see above for blastn's error. Re-run this script to retry" \
             "(already-completed samples will be skipped automatically)." >&2
        rm -f "$out_file"   # don't leave a partial/empty file that looks "done"
    fi

    sleep "$SLEEP_BETWEEN_JOBS"
done

echo ""
echo "All samples processed. Results in: $OUTDIR/"
echo "Next: python3 scripts/summarize_blast_hits.py --blast-dir $OUTDIR --out $OUTDIR/summary.tsv"
