#!/bin/bash
# run_blast_no_itsx_fallback.sh
#
# BLASTs the top-N-longest-contig fallback files (from
# extract_top_contigs_no_itsx.py) for samples where ITSx found no ITS
# content. Uses the same staxids-included format as the current
# run_blast_all_samples.sh, so the output plugs directly into your
# existing lca_assign.py / build_final_report_lca.py with no lookup file
# needed.
#
# Output goes to a SEPARATE directory from your main blast_results/ --
# keep these results clearly labeled as the lower-confidence,
# whole-contig-conserved-region fallback, not mixed in with your
# ITS-region-restricted main results, when you write this up.
#
# Usage:
#   ./run_blast_no_itsx_fallback.sh [fallback_fasta_dir] [output_dir]

set -euo pipefail

FALLBACK_DIR="${1:-no_itsx_fallback}"
OUTDIR="${2:-blast_results_no_itsx}"
DB="core_nt"
SLEEP_BETWEEN_JOBS=30

mkdir -p "$OUTDIR"

shopt -s nullglob
fasta_files=("$FALLBACK_DIR"/*_top_contigs.fasta)

if [[ ${#fasta_files[@]} -eq 0 ]]; then
    echo "No *_top_contigs.fasta files found in $FALLBACK_DIR — run extract_top_contigs_no_itsx.py first." >&2
    exit 1
fi

for fasta in "${fasta_files[@]}"; do
    sample=$(basename "$fasta" "_top_contigs.fasta")
    out_file="$OUTDIR/${sample}_blast.tsv"

    if [[ -s "$out_file" ]]; then
        echo "[skip] $sample already has results: $out_file"
        continue
    fi

    n_seqs=$(grep -c ">" "$fasta" || true)
    echo "[submit] $sample: $n_seqs contig(s) -> $DB (task=blastn, remote, no-ITSx fallback)"

    if blastn -remote -db "$DB" \
        -task blastn \
        -query "$fasta" \
        -outfmt "6 qseqid sacc staxids pident length qcovs evalue bitscore stitle" \
        -max_target_seqs 10 \
        -out "$out_file"
    then
        n_hits=$(wc -l < "$out_file")
        echo "[done]   $sample: $n_hits hit rows -> $out_file"
    else
        echo "[FAILED] $sample — re-run this script to retry." >&2
        rm -f "$out_file"
    fi

    sleep "$SLEEP_BETWEEN_JOBS"
done

echo ""
echo "Done. Results in: $OUTDIR/"
echo "Next:"
echo "  python3 scripts/lca_assign.py --blast-dir $OUTDIR --taxdump-dir taxdump --out $OUTDIR/lca_summary.tsv"
echo "  python3 scripts/build_final_report_lca.py --lca-summary $OUTDIR/lca_summary.tsv --pipeline-out pipeline_out --out final_report_no_itsx.tsv"
