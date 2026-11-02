#!/bin/bash
# run_blast_local.sh
#
# Local BLASTN of every sample against the custom ITS database. Output uses
# the SAME 9-column format as the NCBI-remote runners, so lca_assign.py and
# build_final_report_lca.py work unchanged. No NCBI rate limits: all samples
# run in parallel, in minutes.
#
# Modes (what to use as the query):
#   its      pipeline_out/<sample>/its_extracted/*_its.fasta   (main, ITS-restricted)
#   noitsx   no_itsx_fallback/<sample>_top_contigs.fasta       (lower confidence)
#   reads    reads_blast/<sample>_reads.fasta                  (failed-assembly samples)
#
# Usage:
#   ./scripts/run_blast_local.sh <mode> [db_prefix] [pipeline_out] [out_dir] [jobs] [threads_per_job]
# Example:
#   ./scripts/run_blast_local.sh its custom_db/build/db/its_custom pipeline_out blast_results_custom 6 4
#   python3 scripts/lca_assign.py --blast-dir blast_results_custom --taxdump-dir taxdump \
#           --out blast_results_custom/lca_summary.tsv

set -uo pipefail
MODE="${1:?mode: its | noitsx | reads}"
DB="${2:-custom_db/build/db/its_custom}"
PIPELINE_OUT="${3:-pipeline_out}"
OUTDIR="${4:-blast_results_custom}"
JOBS="${5:-6}"
THREADS="${6:-4}"
QDIR="$OUTDIR/_queries"
mkdir -p "$OUTDIR" "$QDIR"

[[ -f "${DB}.nin" || -f "${DB}.nal" ]] || { echo "ERROR: BLAST db not found at $DB (run build_custom_db.sh)" >&2; exit 1; }

shopt -s nullglob
declare -A QUERIES
case "$MODE" in
  its)
    for d in "$PIPELINE_OUT"/*/its_extracted; do
        s=$(basename "$(dirname "$d")")
        files=("$d"/*_its.fasta)
        (( ${#files[@]} )) || { echo "[skip] $s: no ITS fasta files"; continue; }
        cat "${files[@]}" > "$QDIR/${s}.fasta"; QUERIES[$s]="$QDIR/${s}.fasta"
    done ;;
  noitsx)
    for f in no_itsx_fallback/*_top_contigs.fasta; do
        s=$(basename "$f" _top_contigs.fasta); QUERIES[$s]="$f"; done ;;
  reads)
    for f in reads_blast/*_reads.fasta; do
        s=$(basename "$f" _reads.fasta); QUERIES[$s]="$f"; done ;;
  *) echo "unknown mode: $MODE" >&2; exit 1 ;;
esac
(( ${#QUERIES[@]} )) || { echo "No query files found for mode '$MODE'." >&2; exit 1; }

run_one() {
    local s="$1" q="$2" out="$OUTDIR/${s}_blast.tsv"
    if [[ -s "$out" ]]; then echo "[skip] $s (exists)"; return 0; fi
    # -dust no: short, low-complexity-rich ITS queries must not be masked
    # -max_target_seqs 100: LCA needs ALL hits within 2% of the top score
    if blastn -task blastn -db "$DB" -query "$q" -dust no \
        -evalue 1e-5 -max_target_seqs 100 -num_threads "$THREADS" \
        -outfmt "6 qseqid sacc staxids pident length qcovs evalue bitscore stitle" \
        -out "$out.tmp" 2> "$OUTDIR/${s}.err"; then
        mv "$out.tmp" "$out"; echo "[done] $s: $(wc -l < "$out") hit rows"
    else
        echo "[FAILED] $s -> see $OUTDIR/${s}.err" >&2; rm -f "$out.tmp"
    fi
}
export -f run_one; export DB OUTDIR THREADS

printf '%s\n' "${!QUERIES[@]}" | sort | while read -r s; do echo "$s ${QUERIES[$s]}"; done \
  | xargs -P "$JOBS" -L 1 bash -c 'run_one "$0" "$1"'

echo "Done -> $OUTDIR/. Next: lca_assign.py (see header of this script)."
echo "Samples with no hit file or an empty one had no confident match in the custom DB."
