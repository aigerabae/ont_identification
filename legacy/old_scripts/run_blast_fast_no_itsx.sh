#!/bin/bash
# run_blast_fast_no_itsx.sh
#
# Two independent, NCBI-sanctioned speedups over the original sequential
# run_blast_*.sh scripts:
#
# 1. BOUNDED CONCURRENCY. NCBI's actual published guideline
#    (https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html)
#    restricts SUBMISSION rate (no more than once every 10 seconds) and
#    PER-JOB POLLING rate (no more than once a minute) -- it does NOT
#    prohibit multiple searches running concurrently. This script runs up
#    to --max-concurrent searches at once, staggering each new submission
#    by --stagger-seconds (default 15s, safely above the 10s minimum), so
#    wall-clock time drops roughly in proportion to --max-concurrent
#    without violating the actual rule.
#
# 2. SEARCH-SPACE RESTRICTION via -entrez_query. NCBI's own support
#    documentation explicitly recommends this to speed up searches
#    ("Use organism limit to target a smaller subset in the database" --
#    https://support.nlm.nih.gov/kbArticle/?pn=KA-03285). Default here
#    restricts to Eukaryota, which excludes the large volume of
#    bacterial/viral/archaeal sequence in core_nt while keeping every
#    taxon this project has actually found so far (plants, fungi, an
#    insect). Override with --entrez-query "" to disable if you want the
#    full, unrestricted database.
#
# Operates on any directory of one-fasta-per-sample query files (matches
# no_itsx_fallback/*_top_contigs.fasta, or adapt --suffix for other dirs).
# Output uses the same staxids-included format as run_blast_all_samples.sh
# -- drops straight into lca_assign.py with no lookup file needed.
#
# Usage:
#   ./run_blast_fast_no_itsx.sh <input_dir> <output_dir> [suffix] [max_concurrent] [entrez_query]
#
# Example (the no-ITSx fallback batch):
#   ./run_blast_fast_no_itsx.sh no_itsx_fallback blast_results_no_itsx _top_contigs.fasta 4 "eukaryota[organism]"

set -uo pipefail  # not -e: individual job failures shouldn't kill the scheduler

INPUT_DIR="${1:?Usage: $0 <input_dir> <output_dir> [suffix] [max_concurrent] [entrez_query]}"
OUTDIR="${2:?Usage: $0 <input_dir> <output_dir> [suffix] [max_concurrent] [entrez_query]}"
SUFFIX="${3:-_top_contigs.fasta}"
MAX_CONCURRENT="${4:-4}"
ENTREZ_QUERY="${5-eukaryota[organism]}"   # pass "" explicitly to disable
DB="core_nt"
STAGGER_SECONDS="${STAGGER_SECONDS_OVERRIDE:-15}"

mkdir -p "$OUTDIR"

shopt -s nullglob
fasta_files=("$INPUT_DIR"/*"$SUFFIX")

if [[ ${#fasta_files[@]} -eq 0 ]]; then
    echo "No *${SUFFIX} files found in $INPUT_DIR" >&2
    exit 1
fi

echo "Found ${#fasta_files[@]} sample(s). max_concurrent=$MAX_CONCURRENT, entrez_query='${ENTREZ_QUERY:-<none>}'"

run_one() {
    local fasta="$1"
    local sample
    sample=$(basename "$fasta" "$SUFFIX")
    local out_file="$OUTDIR/${sample}_blast.tsv"

    local entrez_args=()
    if [[ -n "$ENTREZ_QUERY" ]]; then
        entrez_args=(-entrez_query "$ENTREZ_QUERY")
    fi

    if blastn -remote -db "$DB" \
        -task blastn \
        -query "$fasta" \
        "${entrez_args[@]}" \
        -outfmt "6 qseqid sacc staxids pident length qcovs evalue bitscore stitle" \
        -max_target_seqs 10 \
        -out "$out_file"
    then
        echo "[done]   $sample -> $out_file"
    else
        echo "[FAILED] $sample" >&2
        rm -f "$out_file"
    fi
}

n_running=0
for fasta in "${fasta_files[@]}"; do
    sample=$(basename "$fasta" "$SUFFIX")
    out_file="$OUTDIR/${sample}_blast.tsv"

    if [[ -s "$out_file" ]]; then
        echo "[skip] $sample already has results"
        continue
    fi

    echo "[submit] $sample"
    run_one "$fasta" &
    n_running=$((n_running + 1))

    if [[ "$n_running" -ge "$MAX_CONCURRENT" ]]; then
        wait -n   # free up one slot before continuing
        n_running=$((n_running - 1))
    fi

    sleep "$STAGGER_SECONDS"   # keeps NEW SUBMISSIONS >=10s apart regardless of concurrency
done

wait   # let all remaining background jobs finish
echo ""
echo "All submitted. Results in: $OUTDIR/"
