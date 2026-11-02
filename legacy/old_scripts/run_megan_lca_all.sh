#!/bin/bash
# run_megan_lca_all.sh
#
# Runs blast2lca (MEGAN Community Edition, installed via bioconda) on every
# sample's classical-format BLAST output from run_blast_for_megan.sh.
# Only samples that actually have a *_blast_megan.tsv file get processed --
# this naturally excludes any sample where BLAST itself was skipped because
# it had no ITS-extracted contigs to begin with.
#
# Parameters below (-ms, -me, -top, -supp, -sup) are blast2lca's own
# documented defaults, confirmed from MEGAN's source. -mpi (min percent
# identity) is left at 0.0 (no filter) rather than borrowing a threshold
# from an unrelated paper -- adjust if you want a stricter cutoff.
#
# IMPORTANT: check `blast2lca --help` yourself for the exact flag name to
# select "naive" vs "weighted" LCA algorithm before assuming a default --
# this was not independently confirmed and is deliberately left unset here
# rather than guessed.
#
# Usage:
#   ./run_megan_lca_all.sh [megan_blast_dir] [mapping_db_path] [output_dir]
#
# Example:
#   ./run_megan_lca_all.sh blast_results_megan /path/to/megan-map-XXXX.db megan_lca_results

set -euo pipefail

BLAST_DIR="${1:-blast_results_megan}"
MAPPING_DB="${2:?Usage: $0 <blast_dir> <mapping_db_path> [output_dir]}"
OUTDIR="${3:-megan_lca_results}"

mkdir -p "$OUTDIR"

shopt -s nullglob
blast_files=("$BLAST_DIR"/*_blast_megan.tsv)

if [[ ${#blast_files[@]} -eq 0 ]]; then
    echo "No *_blast_megan.tsv files found in $BLAST_DIR — run run_blast_for_megan.sh first." >&2
    exit 1
fi

if [[ ! -f "$MAPPING_DB" ]]; then
    echo "Mapping DB not found at: $MAPPING_DB" >&2
    echo "Download it from https://software-ab.cs.uni-tuebingen.de/download/megan6/welcome.html" >&2
    exit 1
fi

n_done=0
n_failed=0

for blast_file in "${blast_files[@]}"; do
    sample=$(basename "$blast_file" "_blast_megan.tsv")
    out_file="$OUTDIR/${sample}_megan_lca.txt"

    if [[ -s "$out_file" ]]; then
        echo "[skip] $sample already has LCA output: $out_file"
        continue
    fi

    n_seqs=$(cut -f1 "$blast_file" | sort -u | wc -l)
    echo "[running] $sample: $n_seqs query sequence(s) -> blast2lca"

    if blast2lca \
        -i "$blast_file" \
        -f BlastTAB \
        -m BlastN \
        -mdb "$MAPPING_DB" \
        -o "$out_file" \
        -ms 50.0 \
        -me 0.01 \
        -mpi 0.0 \
        -top 10.0 \
        -supp 0.0 \
        -sup 0
    then
        echo "[done] $sample -> $out_file"
        n_done=$((n_done + 1))
    else
        echo "[FAILED] $sample — see blast2lca's error output above." >&2
        n_failed=$((n_failed + 1))
        rm -f "$out_file"
    fi
done

echo ""
echo "Done. $n_done sample(s) succeeded, $n_failed failed."
echo "Results in: $OUTDIR/"
