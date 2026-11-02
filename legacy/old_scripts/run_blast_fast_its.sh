#!/bin/bash
# run_blast_fast_its.sh
#
# Fast version of run_blast_all_samples.sh, for the samples where ITSx
# successfully found ITS-flagged contigs. Combines two things:
#
# 1. Per-sample query preparation: concatenates each sample's multiple
#    its_extracted/*_its.fasta files into one combined query (same as
#    run_blast_all_samples.sh) -- unlike the no-ITSx fallback samples,
#    which already have one pre-combined file per sample.
# 2. The same NCBI-sanctioned speedups as run_blast_fast_no_itsx.sh:
#    bounded concurrency (submissions staggered >=10s apart, several
#    searches running at once -- compliant with NCBI's actual published
#    guideline, see https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html)
#    and an optional -entrez_query search-space restriction (default:
#    Eukaryota, excludes bacterial/archaeal/viral sequence from core_nt).
#
# Automatically skips any sample with no *_its.fasta files (i.e. no ITS
# content) -- same behavior as run_blast_all_samples.sh / run_blast_for_megan.sh.
#
# Output format matches run_blast_all_samples.sh (staxids included) --
# drops straight into lca_assign.py with no lookup file needed.
#
# Usage:
#   ./run_blast_fast_its.sh [pipeline_out_dir] [output_dir] [max_concurrent] [entrez_query]
#
# Example:
#   ./run_blast_fast_its.sh pipeline_out blast_results 4 "eukaryota[organism]"

set -uo pipefail  # not -e: individual job failures shouldn't kill the scheduler

PIPELINE_OUT="${1:-pipeline_out}"
OUTDIR="${2:-blast_results}"
MAX_CONCURRENT="${3:-4}"
ENTREZ_QUERY="${4-eukaryota[organism]}"   # pass "" explicitly to disable
DB="core_nt"
STAGGER_SECONDS="${STAGGER_SECONDS_OVERRIDE:-15}"

mkdir -p "$OUTDIR"

shopt -s nullglob
sample_dirs=("$PIPELINE_OUT"/*/its_extracted)

if [[ ${#sample_dirs[@]} -eq 0 ]]; then
    echo "No */its_extracted directories found under $PIPELINE_OUT — check the path." >&2
    exit 1
fi

echo "Found ${#sample_dirs[@]} sample folder(s) to check. max_concurrent=$MAX_CONCURRENT, entrez_query='${ENTREZ_QUERY:-<none>}'"

run_one() {
    local combined_query="$1"
    local out_file="$2"
    local sample="$3"

    local entrez_args=()
    if [[ -n "$ENTREZ_QUERY" ]]; then
        entrez_args=(-entrez_query "$ENTREZ_QUERY")
    fi

    if blastn -remote -db "$DB" \
        -task blastn \
        -query "$combined_query" \
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
for sample_dir in "${sample_dirs[@]}"; do
    sample=$(basename "$(dirname "$sample_dir")")
    combined_query="$OUTDIR/${sample}_query.fasta"
    out_file="$OUTDIR/${sample}_blast.tsv"

    if [[ -s "$out_file" ]]; then
        n_cols=$(head -1 "$out_file" | awk -F'\t' '{print NF}')
        if [[ "$n_cols" -eq 9 ]]; then
            echo "[skip] $sample already has results (new format, staxids present)"
        elif [[ "$n_cols" -eq 10 ]]; then
            echo "[skip] $sample already has results (OLD format, no staxids -- pass --taxid-lookup to lca_assign.py for this sample)"
        else
            echo "[skip] $sample already has results (unrecognized format: $n_cols columns -- verify manually)"
        fi
        continue
    fi

    its_files=("$sample_dir"/*_its.fasta)
    if [[ ${#its_files[@]} -eq 0 ]]; then
        echo "[skip] $sample: no *_its.fasta files (no ITS content for this sample)"
        continue
    fi

    cat "${its_files[@]}" > "$combined_query"
    n_seqs=$(grep -c ">" "$combined_query" || true)
    echo "[submit] $sample: $n_seqs sequence(s)"

    run_one "$combined_query" "$out_file" "$sample" &
    n_running=$((n_running + 1))

    if [[ "$n_running" -ge "$MAX_CONCURRENT" ]]; then
        wait -n
        n_running=$((n_running - 1))
    fi

    sleep "$STAGGER_SECONDS"
done

wait
echo ""
echo "All submitted. Results in: $OUTDIR/"
