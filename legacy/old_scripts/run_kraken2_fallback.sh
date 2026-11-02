#!/bin/bash
# run_kraken2_fallback.sh
#
# Runs Kraken2 + Bracken (genus-level) directly on RAW READS for the 12
# samples that couldn't produce usable ITS-region results (4 where assembly
# itself crashed, 8 where assembly succeeded but ITSx found no ITS content).
# Kraken2's k-mer classification doesn't need clean assembly or amplicon
# boundaries -- it works on raw reads directly, which is why it can reach
# all 12 samples where the assembly-dependent pipeline could not.
#
# Resolution ceiling: SILVA has almost no species-rank plant nodes, so this
# will cap out at genus/family for the plant fraction -- expected, already
# documented in this project, not a bug.
#
# Resilient by design: one sample's failure doesn't stop the others (same
# philosophy as your original run_kraken2_pipeline_v2.sh).
#
# EDIT THE VARIABLES BELOW to match your actual paths before running.

set -uo pipefail  # NOT -e: we want to continue past per-sample failures

KRAKEN_DB="$HOME/biostar/NCB/hp_august/kraken2_db/SILVA_144_k2db"
RAW_FASTQ_DIR="raw_fastqs"
OUTDIR="kraken2_fallback_results"
THREADS=8
BRACKEN_READ_LEN=150   # adjust per-sample below if a sample used different read length

mkdir -p "$OUTDIR/kraken2_out" "$OUTDIR/bracken_out"

# The 12 samples: 4 where assembly crashed + 8 where ITSx found zero ITS content
SAMPLES=(
    18S-1_S37
    18S-3_S41
    18S-9_S38
    18S-31_S33
    18S-13_S39
    18S-17_S40
    18S-23_S45
    18S-25_S46
    18S-32_S44
    18S-34_S43
    18S-37_S42
    18S-4_S47
)

n_ok=0
n_failed=0
failed_samples=()

for sample in "${SAMPLES[@]}"; do
    r1="$RAW_FASTQ_DIR/${sample}_L001_R1_001.fastq.gz"
    r2="$RAW_FASTQ_DIR/${sample}_L001_R2_001.fastq.gz"
    kreport="$OUTDIR/kraken2_out/${sample}.kreport"
    kout="$OUTDIR/kraken2_out/${sample}.kraken"
    bracken_out="$OUTDIR/bracken_out/${sample}.bracken"
    bracken_kreport="$OUTDIR/bracken_out/${sample}_bracken.kreport"

    if [[ -s "$bracken_out" ]]; then
        echo "[skip] $sample already has Bracken output"
        n_ok=$((n_ok + 1))
        continue
    fi

    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "[FAILED] $sample: raw fastq files not found ($r1 / $r2)" >&2
        failed_samples+=("$sample")
        n_failed=$((n_failed + 1))
        continue
    fi

    echo "[running] $sample: Kraken2 ..."
    if ! kraken2 --db "$KRAKEN_DB" \
        --threads "$THREADS" --paired --confidence 0.1 \
        --report "$kreport" --output "$kout" --use-names \
        "$r1" "$r2"
    then
        echo "[FAILED] $sample: kraken2 itself failed" >&2
        failed_samples+=("$sample")
        n_failed=$((n_failed + 1))
        continue
    fi

    echo "[running] $sample: Bracken (genus-level) ..."
    if ! bracken -d "$KRAKEN_DB" -i "$kreport" \
        -o "$bracken_out" -w "$bracken_kreport" \
        -r "$BRACKEN_READ_LEN" -l G -t 1
    then
        echo "[FAILED] $sample: bracken failed (kraken2 output kept at $kreport)" >&2
        failed_samples+=("$sample")
        n_failed=$((n_failed + 1))
        continue
    fi

    echo "[done] $sample"
    n_ok=$((n_ok + 1))
done

echo ""
echo "Kraken2/Bracken fallback complete: $n_ok succeeded, $n_failed failed."
if [[ ${#failed_samples[@]} -gt 0 ]]; then
    echo "Failed samples:"
    for s in "${failed_samples[@]}"; do
        echo "  - $s"
    done
fi
echo "Results in: $OUTDIR/bracken_out/"
echo ""
echo "Optional: combine all samples into one table with Bracken's own script, e.g.:"
echo "  combine_bracken_outputs.py --files $OUTDIR/bracken_out/*.bracken -o $OUTDIR/combined_abundance.tsv"
