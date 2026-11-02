#!/usr/bin/env bash
# map_reads.sh <bt2_index> <out_sorted_bam> <flagstat_out> <R1> <R2> <merged> <threads> [extra bowtie2 args...]
#
# Maps ALL post-fastp reads (merged + both mates of unmerged pairs, NOT the
# deduplicated copies used for assembly) with bowtie2. Empty read sets are
# skipped instead of crashing bowtie2. If there are no reads at all, a
# zero-byte BAM is written (downstream scripts treat that as "no data").
set -euo pipefail
idx="$1"; bam="$2"; fs="$3"; r1="$4"; r2="$5"; mg="$6"; thr="$7"; shift 7
extra=("$@")

nonempty() { [ -n "$(zcat "$1" | head -c 1 || true)" ]; }

args=()
if nonempty "$r1" && nonempty "$r2"; then args+=(-1 "$r1" -2 "$r2"); fi
if nonempty "$mg"; then args+=(-U "$mg"); fi

if [ ${#args[@]} -eq 0 ]; then
    : > "$bam"
    echo "NO_READS" > "$fs"
    exit 0
fi

bowtie2 -x "$idx" "${args[@]}" --threads "$thr" "${extra[@]}" 2> "${bam}.bt2.log" \
  | samtools sort -@ "$thr" -o "$bam" -
samtools index "$bam"
samtools flagstat "$bam" > "$fs"
