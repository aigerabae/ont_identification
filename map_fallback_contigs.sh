#!/usr/bin/env bash
# map_fallback_contigs.sh
#
# Maps pre-deduplication, post-fastp reads against each sample's no-ITSx
# fallback contigs (no_itsx_fallback/{sample}_top_contigs.fasta), for the
# 8 samples where ITSx found no extractable ITS region.
#
# Consistent with the dedup fix applied to the main ITS-positive samples:
# uses fastp_merge's pre-dedup output (unmerged_r1/r2, merged), not the
# deduplicated reads used for assembly.
#
# Usage: ./map_fallback_contigs.sh
set -euo pipefail

SAMPLES=(18S-13_S39 18S-17_S40 18S-23_S45 18S-25_S46 18S-32_S44 18S-34_S43 18S-37_S42 18S-4_S47)
FALLBACK_DIR="no_itsx_fallback"
PIPELINE_OUT="pipeline_out"
THREADS=4

for sample in "${SAMPLES[@]}"; do
    contigs="${FALLBACK_DIR}/${sample}_top_contigs.fasta"
    idx="${FALLBACK_DIR}/${sample}_fallback_idx"
    r1="${PIPELINE_OUT}/${sample}/qc_reads/unmerged_R1.fastq.gz"
    r2="${PIPELINE_OUT}/${sample}/qc_reads/unmerged_R2.fastq.gz"
    merged="${PIPELINE_OUT}/${sample}/qc_reads/merged.fastq.gz"
    bam="${FALLBACK_DIR}/${sample}_mapped_fallback.sorted.bam"
    flagstat="${FALLBACK_DIR}/${sample}_flagstat_fallback.txt"

    if [ ! -f "$contigs" ]; then
        echo "ERROR [$sample]: $contigs not found, skipping." >&2
        continue
    fi
    if [ ! -f "$r1" ] || [ ! -f "$r2" ] || [ ! -f "$merged" ]; then
        echo "ERROR [$sample]: pre-dedup fastp reads not found under $PIPELINE_OUT/$sample/qc_reads/, skipping." >&2
        continue
    fi

    echo "=== $sample ==="
    bowtie2-build --quiet "$contigs" "$idx"
    bowtie2 -x "$idx" -1 "$r1" -2 "$r2" -U "$merged" \
        --threads "$THREADS" -S "${bam}.sam"
    samtools sort -@ "$THREADS" -o "$bam" "${bam}.sam"
    samtools index "$bam"
    samtools flagstat "$bam" > "$flagstat"
    rm "${bam}.sam"
    echo "Done: $bam"
done

echo "All samples processed. BAMs and flagstats are in $FALLBACK_DIR/"
