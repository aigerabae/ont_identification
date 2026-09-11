"""
Mixed-species plant ITS pipeline — Snakemake version.

Pipeline stages per sample:
  1. fastp: adapter/quality trim + merge overlapping pairs
  2. clumpify: dedupe merged and unmerged read sets separately
  3. Split assembly:
       - unmerged (long-insert) pairs -> metaspades.py (paired mode)
       - merged   (short-insert) reads -> spades.py --only-assembler (single-end mode)
  4. Tag + combine both contig sets, cd-hit-est to remove exact duplicates
  5. vsearch: annotate contig headers with ;size= from SPAdes coverage,
     sort by size, uchime_denovo chimera detection -> nonchimeric contigs
  6. ITSx on the nonchimeric contigs -> region-level ITS1/ITS2/full FASTA + positions.txt
  7. extract_its_regions.py: pick the best available sequence per ITS-flagged
     contig (full > ITS1+ITS2 > ITS1 only > ITS2 only) -> BLAST-ready FASTA files
  8. bowtie2 read-mapping back onto (a) the full contig set for assembly QC,
     and (b) the ITS-flagged contigs only, for a cleaner per-contig depth signal
     to use as a rough relative-abundance proxy.

BLAST itself is NOT automated here (NCBI's web BLAST is what you've been using,
and remote command-line BLAST is slow/rate-limited) — the pipeline's job is to
hand you clean, correctly-extracted, ready-to-upload FASTA files per sample in
pipeline_out/<sample>/its_extracted/, plus the depth numbers to interpret
whatever species calls you get back.

Usage:
    snakemake -j 8 --use-conda   (if you wrap each env; otherwise ensure all
                                   tools below are on PATH in one env)
    snakemake -j 8 --until extract_its   ONplants-1_S40   # single target example

See config.yaml to set your sample list and paths.
"""

configfile: "config.yaml"

import os

SAMPLES = config["samples"]
RAW = config["raw_fastq_dir"]
OUT = config["outdir"]
THREADS = config["threads"]
SIZE_SCALE = config["chimera_size_scale"]


rule all:
    input:
        expand(f"{OUT}/{{sample}}/its_extracted/extraction_log.tsv", sample=SAMPLES),
        expand(f"{OUT}/{{sample}}/qc/flagstat_full.txt", sample=SAMPLES),
        expand(f"{OUT}/{{sample}}/qc/depth_its_mean.tsv", sample=SAMPLES),
        f"{OUT}/skipped_samples_summary.txt",


# ── 1. fastp: trim + merge ──────────────────────────────────────────────
rule fastp_merge:
    input:
        r1=f"{RAW}/{{sample}}_L001_R1_001.fastq.gz",
        r2=f"{RAW}/{{sample}}_L001_R2_001.fastq.gz",
    output:
        merged=f"{OUT}/{{sample}}/qc_reads/merged.fastq.gz",
        unmerged_r1=f"{OUT}/{{sample}}/qc_reads/unmerged_R1.fastq.gz",
        unmerged_r2=f"{OUT}/{{sample}}/qc_reads/unmerged_R2.fastq.gz",
        json=f"{OUT}/{{sample}}/qc_reads/fastp.json",
        html=f"{OUT}/{{sample}}/qc_reads/fastp.html",
    params:
        qqp=config["fastp"]["qualified_quality_phred"],
        minlen=config["fastp"]["length_required"],
    threads: THREADS
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
          --merge --merged_out {output.merged} \
          --out1 {output.unmerged_r1} --out2 {output.unmerged_r2} \
          --qualified_quality_phred {params.qqp} \
          --length_required {params.minlen} \
          --thread {threads} \
          --json {output.json} --html {output.html}
        """


# ── 2. clumpify: dedupe merged and unmerged separately ──────────────────
rule dedup_unmerged:
    input:
        r1=rules.fastp_merge.output.unmerged_r1,
        r2=rules.fastp_merge.output.unmerged_r2,
    output:
        r1=f"{OUT}/{{sample}}/qc_reads/unmerged_R1_dedup.fastq.gz",
        r2=f"{OUT}/{{sample}}/qc_reads/unmerged_R2_dedup.fastq.gz",
    shell:
        """
        clumpify.sh in1={input.r1} in2={input.r2} \
          out1={output.r1} out2={output.r2} dedupe
        """


rule dedup_merged:
    input:
        rules.fastp_merge.output.merged,
    output:
        f"{OUT}/{{sample}}/qc_reads/merged_dedup.fastq.gz",
    shell:
        """
        clumpify.sh in={input} out={output} dedupe
        """


# ── 3. Split assembly ────────────────────────────────────────────────────
rule assemble_longinsert:
    input:
        r1=rules.dedup_unmerged.output.r1,
        r2=rules.dedup_unmerged.output.r2,
    output:
        contigs=f"{OUT}/{{sample}}/assembly/longinsert/contigs.fasta",
    params:
        outdir=f"{OUT}/{{sample}}/assembly/longinsert",
    threads: THREADS
    shell:
        """
        metaspades.py -1 {input.r1} -2 {input.r2} \
          -o {params.outdir} --threads {threads}
        """


rule assemble_shortinsert:
    input:
        rules.dedup_merged.output,
    output:
        contigs=f"{OUT}/{{sample}}/assembly/shortinsert/contigs.fasta",
    params:
        outdir=f"{OUT}/{{sample}}/assembly/shortinsert",
    threads: THREADS
    shell:
        """
        spades.py -s {input} -o {params.outdir} \
          --only-assembler --threads {threads}
        """


# ── 4. Tag, combine, cd-hit dedup ────────────────────────────────────────
rule combine_contigs:
    input:
        longinsert=rules.assemble_longinsert.output.contigs,
        shortinsert=rules.assemble_shortinsert.output.contigs,
    output:
        combined=f"{OUT}/{{sample}}/assembly/combined_contigs_nr_tagged.fasta",
        clstr=f"{OUT}/{{sample}}/assembly/combined_contigs_nr_tagged.fasta.clstr",
    params:
        tagged=f"{OUT}/{{sample}}/assembly/combined_contigs_tagged.fasta",
    threads: THREADS
    shell:
        """
        sed 's/^>/>longinsert_/' {input.longinsert} > {params.tagged}.long
        sed 's/^>/>shortinsert_/' {input.shortinsert} > {params.tagged}.short
        cat {params.tagged}.long {params.tagged}.short > {params.tagged}
        cd-hit-est -i {params.tagged} -o {output.combined} -c 1.0 -T {threads}
        rm {params.tagged}.long {params.tagged}.short
        """


# ── 5. vsearch: size-annotate, sort, chimera detection ───────────────────
rule size_annotate:
    input:
        rules.combine_contigs.output.combined,
    output:
        f"{OUT}/{{sample}}/chimera/combined_contigs_sized.fasta",
    params:
        scale=SIZE_SCALE,
    shell:
        r"""
        perl -pe '
          if (/^>/ and /cov_([\d.]+)/) {{
            my $size = int($1 * {params.scale});
            $size = 1 if $size < 1;
            chomp;
            $_ .= ";size=$size\n";
          }}
        ' {input} > {output}
        """


rule sort_by_size:
    input:
        rules.size_annotate.output,
    output:
        f"{OUT}/{{sample}}/chimera/combined_contigs_sorted.fasta",
    shell:
        "vsearch --sortbysize {input} --output {output}"


rule chimera_detect:
    input:
        rules.sort_by_size.output,
    output:
        nonchimeric=f"{OUT}/{{sample}}/chimera/combined_contigs_nonchimeric.fasta",
        chimeric=f"{OUT}/{{sample}}/chimera/combined_contigs_chimeric.fasta",
        report=f"{OUT}/{{sample}}/chimera/uchime_report.txt",
    shell:
        """
        vsearch --uchime_denovo {input} \
          --nonchimeras {output.nonchimeric} \
          --chimeras {output.chimeric} \
          --uchimeout {output.report} --uchimeout5
        """


# ── 6. ITSx on the nonchimeric contigs ───────────────────────────────────
rule itsx:
    input:
        rules.chimera_detect.output.nonchimeric,
    output:
        summary=f"{OUT}/{{sample}}/itsx/its.summary.txt",
        positions=f"{OUT}/{{sample}}/itsx/its.positions.txt",
    params:
        prefix=f"{OUT}/{{sample}}/itsx/its",
    threads: THREADS
    shell:
        """
        ITSx -i {input} -o {params.prefix} --cpu {threads} --preserve T
        """


# ── 7. Extract the best-available ITS sequence per flagged contig ───────
rule extract_its:
    input:
        rules.itsx.output.positions,
    output:
        log=f"{OUT}/{{sample}}/its_extracted/extraction_log.tsv",
    params:
        prefix=f"{OUT}/{{sample}}/itsx/its",
        outdir=f"{OUT}/{{sample}}/its_extracted",
    shell:
        """
        python3 scripts/extract_its_regions.py \
          --prefix {params.prefix} --outdir {params.outdir}
        """


# ── 8a. Assembly QC: map reads back onto the FULL contig set ────────────
rule map_full_qc:
    input:
        contigs=rules.combine_contigs.output.combined,
        r1=rules.dedup_unmerged.output.r1,
        r2=rules.dedup_unmerged.output.r2,
        merged=rules.dedup_merged.output,
    output:
        bam=f"{OUT}/{{sample}}/qc/mapped_full.sorted.bam",
        flagstat=f"{OUT}/{{sample}}/qc/flagstat_full.txt",
    params:
        idx=f"{OUT}/{{sample}}/qc/full_idx",
    threads: THREADS
    shell:
        """
        bowtie2-build {input.contigs} {params.idx}
        bowtie2 -x {params.idx} -1 {input.r1} -2 {input.r2} -U {input.merged} \
          --threads {threads} -S {output.bam}.sam
        samtools sort -@ {threads} -o {output.bam} {output.bam}.sam
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        rm {output.bam}.sam
        """


# ── 8b. Targeted depth: map onto ITS-flagged contigs only ───────────────
# Uses positions.txt to build the same subset FASTA the extract_its rule
# already identified as containing real ITS signal, so depth numbers are
# calculated only among genuinely comparable, non-junk contigs.
rule its_subset_fasta:
    input:
        contigs=rules.chimera_detect.output.nonchimeric,
        positions=rules.itsx.output.positions,
    output:
        f"{OUT}/{{sample}}/qc/its_contigs_only.fasta",
    shell:
        r"""
        awk 'NR==FNR {{want[$1]=1; next}}
             /^>/ {{hdr=substr($0,2); keep=(hdr in want)}}
             keep' <(cut -f1 {input.positions}) {input.contigs} > {output}
        """


rule map_its_only:
    input:
        contigs=rules.its_subset_fasta.output,
        r1=rules.dedup_unmerged.output.r1,
        r2=rules.dedup_unmerged.output.r2,
        merged=rules.dedup_merged.output,
    output:
        bam=f"{OUT}/{{sample}}/qc/mapped_its.sorted.bam",
        flagstat=f"{OUT}/{{sample}}/qc/flagstat_its.txt",
        no_its_flag=f"{OUT}/{{sample}}/qc/NO_ITS_CONTENT.flag",
    params:
        idx=f"{OUT}/{{sample}}/qc/its_idx",
        sample=lambda wc: wc.sample,
    threads: THREADS
    shell:
        """
        rm -f {output.no_its_flag}
        if ! grep -q ">" {input.contigs}; then
            echo "WARNING [{params.sample}]: no ITS-flagged contigs -- skipping bowtie2 mapping." >&2
            echo "{params.sample}: no ITS-flagged contigs (its_contigs_only.fasta was empty). Mapping and depth were skipped." > {output.no_its_flag}
            touch {output.bam}
            echo "SKIPPED: no ITS-flagged contigs for this sample; mapping not run." > {output.flagstat}
        else
            bowtie2-build {input.contigs} {params.idx}
            bowtie2 -x {params.idx} -1 {input.r1} -2 {input.r2} -U {input.merged} \
              --threads {threads} -S {output.bam}.sam
            samtools sort -@ {threads} -o {output.bam} {output.bam}.sam
            samtools index {output.bam}
            samtools flagstat {output.bam} > {output.flagstat}
            rm {output.bam}.sam
        fi
        """


rule depth_its_mean:
    input:
        bam=rules.map_its_only.output.bam,
        no_its_flag=rules.map_its_only.output.no_its_flag,
    output:
        raw=f"{OUT}/{{sample}}/qc/depth_its_raw.txt",
        mean=f"{OUT}/{{sample}}/qc/depth_its_mean.tsv",
    shell:
        """
        if [ -s {input.no_its_flag} ]; then
            echo "WARNING: no ITS-flagged contigs for this sample -- writing empty depth files." >&2
            touch {output.raw} {output.mean}
        else
            samtools view -b -q 10 {input.bam} > {input.bam}.q10.bam
            samtools index {input.bam}.q10.bam
            samtools depth -a {input.bam}.q10.bam > {output.raw}
            awk '{{sum[$1]+=$3; count[$1]++}} END {{for (c in sum) print c"\t"sum[c]/count[c]}}' \
              {output.raw} | sort -k2 -rn > {output.mean}
            rm {input.bam}.q10.bam {input.bam}.q10.bam.bai
        fi
        """


# ── 9. Aggregate: flag every sample that had no ITS-flagged contigs ─────
# Depends on every sample's map_its_only output, so this always reflects
# the full, current run -- one place to check instead of scrolling through
# Snakemake's log for WARNING lines.
rule flag_low_yield_samples:
    input:
        flags=expand(f"{OUT}/{{sample}}/qc/NO_ITS_CONTENT.flag", sample=SAMPLES),
    output:
        summary=f"{OUT}/skipped_samples_summary.txt",
    run:
        skipped = []
        for flag_path, sample in zip(input.flags, SAMPLES):
            if os.path.getsize(flag_path) > 0:
                skipped.append(sample)
        with open(output.summary, "w") as fh:
            if skipped:
                fh.write(f"{len(skipped)} of {len(SAMPLES)} sample(s) had NO ITS-flagged contigs "
                          f"(mapping/depth steps were skipped for these):\n")
                for s in skipped:
                    fh.write(f"  - {s}\n")
            else:
                fh.write(f"All {len(SAMPLES)} sample(s) had at least one ITS-flagged contig.\n")
