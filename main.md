Downloading database (only SSU):
```
cd /mnt/harddisk/biostar/NCB/hp_august/kraken2_db
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Kraken2/2.1.6/SSU/SILVA_144_k2db.tgz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Kraken2/2.1.6/SSU/SILVA_144_k2db.tgz.md5

# verify integrity
md5sum -c SILVA_144_k2db.tgz.md5

# extract
tar -xzvf SILVA_144_k2db.tgz
```

Installing kraken and bracken and other software:
```
conda create --name kraken2 and bracken and itsx
conda activate kraken2
conda install bioconda::kraken2
conda install -c bioconda -c conda-forge bracken
conda install bioconda::itsx
conda install bioconda::seqtk
conda install bioconda::spades
conda install bioconda::fastp
conda install bioconda::bbmap
conda install bioconda::cd-hit
conda install bioconda::vsearch
conda install bioconda::seqkit
```

Running kraken (on 1 sample, test mode):
```
kraken2 --db ~/biostar/NCB/hp_august/kraken2_db/SILVA_144_k2db \
  --threads 8 --paired --confidence 0.1 \
  --report /tmp/test.kreport --output /tmp/test.kraken --use-names \
  ~/biostar/NCB/hp_august/raw_fastqs/ONplants-1_S40_L001_R1_001.fastq.gz \
  ~/biostar/NCB/hp_august/raw_fastqs/ONplants-1_S40_L001_R2_001.fastq.gz

cat /tmp/test.kreport | head -30
```

Using all samples (genuslevel):
```
./run_kraken2_pipeline_v2.sh
```

Only did some samples. Might want to enrich the specific ITS region first

Verifying if ITS is present:
```
# Step 1 — ITSx directly on raw reads (needs fasta, not fastq, so convert first)
seqtk seq -a raw_fastqs/ONplants-1_S40_L001_R1_001.fastq.gz > ONplants-1_S40_R1.fasta
ITSx -i ONplants-1_S40_R1.fasta -o itsx_test --cpu 8 --preserve T

# check the summary
cat itsx_test.summary.txt
```
This took 1.5 hours and gave me 0 ITS sequences. I might want to run assembly with metaspades first

Pilot Assembly:
```
mkdir -p ~/biostar/NCB/hp_august/pilot_assembly
metaspades.py \
  -1 ~/biostar/NCB/hp_august/raw_fastqs/ONplants-1_S40_L001_R1_001.fastq.gz \
  -2 ~/biostar/NCB/hp_august/raw_fastqs/ONplants-1_S40_L001_R2_001.fastq.gz \
  -o ~/biostar/NCB/hp_august/pilot_assembly/ONplants-1_S40 \
  --threads 30
```

This only made 88 contigs (very little) but from them i was able to get 4 ITS containign regions. Now I will try to improve my assembly quality:
```
mkdir -p qc_reads_v2
fastp \
  -i raw_fastqs/ONplants-1_S40_L001_R1_001.fastq.gz \
  -I raw_fastqs/ONplants-1_S40_L001_R2_001.fastq.gz \
  -o qc_reads_v2/ONplants-1_S40_trimmed_R1.fastq.gz \
  -O qc_reads_v2/ONplants-1_S40_trimmed_R2.fastq.gz \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --thread 8 \
  --json qc_reads_v2/ONplants-1_S40_fastp.json \
  --html qc_reads_v2/ONplants-1_S40_fastp.html
clumpify.sh in1=qc_reads_v2/ONplants-1_S40_trimmed_R1.fastq.gz \
  in2=qc_reads_v2/ONplants-1_S40_trimmed_R2.fastq.gz \
  out1=qc_reads_v2/ONplants-1_S40_R1_dedup.fastq.gz \
  out2=qc_reads_v2/ONplants-1_S40_R2_dedup.fastq.gz \
  dedupe
```

```
# 1. Long-insert (non-overlapping) pairs -> metaSPAdes, paired mode (this is basically your v2, confirmed working)
metaspades.py \
  -1 qc_reads/ONplants-1_S40_unmerged_R1_dedup.fastq.gz \
  -2 qc_reads/ONplants-1_S40_unmerged_R2_dedup.fastq.gz \
  -o pilot_assembly_v4/longinsert_ONplants-1_S40 \
  --threads 30

# 2. Short-insert (merged) reads -> plain SPAdes, single-end mode
#    (metaSPAdes ignores single reads; regular SPAdes does not)
spades.py \
  -s qc_reads/ONplants-1_S40_merged_dedup.fastq.gz \
  -o pilot_assembly_v4/shortinsert_ONplants-1_S40 \
  --only-assembler --threads 30

# 3. Merge both contig sets, remove duplicates
sed 's/^>/>longinsert_/' pilot_assembly_v4/longinsert_ONplants-1_S40/contigs.fasta \
  > pilot_assembly_v4/longinsert_tagged.fasta

sed 's/^>/>shortinsert_/' pilot_assembly_v4/shortinsert_ONplants-1_S40/contigs.fasta \
  > pilot_assembly_v4/shortinsert_tagged.fasta

cat pilot_assembly_v4/longinsert_tagged.fasta pilot_assembly_v4/shortinsert_tagged.fasta \
  > pilot_assembly_v4/combined_contigs_tagged.fasta

cd-hit-est -i pilot_assembly_v4/combined_contigs_tagged.fasta \
  -o pilot_assembly_v4/combined_contigs_nr_tagged.fasta \
  -c 1.0 -T 8

# now this will actually work:
grep "longinsert" pilot_assembly_v4/combined_contigs_nr_tagged.fasta.clstr | wc -l
```

ITSx:
```
ITSx -i pilot_assembly_v4/combined_contigs_nr_tagged.fasta \
  -o pilot_assembly_v4/itsx_v4_test \
  --cpu 30 --preserve T

cat pilot_assembly_v4/itsx_v4_test.summary.txt
```

Still only 86% similarity in BLAST and 5 sequences (2 full length - 1 from short read and 1 from long read)

Bottom line: treat this as "Hippophae rhamnoides and Rhaponticum carthamoides are confidently present" — that claim is well-supported. Do not treat it as "this sample contains only these 2 plant species" — that claim isn't supported yet, and would need either deeper sequencing, a more sensitive extraction strategy, or acceptance that some real components will stay below this pipeline's detection floor.

New plan
# 1 -check quality of assembly:
```
cd pilot_assembly_v4

perl -pe '
  if (/^>/ and /cov_([\d.]+)/) {
    my $size = int($1 * 10);
    $size = 1 if $size < 1;
    chomp;
    $_ .= ";size=$size\n";
  }
' combined_contigs_nr_tagged.fasta > combined_contigs_sized.fasta
```

```
vsearch --sortbysize combined_contigs_sized.fasta \
  --output combined_contigs_sorted.fasta
```

Reading file combined_contigs_sized.fasta 100%  
59490 nt in 250 seqs, min 57, max 2414, avg 238
Getting sizes 100% 
Sorting 100%
Median abundance: 1354
Writing output 100% 

Chimera detection:
```
vsearch --uchime_denovo combined_contigs_sorted.fasta \
  --nonchimeras combined_contigs_nonchimeric.fasta \
  --chimeras combined_contigs_chimeric.fasta \
  --uchimeout uchime_report.txt \
  --uchimeout5
```

Reading file combined_contigs_sorted.fasta 100%  
59490 nt in 250 seqs, min 57, max 2414, avg 238
Masking 100% 
Sorting by abundance 100%
Counting k-mers 100% 
Detecting chimeras 100%  
Found 4 (1.6%) chimeras, 246 (98.4%) non-chimeras,
and 0 (0.0%) borderline sequences in 250 unique sequences.
Taking abundance information into account, this corresponds to
2853 (0.5%) chimeras, 591106 (99.5%) non-chimeras,
and 0 (0.0%) borderline sequences in 593959 total sequences.

```
grep -c ">" combined_contigs_chimeric.fasta   # should show 4
grep ">" combined_contigs_chimeric.fasta       # see which ones
```

4
>shortinsert_NODE_76_length_190_cov_103.634921;size=1036
>shortinsert_NODE_21_length_265_cov_75.927536;size=759
>shortinsert_NODE_16_length_316_cov_60.534392;size=605
>shortinsert_NODE_15_length_330_cov_45.325123;size=453

Checking if my chimeras made it to the after itsx results:
```
cd pilot_assembly_v4
ITSx -i combined_contigs_nonchimeric.fasta \
  -o itsx_v4_nonchimeric --cpu 30 --preserve T
cat itsx_v4_nonchimeric.summary.txt
diff <(sort itsx_v4_test.positions.txt) <(sort itsx_v4_nonchimeric.positions.txt)
```

Checking assembly stats:
```
seqkit stats -a combined_contigs_nr_tagged.fasta
```
file                              format  type  num_seqs  sum_len  min_len  avg_len  max_len   Q1   Q2   Q3  sum_gap  N50  N50_num  Q20(%)  Q30(%)  AvgQual  GC(%)  sum_n
combined_contigs_nr_tagged.fasta  FASTA   DNA        250   59,490       57      238    2,414  143  171  244        0  240       52       0       0        0  50.29      0

Separately for long insert and short insert:
```
seqkit stats -a pilot_assembly_v4/longinsert_ONplants-1_S40/contigs.fasta
seqkit stats -a pilot_assembly_v4/shortinsert_ONplants-1_S40/contigs.fasta
```

Installing quast:
```
conda create --name quast python=3.10
conda activate quast
conda install bioconda::quast
conda install bioconda::bowtie2 bioconda::samtools
```

More detailed stats with QUAST:
```
quast.py pilot_assembly_v4/combined_contigs_nr_tagged.fasta \
  pilot_assembly_v4/longinsert_ONplants-1_S40/contigs.fasta \
  pilot_assembly_v4/shortinsert_ONplants-1_S40/contigs.fasta \
  -o quast_out_v4 --threads 8
```

Read mapping back QC:
```
bowtie2-build combined_contigs_nr_tagged.fasta contigs_idx

bowtie2 -x contigs_idx \
  -1 ../qc_reads/ONplants-1_S40_unmerged_R1_dedup.fastq.gz \
  -2 ../qc_reads/ONplants-1_S40_unmerged_R2_dedup.fastq.gz \
  -U ../qc_reads/ONplants-1_S40_merged_dedup.fastq.gz \
  --threads 30 -S mapped.sam

samtools sort -@8 -o mapped.sorted.bam mapped.sam
samtools index mapped.sorted.bam
samtools flagstat mapped.sorted.bam        # % reads mapping back — sanity check
samtools depth -a mapped.sorted.bam > depth_per_contig.txt
```

Using only those 5 fragments with ITSx for reference:
```
seqkit grep -n -f <(echo -e "longinsert_NODE_1_length_1872_cov_7.430380\nlonginsert_NODE_2_length_1755_cov_2303.875882\nlonginsert_NODE_11_length_389_cov_1.673653\nshortinsert_NODE_1_length_2414_cov_1.659816\nshortinsert_NODE_8_length_441_cov_0.914013") \
  combined_contigs_nr_tagged.fasta > its_contigs_only.fasta

bowtie2-build its_contigs_only.fasta its_contigs_idx

bowtie2 -x its_contigs_idx \
  -1 ../qc_reads/ONplants-1_S40_unmerged_R1_dedup.fastq.gz \
  -2 ../qc_reads/ONplants-1_S40_unmerged_R2_dedup.fastq.gz \
  -U ../qc_reads/ONplants-1_S40_merged_dedup.fastq.gz \
  --threads 30 -S mapped_its_only.sam

samtools sort -@8 -o mapped_its_only.sorted.bam mapped_its_only.sam
samtools index mapped_its_only.sorted.bam
samtools flagstat mapped_its_only.sorted.bam

awk '/^>shortinsert_NODE_1_length_2414/{flag=1; print; next} /^>/{flag=0} flag'   combined_contigs_nr_tagged.fasta > node1_full.fasta
awk '/^>longinsert_NODE_2_length_1755/{flag=1; print; next} /^>/{flag=0} flag'   combined_contigs_nr_tagged.fasta > node2_full.fasta
```

I then made a custom script to extract specific fasta files that i want to blast. For my pilot sample it resulted in 5 sequences 

Hippophae rhamnoides was detected via two distinct ITS sequence variants: (1) a dominant, high-identity variant (99.27% to reference, longinsert_NODE_2, depth 5,485x) and (2) a divergent variant with a real reference-coverage gap (86.75% best identity, depth ~110x combined), independently reconstructed by both assembly branches from different read populations, confirming it's a genuine sequence rather than an assembly artifact.

I used this pipeline to scale it to all 18 samples:
```
snakemake -j 30 --scheduler greedy
# snakemake -j 30 --forcerun its_subset_fasta --scheduler greedy              # when my run crashed i re-ran the specific bit that wasn't working properly; i fixed the code first tho
chmod +x scripts/run_blast_all_samples.sh
./scripts/run_blast_all_samples.sh pipeline_out blast_results
python3 scripts/summarize_blast_hits.py --blast-dir blast_results  --out blast_results/summary.tsv --top-n 5
python3 scripts/build_final_report.py  --blast-summary blast_results/summary.tsv  --pipeline-out pipeline_out --out final_report.tsv
python3 scripts/genus_rollup.py --report final_report.tsv --out genus_summary.tsv
python3 scripts/pivot_genus_matrix.py --genus-summary genus_summary.tsv --out final_matrix.tsv
```
