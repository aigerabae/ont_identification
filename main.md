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
