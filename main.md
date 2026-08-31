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

Installing kraken and bracken:
```
conda create --name kraken2 and bracken
conda activate kraken2
conda install bioconda::kraken2
conda install -c bioconda -c conda-forge bracken
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

Only did some samples. Might want to enrich the specific regions first
