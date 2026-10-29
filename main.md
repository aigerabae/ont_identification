Test run:
```
cd ~/biostar/NCB/hp_august2
conda env create -f environment.yml && conda activate plant_its
bash tests/run_test.sh        # optional: about 5 minutes, no internet; should print the test table
```

Building database (once):
```
cd ~/biostar/NCB/hp_august2
conda activate plant_its
chmod +x db/*.sh db/*.py scripts/*.sh

bash tests/run_test.sh          # optional sanity check, ~5 min; wait for the final table
bash db/get_taxdump.sh          # must finish first: creates db/taxdump/nodes.dmp
bash db/build_custom_db.sh      # hours; run it in tmux or screen so it survives a closed terminal
REUSE_ITSX=1 bash db/build_custom_db.sh 2>&1 | tee db/build/rerun.log                   # download completed but test didn't work; I re wrote the build custom db script and reran without re-downloading the complete database

awk '/^>/{p=($0==">EF104182.1_its12")} p' db/build/refs.fasta > /tmp/q.fa
grep -vc '>' /tmp/q.fa          # sequence lines, should be more than 1
blastn -task blastn -db db/build/db/its_custom -query /tmp/q.fa -dust no -max_target_seqs 500 \
  -outfmt "6 qseqid sacc staxids pident qcovs length" | sort -k4,4nr -k5,5nr | head -5

rm -f db/taxdump/taxdump.tar.gz         # previous download didn't complete so i am deleting previous run and re-downloading
bash db/get_taxdump.sh
ls db/taxdump/        # must list nodes.dmp names.dmp merged.dmp
```

Downloading human for contamination check:
```
conda activate plant_its
bash db/build_human_index.sh 2>&1 | tee db/human_index.log
```

Validating the database (without humans):
```
cd ~/biostar/NCB/hp_august2 && conda activate plant_its
python3 db/validate_db.py --fasta db/build/refs.fasta --taxid-map db/build/refs.taxid.map \
  --blast-db db/build/db/its_custom --taxdump-dir db/taxdump --out db/build/validation \
  --max-queries 3000 --threads 8
cat db/build/validation/validation_summary.txt
```


Running anew:
```
cd ~/biostar/NCB/hp_august2 && conda activate plant_its
snakemake --profile profiles/default --cores 16 --config 'samples=[ONplants-1_S40,18S-1_S37]' 2>&1 | tee trial.log
```

After 20 hours it's still not done; The control for ONplants-1_S40’s five ITS contigs took 5 h 53 min (17:03 → 22:56) for five sequences. NCBI queues each web BLAST job. The pipeline allows only one NCBI job at a time, as NCBI’s usage policy asks, so everything else waits (“Waiting for more resources” is normal). It is now working on the residual-reads control, which can mean up to 8 batches of 25 sequences, each possibly hours, and then the second sample’s jobs. So the remote control would take days for two samples and weeks for 32. It isn’t usable at full scale.



I decided to revert back on the logic i used before and just add primers cutting using my primers + changing the database to local
```
conda activate plant_its
snakemake --profile profiles/default --cores 16 --forcerun select_queries \
  --config 'samples=[ONplants-1_S40,18S-1_S37]' 2>&1 | tee trial3.log
cat results/final/final_table.tsv
```
