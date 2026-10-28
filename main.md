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
bash db/build_human_index.sh
ls db/taxdump/        # must list nodes.dmp names.dmp merged.dmp
```
