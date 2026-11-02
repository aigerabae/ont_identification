# How I ran the pipeline (run log)

All commands are run from the repository root in the `plant_its` conda environment:
```
cd ~/biostar/NCB/hp_august2 && conda activate plant_its
```

## 1. Setup and test run
```
conda env create -f environment.yml && conda activate plant_its
chmod +x db/*.sh db/*.py scripts/*.sh
bash tests/run_test.sh        # optional: about 5 minutes, no internet; should print the test table
```

## 2. ITS database (once)
```
bash db/get_taxdump.sh          # must finish first: creates db/taxdump/nodes.dmp
ls db/taxdump/                  # must list nodes.dmp names.dmp merged.dmp
bash db/build_custom_db.sh      # hours; run it in tmux/screen or with nohup so it survives a closed terminal
```
The download completed but the self-test didn't work; I rewrote the build script and reran it without
re-downloading the database:
```
REUSE_ITSX=1 bash db/build_custom_db.sh 2>&1 | tee db/build/rerun.log
```
Spot check that a reference finds itself:
```
awk '/^>/{p=($0==">EF104182.1_its12")} p' db/build/refs.fasta > /tmp/q.fa
grep -vc '>' /tmp/q.fa          # sequence lines, should be more than 1
blastn -task blastn -db db/build/db/its_custom -query /tmp/q.fa -dust no -max_target_seqs 500 \
  -outfmt "6 qseqid sacc staxids pident qcovs length" | sort -k4,4nr -k5,5nr | head -5
```
An earlier taxdump download hadn't completed, so I deleted it and downloaded it again:
```
rm -f db/taxdump/taxdump.tar.gz
bash db/get_taxdump.sh
ls db/taxdump/
```

## 3. Human reference (contamination screen)
```
bash db/build_human_index.sh 2>&1 | tee db/human_index.log
```

## 4. Validating the ITS database
```
python3 db/validate_db.py --fasta db/build/refs.fasta --taxid-map db/build/refs.taxid.map \
  --blast-db db/build/db/its_custom --taxdump-dir db/taxdump --out db/build/validation \
  --max-queries 3000 --threads 8
cat db/build/validation/validation_summary.txt
```

## 5. First runs: the remote core_nt control does not scale
Trial with two samples:
```
snakemake --profile profiles/default --cores 16 --config 'samples=[ONplants-1_S40,18S-1_S37]' 2>&1 | tee trial.log
```
After 20 hours it was still not done. The control for ONplants-1_S40's five ITS contigs took 5 h 53 min
(17:03 → 22:56) for five sequences. NCBI queues each web BLAST job, and the pipeline allows only one NCBI job at a
time, as NCBI's usage policy asks, so everything else waits ("Waiting for more resources" is normal). The
residual-read control could mean up to 8 more batches of 25 sequences, each possibly hours, then the second
sample. So the remote control would take days for two samples and weeks for 32: not usable at full scale.

I went back to the logic I used before, added primer trimming with the lab's primers, and switched the
classification to the local database only (`core_nt.mode: skip` in `config/config.yaml`):
```
snakemake --profile profiles/default --cores 16 --forcerun select_queries \
  --config 'samples=[ONplants-1_S40,18S-1_S37]' 2>&1 | tee trial3.log
cat results/final/final_table.tsv
```
The two samples ran as expected. Then the full 32-sample run:
```
grep -n "raw_fastq_dir" config/config.yaml      # must point to the real FASTQ folder, e.g.:
# sed -i 's#/path/to/raw_fastqs#/mnt/harddisk/biostar/NCB/hp_august/raw_fastqs#' config/config.yaml
snakemake --profile profiles/default --cores 32 2>&1 | tee full_run.log
```
Rerun of the LCA step after adding the identity floor for fallback contigs
(`classification.fallback_contigs_min_pident`):
```
grep -A3 "^core_nt" config/config.yaml        # must show: mode: skip
snakemake --profile profiles/default --cores 4 --forcerun lca_custom 2>&1 | tee rerun5.log
```
Note: `lca_custom` runs at `--cores 4` because each LCA job used to load all BLAST hits into memory; at 32 cores
it froze the machine (fixed later, section 9).

## 6. Removing the false fungal calls
Diagnostics showed that a few ITS references still carried 18S/28S flanks (e.g. AF250821.1_its12, 2108 bp) and
attracted plant 18S/28S reads as "fungi", and that fallback-contig calls rested on matches covering only 5-11 %
of the contig. Fixes:

Clean ITS database (`db/build_v2`): drop records whose length is impossible for the region they claim
(single spacer > 450 bp, ITS1-5.8S-ITS2 > 1000 bp; animal decoys and gymnosperms exempt). 68 of 108 902 dropped.
The same filter is now step 4b of `db/build_custom_db.sh`.
```
mkdir -p db/build_v2/db
python3 db/filter_long_refs.py --in-fasta db/build/refs.fasta --in-map db/build/refs.taxid.map \
  --taxdump-dir db/taxdump --out-fasta db/build_v2/refs.fasta --out-map db/build_v2/refs.taxid.map \
  --out-dropped db/build_v2/dropped_records.tsv
makeblastdb -in db/build_v2/refs.fasta -dbtype nucl -parse_seqids \
  -taxid_map db/build_v2/refs.taxid.map -out db/build_v2/db/its_custom -title its_custom_v2
md5sum db/build_v2/refs.fasta db/build_v2/refs.taxid.map | tee db/build_v2/MD5.txt
# expected: 0346d6202a0d84e27ed59632c5b00ef7 / 433690d73f1f3cdb22711e6aeac462cf
```
(I also wrote `db/build_v2/MANIFEST.txt` by hand with the build date, filter settings, counts and md5, so
`versions.txt` picks it up.)

Config changes: point the pipeline to the cleaned database and require ≥50 % coverage for contig hits:
```
cp config/config.yaml config/config.yaml.before_v2
sed -i 's#db/build/db/its_custom#db/build_v2/db/its_custom#' config/config.yaml
# classification.contigs.min_qcov: 0 -> 50 (edited by hand)
snakemake --profile profiles/default --cores 4 --forcerun blast_custom lca_custom
```
Result: fungal reads 117 064 → 5 038, plant ITS calls unchanged.

## 7. 18S/28S (nuclear rDNA) database for the samples without ITS
Nine ONplants samples had no ITS contigs and the 18S samples had no ITS at all; their plant DNA sits on 18S/28S
contigs (manual NCBI BLAST of ONplants-2 and 18S-34 showed Elaeagnaceae and Cardueae). So I built a second local
database of 18S/28S records, used only where ITS gives no call.

I use my NCBI API key, entered once into a private file (the build scripts read it and never print it):
```
read -s -p "NCBI API key: " K && echo "$K" > ~/.ncbi_api_key && chmod 600 ~/.ncbi_api_key && unset K && echo " saved"
```
Build (~1.5 h; resumable, finished groups are skipped on a rerun):
```
nohup db/build_rdna_db.sh > db/build_rdna.log 2>&1 &
tail -40 db/build_rdna.log        # ends with "self-test OK" and "Done. DB prefix: db/build_rdna/db/rdna_custom"
# 171 592 records downloaded -> 51 526 references; md5 refs.fasta 2b34be0085ffb515f2a43f0764a8153d
```

## 8. Calibrating the 18S/28S cut-offs (leave-one-out)
Each test fragment (1000 bp and 150 bp from 2 293 plant references) is classified with its own species hidden,
and with its whole genus hidden:
```
nohup python3 db/validate_rdna_db.py --fasta db/build_rdna/refs.fasta --taxid-map db/build_rdna/refs.taxid.map \
  --blast-db db/build_rdna/db/rdna_custom --taxdump-dir db/taxdump --out db/build_rdna/validation \
  --threads 16 > db/build_rdna/validation.log 2>&1 &
```
Rerun of the classification part only (reuses the BLAST results), after adding the tribe option:
```
mkdir -p db/build_rdna/validation_tribe && cp db/build_rdna/validation/loo.blast.tsv db/build_rdna/validation_tribe/
python3 db/validate_rdna_db.py --fasta db/build_rdna/refs.fasta --taxid-map db/build_rdna/refs.taxid.map \
  --blast-db db/build_rdna/db/rdna_custom --taxdump-dir db/taxdump --out db/build_rdna/validation_tribe \
  --reuse-blast > db/build_rdna/validation_tribe.log 2>&1
cat db/build_rdna/validation_tribe/calibration_summary.txt
```
Decision: never call genus from 18S/28S (genus calls were wrong 4.5-5.3 % of the time when the true genus is
missing); tribe and family need ≥95 % identity. Expected error ~1.5 % (genus in database) to ~4 % (genus missing).
Final `rdna:` settings in `config/config.yaml`:
```
sed -i 's/min_pident: 97/min_pident: 95/' config/config.yaml
sed -i 's/species:101,genus:101,family:9[57]/species:101,genus:101,tribe:95,family:95/; s/max_unique: 20000 /max_unique: 50000 /' config/config.yaml
sed -n '/^rdna:/,/max_unique/p' config/config.yaml | grep -E "contigs:|residual:|min_pident|rank_caps|max_unique"
```

## 9. Testing the 18S/28S pass on 4 samples, and the memory fix
One sample from each group, asking only for their per-sample reports (so `results/final/` is not overwritten):
```
T="results/ONplants-3_S42/report/rows.tsv results/ONplants-17_S56/report/rows.tsv results/18S-34_S43/report/rows.tsv results/18S-1_S37/report/rows.tsv"
snakemake --profile profiles/default --cores 16 -n $T | sed -n '/^job/,/^total/p'   # dry run first
nohup snakemake --profile profiles/default --cores 16 $T > preview.log 2>&1 &
```
ONplants-17's LCA step was killed by the system (out of memory: 20 million BLAST rows, 15.8 GB). I updated
`scripts/lca_assign.py` to classify one read at a time (identical output, 1.3 GB instead of 15.8 GB) and added
the optional `tribe` cut-off. Check:
```
/usr/bin/time -v python3 scripts/lca_assign.py --blast results/ONplants-17_S56/classification/rdna_residual.blast.tsv \
  --sample ONplants-17_S56 --taxdump-dir db/taxdump --out /tmp/test17.lca.tsv --tie-fraction 0.02 \
  --min-bitscore 100 --max-evalue 1e-10 --min-pident 95 --min-qcov 80 --rank-caps "species:101,genus:101,family:95" \
  2>&1 | grep -E "LCA calls|Maximum resident|Exit status"
```

## 10. Final full run (all 32 samples)
`run_all.sh` checks the `rdna:` settings, runs the 18S/28S BLASTs at 32 cores (`--until blast_rdna`, eight at a
time, ~180 MB each), then the LCA, reports and final tables at 16 cores (~1.3 GB per LCA job). Start it once only:
```
chmod +x run_all.sh
nohup ./run_all.sh > run_all.log 2>&1 &
cat run_all.log                                   # stage
pgrep -c -f run_all.sh                            # must be 1
grep -c "Finished job" rdna_blast.log; pgrep -c -x blastn     # BLAST progress
grep "steps (" rdna_rest.log | tail -1                         # second stage
```
It took ~3-4 h. `run_all.log` ends with `done` and `samples in final_table.tsv: 32 (expected 32)`.

Optional: upload the results to Google Drive when the run finishes (`rclone config` once, remote `gdrive`):
```
nohup ./upload_when_done.sh > upload.log 2>&1 &
```

Results are in `results/final/`: `final_table.tsv` (the 5 requested columns), `final_table_detailed.tsv`
(per evidence line: ITS / 18S/28S / human screen, with divergence class), `qc_summary.tsv`, `versions.txt`
(tool versions and both database manifests).

## Useful
```
snakemake --profile profiles/default --unlock      # after a crash ("Directory cannot be locked")
pgrep -af "snakemake|blastn"                        # what is running
snakemake --profile profiles/default -n | sed -n '/^job/,/^total/p'   # what a run would redo
```
