# Deploying the pipeline and the local databases

## Requirements
Linux, 32 GB RAM, ≥50 GB disk, Conda/mamba. Internet only for the one-off database builds; the pipeline itself
runs offline. Optional NCBI API key (free, from an NCBI account) for faster downloads.

## Option A: use the frozen databases (recommended, no download)
Copy `db/build_v2/` and `db/build_rdna/` (each with `refs.fasta`, `refs.taxid.map`, `MANIFEST.txt`) and
`db/taxdump/`, then rebuild the BLAST files:
```bash
makeblastdb -in db/build_v2/refs.fasta   -dbtype nucl -parse_seqids -taxid_map db/build_v2/refs.taxid.map   -out db/build_v2/db/its_custom
makeblastdb -in db/build_rdna/refs.fasta -dbtype nucl -parse_seqids -taxid_map db/build_rdna/refs.taxid.map -out db/build_rdna/db/rdna_custom
md5sum db/build_v2/refs.fasta db/build_rdna/refs.fasta     # compare with DATABASE.md
```

## Option B: build from GenBank
```bash
git clone <repo> && cd plant_its_pipeline
mamba env create -f environment.yml && conda activate plant_its
bash tests/run_test.sh                # must finish and print the test table (see tests/EXPECTED.md)
db/get_taxdump.sh                     # -> db/taxdump/ (also needed by the length filter)

# optional API key: saved once in a private file, read by both build scripts, never written to logs
read -s -p "NCBI API key: " K && echo "$K" > ~/.ncbi_api_key && chmod 600 ~/.ncbi_api_key && unset K

db/build_custom_db.sh                 # ITS database -> db/build/ (download tens of minutes; ITSx hours)
nohup db/build_rdna_db.sh > db/build_rdna.log 2>&1 &   # 18S/28S database -> db/build_rdna/ (~1.5 h)
db/build_human_index.sh               # -> db/human/GRCh38.* (~3 GB download)

# validation / calibration (optional, recommended after every rebuild)
python3 db/validate_db.py --fasta db/build/refs.fasta --taxid-map db/build/refs.taxid.map \
        --blast-db db/build/db/its_custom --taxdump-dir db/taxdump --out db/build/validation
python3 db/validate_rdna_db.py --fasta db/build_rdna/refs.fasta --taxid-map db/build_rdna/refs.taxid.map \
        --blast-db db/build_rdna/db/rdna_custom --taxdump-dir db/taxdump --out db/build_rdna/validation --threads 16
```
Both build scripts resume after interruption (finished groups are skipped). If downloads fail with
`curl: (56) ... unexpected eof`, EDirect retries automatically; the build still completes.
Point `database.blast_db` and `rdna.blast_db` in `config/config.yaml` to the databases you built.

## Running
Edit `config/config.yaml` (`raw_fastq_dir`, sample names, `primers.cutadapt_args`); keep `core_nt.mode: skip`.
```bash
snakemake --profile profiles/default -n          # dry run
nohup ./run_all.sh > run_all.log 2>&1 &          # full run, unattended
```
`run_all.sh` checks the 18S/28S settings, runs the 18S/28S BLASTs at 32 cores
(`--until blast_rdna`; ~180 MB memory per BLAST), then everything else at 16 cores. Each LCA job needs ~1.3 GB
(the NCBI taxonomy); do not run many more than 16 at once on 32 GB. Re-running only repeats what changed.
Typical times on 32 cores: ITS part ~30 min; 18S/28S BLAST ~3-4 h; reports ~15 min.
If clumpify aborts with a JVM memory error, set `clumpify_memory_gb`.
To switch the 18S/28S pass off: `rdna.enabled: false` (all its rules then write empty outputs).
