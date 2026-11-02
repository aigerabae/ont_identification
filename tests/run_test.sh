#!/bin/bash
# End-to-end test on synthetic data (no internet needed). Run from the repository root:
#   bash tests/run_test.sh [snakemake args, e.g. -n]
# Requires the conda env of environment.yml to be active.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
T="${TEST_WORK:-tests/work}"      # set TEST_WORK=/some/fast/local/dir if the repo is on a network drive
rm -rf "$T/results" "$T/.snakemake"; mkdir -p "$T/db" "$T/human"
sed "s#__WORK__#$T#g" tests/test_config.yaml > "$T/config.yaml"
python3 tests/test_taxdump.py
python3 tests/make_test_data.py
bash db/make_blast_db_from_frozen.sh tests/data/refs.fasta tests/data/refs.taxid.map "$T/db/its_custom" "test ITS db"
bash db/make_blast_db_from_frozen.sh tests/data/core_nt_mock.fasta tests/data/core_nt_mock.taxid.map "$T/db/core_nt_mock" "mock core_nt"
HUMAN_FASTA=tests/data/human/human_ref.fasta bash db/build_human_index.sh "$T/human" 4
echo "test db" > "$T/MANIFEST.txt"
BLASTDB="$PWD/$T/db" snakemake --directory "$PWD" --configfile "$T/config.yaml" --cores 4 --keep-going \
    --rerun-incomplete --resources ncbi=1 --printshellcmds "$@"
echo; echo "=== final_table.tsv ==="; cat "$T/results/final/final_table.tsv" 2>/dev/null || true
