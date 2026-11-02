#!/bin/bash
# Full run, all samples: 18S/28S BLAST at 32 cores, then classification + reports at 16 cores.
set -uo pipefail
cd ~/biostar/NCB/hp_august2
echo "== $(date '+%F %T') start"

# 1. config check: calibrated 18S/28S cut-offs and read cap must be in place
CFG=$(sed -n '/^rdna:/,/max_unique/p' config/config.yaml)
n_caps=$(grep -c 'rank_caps: "species:101,genus:101,tribe:95,family:95"' <<< "$CFG")
n_pid=$(grep -c 'min_pident: 95' <<< "$CFG")
grep -q 'max_unique: 50000' <<< "$CFG" && cap_ok=1 || cap_ok=0
if [ "$n_caps" -ne 2 ] || [ "$n_pid" -ne 2 ] || [ "$cap_ok" -ne 1 ]; then
    echo "STOP: rdna section of config/config.yaml is not as expected:"; echo "$CFG"; exit 1
fi
echo "config OK: tribe:95,family:95 / min_pident 95 / max_unique 50000"

snakemake --profile profiles/default --unlock > /dev/null 2>&1

# 2. 18S/28S BLASTs (low memory; 8 parallel at 4 threads each)
echo "== $(date '+%F %T') step 1/2: 18S/28S BLAST (32 cores)"
if ! snakemake --profile profiles/default --cores 32 --until blast_rdna > rdna_blast.log 2>&1; then
    echo "STOP: BLAST stage failed -- see: grep -n -A25 'Error in rule' rdna_blast.log"; exit 1
fi
echo "   $(grep 'steps (' rdna_blast.log | tail -1)"

# 3. everything else: 18S/28S classification, sample reports, final tables, versions
echo "== $(date '+%F %T') step 2/2: classification + reports (16 cores)"
if ! snakemake --profile profiles/default --cores 16 > rdna_rest.log 2>&1; then
    echo "STOP: second stage failed -- see: grep -n -A25 'Error in rule' rdna_rest.log"; exit 1
fi
echo "   $(grep 'steps (' rdna_rest.log | tail -1)"

# 4. quick check of the result
echo "== $(date '+%F %T') done"
echo "samples in final_table.tsv: $(tail -n +2 results/final/final_table.tsv | cut -f1 | sort -u | wc -l) (expected 32)"
grep -c "rDNA DB" results/final/final_table_detailed.tsv | xargs echo "rows from the 18S/28S pass:"
awk -F'\t' 'NR>1 && $NF!=""{print "warning:", $1, $NF}' results/final/qc_summary.tsv
