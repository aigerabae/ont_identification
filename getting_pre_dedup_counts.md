```
snakemake -j 4 \
  $(ls pipeline_out/ | grep -E "^ONplants-|^18S-30_S32$|^18S-36_S36$" | sed 's#^#/home/aygera/biostar/NCB/hp_august/pipeline/pipeline_out/#; s#$#/qc/mapped_its_prededup.sorted.bam#') \
  -k
  
python3 scripts/recount_reads_prededup.py   --lca-summary blast_results/lca_summary.tsv   --pipeline-out /home/aygera/biostar/NCB/hp_august/pipeline/pipeline_out   --build-script scripts/build_final_report_lca.py   --mapq 10   --out read_counts_report_prededup.tsv
```
