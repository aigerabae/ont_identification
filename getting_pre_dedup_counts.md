```
snakemake -j 4 \
  $(ls pipeline_out/ | grep -E "^ONplants-|^18S-30_S32$|^18S-36_S36$" | sed 's#^#/home/aygera/biostar/NCB/hp_august/pipeline/pipeline_out/#; s#$#/qc/mapped_its_prededup.sorted.bam#') \
  -k
  
python3 scripts/recount_reads_prededup.py   --lca-summary blast_results/lca_summary.tsv   --pipeline-out /home/aygera/biostar/NCB/hp_august/pipeline/pipeline_out   --build-script scripts/build_final_report_lca.py   --mapq 10   --out read_counts_report_prededup.tsv
```


Getting species id for non assembled reads:
```
# 1. Generate FASTA files, one per sample, all unique reads (not just top 20)
python3 scripts/export_unique_reads_fasta.py \
  --pipeline-out pipeline_out \
  --samples 18S-1_S37 18S-3_S41 18S-9_S38 18S-31_S33 \
  --out-dir failed_samples_diagnostic

# 2. For each sample, open the resulting *_all_unique.fasta, paste its full
#    contents into https://blast.ncbi.nlm.nih.gov/Blast.cgi (blastn, core_nt),
#    run it, then click "Download All" -> "Text" and save as e.g.
#    18S-1_S37_results.txt (repeat for all 4 samples)

# 3. Parse all four downloaded files together into one report
python3 scripts/parse_ncbi_text_results.py \
  failed_samples_diagnostic/18S-1_S37_results.txt \
  failed_samples_diagnostic/18S-3_S41_results.txt \
  failed_samples_diagnostic/18S-9_S38_results.txt \
  failed_samples_diagnostic/18S-31_S33_results.txt \
  --out species_id_report.tsv

cat species_id_report.tsv
```
