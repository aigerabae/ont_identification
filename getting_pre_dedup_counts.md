snakemake -j 4 \
  $(ls pipeline_out/ | grep -E "^ONplants-|^18S-30_S32$|^18S-36_S36$" | sed 's#^#/home/aygera/biostar/NCB/hp_august/pipeline/pipeline_out/#; s#$#/qc/mapped_its_prededup.sorted.bam#') \
  -k
