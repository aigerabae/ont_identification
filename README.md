# plant_its_pipeline

Species composition of mixed herbal/plant material from Illumina MiSeq reads of the rRNA cistron
(18S–ITS1–5.8S–ITS2–28S). Snakemake workflow; **primary classification = local BLAST against a
custom ITS reference database**; an NCBI core_nt control is optional and off by default.

## Final output
`results/final/final_table.tsv` – one row per taxon per sample:

| Образец | Таксон | Таксономический ранг | Количество прочтений | Доля прочтений, % |
|---|---|---|---|---|

(`report_language: en` in the config gives English headers.) Every read of a sample lands in exactly one row
(taxon, human, other non-plant taxon, or "Unclassified"); fractions sum to 100 %.
When several species tie (hits within 2 % of the best bit score) the **genus** is reported.
More columns (taxid, evidence, divergence from the database, % of plant reads): `final_table_detailed.tsv`;
per-sample QC: `qc_summary.tsv`; custom-DB vs core_nt agreement per ITS contig: `concordance.tsv`;
tool versions + DB manifest: `versions.txt`.

## Quick start
```bash
mamba env create -f environment.yml && conda activate plant_its
bash tests/run_test.sh                    # 5-minute synthetic end-to-end test, no internet
db/get_taxdump.sh                         # NCBI taxonomy
db/build_custom_db.sh                     # download + ITSx-trim + makeblastdb (hours; see docs/DATABASE.md)
db/build_human_index.sh                   # human reference for the contamination screen
# edit config/config.yaml: raw_fastq_dir, primers.cutadapt_args, core_nt.mode
snakemake --profile profiles/default -n   # dry run
snakemake --profile profiles/default
```
Documentation: `docs/METHODS.md` (methods), `docs/DEPLOYMENT.md` (deploying the DB locally),
`docs/DATABASE.md` (contents, seed lists, validation), `docs/LIMITATIONS.md`,
`docs/CHANGES_FROM_PREVIOUS_VERSION.md`, `docs/TIMELINE.md`.

## Layout
`workflow/Snakefile` · `scripts/` (python helpers) · `db/` (database build + validation) ·
`config/` · `profiles/default` · `tests/` · `legacy/` (old scripts, not used).
