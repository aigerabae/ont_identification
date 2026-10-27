What changed from the old version

Primary method: classification is now a local BLAST against a custom ITS database, with taxonomic annotation. The database covers whole families and genera plus decoys (human, insects, fungi).
core_nt role: it is only a control. It checks every ITS contig and anything the custom database couldn’t explain, and writes the agreement to concordance.tsv.
Ties: when several species score within 2% of the best hit, the table reports the genus. If the best identity is too low for a species or genus call, the call moves up a rank.
Output: final_table.tsv has his five columns, and report_language: en switches the headers to English. Every read of a sample lands in exactly one row, and the percentages sum to 100.
Failed samples: the four tiny samples go through read-level classification plus a human-reference screen, so no manual web BLAST is needed.
Fixes: it counts reads before deduplication, and the denominator is now merged + R1 + R2.
Docs and tests: the repo includes methods, deployment, database and limitations docs, a changes list and a timeline. There is also a leave-one-out validation script for the database.

What the synthetic test showed

Test data: five samples covering a clean mixture, near-identical congeners, a species missing from the database, a failed assembly with human reads, and a no-ITS sample.
Results: the ~60:40 mixture came out as 60.3:39.7. Near-identical congeners gave a genus call, not a false species. The failed-assembly sample came out 87% human. The no-ITS sample was named by the core_nt control.
Sample rule: in the mixed sample, Alpha is reported as a genus because the database holds a species one SNP away. That is by design.
Bugs fixed on the way: a symlink failure, a clumpify memory error (new clumpify_memory_gb setting) and a SPAdes cleanup error.

Not tested here

NCBI downloads (database build, taxdump, human genome), real core_nt, real primers, and the real 32 samples, because the sandbox blocks NCBI.
Before you send this to Vlad, run db/build_custom_db.sh and db/validate_db.py, and put your lab’s primers into primers.cutadapt_args in config/config.yaml.
The tests/work and tests/data folders were deleted before zipping; bash tests/run_test.sh regenerates them.
