# Methods

Based on the "NGS data analysis" section of PLoS ONE 0156426 (read processing, taxonomic assignment against a
specialised reference database, divergence from the reference), adapted to long-range-PCR amplicons of the
whole nuclear rRNA cistron (18S-ITS1-5.8S-ITS2-28S). All classification is local and offline; no public
general database (core_nt) is used.

## 1. Read processing
1. Primer removal (Cutadapt; the lab's 18S-cistron primer table in `config/primers.fasta`, 33 oligos in sets 1-3,
   anchored at the read start on R1 and R2, error rate 0.1; arguments in `primers.cutadapt_args`).
2. fastp: adapter/quality trimming and merging of overlapping pairs → merged reads + unmerged R1/R2.
   Total reads of a sample = merged + R1 + R2 (post-fastp). All fractions use this denominator.
3. clumpify (BBTools) removes PCR duplicates **for assembly only**; quantification uses non-deduplicated reads
   (duplicate collision rises with depth and would flatten real abundance differences).

## 2. Assembly and ITS extraction
4. Split assembly: unmerged pairs → metaSPAdes; merged reads → SPAdes single-end (`--only-assembler`);
   contigs pooled and made non-redundant with cd-hit-est (100 %).
5. Chimera filter: vsearch `--uchime_denovo` (abundance from SPAdes coverage).
6. ITSx locates ITS1/5.8S/ITS2; one region per contig (full ITS1-5.8S-ITS2 > ITS1+ITS2 > single spacer).
   Only contigs with an extracted ITS region are used as ITS contigs (contigs on which ITSx found only
   18S/28S anchors are not; see LIMITATIONS and CHANGES).
7. Query set per sample (automatic, `mode` in `qc_summary.tsv`): `its` = ITS regions; `contigs` = assembly OK
   but no ITS: the 10 longest contigs (≥100 bp); `none` = assembly failed.

## 3. Quantification (every read counted once)
8. All post-fastp reads are mapped (bowtie2) to the query contigs; counts per contig (samtools idxstats, no MAPQ filter).
9. Human screen: reads not on a contig are mapped to GRCh38; those reads are assigned to *Homo sapiens*.
10. Remaining reads ("residual") are dereplicated (exact sequence). The 50 000 most abundant sequences are
    classified (`residual.max_unique_sequences`); rarer sequences are counted as Unclassified.
11. Per sample, the rows add up to the post-fastp read total (checked for all 32 samples; a mismatch is a QC warning).

## 4. Primary classification: custom ITS database
12. `blastn -task blastn` against the local ITS database (`db/build_v2`, 108 834 references; see DATABASE.md).
    Contig queries need ≥50 % query coverage; read queries need ≥80 % coverage and ≥90 % identity. Fallback
    contigs (mode `contigs`) need ≥90 % identity for any call.
13. LCA: all hits within 2 % of the best bit score form the tie group; the call is their lowest common ancestor
    in the NCBI taxonomy (sub-specific nodes reported as species). Tie between congeners → **genus**.
    Identity caps: best identity < 97 % → not species; < 90 % → not genus; < 80 % → not family. Calls broader
    than family are reported as unassigned.

## 5. Secondary classification: 18S/28S (nuclear rDNA) database -- only where ITS gives no call
14. Queries: (a) the fallback contigs of samples in mode `contigs`; (b) residual read sequences that have no
    informative ITS call. A read or contig with an ITS call keeps it; the rDNA database never overrides ITS.
15. `blastn -task megablast` against the local rDNA database (`db/build_rdna`, 51 526 references; see DATABASE.md),
    up to 1000 hits per query (18S/28S produce many near-identical hits; too few would make the LCA too narrow).
16. Same LCA, with cut-offs set by leave-one-out calibration (section 7): **no genus or species calls**;
    a capped genus call stops at its tribe; tribe and family/subfamily calls need **≥95 % identity**; hits below
    95 % identity are ignored. Settings: `rank_caps: species:101,genus:101,tribe:95,family:95`, `min_pident: 95`.
17. Rows from this step are labelled `rDNA DB (18S/28S) / contig, no ITS` or `rDNA DB (18S/28S) / read-level`
    in `final_table_detailed.tsv`.

## 6. Divergence from the reference
18. Every assigned taxon carries a divergence class from the best hit's identity: identical (≥99.5 %),
    low (97-99.5 %), moderate (93-97 %), high (<93 %). Adapted from the paper's 0 / 0-2 / 2-5 / >5 % distance bins
    to the identity cut-offs used here and to Illumina error rates.

## 7. Calibration of the 18S/28S cut-offs (`db/validate_rdna_db.py`)
19. 2 293 plant references (random sample plus ≥60 per target family) were cut into 6 162 query fragments
    of 1000 bp (contig-like) and 150 bp (read-like); whole 18S-ITS-28S records contribute only their 18S and 28S
    ends, so the ITS cannot help. Each fragment was classified with its own **species** hidden, and again with
    its whole **genus** hidden (the case of a genus without 18S/28S in GenBank, e.g. *Leuzea*).
20. Result for the chosen setting (`tribe:95,family:95`, genus never called):

    | | correct tribe | correct family | no call | wrong within family | wrong family |
    |---|---|---|---|---|---|
    | species hidden, 1000 bp | 35.6 % | 43.8 % | 19.1 % | 0.5 % | 1.1 % |
    | species hidden, 150 bp  | 27.8 % | 30.7 % | 40.2 % | 0.3 % | 1.0 % |
    | genus hidden, 1000 bp   | 18.8 % | 30.8 % | 46.4 % | 1.3 % | 2.6 % |
    | genus hidden, 150 bp    | 12.9 % | 20.1 % | 62.6 % | 1.5 % | 2.9 % |

    Expected error: ~1.5 % of calls when the genus is in the database, ~4 % when it is not.
    Why no genus calls: with the genus hidden, genus calls at ≥99 % identity were wrong in 4.5-5.3 % of
    fragments (9.3 % of 1000-bp calls at ≥99.5 % identity). Why 95 % and not 97 %: 97 % lost ~7 points of correct
    calls for ~0.2 points less error. Full tables: `db/build_rdna/validation_tribe/`.

## 8. Reporting
21. `final_table.tsv`: sample, taxon, rank, reads, % of all reads (the 5 requested columns; reads of the same
    taxon from different evidence lines are summed). `final_table_detailed.tsv`: one row per taxon × evidence ×
    divergence class, with taxid, group and % of plant reads. `qc_summary.tsv`: per-sample read accounting,
    including `rdna_assigned_reads` and `rdna_plant_reads`. `versions.txt`: tool versions and both database manifests.

## 9. Not used: core_nt control (`core_nt.mode: skip`)
The code for an optional control search against NCBI core_nt is kept but was not run; the specialised local
databases are the only classification method.
