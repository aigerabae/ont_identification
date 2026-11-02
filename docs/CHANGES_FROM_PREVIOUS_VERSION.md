# Changes relative to the previous version (and why)

## A. Previous report (core_nt) → custom-database pipeline
| Supervisor comment | Change |
|---|---|
| Custom local reference DB as primary method | `db/` builds it (ITSx-trimmed, taxonomy-annotated, decoys, frozen files); `blast_custom` + LCA is the classification; core_nt is not used |
| Include related species / other taxa | whole families/genera in the seed list; decoys for human/insects/fungi |
| Genus when species tie | LCA over hits within 2 % of the best score; identity caps |
| Output table format | `final_table.tsv` with the 5 requested columns |
| Reproducible, finished, documented | one Snakemake workflow for all 32 samples, config, env, tests, docs |

Other fixes: quantification on non-deduplicated reads; one denominator (merged+R1+R2); every read in exactly one
row; human reads counted from a human reference screen; primer trimming (Cutadapt); divergence class per call;
failed-assembly samples handled by read-level classification; no manual sample lists; uninformative ranks
removed from the table.

## B. Revision of 30 Oct 2026 (after review of the first custom-database run)
| Finding | Change | Effect |
|---|---|---|
| 3 ITS references (2.1, 1.6, 1.6 kb) still carried 18S/28S flanks; plant 18S/28S reads matched them at 78-90 % and were reported as fungi (e.g. 12 % "*Paraphaeosphaeria*" in ONplants-16) | length filter in the database build (`filter_long_refs.py`); 68 records dropped → `db/build_v2` | fungal reads 117 064 → 5 038; plant ITS calls unchanged |
| Fallback-contig calls (e.g. *Scutellaria* in 18S-17/-23/-32) rested on 65-bp matches covering 5-11 % of the contig | `classification.contigs.min_qcov` 0 → 50 | these calls removed |
| ITS-less samples and 9 ONplants samples without ITS contigs had no plant composition, though manual NCBI BLAST showed plant 18S/28S contigs (Elaeagnaceae, Cardueae) | second local database of nuclear 18S/28S (`db/build_rdna`), used only where ITS gives no call; cut-offs from leave-one-out calibration (no genus; tribe/family ≥95 %) | samples with a plant composition 11 → 26; plant reads 1.30 M → 2.20 M |
| `lca_assign.py` held all hits in memory (15.8 GB for one sample, killed) | streaming, one query at a time (identical output) | 1.3 GB per job |
| A capped genus call skipped the tribe (single-genus hit → subfamily) | optional `tribe` cap, used only by the 18S/28S settings | Cardueae/Anthemideae reported as tribes where supported |

Note on comparing with the previous report: its "mapped to contig, no confident taxon" (up to 77 % in ONplants-7)
were reads on 18S/28S contigs, because it mapped reads to every contig on which ITSx detected anything.
This version maps only to contigs with an extracted ITS region and classifies the other reads with 18S/28S.

## C. Agreement with the previous report for samples without ITS (18S/28S vs core_nt)
| Sample | Previous (core_nt) | This version (18S/28S database) |
|---|---|---|
| 18S-17 | Mentheae 90 % | Nepetoideae 90 % + Mentheae 4 % |
| 18S-23 | *Humulus lupulus* 88 % | Cannabaceae 51 % |
| 18S-25 | Nepetoideae 49 %, *Phlomoides* 6 % | Nepetoideae 49 % + Mentheae 15 % |
| 18S-34 | *Cousinia* / Cardueae 99 % | Asteraceae 66 % |
| 18S-37 | *Zea mays* 98 % | Poaceae 43 % + Andropogoneae 24 % |
| 18S-4 | Anthemideae 20 % + fungi | Anthemideae 15 % + fungi 25 % |
| 18S-13 | asterids / *Lactuca* 38 % | plants 1.8 % (not reproduced) |
| 18S-32 | Asteraceae 96 % | plants 3.2 % (not reproduced) |
