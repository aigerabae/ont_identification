# Timeline (indicative; machine time measured on the 32-core workstation, Oct 2026)

| Step | Time |
|---|---|
| Build + freeze ITS database (download + ITSx), length filter | 0.5-1 day (ITSx hours); filter minutes |
| Build + freeze 18S/28S database (download, with NCBI API key) | ~1.5 h |
| Validation (ITS) and calibration (18S/28S) | ~30 min each |
| Human index, taxdump, config (primers) | 0.5 day |
| Run 32 samples: ITS part | ~30 min |
| Run 32 samples: 18S/28S BLAST (32 cores) + reports | ~3-4 h |
| Review results, adjust seed list / cut-offs, rerun if needed | 1-2 days |
| Final tables + report text | 1 day |
| Validation on the finished pipeline (supervisor) | as planned by the lab |

Reruns after a settings change repeat only what changed: a new cut-off reruns the LCA and the reports
(minutes); a new database or read cap reruns the corresponding BLASTs (hours).
