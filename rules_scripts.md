Here's the complete rule set across all four scripts, in the order data flows through them.

## `summarize_blast_hits.py`
*Input: raw `blastn -remote` tabular output per sample → Output: ranked, deduplicated hits per contig*

1. **Column contract**: expects exactly 10 tab-separated fields per row (`qseqid, sacc, sscinames, scomnames, pident, length, qcovs, evalue, bitscore, stitle`). If a row has *more* than 10 fields, everything past the 9th tab is rejoined into `stitle` (guards against stray tabs in the title text). If it has *fewer*, the row is silently skipped.
2. **Numeric fields must parse**: `pident`, `qcovs`, `evalue`, `bitscore` must convert to `float`; any row that fails is dropped.
3. **Ranking**: within each contig (`qseqid`), hits are sorted by `bitscore` descending — highest-scoring alignment first.
4. **Species-name resolution**: if `sscinames` is missing (empty string, `"N/A"`, `"NA"`, or `"-"` — case-insensitive), the organism name is instead extracted from `stitle` via regex: an optional `PREDICTED:` prefix is stripped, then the first `Genus species` (optionally followed by `subsp. X` or `var. X`) is captured. If that regex fails to match anything at all, the first 60 characters of the title are used verbatim as a fallback label.
5. **Deduplication**: within a contig's ranked hit list, only the *first* occurrence of each distinct species-name (resolved per rule 4) is kept — later hits to the same species (e.g. 10 near-identical GenBank accessions) are dropped. **Critically, this same resolved name is what's compared for dedup**, not the raw `sscinames` string — otherwise every hit with literal `"N/A"` would collapse into one, deleting real alternative species (this was the bug found and fixed mid-conversation).
6. **Cap**: at most `--top-n` (default 5) distinct-species hits are kept per contig.
7. **Output**: one row per (sample, contig, rank).

## `build_final_report.py`
*Input: `summary.tsv` + per-sample `depth_its_mean.tsv` → Output: species abundance per sample*

1. **Only rank 1 is used** — the single best BLAST hit per contig is what feeds the report; ranks 2+ are informational only (visible in `summary.tsv` for manual review) but don't contribute to abundance math.
2. **Quality filter**: a hit is excluded unless `bitscore >= --min-bitscore` (default 100) **and** `evalue <= --max-evalue` (default 1e-15). Excluded hits are written to a `.excluded_low_confidence.tsv` sidecar file rather than silently dropped.
3. **Contig-ID suffix stripping for the join**: `extract_its_regions.py` appends `_ITS1_ITS2_combined` to any contig whose ITS1+ITS2 were concatenated, but the depth file's contig IDs come from the original assembly (no suffix). This suffix is stripped via regex before looking up that contig's depth, so those rows aren't silently treated as "no depth found."
4. **Missing depth defaults to 0.0** — if a contig genuinely isn't in the depth file, it contributes zero rather than crashing or being dropped.
5. **Multiple contigs → one species, depth summed**: all contigs in a sample whose top hit resolved to the same species have their depths **added together** before computing percentage — this is what correctly reproduced your pilot's *Hippophae* case (2-3 separate contigs, one combined abundance).
6. **Representative hit per species**: when multiple contigs support the same species, the one with the **highest bitscore** is kept as the "best_pident / best_bitscore / best_accession" record shown for that species row.
7. **Percentage basis**: `% = species_depth / sum(all passing species' depths in that sample) × 100` — i.e. percentages are relative to the sample's total *passing* ITS-region depth, not to total sequencing depth.
8. **Missing depth file for a whole sample**: if a sample's `depth_its_mean.tsv` doesn't exist or is empty, its species calls are still included (with depth 0), and the sample name is listed in a printed warning rather than causing a crash.

## `genus_rollup.py`
*Input: `final_report.tsv` → Output: genus-level totals per sample*

1. **Genus extraction**: strip the `" (from title, taxdb not installed)"` marker if present, then take the first capitalized word matching `^[A-Z][a-z]+` as the genus. If nothing matches, the full (cleaned) name is used as its own bucket.
2. **"Uncultured" special case**: any name starting with `uncultured` (case-insensitive) is bucketed as `"Uncultured/environmental (unresolved)"` instead of trying to extract a genus — there isn't a real genus to roll up to, and this is a deliberately distinct interpretive category (background signal, not a named organism).
3. **Synonym merging**: a small hardcoded table (`{"Leuzea": "Rhaponticum"}`) merges genera known to be **true taxonomic synonyms** of the same species — applied *before* bucketing, so `Leuzea carthamoides` and `Rhaponticum carthamoides` land in the same genus bucket.
4. **Deliberately not merged**: `Serratula`, `Stemmacantha`, `Acroptilon` — close relatives of *Rhaponticum* in the same tribe (Cardueae), but *not* confirmed synonyms, so they're kept as separate genus buckets rather than assumed to be the same species.
5. **Per-sample, per-genus depth is summed** across every species label that collapsed into that genus (mirrors rule 5 in `build_final_report.py`, one level up).
6. **Percentage recomputed at genus level**: `% = genus_depth / sample's total depth (across all genera present)`, not just re-summed from the species-level percentages (avoids compounding rounding error).
7. **Species labels preserved for audit**: every genus row keeps a semicolon-joined list of the original species-level labels that fed into it (`species_labels_seen`), so you can trace a genus total back to what actually produced it.
8. **Watch-list warning**: if any of the deliberately-unmerged tribe-mate genera (rule 4) appear in the output, a printed note flags them for manual review rather than letting them pass silently as if they were fully resolved.

## `pivot_genus_matrix.py`
*Input: `genus_summary.tsv` → Output: sample × genus wide matrix*

1. **Sparse-to-dense conversion**: every (sample, genus) pair not present in the input is filled with `0.0` in the output — no blank cells, so every value is directly comparable.
2. **Column order = prevalence-first**: genera are ordered by how many samples they appear in (descending), ties broken alphabetically — so the most biologically important columns (e.g. *Hippophae*, *Rhaponticum*, *Rhodiola*) appear first, and rare contamination hits are pushed to the right.
3. **Row order = first-seen order** in the input file (not re-sorted), preserved via a list+set combo for both correctness and lookup speed.
4. **`row_total_pct` sanity column**: computed by re-summing that row's genus percentages — should land near 100 given how upstream percentages were built; a lower total is a signal (not a guarantee of a bug) that some of that sample's depth came from hits the quality filter in `build_final_report.py` excluded.
