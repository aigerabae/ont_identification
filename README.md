# ont_identification

Studies of similar scope:  
https://pmc.ncbi.nlm.nih.gov/articles/PMC10317219/   
Validation of the shotgun metabarcoding approach for comprehensively identifying herbal products containing plant, fungal, and animal ingredients (Zhang et al., 2023)  

https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2024.1358136/full#h1  
Species identification of biological ingredients in herbal product, Gurigumu-7, based on DNA barcoding and shotgun metagenomics (Wei et al., 2024)  

https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2021.607200/full#s2    
The Species Identification in Traditional Herbal Patent Medicine, Wuhu San, Based on Shotgun Metabarcoding (Liu et al., 2021)  

https://pmc.ncbi.nlm.nih.gov/articles/PMC9187672/#Sec9  
Development of a DNA barcode library of plants in the Thai Herbal Pharmacopoeia and Monographs for authentication of herbal products (Urumarudappa, et al., 2022)  

Pipeline overview

Goal: identify plant species present in mixed-origin herbal medicine samples from long-range PCR-tiled rRNA cistron sequencing (18S–ITS1–5.8S–ITS2–28S), and estimate their relative proportions — implemented as a Snakemake workflow, run per sample, scalable across all 18.

Key steps
1) QC + merge (fastp) — adapter/quality trim, then merge overlapping R1/R2 pairs into single longer reads where fragment size allows.  
2) Deduplication (clumpify) — remove PCR duplicates, run separately on merged vs. unmerged reads.  
3) Split assembly — long-insert (unmerged) pairs → metaspades.py (paired mode); short-insert (merged) reads → spades.py --only-assembler (single-end mode). Combined contig sets afterward via cd-hit-est.  
4) Chimera detection (vsearch/uchime_denovo) — contigs annotated with abundance (;size= from SPAdes coverage), sorted, screened for chimeric assemblies before trusting any species call.  
5) ITS extraction (ITSx) — HMM-based identification of ITS1/5.8S/ITS2/flanking SSU-LSU boundaries per contig, run on assembled contigs (not raw reads).  
6) Region selection script — per ITS-flagged contig, picks the best available sequence for BLAST: full ITS1–5.8S–ITS2 span if available, else ITS1+ITS2 concatenated, else whichever single region ITSx recovered.  
7) BLAST (manual, NCBI web) — species ID per extracted sequence; deliberately not automated, since accuracy on ambiguous/divergent hits mattered more than throughput.  
8) Targeted depth mapping (bowtie2 + samtools) — reads mapped back onto ITS-flagged contigs only (not the full 250-contig assembly), MAPQ-filtered, to get a per-contig depth number as a relative-abundance proxy between species.  

Key decisions and why  
1) Split assembly by fragment type, not pooled. Mixing merged (short-insert) and unmerged (long-insert) reads in one metaSPAdes run broke insert-size estimation and repeat resolution; assembling separately and merging contigs afterward avoided this and let each branch recover sequence the other missed.  
2) BLAST the ITSx-extracted region, never the whole contig. Full contigs include long, highly conserved 18S/28S flanking sequence; whole-contig BLAST lets that conserved region dominate the score and return unrelated genera. Only the isolated ITS1/ITS2 window is species-diagnostic.  
3) Target bowtie2 mapping at ITS-flagged contigs only, not the full assembly. Mapping against all 250 contigs caused ~85–90% multi-mapping (near-identical short/junk contigs competing for reads), making depth numbers meaningless. Restricting to the handful of real ITS-bearing contigs collapsed multi-mapping to a manageable level and made depth differences between species interpretable.  
4) Chimera-check before, not after, trusting species calls. Cheap insurance against a false species ID propagating silently across all 18 samples.  
5) Depth ≠ abundance, stated explicitly. rRNA copy number varies between species and PCR efficiency varies between templates, so per-contig read depth is treated as a relative signal, not a precise proportion — this caveat is carried through rather than presented as a clean percentage.  





Comparing my pipeline to these papers:
## Quick comparison table

| | Your pipeline | Zhang 2023 | Wei 2024 (Gurigumu-7) | Liu 2021 (Wuhu San) | Urumarudappa 2022 (Thai) |
|---|---|---|---|---|---|
| Library prep | PCR-tiled long amplicon + shear | True shotgun, no PCR | True shotgun, no PCR | True shotgun, no PCR | Single-locus PCR |
| Sequencing | Illumina 2×150bp | Illumina, ~6.8Gb/sample | Illumina, ~7Gb/sample | Illumina | Sanger |
| Target locus | Whole 18S-28S cistron | ITS2, psbA-trnH, rbcL, matK, COI | ITS2, matK, rbcL | ITS2, rbcL, matK, psbA-trnH | ITS2, matK, rbcL, trnH-psbA |
| Assembler | Split metaSPAdes/SPAdes | Local pipeline (assembly-based) | Assembly-based | MEGAHIT + metaSPAdes | None (Sanger, no assembly) |
| Chimera check | vsearch uchime_denovo | Addressed as known PCR-free advantage | Chimeric detection step named | UCHIME | N/A |
| Reference DB | GenBank core_nt | TCM-BOL, BOLD, GenBank | TCM-BOL + own Sanger references | GenBank, TCM-BOL, BOLD | Custom-built barcode library |
| Adulterant framework | Not yet applied | Detects unlabeled species | Detects known toxic substitute | ASAC framework (Authentic/Substitution/Adulterant/Contaminant) | Label-vs-BLAST match |
| Quantification | Read-depth proxy (yours, in progress) | Not primary focus | Not primary focus | Read-depth/coverage via bowtie2 | N/A (presence/absence only) |

## Zhang et al. 2023 (Validation of shotgun metabarcoding)

**Similar**: Same overall logic as yours — recover marker regions from complex mixtures without relying on primer-anchored amplicon boundaries, then BLAST against curated databases (TCM-BOL/BOLD/GenBank, same idea as your core_nt use). They also validated against orthogonal traditional methods (microscopy/TLC/HPLC), the wet-lab equivalent of your Kraken2-vs-assembly cross-check.

**Different**: Their libraries are **true shotgun** (no PCR at all) — this is the key architectural difference from your data. Because there's no PCR step, they don't inherit the primer-efficiency bias or PCR-chimera risk your README flagged as a caveat; their chimera risk instead comes from library-prep artifacts, not amplification. They also extracted four barcode types (ITS2, psbA-trnH, rbcL, matK) plus COI for the animal component, giving them cross-marker corroboration you don't currently have (you're ITS-only). Why: they were mixture-composition-agnostic from the start (plant+fungal+animal), so needed a marker per kingdom; you're plant-focused with one long, opportunistically-tiled amplicon.

## Wei et al. 2024 (Gurigumu-7)

**Similar**: Same research lineage as Zhang 2023 (overlapping authors), same three-marker ITS2/matK/rbcL approach, same assembly-based extraction logic.

**Different, and worth adopting**: They generated **Sanger sequences from the individual known ingredients first**, then used those as curated references to disambiguate the shotgun-derived OTUs — a targeted reference-refinement step you haven't done. This is likely why they could confidently detect a specific toxic adulteration event (*Akebia quinata* substituted by *Aristolochia manshuriensis* in one sample) rather than stopping at "species present." Why this matters for you: your *Hippophae* calls plateaued at ~86-87% identity due to a real GenBank reference-coverage gap — Wei's approach (sequencing your own known/expected reference material via Sanger first) is a direct, cheap fix for exactly that kind of gap if you have access to authenticated reference specimens.

## Liu et al. 2021 (Wuhu San)

This is the one your own README already identified as closest — confirmed. **Similar**: same core logic (assemble target loci from fragmented reads rather than relying on clean amplicon boundaries), MEGAHIT+metaSPAdes assembly, cd-hit dedup, UCHIME chimera detection (your vsearch/uchime_denovo step is the modern equivalent), bowtie2 read-mapping back for depth/coverage, BLAST + database cross-referencing.

**Different**: Two things you don't have yet. First, they run **both MEGAHIT and metaSPAdes** in parallel for extra assembly diversity — you've noted this as a planned next step, not yet done. Second, and more consequential: they apply an explicit **Authentic/Substitution/Adulterant/Contaminant (ASAC) classification framework**, which requires a known/expected ingredient list to compare findings against — something your README already flagged as "only relevant if there's a known/expected ingredient list," which you don't currently have for your 18 samples. Why the difference: Wuhu San is a single, well-defined patent medicine with a known formula; your samples appear to be more exploratory/unlabeled, so this framework isn't directly applicable until (if) you get ingredient-list ground truth.

## Urumarudappa et al. 2022 (Thai Pharmacopoeia)

**Fundamentally different approach**, not just different parameters. This is **not** a shotgun/mixture-deconvolution study — it's a **reference-library-building study**: 101 single, authenticated plant species were PCR-amplified and Sanger-sequenced across ITS2, matK, rbcL, and trnH-psbA to build a curated barcode database, then twenty single-herb commercial products were tested by direct amplification, Sanger sequencing, and BLAST against that reference library — no shotgun sequencing, no assembly, no mixture deconvolution at all. Why: their goal was building **regulatory reference infrastructure** for single-ingredient products, not resolving what's inside an already-mixed sample the way yours (and the other three) are. This is the "clean assumption" case your own literature review already noted sets realistic expectations for what markers can resolve species-wise even under ideal, non-mixed, Sanger-quality conditions — worth keeping as your resolution benchmark, not as a pipeline template.

## What this means for your pipeline specifically

Two concrete, actionable gaps stand out against this comparison set:

1. **Reference curation (Wei's approach)** — if you have or can obtain authenticated *Hippophae*/*Rhaponticum* reference specimens, Sanger-sequencing them directly would resolve your 86-87% identity ceiling far more definitively than deeper shotgun sequencing would.
2. **Dual assembler (Liu's approach)** — adding MEGAHIT alongside your split metaSPAdes/SPAdes runs, exactly as your README's next-steps list already planned, is validated by two of these four papers as standard practice for this problem, not just a nice-to-have.

The ASAC framework and multi-kingdom marker panel (Zhang's COI-for-animal, psbA-trnH-for-plants breadth) are both good ideas but contingent on things you don't have yet (a known ingredient list; evidence you need non-plant kingdom coverage) — I'd treat those as "keep in mind" rather than immediate next steps.
