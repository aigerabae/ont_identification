# Limitations (read before interpreting results)

## What a read fraction means
* Read fractions are *sequence* fractions, not biomass: rDNA copy number differs between species and with PCR
  amplification of long targets.
* ITS rows and 18S/28S rows count different parts of the cistron. Their proportions are not directly
  comparable: e.g. in samples without ITS contigs, Crassulaceae (18S/28S) exceeds Elaeagnaceae, whereas ITS
  in the other samples shows *Hippophae* far above *Rhodiola*.
* The same plant can appear at two ranks in one sample: e.g. *Hippophae* (genus, from ITS) and Elaeagnaceae
  (family, from ITS reads that could not be resolved further and from 18S/28S). They are not two different
  plants; the family row holds the reads that cannot be resolved to genus.

## Resolution
* ITS often cannot separate close congeners; the pipeline then reports the genus on purpose.
* A species absent from GenBank cannot be named; identity caps move such calls up to genus or family.
* 18S/28S calls are never below tribe (calibration: genus calls were wrong in 4.5-5.3 % of cases when the true
  genus was missing). Elaeagnaceae and Crassulaceae have no tribes in the NCBI taxonomy, so their 18S/28S calls
  stay at family level (*Hippophae* / *Rhodiola* cannot be named from 18S/28S).
* Expected error of 18S/28S calls (leave-one-out): ~1.5 % when the genus is in the database, ~4 % when not;
  higher for Cannabaceae (19 % wrong family for 1000-bp queries when the genus is missing) and Rosaceae (7.5 %).
  Calls of unexpected families at ≤0.2 % of a sample (e.g. Viscaceae, Chloranthaceae) are within this error.
* Read-level calls (150-bp reads) are less resolved than contig calls; conserved 28S reads that tie across
  tribes or families get no call (e.g. ONplants-17: a 28S contig matched Cardueae in a manual NCBI BLAST, but
  most of its reads tie above family and stay unassigned).

## Coverage
* Only the 50 000 most abundant residual read sequences per sample are classified. In 8 samples
  (ONplants-2, -5, -7, -8, -11, -13, -14, -17) the rarer sequences make up 12-26 % of all reads; they are counted
  as Unclassified. Raising `residual.max_unique_sequences` would classify them (longer run).
* The ITS database contains only the seed-list taxa (chosen by the analyst); plants outside it cannot get an
  ITS call. The 18S/28S database covers all flowering-plant families, at family/tribe resolution.
* Contigs on which ITSx found only 18S/28S (or mislabelled 28S as a partial ITS2) are not ITS contigs; their
  reads are classified read by read, with 18S/28S.
* Samples 18S-1, -3, -9, -31 have <1000 raw reads, 78-94 % human: no plant composition can be derived from them.
  18S-13 and 18S-32 give <4 % plant reads.

## Other
* Remaining fungal ITS read-level calls (~5 000 reads, mostly Pleosporaceae/*Alternaria*) are at 90-93 % identity
  and tentative.
* The previous version's "mapped to contig, no confident taxon" reads were reads on 18S/28S contigs (that version
  mapped to every contig ITSx flagged); they are not comparable to this version's "Unclassified".
* Validation so far: leave-one-out on the databases and agreement with the previous (core_nt) report for the
  samples without ITS. Validation against samples of known composition has not been done.
