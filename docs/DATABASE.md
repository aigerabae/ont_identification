# Local reference databases

Two local BLAST databases, both built from GenBank by scripts in `db/`, frozen and checksummed.
The ITS database is the primary classifier; the 18S/28S database is used only where ITS gives no call.

| | ITS database (primary) | 18S/28S database (secondary) |
|---|---|---|
| Build script | `db/build_custom_db.sh` | `db/build_rdna_db.sh` |
| Location used | `db/build_v2/db/its_custom` | `db/build_rdna/db/rdna_custom` |
| References | 108 834 | 51 526 (29 628 taxids) |
| md5 `refs.fasta` | 0346d6202a0d84e27ed59632c5b00ef7 | 2b34be0085ffb515f2a43f0764a8153d |
| md5 `refs.taxid.map` | 433690d73f1f3cdb22711e6aeac462cf | a307d7ba1a5cc2c8e7c8c36dcd2c389b |
| Built | 2026-10-28 (download), filtered 2026-10-30 | 2026-10-30 |
| Validation | `db/validate_db.py` (leave-one-out) | `db/validate_rdna_db.py` (leave-one-out calibration) |

## ITS database
* **Seed list** `db/taxa_seed.txt`: families, tribes and genera expanded to *all* GenBank ITS records
  (Elaeagnaceae, Crassulaceae, Cardueae, Cichorieae, Anthemideae, *Lactuca*, Rosaceae, Lamiaceae, Apiaceae,
  Ranunculaceae, Fabaceae (≤20 000 records), Cannabaceae, Zingiberaceae, *Hypericum*, *Ginkgo*, *Panax*,
  *Glycyrrhiza*, *Valeriana*; fungi Cladosporiaceae, Pleosporaceae, Nectriaceae). The list was chosen by the
  analyst from the expected components and the previous results. Extra single accessions: `db/extra_accessions.txt`.
* **Decoys** `db/decoys.txt`: human rDNA, insect ITS, general fungal ITS, so off-target reads get their real
  identity instead of being forced onto a plant.
* **Processing**: records trimmed to the ITS region with ITSx; one region per accession
  (full > ITS1+ITS2 > ITS1 > ITS2); exact duplicates removed only *within* a taxid, so identical sequences
  shared by different species remain and produce a genus-level call.
* **Length filter** (`db/filter_long_refs.py`, step 4b of the build): records whose length is impossible for
  the region they claim are dropped (single spacer > 450 bp, ITS1-5.8S-ITS2 > 1000 bp). Such records still carry
  18S/28S flanks; in the first build three of them (AF250821, AF090879, AF090882) attracted ~75 000 plant
  18S/28S reads and reported them as fungi. Animal decoys (Metazoa, 33208) and gymnosperms
  (Acrogymnospermae, 1437180; *Ginkgo*) are exempt, their ITS is genuinely long. 68 of 108 902 records
  dropped (36 fungal, 32 flowering-plant); every one is listed in `dropped_records.tsv`.
* **Validation**: `db/validate_db.py`, leave-one-out with the species hidden. Note: it counts a wrong genus
  inside the right family as `correct_family`, which flatters the numbers; `validate_rdna_db.py` counts it as an error.

## 18S/28S database
* **Content** (`queries.tsv`, `download_counts.tsv`): nuclear 18S/25-28S records of 1-10 kb for all flowering
  plants (Magnoliopsida, 117 854 records) and gymnosperms (≤5 000); 300 bp-10 kb for thinly covered target taxa
  (Elaeagnaceae, Crassulaceae, Cardueae, Anthemideae, *Hypericum*, *Valeriana*, *Glycyrrhiza*, *Ginkgo*, *Panax*);
  controls: fungal seed families, a fungal backbone (≤20 000), insects (≤10 000), human. Chloroplast,
  mitochondrial and genome/chromosome/scaffold records excluded. Coverage across *families* matters more than
  depth: 18S/28S resolve family or tribe, and a missing family would attract reads to its nearest neighbour.
* **Processing**: no ITSx (records kept whole); one copy per accession; records with >5 % N dropped;
  identical sequences collapsed within a taxid; at most 10 records per taxid (longest first), which removes
  the many near-identical "PREDICTED" rRNA copies from genome annotations. 171 592 downloaded accessions →
  51 526 references (removed: 110 397 over the per-taxid cap, 8 121 identical, 1 040 N-rich, 508 without taxid).
* **Known gaps**: *Ginkgo* short records 28 of 29 (one record not returned by NCBI). Some genera have no
  18S/28S in GenBank (e.g. *Leuzea* 28S); their reads are called at tribe/family level of the nearest relatives.
* **Calibration**: see METHODS section 7. Results in `db/build_rdna/validation_tribe/`.

## Reproducibility
* The frozen `refs.fasta` + `refs.taxid.map` of each database rebuild it exactly, offline
  (`makeblastdb -parse_seqids -taxid_map`); MANIFEST.txt records date, settings, tool versions, counts and md5.
* A fresh download will not be identical (GenBank grows), but follows the same rules: the length filter is part
  of `build_custom_db.sh`; caps and filters are part of `build_rdna_db.sh`.
