# legacy/ -- previous version (NOT part of the new method)

`old_scripts/` holds the scripts of the earlier workflow, kept only for traceability.
That workflow classified by **manual NCBI web BLAST (core_nt) of ITS regions**, top-hit based,
with read counts from merged reads only. The supervisor asked for a **local custom reference
database as the primary method**; the new pipeline (`workflow/`, `scripts/`, `db/`) replaces it.
Do not mix outputs of the two. See `docs/CHANGES_FROM_PREVIOUS_VERSION.md`.
