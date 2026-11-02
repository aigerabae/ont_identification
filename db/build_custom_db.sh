#!/bin/bash
# build_custom_db.sh
#
# Builds the local, reproducible ITS reference BLAST database:
#   1. downloads every ITS record for each taxon in taxa_seed.txt (whole
#      families/genera -- NOT only the species already seen in core_nt)
#   2. downloads the decoy (non-target) records in decoys.txt
#   3. trims every record to its ITS region with ITSx (so 18S/28S flanks of
#      long GenBank records cannot dominate alignment scores)
#   4. keeps one ITS reference per accession (full > ITS1+ITS2 > ITS1 > ITS2),
#      de-duplicates only within the same taxid
#   4b. length filter (filter_long_refs.py): drops records whose length is impossible for the
#      ITS region they claim -- single spacer (ITS1 or ITS2) > 450 bp, ITS1-5.8S-ITS2 > 1000 bp.
#      Such records still carry 18S/28S flanks and make plant 18S/28S reads look like fungi.
#      Animal decoys (Metazoa) and gymnosperms (Ginkgo) are exempt: their ITS is genuinely long.
#      Every dropped record is listed in dropped_records.tsv.
#   5. makeblastdb with NCBI taxids (so lca_assign.py works unchanged)
#   6. self-test + MANIFEST (date, tool versions, md5, counts) for the methods
#
# Requirements: the conda env of environment.yml (entrez-direct, blast, itsx, python)
#               NCBI taxonomy dump (nodes.dmp, names.dmp) in TAXDUMP_DIR -- db/get_taxdump.sh
# Optional but recommended: an NCBI API key (faster downloads).
#   Either save it once in a private file:   echo 'YOUR_KEY' > ~/.ncbi_api_key && chmod 600 ~/.ncbi_api_key
#   or export it in the shell that starts the build:   export NCBI_API_KEY=YOUR_KEY
#
# Usage:
#   db/build_custom_db.sh [taxa_file] [decoys_file] [outdir]
#   env: MAX_PER_TAXON=20000  THREADS=16  DB_NAME=its_custom  MIN_LEN=100
#        LENGTH_FILTER=1  FILTER_MAX_SPACER=450  FILTER_MAX_ITS=1000
#        FILTER_EXEMPT_TAXIDS=33208,1437180  TAXDUMP_DIR=db/taxdump
#
# Expect: download = tens of minutes; ITSx on ~100k records = hours (use THREADS).

set -euo pipefail

TAXA_FILE="${1:-db/taxa_seed.txt}"
DECOY_FILE="${2:-db/decoys.txt}"
OUTDIR="${3:-db/build}"
EXTRA_FILE="${EXTRA_FILE:-db/extra_accessions.txt}"
MAX_PER_TAXON="${MAX_PER_TAXON:-20000}"
THREADS="${THREADS:-16}"
DB_NAME="${DB_NAME:-its_custom}"
MIN_LEN="${MIN_LEN:-100}"
LENGTH_FILTER="${LENGTH_FILTER:-1}"                       # 1 = apply step 4b (default), 0 = skip
FILTER_MAX_SPACER="${FILTER_MAX_SPACER:-450}"             # max bp for an ITS1-only or ITS2-only record
FILTER_MAX_ITS="${FILTER_MAX_ITS:-1000}"                  # max bp for an ITS1-5.8S-ITS2 record
FILTER_EXEMPT_TAXIDS="${FILTER_EXEMPT_TAXIDS:-33208,1437180}"   # Metazoa, Acrogymnospermae: never dropped
TAXDUMP_DIR="${TAXDUMP_DIR:-db/taxdump}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # = db/

for tool in esearch efetch xtract ITSx makeblastdb blastn python3; do
    command -v "$tool" >/dev/null || { echo "ERROR: '$tool' not found in PATH" >&2; exit 1; }
done
[[ -f "$TAXA_FILE" ]] || { echo "ERROR: taxa file not found: $TAXA_FILE" >&2; exit 1; }

# NCBI API key: EDirect reads NCBI_API_KEY from the environment. If it is not set, read it from a
# private file (default ~/.ncbi_api_key). The key itself is never printed or written to logs / MANIFEST.
NCBI_API_KEY_FILE="${NCBI_API_KEY_FILE:-$HOME/.ncbi_api_key}"
if [[ -z "${NCBI_API_KEY:-}" && -s "$NCBI_API_KEY_FILE" ]]; then
    NCBI_API_KEY="$(tr -d ' \t\r\n' < "$NCBI_API_KEY_FILE")"
    export NCBI_API_KEY
fi
if [[ -n "${NCBI_API_KEY:-}" ]]; then
    echo "NCBI API key: in use (${#NCBI_API_KEY} characters, ending ...${NCBI_API_KEY: -4}); up to 10 requests/s"
else
    echo "NCBI API key: not set -> slow downloads. Save it in ~/.ncbi_api_key (chmod 600) or export NCBI_API_KEY." >&2
fi
# Checked now, not at step 4b, so a missing taxonomy does not fail the build after hours of ITSx.
if [[ "$LENGTH_FILTER" == "1" ]]; then
    [[ -f "$TAXDUMP_DIR/nodes.dmp" && -f "$TAXDUMP_DIR/names.dmp" ]] || {
        echo "ERROR: NCBI taxonomy not found in $TAXDUMP_DIR (needed by the length filter)." >&2
        echo "       Run db/get_taxdump.sh, set TAXDUMP_DIR=..., or LENGTH_FILTER=0 to skip the filter." >&2
        exit 1; }
    [[ -f "$SCRIPT_DIR/filter_long_refs.py" ]] || {
        echo "ERROR: $SCRIPT_DIR/filter_long_refs.py not found" >&2; exit 1; }
fi

RAW="$OUTDIR/raw"; DBDIR="$OUTDIR/db"
mkdir -p "$RAW" "$DBDIR"
# EDirect prints the full request URL (including api_key=...) when a request fails; those messages go
# to edirect.log. Mask the key there whenever the script exits (also after an error or Ctrl-C).
scrub_logs() { [[ -f "$OUTDIR/edirect.log" ]] && sed -i -E 's/api_key=[A-Za-z0-9]+/api_key=REDACTED/g' "$OUTDIR/edirect.log" || true; }
trap scrub_logs EXIT
: > "$OUTDIR/all_raw.fasta"
: > "$OUTDIR/accession_taxid_group.tsv"
: > "$OUTDIR/download_counts.tsv"
echo -e "group\tn_available\tn_downloaded" >> "$OUTDIR/download_counts.tsv"

ITS_TERMS='(internal transcribed spacer[Title] OR ITS1[Title] OR ITS2[Title] OR "5.8S"[Title])'
EXCL='NOT (genome[Title] OR chromosome[Title] OR scaffold[Title] OR contig[Title])'

# A group counts as complete when >=99 % of the requested records arrived (NCBI sometimes cannot
# return 1-3 records out of thousands: withdrawn/suppressed entries). Complete groups are never
# downloaded again; incomplete ones are kept (so the database is not silently missing a group),
# listed in INCOMPLETE_GROUPS.txt, and retried the next time the script runs.
: > "$OUTDIR/INCOMPLETE_GROUPS.txt"
fetch_group() {   # label  query  cap
    local label="$1" query="$2" cap="$3"
    local safe; safe=$(echo "$label" | tr -c 'A-Za-z0-9_.-' '_')
    local fa="$RAW/${safe}.fasta" at="$RAW/${safe}.acc_taxid.tsv" sx="$RAW/${safe}.search.xml"
    local reuse=0
    if [[ -s "$RAW/${safe}.done" ]]; then
        reuse=1
    elif [[ -s "$fa" && -s "$at" ]]; then      # files from an earlier run without a .done marker
        local g0 n0; g0=$(grep -c '^>' "$fa" || true); n0=$(wc -l < "$at")
        if (( n0 > 0 && g0 * 100 >= n0 * 99 )); then reuse=1; touch "$RAW/${safe}.done"; fi
    fi
    if (( reuse )); then
        echo "[skip] $label (already downloaded: $(grep -c '^>' "$fa") sequences)"
    else
        esearch -db nuccore -query "$query" > "$sx" 2>> "$OUTDIR/edirect.log"
        local count; count=$(xtract -input "$sx" -pattern ENTREZ_DIRECT -element Count)
        count=${count:-0}
        local stop=$(( count < cap ? count : cap ))
        if (( stop == 0 )); then
            echo "[warn] $label: 0 records for query: $query" >&2
            echo -e "${label}\t0\t0" >> "$OUTDIR/download_counts.tsv"
            return 0
        fi
        echo "[fetch] $label: $stop of $count records"
        efetch -format fasta -start 1 -stop "$stop" < "$sx" > "$fa" 2>> "$OUTDIR/edirect.log"
        efetch -format docsum -start 1 -stop "$stop" < "$sx" 2>> "$OUTDIR/edirect.log" \
            | xtract -pattern DocumentSummary -element AccessionVersion TaxId > "$at"
        echo -e "${label}\t${count}\t${stop}" >> "$OUTDIR/download_counts.tsv"
        local got nacc; got=$(grep -c '^>' "$fa" || true); nacc=$(wc -l < "$at")
        if (( got * 100 >= stop * 99 && nacc * 100 >= stop * 99 )); then
            echo "[ok]   $label: $got of $stop sequences (retries are logged in edirect.log)"
            touch "$RAW/${safe}.done"
        else
            echo "[WARN] $label: only $got sequences / $nacc taxid rows of $stop requested -> kept, but re-run the script to retry" >&2
            echo -e "${label}\t${stop}\t${got}" >> "$OUTDIR/INCOMPLETE_GROUPS.txt"
        fi
    fi
    cat "$fa" >> "$OUTDIR/all_raw.fasta"
    awk -v g="$label" 'BEGIN{OFS="\t"}{print $1,$2,g}' "$at" >> "$OUTDIR/accession_taxid_group.tsv"
}

echo "== Step 1: target taxa =="
while IFS= read -r -u 3 taxon; do
    [[ -z "${taxon// }" || "$taxon" =~ ^# ]] && continue
    fetch_group "target:${taxon}" "\"${taxon}\"[Organism] AND ${ITS_TERMS} AND 200:10000[SLEN] ${EXCL}" "$MAX_PER_TAXON"
done 3< "$TAXA_FILE"

echo "== Step 2: decoys =="
if [[ -f "$DECOY_FILE" ]]; then
    while IFS=$'\t' read -r -u 3 label query cap; do
        [[ -z "${label// }" || "$label" =~ ^# ]] && continue
        fetch_group "decoy:${label}" "$query" "${cap:-1000}"
    done 3< "$DECOY_FILE"
fi

echo "== Step 2b: extra accessions =="
if [[ -f "$EXTRA_FILE" ]]; then
    while IFS=$'\t' read -r -u 3 acc label; do
        [[ -z "${acc// }" || "$acc" =~ ^# ]] && continue
        fetch_group "extra:${acc}" "${acc}[Accession]" 1
    done 3< "$EXTRA_FILE"
fi

if [[ -s "$OUTDIR/INCOMPLETE_GROUPS.txt" ]]; then
    echo "!! incomplete groups (group, requested, received) -- re-run this script later to retry them:"
    cat "$OUTDIR/INCOMPLETE_GROUPS.txt"
fi
echo "== Step 3: ITSx trimming =="
# The same accession can be returned by overlapping groups (e.g. Lactuca inside Cichorieae):
# keep one copy per accession (ITSx warns about, and may overwrite, identical identifiers).
awk '/^>/{split(substr($0,2),a," "); keep=!(a[1] in seen); seen[a[1]]=1} keep' \
    "$OUTDIR/all_raw.fasta" > "$OUTDIR/all_raw.unique.fasta"
echo "   $(grep -c '>' "$OUTDIR/all_raw.fasta") downloaded records -> $(grep -c '>' "$OUTDIR/all_raw.unique.fasta") unique accessions"
RAW_MD5=$(md5sum "$OUTDIR/all_raw.unique.fasta" | cut -d' ' -f1)
if [[ "${REUSE_ITSX:-0}" == "1" && -s "$OUTDIR/itsx_out.positions.txt" ]]; then
    echo "[skip] REUSE_ITSX=1: keeping existing ITSx output"
elif [[ -f "$OUTDIR/itsx_out.input.md5" && "$(cat "$OUTDIR/itsx_out.input.md5")" == "$RAW_MD5" ]]; then
    echo "[skip] ITSx output already present for this exact input"
else
    rm -f "$OUTDIR"/itsx_out.*
    echo "   running ITSx (hours for >100k records); detailed progress: $OUTDIR/itsx.log"
    ITSx -i "$OUTDIR/all_raw.unique.fasta" -o "$OUTDIR/itsx_out" --cpu "$THREADS" \
         --preserve T --graphical F -t all > "$OUTDIR/itsx.log" 2>&1
    echo "$RAW_MD5" > "$OUTDIR/itsx_out.input.md5"
fi

echo "== Step 4: reference FASTA =="
python3 "$SCRIPT_DIR/build_ref_fasta.py" \
    --itsx-prefix "$OUTDIR/itsx_out" \
    --accession-map "$OUTDIR/accession_taxid_group.tsv" \
    --out-fasta "$OUTDIR/refs.fasta" \
    --out-taxid-map "$OUTDIR/refs.taxid.map" \
    --out-manifest "$OUTDIR/refs_manifest.tsv" \
    --min-len "$MIN_LEN"

echo "== Step 4b: length filter =="
N_BEFORE_FILTER=$(grep -c '>' "$OUTDIR/refs.fasta")
if [[ "$LENGTH_FILTER" == "1" ]]; then
    mv "$OUTDIR/refs.fasta" "$OUTDIR/refs.unfiltered.fasta"
    mv "$OUTDIR/refs.taxid.map" "$OUTDIR/refs.unfiltered.taxid.map"
    python3 "$SCRIPT_DIR/filter_long_refs.py" \
        --in-fasta "$OUTDIR/refs.unfiltered.fasta" \
        --in-map "$OUTDIR/refs.unfiltered.taxid.map" \
        --taxdump-dir "$TAXDUMP_DIR" \
        --out-fasta "$OUTDIR/refs.fasta" \
        --out-map "$OUTDIR/refs.taxid.map" \
        --out-dropped "$OUTDIR/dropped_records.tsv" \
        --max-spacer "$FILTER_MAX_SPACER" --max-its "$FILTER_MAX_ITS" \
        --exempt-taxids "$FILTER_EXEMPT_TAXIDS"
    N_DROPPED=$(( $(wc -l < "$OUTDIR/dropped_records.tsv") - 1 ))
    echo "   $N_BEFORE_FILTER references -> $(grep -c '>' "$OUTDIR/refs.fasta") ($N_DROPPED dropped, see dropped_records.tsv)"
    echo "   note: refs_manifest.tsv still lists all $N_BEFORE_FILTER pre-filter references"
else
    N_DROPPED=0
    echo "[skip] LENGTH_FILTER=0: all $N_BEFORE_FILTER references kept"
fi

echo "== Step 5: makeblastdb =="
makeblastdb -in "$OUTDIR/refs.fasta" -dbtype nucl -parse_seqids \
    -taxid_map "$OUTDIR/refs.taxid.map" \
    -title "Custom ITS reference ($(date +%F))" -out "$DBDIR/$DB_NAME"

echo "== Step 6: self-test =="
# 20 random references must find themselves at >=99.9 % identity, carrying a taxid
python3 - "$OUTDIR/refs.fasta" "$OUTDIR/selftest_query.fasta" <<'PY'
import sys
recs=[];h=None
for l in open(sys.argv[1]):
    if l.startswith(">"): recs.append([l,""]); 
    else: recs[-1][1]+=l
step=max(len(recs)//50,1)
open(sys.argv[2],"w").write("".join(h+s for h,s in recs[::step][:50]))
PY
blastn -task blastn -db "$DBDIR/$DB_NAME" -query "$OUTDIR/selftest_query.fasta" -dust no \
    -outfmt "6 qseqid sacc staxids pident qcovs" -max_target_seqs 500 > "$OUTDIR/selftest_result.tsv"
python3 - "$OUTDIR/selftest_query.fasta" "$OUTDIR/selftest_result.tsv" <<'PY' || { echo "SELF-TEST FAILED: check taxid map / makeblastdb" >&2; exit 1; }
import sys
q=[l[1:].split()[0] for l in open(sys.argv[1]) if l.startswith(">")]
own, any_ok, nhits = set(), set(), {}
for l in open(sys.argv[2]):
    f=l.rstrip("\n").split("\t")
    nhits[f[0]]=nhits.get(f[0],0)+1
    good = float(f[3])>=99.9 and float(f[4])>=99 and f[2].strip() not in ("0","N/A","")
    if good: any_ok.add(f[0])
    if good and f[0]==f[1]: own.add(f[0])
print(f"self-test: {len(any_ok)}/{len(q)} references have a >=99.9% hit with a taxid; "
      f"{len(own)}/{len(q)} find their own accession")
# A reference that is shared by >500 identical records can lose its own accession to BLAST's
# arbitrary tie truncation; that is harmless. Failing = taxid/db broken = nearly all queries lack a good hit.
for x in q:
    if x not in any_ok: print(f"   no good hit: {x} (hits: {nhits.get(x,0)})")
sys.exit(0 if (len(any_ok) >= 0.95*len(q) and len(own) >= 0.9*len(q)) else 1)
PY
echo "self-test OK"

echo "== Manifest =="
{
    echo "build_date: $(date -u +%FT%TZ)"
    echo "taxa_file: $TAXA_FILE"; echo "decoy_file: $DECOY_FILE"
    echo "max_per_taxon: $MAX_PER_TAXON  min_len: $MIN_LEN"
    if [[ "$LENGTH_FILTER" == "1" ]]; then
        echo "length_filter: on  max_spacer: $FILTER_MAX_SPACER  max_its: $FILTER_MAX_ITS  exempt_taxids: $FILTER_EXEMPT_TAXIDS  taxdump: $TAXDUMP_DIR"
    else
        echo "length_filter: off"
    fi
    echo "n_references_before_filter: $N_BEFORE_FILTER  n_dropped_by_filter: $N_DROPPED"
    echo "blast: $(blastn -version | head -1)"
    echo "itsx: $(ITSx --version 2>&1 | head -2 | tr '\n' ' ')"
    echo "edirect: $(esearch -version 2>&1 | head -1)"
    echo "python: $(python3 --version)"
    echo "n_reference_sequences: $(grep -c '>' "$OUTDIR/refs.fasta")"
    echo "md5_refs_fasta: $(md5sum "$OUTDIR/refs.fasta" | cut -d' ' -f1)"
    echo "md5_taxid_map: $(md5sum "$OUTDIR/refs.taxid.map" | cut -d' ' -f1)"
} > "$OUTDIR/MANIFEST.txt"
cat "$OUTDIR/MANIFEST.txt"
echo
echo "Done. DB prefix: $DBDIR/$DB_NAME"
echo "Archive for reproducibility: refs.fasta, refs.taxid.map, refs_manifest.tsv, dropped_records.tsv, MANIFEST.txt, taxa_seed.txt, decoys.txt"
echo "(anyone can rebuild the BLAST db from refs.fasta + refs.taxid.map with the makeblastdb command above, without re-downloading)"
