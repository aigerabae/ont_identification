#!/bin/bash
# build_rdna_db.sh
#
# Builds the SECONDARY local reference BLAST database of nuclear rDNA (18S / 25-28S, incl. whole
# 18S-ITS-28S records). It is used only for contigs and reads that carry no ITS; the custom ITS
# database (build_custom_db.sh) stays the primary classifier. 18S/28S resolve family or tribe,
# rarely genus, so coverage across families matters more than depth within a genus:
#   1. backbone: every nuclear rDNA record of 1-10 kb for all flowering plants (Magnoliopsida)
#      and gymnosperms (Acrogymnospermae, for Ginkgo)
#   2. thin target taxa (few long records in GenBank): all records of 300 bp-10 kb
#   3. controls: fungal seed families, a capped fungal backbone, insects, human
#   4. one copy per accession; drop records with >5 % N; identical sequences collapsed within a
#      taxid; at most PER_TAXID_CAP records per taxid (longest first) -- removes the many
#      near-identical RefSeq "PREDICTED" rRNA copies from genome annotations
#   5. makeblastdb with NCBI taxids (so lca_assign.py works unchanged)
#   6. self-test + MANIFEST (date, queries, caps, tool versions, md5, counts)
#
# No ITSx: records are kept whole (ITS + flanks), because the queries are 18S/28S contigs/reads.
# Chloroplast / mitochondrial / plastid rRNA and genome / chromosome / scaffold records are excluded.
#
# Requirements: the conda env of environment.yml (entrez-direct, blast, python)
# Optional but recommended: an NCBI API key (much faster; without one NCBI answered ~1 query/20 s).
#   Either save it once in a private file:   echo 'YOUR_KEY' > ~/.ncbi_api_key && chmod 600 ~/.ncbi_api_key
#   or export it in the shell that starts the build:   export NCBI_API_KEY=YOUR_KEY
#   (NCBI_API_KEY_FILE=... points to a different file)
#
# Usage:
#   db/build_rdna_db.sh [outdir]          (default outdir: db/build_rdna)
#   env: DB_NAME=rdna_custom  PER_TAXID_CAP=10  MIN_LEN=300  MAX_LEN=10000  MAX_N_FRAC=0.05
#        THIN_TAXA="Elaeagnaceae Crassulaceae ..."   (300 bp-10 kb instead of 1-10 kb)
#        CAP_FUNGI=20000  CAP_INSECTA=10000  CAP_GYMNO=5000  CAP_DEFAULT=200000
#   Run in the background, e.g.:
#     nohup db/build_rdna_db.sh > db/build_rdna.log 2>&1 &
#
# Expect: download 1-3 h without an API key (~170k records, ~350 MB); processing minutes.

set -euo pipefail

OUTDIR="${1:-db/build_rdna}"
DB_NAME="${DB_NAME:-rdna_custom}"
PER_TAXID_CAP="${PER_TAXID_CAP:-10}"
MIN_LEN="${MIN_LEN:-300}"
MAX_LEN="${MAX_LEN:-10000}"
MAX_N_FRAC="${MAX_N_FRAC:-0.05}"
THIN_TAXA="${THIN_TAXA:-Elaeagnaceae Crassulaceae Cardueae Anthemideae Hypericum Valeriana Glycyrrhiza Ginkgo Panax}"
CAP_FUNGI="${CAP_FUNGI:-20000}"
CAP_INSECTA="${CAP_INSECTA:-10000}"
CAP_GYMNO="${CAP_GYMNO:-5000}"
CAP_DEFAULT="${CAP_DEFAULT:-200000}"

for tool in esearch efetch xtract makeblastdb blastn python3 md5sum; do
    command -v "$tool" >/dev/null || { echo "ERROR: '$tool' not found in PATH" >&2; exit 1; }
done

# NCBI API key: EDirect reads NCBI_API_KEY from the environment. If it is not set, read it from a
# private file (default ~/.ncbi_api_key), so the key never has to be typed on the command line.
# The key itself is never printed or written to logs / MANIFEST.
NCBI_API_KEY_FILE="${NCBI_API_KEY_FILE:-$HOME/.ncbi_api_key}"
if [[ -z "${NCBI_API_KEY:-}" && -s "$NCBI_API_KEY_FILE" ]]; then
    NCBI_API_KEY="$(tr -d ' \t\r\n' < "$NCBI_API_KEY_FILE")"
    export NCBI_API_KEY
fi
if [[ -n "${NCBI_API_KEY:-}" ]]; then
    echo "NCBI API key: in use (${#NCBI_API_KEY} characters, ending ...${NCBI_API_KEY: -4}); up to 10 requests/s"
    API_KEY_STATUS="yes"
else
    echo "NCBI API key: not set -> slow (about 3 requests/s, often much slower). See the header of this script." >&2
    API_KEY_STATUS="no"
fi

RAW="$OUTDIR/raw"; DBDIR="$OUTDIR/db"
mkdir -p "$RAW" "$DBDIR"
# EDirect prints the full request URL (including api_key=...) when a request fails; those messages go
# to edirect.log. Mask the key there whenever the script exits (also after an error or Ctrl-C).
scrub_logs() { [[ -f "$OUTDIR/edirect.log" ]] && sed -i -E 's/api_key=[A-Za-z0-9]+/api_key=REDACTED/g' "$OUTDIR/edirect.log" || true; }
trap scrub_logs EXIT
: > "$OUTDIR/all_raw.fasta"
: > "$OUTDIR/accession_taxid_group.tsv"
: > "$OUTDIR/queries.tsv"
echo -e "group\tn_available\tn_downloaded" > "$OUTDIR/download_counts.tsv"

RRNA='("18S"[Title] OR "28S"[Title] OR "25S"[Title] OR "26S"[Title] OR "small subunit ribosomal"[Title] OR "large subunit ribosomal"[Title])'
EXCL='NOT (chloroplast[Title] OR mitochondrial[Title] OR mitochondrion[Title] OR plastid[Title] OR "16S"[Title] OR "23S"[Title] OR genome[Title] OR chromosome[Title] OR scaffold[Title] OR contig[Title])'
LONG="1000:${MAX_LEN}[SLEN]"
SHORT="${MIN_LEN}:${MAX_LEN}[SLEN]"

# Same resumable download logic as build_custom_db.sh: a group is complete when >=99 % of the
# requested records arrived; complete groups are never downloaded again, incomplete ones are kept,
# listed in INCOMPLETE_GROUPS.txt and retried on the next run.
: > "$OUTDIR/INCOMPLETE_GROUPS.txt"
fetch_group() {   # label  query  cap
    local label="$1" query="$2" cap="$3"
    local safe; safe=$(echo "$label" | tr -c 'A-Za-z0-9_.-' '_')
    local fa="$RAW/${safe}.fasta" at="$RAW/${safe}.acc_taxid.tsv" sx="$RAW/${safe}.search.xml"
    echo -e "${label}\t${cap}\t${query}" >> "$OUTDIR/queries.tsv"
    if [[ -s "$RAW/${safe}.done" ]]; then
        echo "[skip] $label (already downloaded: $(grep -c '^>' "$fa") sequences)"
    else
        esearch -db nuccore -query "$query" > "$sx" 2>> "$OUTDIR/edirect.log"
        local count; count=$(xtract -input "$sx" -pattern ENTREZ_DIRECT -element Count)
        count=${count:-0}
        local stop=$(( count < cap ? count : cap ))
        if (( stop == 0 )); then
            echo "[warn] $label: 0 records" >&2
            echo -e "${label}\t0\t0" >> "$OUTDIR/download_counts.tsv"
            return 0
        fi
        echo "[fetch] $label: $stop of $count records ($(date +%T))"
        efetch -format fasta -start 1 -stop "$stop" < "$sx" > "$fa" 2>> "$OUTDIR/edirect.log"
        efetch -format docsum -start 1 -stop "$stop" < "$sx" 2>> "$OUTDIR/edirect.log" \
            | xtract -pattern DocumentSummary -element AccessionVersion TaxId > "$at"
        echo -e "${label}\t${count}\t${stop}" >> "$OUTDIR/download_counts.tsv"
        local got nacc; got=$(grep -c '^>' "$fa" || true); nacc=$(wc -l < "$at")
        if (( got * 100 >= stop * 99 && nacc * 100 >= stop * 99 )); then
            echo "[ok]   $label: $got of $stop sequences ($(date +%T))"
            date -u +%FT%TZ > "$RAW/${safe}.done"
        else
            echo "[WARN] $label: only $got sequences / $nacc taxid rows of $stop -> kept; re-run to retry" >&2
            echo -e "${label}\t${stop}\t${got}" >> "$OUTDIR/INCOMPLETE_GROUPS.txt"
        fi
    fi
    cat "$fa" >> "$OUTDIR/all_raw.fasta"
    awk -v g="$label" 'BEGIN{OFS="\t"}{print $1,$2,g}' "$at" >> "$OUTDIR/accession_taxid_group.tsv"
}

echo "== Step 1: backbone (1-10 kb) =="
fetch_group "backbone:Magnoliopsida"   "\"Magnoliopsida\"[Organism] AND $RRNA AND $LONG $EXCL"    "$CAP_DEFAULT"
fetch_group "backbone:Acrogymnospermae" "\"Acrogymnospermae\"[Organism] AND $RRNA AND $LONG $EXCL" "$CAP_GYMNO"

echo "== Step 2: thin target taxa (${MIN_LEN} bp-10 kb) =="
for t in $THIN_TAXA; do
    fetch_group "thin:$t" "\"$t\"[Organism] AND $RRNA AND $SHORT $EXCL" "$CAP_DEFAULT"
done

echo "== Step 3: controls =="
for t in Cladosporiaceae Pleosporaceae Nectriaceae; do
    fetch_group "control:$t" "\"$t\"[Organism] AND $RRNA AND $LONG $EXCL" "$CAP_DEFAULT"
done
fetch_group "control:Fungi"        "\"Fungi\"[Organism] AND $RRNA AND $LONG $EXCL"        "$CAP_FUNGI"
fetch_group "control:Insecta"      "\"Insecta\"[Organism] AND $RRNA AND $LONG $EXCL"      "$CAP_INSECTA"
fetch_group "control:Homo_sapiens" "\"Homo sapiens\"[Organism] AND $RRNA AND $LONG $EXCL" "$CAP_DEFAULT"

if [[ -s "$OUTDIR/INCOMPLETE_GROUPS.txt" ]]; then
    echo "!! incomplete groups (group, requested, received) -- re-run this script later to retry:"
    cat "$OUTDIR/INCOMPLETE_GROUPS.txt"
fi

echo "== Step 4: reference FASTA =="
python3 - "$OUTDIR" "$MIN_LEN" "$MAX_LEN" "$MAX_N_FRAC" "$PER_TAXID_CAP" <<'PY'
import sys, collections, hashlib
outd, min_len, max_len, max_n, cap = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
tax, grp = {}, {}
for line in open(f"{outd}/accession_taxid_group.tsv"):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 2 and p[1].strip().isdigit():
        tax.setdefault(p[0], p[1].strip()); grp.setdefault(p[0], p[2] if len(p) > 2 else "")
seqs, order, cur = {}, [], None
for line in open(f"{outd}/all_raw.fasta"):
    if line.startswith(">"):
        acc = line[1:].split()[0]
        cur = None if acc in seqs else acc        # one copy per accession (groups overlap)
        if cur: seqs[cur] = []; order.append(cur)
    elif cur:
        seqs[cur].append(line.strip().upper())
stats = collections.Counter()
by_taxid = collections.defaultdict(list)
for acc in order:
    s = "".join(seqs[acc]); n = len(s)
    if acc not in tax:              stats["no_taxid"] += 1; continue
    if not (min_len <= n <= max_len): stats["length"] += 1; continue
    if s.count("N") > max_n * n:    stats["too_many_N"] += 1; continue
    by_taxid[tax[acc]].append((n, acc, s))
kept = []
for t, recs in by_taxid.items():
    recs.sort(key=lambda r: (-r[0], r[1]))        # longest first, accession breaks ties (deterministic)
    seen, chosen = set(), []
    for n, acc, s in recs:
        h = hashlib.md5(s.encode()).hexdigest()
        if h in seen: stats["identical_within_taxid"] += 1; continue
        if len(chosen) >= cap: stats["over_taxid_cap"] += 1; continue
        seen.add(h); chosen.append((acc, s))
    kept.extend(chosen)
kept.sort(key=lambda r: r[0])                     # sorted by accession: stable output, stable md5
with open(f"{outd}/refs.fasta", "w") as fa, open(f"{outd}/refs.taxid.map", "w") as fm, \
     open(f"{outd}/refs_manifest.tsv", "w") as mf:
    mf.write("accession\ttaxid\tlength_bp\tgroup\n")
    for acc, s in kept:
        fa.write(f">{acc}\n")
        for i in range(0, len(s), 80): fa.write(s[i:i+80] + "\n")
        fm.write(f"{acc}\t{tax[acc]}\n")
        mf.write(f"{acc}\t{tax[acc]}\t{len(s)}\t{grp.get(acc, '')}\n")
print(f"   {len(order)} unique accessions -> {len(kept)} references in {len(by_taxid)} taxids; removed: {dict(stats)}")
open(f"{outd}/filter_stats.txt", "w").write(
    "\n".join(f"{k}\t{v}" for k, v in [("unique_accessions", len(order)), ("references", len(kept)),
                                        ("taxids", len(by_taxid))] + sorted(stats.items())) + "\n")
PY

echo "== Step 5: makeblastdb =="
makeblastdb -in "$OUTDIR/refs.fasta" -dbtype nucl -parse_seqids \
    -taxid_map "$OUTDIR/refs.taxid.map" \
    -title "Custom nuclear rDNA reference ($(date +%F))" -out "$DBDIR/$DB_NAME"

echo "== Step 6: self-test =="
# ~50 evenly spaced references must find a >=99.9 % hit with a taxid (usually themselves)
python3 - "$OUTDIR/refs.fasta" "$OUTDIR/selftest_query.fasta" <<'PY'
import sys
recs = []
for l in open(sys.argv[1]):
    if l.startswith(">"): recs.append([l, ""])
    else: recs[-1][1] += l
step = max(len(recs) // 50, 1)
open(sys.argv[2], "w").write("".join(h + s for h, s in recs[::step][:50]))
PY
blastn -task megablast -db "$DBDIR/$DB_NAME" -query "$OUTDIR/selftest_query.fasta" -dust no \
    -outfmt "6 qseqid sacc staxids pident qcovs" -max_target_seqs 50 > "$OUTDIR/selftest_result.tsv"
python3 - "$OUTDIR/selftest_query.fasta" "$OUTDIR/selftest_result.tsv" <<'PY' || { echo "SELF-TEST FAILED: check taxid map / makeblastdb" >&2; exit 1; }
import sys
q = [l[1:].split()[0] for l in open(sys.argv[1]) if l.startswith(">")]
ok = set()
for l in open(sys.argv[2]):
    f = l.rstrip("\n").split("\t")
    if float(f[3]) >= 99.9 and float(f[4]) >= 99 and f[2].strip() not in ("0", "N/A", ""):
        ok.add(f[0])
print(f"self-test: {len(ok)}/{len(q)} references have a >=99.9% full-length hit with a taxid")
for x in q:
    if x not in ok: print(f"   no good hit: {x}")
sys.exit(0 if len(ok) >= 0.95 * len(q) else 1)
PY
echo "self-test OK"

echo "== Manifest =="
{
    echo "build_date: $(date -u +%FT%TZ)"
    echo "purpose: secondary database for contigs/reads without ITS (18S/28S); primary = custom ITS database"
    echo "length: ${MIN_LEN}-${MAX_LEN} bp (thin taxa), 1000-${MAX_LEN} bp (all other groups)  max_N_fraction: $MAX_N_FRAC"
    echo "per_taxid_cap: $PER_TAXID_CAP (longest first, identical sequences collapsed within taxid)"
    echo "thin_taxa: $THIN_TAXA"
    echo "caps: fungi=$CAP_FUNGI insecta=$CAP_INSECTA gymnosperms=$CAP_GYMNO other=$CAP_DEFAULT"
    echo "queries: see queries.tsv; downloaded counts: see download_counts.tsv"
    echo "blast: $(blastn -version | head -1)"
    echo "edirect: $(esearch -version 2>&1 | head -1)  ncbi_api_key_used: $API_KEY_STATUS"
    echo "python: $(python3 --version)"
    echo "n_reference_sequences: $(grep -c '>' "$OUTDIR/refs.fasta")"
    echo "md5_refs_fasta: $(md5sum "$OUTDIR/refs.fasta" | cut -d' ' -f1)"
    echo "md5_taxid_map: $(md5sum "$OUTDIR/refs.taxid.map" | cut -d' ' -f1)"
} > "$OUTDIR/MANIFEST.txt"
cat "$OUTDIR/MANIFEST.txt"
echo
echo "Done. DB prefix: $DBDIR/$DB_NAME"
echo "Archive for reproducibility: refs.fasta, refs.taxid.map, refs_manifest.tsv, queries.tsv, download_counts.tsv, MANIFEST.txt"
