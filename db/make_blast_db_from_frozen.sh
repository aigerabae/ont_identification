#!/bin/bash
# Rebuild the BLAST database from the frozen files (no downloads, fully reproducible):
#   refs.fasta + refs.taxid.map   (produced by build_custom_db.sh, archived with the project)
# Usage: db/make_blast_db_from_frozen.sh <refs.fasta> <refs.taxid.map> <outprefix> [title]
set -euo pipefail
FA="$1"; MAP="$2"; PREFIX="$3"; TITLE="${4:-Custom ITS reference}"
mkdir -p "$(dirname "$PREFIX")"
makeblastdb -in "$FA" -dbtype nucl -parse_seqids -taxid_map "$MAP" -title "$TITLE" -out "$PREFIX" > /dev/null
blastdbcmd -db "$PREFIX" -info > "$PREFIX.info.txt"; sed -n 1,3p "$PREFIX.info.txt"
