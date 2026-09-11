#!/usr/bin/env python3
"""
parse_ncbi_text_results.py

Parses NCBI web BLAST results downloaded as plain text ("Download All" ->
"Text" from the results page) and computes exact, read-weighted per-species
counts and percentages per sample.

Expects query IDs to end in '_count{N}' or '__count{N}' (matches both
export_unique_reads_fasta.py's and diagnose_failed_samples.py's ID formats),
so each query's real read count can be recovered directly -- no separate
lookup file needed. The sample name for each input file is given explicitly
on the command line (SAMPLE=path), not inferred from the query IDs, so
files produced by either script can be freely mixed.

For each query, only the TOP hit is used (NCBI always lists hits sorted by
score, best first). Queries with no significant hit are tallied separately
rather than silently dropped.

Usage:
    python3 parse_ncbi_text_results.py \
        18S-1_S37=18S-1_S37_results.txt \
        18S-3_S41=18S-3_S41_results.txt \
        18S-9_S38=18S-9_S38_results.txt \
        18S-31_S33=18S-31_S33_results.txt \
        --out species_id_report.tsv
"""
import argparse
import csv
import re
import sys
from collections import defaultdict

QUERY_HEADER_RE = re.compile(r"Query #\d+:\s*(\S+)")
# Matches trailing "_count123" or "__count123" regardless of how many
# underscores or what prefix precedes it -- robust to both
# export_unique_reads_fasta.py's "{sample}__u{N}__count{M}" and
# diagnose_failed_samples.py's "{sample}_{source}_rank{N}_count{M}" formats.
COUNT_RE = re.compile(r"_count(\d+)$")
NO_HITS_RE = re.compile(r"No significant similarity found|No hits found", re.IGNORECASE)

# Anchors on the reliable, fixed-format trailing columns (taxid, scores,
# query cover %, e-value, % identity, accession length, accession) rather
# than assuming a fixed number of spaces between columns -- NCBI's text
# export spacing is not perfectly consistent, especially after copy/paste,
# so \s+ (any run of whitespace) is used instead of a hard-coded \s{2,}.
HIT_LINE_RE = re.compile(
    r"^(?P<rest>.+?)\s+"
    r"(?P<taxid>\d+)\s+"
    r"(?P<maxscore>[\d.]+)\s+"
    r"(?P<totalscore>[\d.]+)\s+"
    r"(?P<cover>\d+)%\s+"
    r"(?P<evalue>\S+)\s+"
    r"(?P<ident>[\d.]+)\s+"
    r"(?P<acclen>\d+)\s+"
    r"(?P<accession>\S+)\s*$"
)


def parse_count(qid):
    """Returns read count parsed from the trailing '_count{N}' or
    '__count{N}' in the query ID, or None if not found."""
    m = COUNT_RE.search(qid)
    return int(m.group(1)) if m else None


def extract_label_from_rest(rest):
    """'rest' is everything before the taxid anchor: Description, Scientific
    Name, and Common Name concatenated with inconsistent spacing. Taxid
    itself is the reliable grouping key (used by the caller); this just
    tries to recover a human-readable label for display purposes, with a
    safe fallback if the split doesn't come out cleanly."""
    parts = re.split(r"\s{2,}", rest.strip())
    parts = [p.strip() for p in parts if p.strip()]
    if len(parts) >= 3:
        # last 2 parts should be Scientific Name, Common Name;
        # everything before that is Description.
        return parts[-2]
    if len(parts) == 2:
        return parts[-1]
    # Spacing collapsed entirely -- fall back to a genus-species pattern
    # search (capitalized word + lowercase word) as a best-effort label.
    m = re.search(r"[A-Z][a-z]+ [a-z]+(?: [a-z]+)?", rest)
    if m:
        return m.group(0)
    return rest.strip()[:60]  # last resort: truncated raw text


def parse_blast_text(text):
    """Returns list of (query_id, top_hit_or_None) tuples.
    top_hit_or_None is a dict with 'description', 'sciname', 'taxid', etc.,
    or None if that query had no significant hit."""
    # Defensive: if line breaks got collapsed somewhere (has happened with
    # copy-pasted browser text), force a break before every 'Query #N:'
    # marker so block-splitting below still works on a downloaded file
    # that's mostly well-formed.
    text = re.sub(r"(?<!\n)(Query #\d+:)", r"\n\1", text)

    blocks = re.split(r"(?=Query #\d+:)", text)
    results = []

    for block in blocks:
        header_match = QUERY_HEADER_RE.search(block)
        if not header_match:
            continue
        qid = header_match.group(1)

        if NO_HITS_RE.search(block):
            results.append((qid, None))
            continue

        lines = block.splitlines()
        # Find the column-header line (starts with 'Description'), the
        # actual top hit is the next non-empty line after it.
        top_hit = None
        for i, line in enumerate(lines):
            if line.strip().startswith("Description") and "Scientific" not in line:
                # this is the second header line; the real data starts after it
                for data_line in lines[i + 1:]:
                    if data_line.strip() == "":
                        continue
                    if data_line.strip().startswith("Query #"):
                        break
                    m = HIT_LINE_RE.match(data_line.strip())
                    if m:
                        top_hit = {
                            "taxid": m.group("taxid"),
                            "label": extract_label_from_rest(m.group("rest")),
                        }
                    break
                break
        results.append((qid, top_hit))

    return results


def main():
    ap = argparse.ArgumentParser(
        description="Each input is 'SAMPLE=path/to/results.txt' -- the sample "
                     "name is taken from this argument directly, not parsed "
                     "from query IDs, so files from either "
                     "export_unique_reads_fasta.py or diagnose_failed_samples.py "
                     "(different ID formats) can be mixed freely."
    )
    ap.add_argument("sample_files", nargs="+",
                     help="One or more 'SAMPLE=path/to/results.txt' pairs")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # sample -> taxid -> {"count": N, "label": display_name}
    tally = defaultdict(lambda: defaultdict(lambda: {"count": 0, "label": None}))
    unparsed_queries = []

    for entry in args.sample_files:
        if "=" not in entry:
            sys.exit(f"ERROR: '{entry}' is not in 'SAMPLE=path' format.")
        sample, path = entry.split("=", 1)

        with open(path, encoding="utf-8", errors="replace") as f:
            text = f.read()

        results = parse_blast_text(text)
        if not results:
            print(f"WARNING: no 'Query #' blocks found in {path} -- "
                  f"check this is the right file / format.", file=sys.stderr)
            continue

        for qid, top_hit in results:
            count = parse_count(qid)
            if count is None:
                unparsed_queries.append(qid)
                continue

            if top_hit is None:
                key = "No significant hit"
                tally[sample][key]["count"] += count
                tally[sample][key]["label"] = key
            else:
                key = top_hit["taxid"]  # reliable grouping key
                tally[sample][key]["count"] += count
                tally[sample][key]["label"] = top_hit["label"]  # human-readable

    if unparsed_queries:
        print(f"WARNING: {len(unparsed_queries)} query ID(s) had no parseable "
              f"'_count{{N}}' suffix and were skipped: "
              f"{unparsed_queries[:5]}{'...' if len(unparsed_queries) > 5 else ''}",
              file=sys.stderr)

    rows_out = []
    for sample, taxid_counts in tally.items():
        total = sum(v["count"] for v in taxid_counts.values())
        for taxid, v in sorted(taxid_counts.items(), key=lambda kv: -kv[1]["count"]):
            pct = v["count"] / total * 100 if total > 0 else 0.0
            rows_out.append({
                "sample": sample,
                "taxid": taxid if taxid != "No significant hit" else "-",
                "species": v["label"],
                "read_count": v["count"],
                "pct_of_total_reads": round(pct, 2),
            })

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "taxid", "species", "read_count", "pct_of_total_reads"],
                                 delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"\nWrote {len(rows_out)} rows to {args.out}\n")
    print(f"{'Sample':<16}{'Total reads':<14}{'Top species (by reads)'}")
    for sample, taxid_counts in tally.items():
        total = sum(v["count"] for v in taxid_counts.values())
        top_taxid, top_v = max(taxid_counts.items(), key=lambda kv: kv[1]["count"])
        pct = top_v["count"] / total * 100 if total > 0 else 0.0
        print(f"{sample:<16}{total:<14}{top_v['label']} ({pct:.1f}%)")


if __name__ == "__main__":
    main()
