#!/usr/bin/env python3
"""
fetch_staxids.py

Replaces the earlier fetch_staxids.sh attempt, which incorrectly assumed
blastdbcmd supports a -remote flag (it doesn't -- that error was mine).

This resolves accession -> taxid using NCBI's E-utilities esummary
endpoint instead, which needs no local database and no BLAST tools at
all -- just plain HTTP requests, batched to stay within NCBI's usage
guidelines (https://www.ncbi.nlm.nih.gov/books/NBK25497/): max 3
requests/second without an API key (10/sec with one).

Usage:
    python3 fetch_staxids.py --blast-dir blast_results --out blast_results/accession_taxid_lookup.tsv
    python3 fetch_staxids.py --blast-dir blast_results --out blast_results/accession_taxid_lookup.tsv --api-key YOUR_KEY

Resumable: already-resolved accessions (already in --out) are skipped on
re-run, so if this is interrupted or a batch fails, just run it again.
"""
import argparse
import glob
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
BATCH_SIZE = 150          # accessions per request; well under URL length limits
SLEEP_NO_KEY = 0.4        # ~3 req/sec ceiling without an API key
SLEEP_WITH_KEY = 0.11     # ~10 req/sec ceiling with an API key


def collect_accessions(blast_dir: str) -> set[str]:
    accs: set[str] = set()
    for path in glob.glob(os.path.join(blast_dir, "*_blast.tsv")):
        with open(path) as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2:
                    accs.add(fields[1])  # sacc is column 2
    return accs


def load_already_done(out_path: str) -> set[str]:
    done: set[str] = set()
    if os.path.exists(out_path):
        with open(out_path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if parts:
                    done.add(parts[0])
    return done


def fetch_batch(accessions: list[str], api_key: str | None) -> dict[str, str]:
    """Returns {accession: taxid_str} for whichever of the given accessions
    NCBI successfully resolved."""
    params = {
        "db": "nuccore",
        "id": ",".join(accessions),
        "retmode": "json",
    }
    if api_key:
        params["api_key"] = api_key
    url = f"{ESUMMARY_URL}?{urllib.parse.urlencode(params)}"

    with urllib.request.urlopen(url, timeout=60) as resp:
        data = json.loads(resp.read().decode("utf-8"))

    result = data.get("result", {})
    resolved: dict[str, str] = {}
    for uid in result.get("uids", []):
        record = result.get(uid, {})
        taxid = record.get("taxid")
        # accessionversion (e.g. "KF620607.1") -- strip version to match our
        # accession-without-version sacc values; caption is a fallback
        acc_with_version = record.get("accessionversion", "") or record.get("caption", "")
        acc_no_version = acc_with_version.split(".")[0]
        if taxid and acc_no_version:
            resolved[acc_no_version] = str(taxid)
    return resolved


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--blast-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--api-key", default=None, help="Optional NCBI API key for a higher rate limit")
    args = ap.parse_args()

    all_accessions = collect_accessions(args.blast_dir)
    already_done = load_already_done(args.out)
    todo = sorted(all_accessions - already_done)

    print(f"Total unique accessions: {len(all_accessions)}")
    print(f"Already resolved: {len(already_done)}")
    print(f"Remaining to look up: {len(todo)}")

    if not todo:
        print("Nothing to do.")
        return

    sleep_time = SLEEP_WITH_KEY if args.api_key else SLEEP_NO_KEY
    n_resolved_this_run = 0
    unresolved: list[str] = []

    with open(args.out, "a") as out_fh:
        for i in range(0, len(todo), BATCH_SIZE):
            batch = todo[i : i + BATCH_SIZE]
            print(f"[batch {i // BATCH_SIZE + 1}] resolving {len(batch)} accessions ...")
            try:
                resolved = fetch_batch(batch, args.api_key)
            except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as e:
                print(f"  ERROR on this batch: {e} -- will remain unresolved, re-run script to retry")
                unresolved.extend(batch)
                time.sleep(sleep_time)
                continue

            for acc in batch:
                if acc in resolved:
                    out_fh.write(f"{acc}\t{resolved[acc]}\n")
                    n_resolved_this_run += 1
                else:
                    unresolved.append(acc)
            out_fh.flush()
            time.sleep(sleep_time)

    print(f"\nResolved {n_resolved_this_run}/{len(todo)} accessions this run -> {args.out}")
    if unresolved:
        print(f"{len(unresolved)} accession(s) could not be resolved. First few: {unresolved[:10]}")
        print("Re-running this script will retry only the still-unresolved ones.")


if __name__ == "__main__":
    main()
