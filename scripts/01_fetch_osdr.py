#!/usr/bin/env python3
"""Acquire the cohort studies from NASA OSDR.

For each accession in `cohorts.py` this pulls three things:

  data/counts/<acc>_counts.csv     the unnormalized count matrix (RSEM preferred)
  data/runsheets/<acc>_runsheet.csv  the GeneLab runsheet — sample IDs that match the
                                     count-matrix columns, with clean Factor Value columns
  data/osdr/study_catalog.json     title, publication, factors and MEASURED file counts

Every count in the catalogue is measured from the live file listing. Nothing about a
study is asserted from memory: if OSDR revises a deposit, the catalogue changes and
`02_harmonise_metadata.py` fails loudly rather than analysing stale assumptions.

  python3 scripts/01_fetch_osdr.py                     # both cohorts
  python3 scripts/01_fetch_osdr.py --cohort A          # one cohort
  python3 scripts/01_fetch_osdr.py --dry-run           # resolve URLs, download nothing
  python3 scripts/01_fetch_osdr.py --include-microarray
"""
from __future__ import annotations

import argparse
import datetime as dt
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cohorts import COHORTS, EXCLUSIONS, OPTIONAL_MICROARRAY, cohort_of  # noqa: E402
from lib_sources import (DATA, OSDR_STUDY_URL, download, log, osdr_file_index,  # noqa: E402
                         osdr_meta, osdr_publication, osdr_url, write_json)

COUNTS = DATA / "counts"
RUNSHEETS = DATA / "runsheets"
OSDR_DIR = DATA / "osdr"

# Preference order for the expression matrix. DREM models log-ratios between
# timepoints, so we start from *unnormalized* counts and do our own normalisation with
# the batch structure in hand — GeneLab's own normalised tables are normalised per
# study, which would bake each study's own reference into the cross-study series.
# rRNA-removed variants are skipped: only some studies have them, and mixing the two
# would change the library-size denominator between studies.
COUNT_PREFS = [
    lambda n: n.endswith("RSEM_Unnormalized_Counts_GLbulkRNAseq.csv") and "rRNArm" not in n,
    lambda n: n.endswith("STAR_Unnormalized_Counts_GLbulkRNAseq.csv") and "rRNArm" not in n,
    lambda n: bool(re.search(r"Unnormalized_Counts.*\.csv$", n)) and "rRNArm" not in n,
    lambda n: bool(re.search(r"Unnormalized_Counts.*\.csv$", n)),
]

RUNSHEET_PREFS = [
    lambda n: n.endswith("_runsheet.csv"),
    lambda n: n.endswith("SampleTable_GLbulkRNAseq.csv"),
    lambda n: bool(re.search(r"SampleTable.*\.csv$", n)),
]


def pick(pairs: list[tuple[str, str]], prefs) -> tuple[str | None, str | None]:
    """First (name, url) matching the highest-priority predicate that matches anything."""
    for pred in prefs:
        hits = [(n, u) for n, u in pairs if pred(n)]
        if hits:
            return sorted(hits)[0]
    return None, None


def fetch_study(acc: str, dry_run: bool) -> dict:
    pairs = osdr_file_index(acc)
    counts_name, counts_url = pick(pairs, COUNT_PREFS)
    run_name, run_url = pick(pairs, RUNSHEET_PREFS)

    try:
        meta = osdr_meta(acc)
    except RuntimeError as e:
        log(f"  {acc}: metadata unavailable ({e})")
        meta = {}
    try:
        pub = osdr_publication(acc)
    except RuntimeError:
        pub = None

    card = {
        "accession": acc,
        "cohort": cohort_of(acc),
        "url": OSDR_STUDY_URL.format(acc=acc),
        "title": meta.get("study title") or meta.get("title"),
        "assay_technology": meta.get("study assay technology type"),
        "publication": pub,
        # measured, not asserted
        "n_files_total": len(pairs),
        "counts_file": counts_name,
        "runsheet_file": run_name,
        "retrieved": dt.date.today().isoformat(),
    }

    if not counts_name:
        card["warning"] = "no unnormalized count matrix in the file listing"
        log(f"  {acc}: NO count matrix ({len(pairs)} files)")
    else:
        log(f"  {acc}: {counts_name}  ({len(pairs)} files total)")

    if dry_run:
        card["counts_url"] = counts_url
        card["runsheet_url"] = run_url
        return card

    if counts_name:
        dest = COUNTS / f"{acc}_counts.csv"
        download(osdr_url(counts_url), dest)
        card["counts_path"] = str(dest.relative_to(DATA.parent))
        card["counts_bytes"] = dest.stat().st_size
    if run_name:
        dest = RUNSHEETS / f"{acc}_runsheet.csv"
        download(osdr_url(run_url), dest)
        card["runsheet_path"] = str(dest.relative_to(DATA.parent))
        card["runsheet_bytes"] = dest.stat().st_size
    return card


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", choices=["A", "B"], help="restrict to one cohort")
    ap.add_argument("--dry-run", action="store_true",
                    help="resolve every URL and report measured file counts, download nothing")
    ap.add_argument("--include-microarray", action="store_true",
                    help="also fetch the optional cross-platform microarray studies")
    args = ap.parse_args()

    accs: list[str] = []
    for key, cohort in COHORTS.items():
        if args.cohort and key != args.cohort:
            continue
        accs += list(cohort["studies"])
    if args.include_microarray:
        accs += OPTIONAL_MICROARRAY
    accs = sorted(set(accs), key=lambda a: int(a.split("-")[1]))

    log(f"{'DRY RUN: ' if args.dry_run else ''}fetching {len(accs)} studies from OSDR")
    cards, failed = [], []
    for acc in accs:
        try:
            cards.append(fetch_study(acc, args.dry_run))
        except RuntimeError as e:
            log(f"  {acc}: FAILED {e}")
            failed.append(acc)

    catalog = {
        "generated": dt.date.today().isoformat(),
        "cohorts": {k: {kk: vv for kk, vv in c.items() if kk != "studies"}
                    for k, c in COHORTS.items()},
        "exclusions": EXCLUSIONS,
        "studies": cards,
        "failed": failed,
    }
    write_json(OSDR_DIR / "study_catalog.json", catalog)
    log(f"\nwrote {OSDR_DIR / 'study_catalog.json'}  "
        f"({len(cards)} studies, {len(failed)} failed)")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
