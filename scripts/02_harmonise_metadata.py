#!/usr/bin/env python3
"""Harmonise the cohort sample metadata into one table with a real time axis.

Reads the OSDR biodata factor-value table (sample-level, every Arabidopsis study) plus
the per-study GeneLab runsheets, and writes:

  data/metadata_master.csv    one row per sample: accession, sample id, cohort, arm,
                              genotype, irradiated/flight flag, time in minutes (A) or
                              days (B), and the count-matrix column it maps to
  results/qc/metadata_qc.json the checks below, pass or fail

Five checks decide whether the pseudo-time-series is safe to build. Each *fails the
run* rather than warning, because every one of them silently corrupts the model:

  1. Duplicate deposits. OSD-219 re-deposits OSD-218's 32 BioSamples verbatim. This is
     re-detected from sample-name overlap, not hardcoded — if OSDR ever de-duplicates
     the pair, the check notices the exclusion is no longer needed.
  2. Timepoint drift. The timepoints recovered from the live API must match
     `cohorts.py`'s expectations. OSDR revises deposits.
  3. Sample-ID join. Every metadata row must map to a real count-matrix column, or the
     expression matrix silently loses samples.
  4. Cohort shape. The assembled series must equal the union of the per-study
     expectations — catches a study contributing an arm nobody declared.
  5. Arm balance. Each (arm, timepoint) cell must have >= 2 replicates.

As of 2026-08-28 these pass with cohort A at 108 timed samples over 10 timepoints
(10 min - 72 h) plus 36 untreated anchor samples, and cohort B at 115 samples over
3 timepoints (4/6/8 d). Those numbers are reported by the checks, never assumed by them.

  python3 scripts/02_harmonise_metadata.py
  python3 scripts/02_harmonise_metadata.py --refresh   # re-download the factor table
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import io
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cohorts import COHORTS, EXCLUSIONS  # noqa: E402
from lib_sources import DATA, RESULTS, log, osdr_factor_table, write_json  # noqa: E402

FACTOR_CACHE = DATA / "osdr" / "arabidopsis_factor_values.csv"
MASTER = DATA / "metadata_master.csv"
QC = RESULTS / "qc" / "metadata_qc.json"

# The OSDR factor-value columns that can carry a time axis, in priority order. Studies
# disagree about which one they use for the same concept: the Bourbousse accessions put
# post-irradiation time in "time of sample collection after treatment", while the
# spaceflight studies split between "age" and "time" for seedling age.
TIME_COLUMNS = [
    "study.factor value.time of sample collection after treatment",
    "study.factor value.time",
    "study.factor value.age",
    "study.factor value.exposure duration",
    "study.factor value.growth time",
]

UNIT_MINUTES = {"minute": 1.0, "hour": 60.0, "day": 1440.0, "week": 10080.0}

NA = {"NaN", "nan", "", "{Not Applicable}", "Not Applicable", None}


def nz(v) -> bool:
    return v not in NA


def parse_time(raw: str) -> float | None:
    """'90 {minute}' -> 90.0 minutes. Returns None for free-text time labels.

    Several parabolic-flight studies use prose in the time field ('Fixed at the end of
    the first parabola'); those are not on the ordered numeric axis DREM needs, so they
    are dropped rather than guessed at.
    """
    if not nz(raw):
        return None
    m = re.match(r"\s*([\d.]+)\s*\{(\w+)\}\s*$", raw)
    if not m:
        return None
    value, unit = float(m.group(1)), m.group(2).lower()
    if unit not in UNIT_MINUTES:
        return None
    return value * UNIT_MINUTES[unit]


def load_factor_table(refresh: bool) -> list[dict]:
    if refresh or not FACTOR_CACHE.exists():
        log("fetching the Arabidopsis factor-value table from the biodata API ...")
        text = osdr_factor_table("Arabidopsis thaliana")
        FACTOR_CACHE.parent.mkdir(parents=True, exist_ok=True)
        FACTOR_CACHE.write_text(text)
    return list(csv.DictReader(io.StringIO(FACTOR_CACHE.read_text())))


def genotype_of(row: dict) -> str:
    """A single genotype label from whichever column the study happened to use."""
    for col in ("study.factor value.genotype", "study.factor value.cultivar",
                "study.factor value.ecotype"):
        if nz(row.get(col)):
            return row[col]
    return "Wild Type"


def treated_of(row: dict, cohort_key: str) -> str:
    """'treated' vs 'control' — the arm split, per cohort."""
    if cohort_key == "A":
        ir = row.get("study.factor value.ionizing radiation")
        if not nz(ir) or "non-irrad" in ir.lower():
            return "control"
        return "treated"
    sf = row.get("study.factor value.spaceflight")
    if nz(sf) and "flight" in sf.lower() and "ground" not in sf.lower():
        return "treated"
    return "control"


def count_matrix_columns(acc: str) -> list[str]:
    path = DATA / "counts" / f"{acc}_counts.csv"
    if not path.exists():
        return []
    with open(path) as fh:
        header = next(csv.reader(fh))
    return [c for c in header[1:] if c]


# --------------------------------------------------------------------------- checks

def check_duplicates(rows: list[dict]) -> dict:
    """Re-derive duplicate deposits from sample-name overlap rather than trusting a list."""
    by_acc = defaultdict(set)
    for r in rows:
        by_acc[r["id.accession"]].add(r["id.sample name"])

    accs = sorted(by_acc, key=lambda a: int(a.split("-")[1]))
    dupes = []
    for i, a in enumerate(accs):
        for b in accs[i + 1:]:
            if by_acc[a] and by_acc[a] == by_acc[b]:
                dupes.append({"kept": a, "dropped": b, "n_samples": len(by_acc[a])})

    declared = {a for a in EXCLUSIONS if "Duplicate" in EXCLUSIONS[a]}
    detected = {d["dropped"] for d in dupes}
    return {
        "detected_duplicate_pairs": dupes,
        "declared_duplicate_exclusions": sorted(declared),
        # A declared exclusion that is no longer a real duplicate would silently drop
        # good data, so surface that direction too.
        "declared_but_not_detected": sorted(declared - detected),
        "detected_but_not_declared": sorted(detected - declared),
        "pass": detected == declared,
    }


def check_timepoints(observed: dict[str, list[float]]) -> dict:
    per_study, ok = {}, True
    for key, cohort in COHORTS.items():
        scale = 1.0 if cohort["time_unit"] == "minute" else 1440.0
        for acc, spec in cohort["studies"].items():
            expected = sorted(float(t) * scale for t in spec["expected_timepoints"])
            got = sorted(set(observed.get(acc, [])))
            per_study[acc] = {"expected_minutes": expected, "observed_minutes": got,
                              "match": expected == got}
            ok &= expected == got
    return {"per_study": per_study, "pass": ok}


def check_joins(rows: list[dict]) -> dict:
    per_study, ok = {}, True
    for acc in sorted({r["id.accession"] for r in rows}, key=lambda a: int(a.split("-")[1])):
        cols = set(count_matrix_columns(acc))
        names = {r["id.sample name"] for r in rows if r["id.accession"] == acc}
        if not cols:
            per_study[acc] = {"status": "count matrix not downloaded", "matched": 0,
                              "n_metadata": len(names)}
            ok = False
            continue
        matched = names & cols
        per_study[acc] = {
            "n_metadata": len(names), "n_matrix_columns": len(cols),
            "matched": len(matched),
            "unmatched_metadata": sorted(names - cols)[:8],
            "unmatched_columns": sorted(cols - names)[:8],
        }
        ok &= len(matched) == len(names)
    return {"per_study": per_study, "pass": ok}


def check_cohort_shape(master: list[dict]) -> dict:
    """The assembled series must match the union of the per-study expectations.

    Derived from `cohorts.py` rather than compared against a remembered sample count,
    so the check stays true if a study is added or dropped. What it catches is a study
    contributing timepoints nobody declared — e.g. a revised deposit gaining an arm.
    """
    per_cohort, ok = {}, True
    for key, cohort in COHORTS.items():
        scale = 1.0 if cohort["time_unit"] == "minute" else 1440.0
        excluded = set(EXCLUSIONS)
        expected = sorted({float(t) * scale
                           for acc, spec in cohort["studies"].items() if acc not in excluded
                           for t in spec["expected_timepoints"]})
        rows = [r for r in master if r["cohort"] == key]
        timed = [r for r in rows if r["time_minutes"] != ""]
        got = sorted({float(r["time_minutes"]) for r in timed})
        per_cohort[key] = {
            "name": cohort["name"],
            "n_samples": len(rows),
            "n_timed_samples": len(timed),
            "n_untimed_anchor_samples": len(rows) - len(timed),
            "expected_timepoints": expected,
            "observed_timepoints": got,
            "match": expected == got,
        }
        ok &= expected == got
    return {"per_cohort": per_cohort, "pass": ok}


def check_balance(master: list[dict], min_reps: int = 2) -> dict:
    cells = Counter((r["cohort"], r["arm"], r["time_minutes"]) for r in master
                    if r["time_minutes"] != "")
    thin = {f"{c}|{a}|{t}": n for (c, a, t), n in sorted(cells.items()) if n < min_reps}
    return {"min_replicates": min_reps, "n_cells": len(cells),
            "under_replicated": thin, "pass": not thin}


# --------------------------------------------------------------------------- main

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refresh", action="store_true", help="re-download the factor table")
    ap.add_argument("--allow-fail", action="store_true",
                    help="write outputs and report failures without exiting non-zero")
    args = ap.parse_args()

    rows = load_factor_table(args.refresh)
    log(f"factor table: {len(rows)} samples across "
        f"{len({r['id.accession'] for r in rows})} Arabidopsis studies")

    dup = check_duplicates(rows)
    excluded = {d["dropped"] for d in dup["detected_duplicate_pairs"]} | set(EXCLUSIONS)

    master: list[dict] = []
    observed: dict[str, list[float]] = defaultdict(list)
    for key, cohort in COHORTS.items():
        for acc in cohort["studies"]:
            if acc in excluded:
                log(f"  {acc}: excluded — {EXCLUSIONS.get(acc, 'duplicate deposit')}")
                continue
            for r in (x for x in rows if x["id.accession"] == acc):
                minutes = next((t for c in TIME_COLUMNS
                                if (t := parse_time(r.get(c, ""))) is not None), None)
                if minutes is not None:
                    observed[acc].append(minutes)
                master.append({
                    "accession": acc,
                    "sample_name": r["id.sample name"],
                    "cohort": key,
                    "cohort_name": cohort["name"],
                    "role": cohort["studies"][acc]["role"],
                    "genotype": genotype_of(r),
                    "arm": treated_of(r, key),
                    "time_minutes": "" if minutes is None else f"{minutes:g}",
                    "time_days": "" if minutes is None else f"{minutes / 1440:g}",
                    "dose": r.get("study.factor value.absorbed radiation dose", ""),
                    "radiation": r.get("study.factor value.ionizing radiation", ""),
                    "spaceflight": r.get("study.factor value.spaceflight", ""),
                    "altered_gravity": r.get("study.factor value.altered gravity", ""),
                    "organism_part": r.get("study.factor value.organism part", ""),
                })

    MASTER.parent.mkdir(parents=True, exist_ok=True)
    with open(MASTER, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(master[0]))
        w.writeheader()
        w.writerows(master)
    log(f"wrote {MASTER}  ({len(master)} samples)")

    checks = {
        "generated": dt.date.today().isoformat(),
        "n_samples": len(master),
        "duplicates": dup,
        "timepoints": check_timepoints(observed),
        "sample_joins": check_joins([r for r in rows
                                     if r["id.accession"] in {m['accession'] for m in master}]),
        "cohort_shape": check_cohort_shape(master),
        "balance": check_balance(master),
    }
    checks["pass"] = all(checks[k]["pass"] for k in
                         ("duplicates", "timepoints", "sample_joins", "cohort_shape", "balance"))
    write_json(QC, checks)

    for name in ("duplicates", "timepoints", "sample_joins", "cohort_shape", "balance"):
        log(f"  {'PASS' if checks[name]['pass'] else 'FAIL'}  {name}")
    if not checks["pass"]:
        log(f"\nQC failures detailed in {QC}")
        log(json.dumps({k: v for k, v in checks.items()
                        if isinstance(v, dict) and not v.get("pass")}, indent=1)[:3000])
    return 0 if (checks["pass"] or args.allow_fail) else 1


if __name__ == "__main__":
    raise SystemExit(main())
