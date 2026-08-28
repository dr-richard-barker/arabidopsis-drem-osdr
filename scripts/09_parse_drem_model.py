#!/usr/bin/env python3
"""Parse DREM's batch output into tidy tables.

DREM writes three artefacts per run and none of them are analysis-ready:

  model.txt        one section per split: an INTERCEPT line followed by a regression
                   coefficient for every TF in the prior
  geneassign.txt   gene, assigned node id, then the gene's trajectory values
  tfscores/*.txt   TF activity scores, in THREE different file types with THREE
                   different column schemas

The tfscores directory is the part that punishes assumptions:

  split_<n>_<x>_<y>_to_<a>_<b>_vsOthers.txt
      One branch against the others. Columns are
      `TF, Coeff, <SIDE> 0, <SIDE> 1, Other 0, Other 1, Avg. <SIDE>, Avg. Other, Diff, Score`
      where SIDE is *either* `High` or `Low` depending on which way the branch went.
      Keying on "Avg. High" therefore drops every low-branch file on the floor.
  split_<n>_<x>_<y>_2way.txt
      The two-way split test: `Avg. Low` and `Avg. High` together.
  path_<n>_<x>_<y>_to_<a>_<b>.txt
      Per-path enrichment, an entirely different schema:
      `TF, Num Total, Num Parent, Num Path, Expect Overall, Diff. Overall, Score Overall,
       Expect Split, Diff. Split, Score Split, % Split`.
      Its TF field carries a trailing token (`AT1G01060 1`) that must be stripped.

`Score` is DREM's significance statistic, where SMALLER means more significant. That
direction is easy to invert by accident, so it is recorded in the output and every
ranking here sorts ascending.

Outputs, per run:
  results/drem/parsed/<run>_tfscores.tsv   long table over all three file types
  results/drem/parsed/<run>_genes.tsv      gene -> node assignment
  results/drem/parsed/<run>_summary.json   counts and the top TF per split

  python3 scripts/09_parse_drem_model.py
  python3 scripts/09_parse_drem_model.py --top 10
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import RESULTS, log, write_json  # noqa: E402

RUNS = RESULTS / "drem" / "runs"
OUT = RESULTS / "drem" / "parsed"

NUM = r"-?[\d.]+"
SPLIT_RE = re.compile(
    rf"^split_(\d+)_({NUM})_({NUM})(?:_to_({NUM})_({NUM}))?_(vsOthers|2way)\.txt$")
PATH_RE = re.compile(rf"^path_(\d+)_({NUM})_({NUM})_to_({NUM})_({NUM})\.txt$")


def _f(row: dict, *names: str) -> float | None:
    """First present, parsable column among `names`."""
    for n in names:
        if n in row and row[n] not in (None, "", "NaN"):
            try:
                return float(row[n])
            except ValueError:
                return None
    return None


def _tf(raw: str) -> str:
    """`AT1G01060` or `AT1G01060 1` -> `AT1G01060`."""
    return str(raw).split()[0].strip().upper() if raw else ""


def parse_tfscores(run_dir: Path) -> tuple[pd.DataFrame, dict]:
    rows: list[dict] = []
    stats = {"vsOthers": 0, "2way": 0, "path": 0, "empty": 0, "unrecognised": []}
    tfdir = run_dir / "tfscores"
    if not tfdir.exists():
        return pd.DataFrame(), stats

    for path in sorted(tfdir.glob("*.txt")):
        text = path.read_text().strip()
        if not text:
            # DREM writes an empty file where a test does not apply. Normal, not a
            # parse failure — but counted, so "no data" never hides behind "no file".
            stats["empty"] += 1
            continue

        m_split = SPLIT_RE.match(path.name)
        m_path = PATH_RE.match(path.name)
        if not (m_split or m_path):
            stats["unrecognised"].append(path.name)
            continue

        reader = csv.DictReader(text.splitlines(), delimiter="\t")
        if m_split:
            node, fx, fy, tx, ty, kind = m_split.groups()
            stats[kind] += 1
            for r in reader:
                tf = _tf(r.get("TF", ""))
                if not tf:
                    continue
                # The branch side is High for an up-branch and Low for a down-branch.
                avg_branch = _f(r, "Avg. High", "Avg. Low")
                rows.append({
                    "kind": kind, "split": int(node),
                    "from_time": float(fx), "from_value": float(fy),
                    "to_time": float(tx) if tx else None,
                    "to_value": float(ty) if ty else None,
                    "TF": tf,
                    "coeff": _f(r, "Coeff"),
                    "avg_branch": avg_branch,
                    "avg_other": _f(r, "Avg. Other", "Avg. Low"),
                    "diff": _f(r, "Diff"),
                    "score": _f(r, "Score"),
                })
        else:
            node, fx, fy, tx, ty = m_path.groups()
            stats["path"] += 1
            for r in reader:
                tf = _tf(r.get("TF", ""))
                if not tf:
                    continue
                rows.append({
                    "kind": "path", "split": int(node),
                    "from_time": float(fx), "from_value": float(fy),
                    "to_time": float(tx), "to_value": float(ty),
                    "TF": tf,
                    "coeff": _f(r, "Diff. Overall"),
                    "avg_branch": _f(r, "Num Path"),
                    "avg_other": _f(r, "Expect Overall"),
                    "diff": _f(r, "Diff. Split", "Diff. Overall"),
                    # Split-level significance where DREM gives it, else path-level.
                    "score": _f(r, "Score Split", "Score Overall"),
                })

    df = pd.DataFrame(rows)
    return df[df["score"].notna()] if not df.empty else df, stats


def parse_geneassign(run_dir: Path) -> pd.DataFrame:
    path = run_dir / "geneassign.txt"
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    rows = []
    for line in path.read_text().splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        try:
            node = int(f[1])
        except ValueError:
            continue
        traj = []
        for v in f[2:]:
            try:
                traj.append(float(v))
            except ValueError:
                pass
        rows.append({"gene": f[0].strip().upper(), "node": node,
                     "trajectory": ",".join(f"{v:g}" for v in traj)})
    return pd.DataFrame(rows)


def parse_model_sections(run_dir: Path) -> int:
    path = run_dir / "model.txt"
    if not path.exists():
        return 0
    return sum(1 for line in path.read_text().splitlines()
               if line.startswith("INTERCEPT"))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top", type=int, default=5,
                    help="how many top TFs per split to record in the summary")
    args = ap.parse_args()

    if not RUNS.exists():
        raise SystemExit("no DREM runs — run 08_run_drem.py first")
    OUT.mkdir(parents=True, exist_ok=True)

    report = {"generated": dt.date.today().isoformat(),
              "score_direction": "lower Score = more significant (DREM convention)",
              "runs": {}}

    for run_dir in sorted(d for d in RUNS.iterdir() if d.is_dir()):
        tf, stats = parse_tfscores(run_dir)
        genes = parse_geneassign(run_dir)
        if tf.empty and genes.empty:
            log(f"  {run_dir.name}: no parsable output — skipped")
            continue

        if not tf.empty:
            tf.sort_values(["kind", "split", "score"]).to_csv(
                OUT / f"{run_dir.name}_tfscores.tsv", sep="\t", index=False)
        if not genes.empty:
            genes.to_csv(OUT / f"{run_dir.name}_genes.tsv", sep="\t", index=False)

        # Prefer the one-vs-others test; fall back to the two-way test, then paths, so
        # a run is summarised on the most specific evidence it actually produced.
        top_per_split = {}
        basis = None
        if not tf.empty:
            for kind in ("vsOthers", "2way", "path"):
                block = tf[tf["kind"] == kind]
                if block.empty:
                    continue
                basis = kind
                for split, sub in block.groupby("split"):
                    best = sub.nsmallest(args.top, "score")
                    top_per_split[str(split)] = [
                        {"TF": r.TF, "score": round(r.score, 6),
                         "diff": None if r.diff is None else round(r.diff, 4)}
                        for r in best.itertuples()]
                break

        report["runs"][run_dir.name] = {
            "n_model_sections": parse_model_sections(run_dir),
            "tfscore_files": {k: v for k, v in stats.items() if k != "unrecognised"},
            "unrecognised_files": stats["unrecognised"],
            "attribution_basis": basis,
            "n_splits": int(tf["split"].nunique()) if not tf.empty else 0,
            "n_tfs_scored": int(tf["TF"].nunique()) if not tf.empty else 0,
            "n_genes_assigned": int(len(genes)),
            "n_nodes": int(genes["node"].nunique()) if not genes.empty else 0,
            "top_tfs_per_split": top_per_split,
        }
        r = report["runs"][run_dir.name]
        log(f"  {run_dir.name}: {r['n_nodes']} nodes, {r['n_splits']} splits "
            f"({basis}), {r['n_tfs_scored']} TFs scored, {len(genes)} genes")

    write_json(RESULTS / "drem" / "parsed_summary.json", report)
    log(f"\nparsed {len(report['runs'])} runs -> {OUT}")
    return 0 if report["runs"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
