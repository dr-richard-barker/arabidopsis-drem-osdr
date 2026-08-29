#!/usr/bin/env python3
"""Build the TF-activity matrix: 479 transcription factors x every OSDR plant contrast.

Each feature is one DREM transcription factor's activity in one contrast: its target set
from `data/tf_prior/tf_gene_edges.tsv` scored against that contrast's log2 fold-change
ranking. The result is what scripts 20-22 do statistics, machine learning and trajectory
projection on.

**The null is analytic, not permuted.** For a size-n subset of N ranks with mean mu and
variance sigma^2, the subset mean has expectation mu and variance
(sigma^2 / n) * (N - n) / (N - 1) --- the finite-population correction for sampling
without replacement. That is exact. Permuting instead would cost 479 TFs x ~80 contrasts
x 2,000 draws, roughly 77 million resamples, to estimate a quantity with a closed form.
`--verify-null` checks the closed form against permutation on a few TFs, and that check is
part of the pipeline's verification rather than a one-off.

Each contrast is annotated with the mission it flew on, because studies sharing a mission
share hardware, launch and ground-control handling and are therefore not independent:
SpaceX-5 alone contributes five studies. Every downstream cross-validation groups on this.
Platform (RNA-seq or microarray) is recorded for the same reason --- it is a confounder that
has to be testable, not assumed away.

  results/tf_activity/tf_activity.tsv      contrasts x TFs, z-scores
  results/tf_activity/contrast_meta.tsv    mission, platform, factor, level, n samples
  results/qc/tf_activity_qc.json           coverage, null verification, positive control

  python3 scripts/19_tf_activity_matrix.py
  python3 scripts/19_tf_activity_matrix.py --verify-null
"""
from __future__ import annotations

import argparse
import datetime as dt
import importlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import (DATA, RESULTS, get_json, log,  # noqa: E402
                         quiet_accelerate_blas_warnings, write_json)

quiet_accelerate_blas_warnings()

scan = importlib.import_module("16_scan_osdr_plants")

OUT = RESULTS / "tf_activity"
PRIOR = DATA / "tf_prior" / "tf_gene_edges.tsv"
BIODATA = "https://visualization.osdr.nasa.gov/biodata/api/v2"

SOG1 = "AT1G25580"


def load_tf_targets(min_score: float) -> dict[str, list[str]]:
    """TF -> target loci, keeping only edges above `min_score`."""
    targets: dict[str, list[str]] = defaultdict(list)
    with open(PRIOR) as fh:
        next(fh)
        for line in fh:
            tf, gene, score = line.rstrip("\n").split("\t")
            if float(score) >= min_score:
                targets[tf].append(gene)
    return dict(targets)


def mission_of(acc: str) -> tuple[str, str]:
    """(mission name, space program) for a study, from the biodata metadata endpoint."""
    try:
        d = get_json(f"{BIODATA}/dataset/{acc}/metadata/")
        m = (d.get(acc) or {}).get("metadata", {})
    except RuntimeError:
        return "", ""

    def flat(v):
        if isinstance(v, list):
            v = v[0] if v else ""
        if isinstance(v, dict):
            v = v.get("name", "")
        return str(v or "").strip()

    mis = m.get("mission")
    if isinstance(mis, list):
        mis = mis[0] if mis else {}
    name = flat(mis.get("name") if isinstance(mis, dict) else mis)
    return name, flat(m.get("space program"))


def analytic_z(ranks: np.ndarray, idx: np.ndarray) -> float:
    """Exact z for the mean of a size-n subset drawn without replacement.

    Var of a subset mean = (sigma^2 / n) * (N - n) / (N - 1). With N in the tens of
    thousands the correction is small, but it is free and it is what makes the result
    exact rather than approximate.
    """
    N = len(ranks)
    n = len(idx)
    if n < 5 or n >= N:
        return float("nan")
    mu = ranks.mean()
    var = ranks.var(ddof=1) / n * (N - n) / (N - 1)
    if var <= 0:
        return float("nan")
    return float((ranks[idx].mean() - mu) / np.sqrt(var))


def verify_null(ranks: np.ndarray, sizes: list[int], n_perm: int,
                rng: np.random.Generator) -> list[dict]:
    """Compare the closed form against permutation. The speed-up is only worth it if exact."""
    out = []
    for n in sizes:
        idx = rng.choice(len(ranks), n, replace=False)
        obs = ranks[idx].mean()
        null = np.array([ranks[rng.choice(len(ranks), n, replace=False)].mean()
                         for _ in range(n_perm)])
        z_perm = (obs - null.mean()) / null.std()
        out.append({"set_size": int(n),
                    "z_analytic": round(analytic_z(ranks, idx), 4),
                    "z_permutation": round(float(z_perm), 4),
                    "abs_difference": round(abs(analytic_z(ranks, idx) - z_perm), 4)})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organism", default="Arabidopsis thaliana")
    ap.add_argument("--min-edge-score", type=float, default=0.5,
                    help="keep only the stronger half of each TF's edges (default 0.5); "
                         "DAP-seq is promiscuous and the weakest edges are mostly noise")
    ap.add_argument("--min-targets", type=int, default=20)
    ap.add_argument("--verify-null", action="store_true",
                    help="check the analytic null against permutation and stop")
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1260)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    targets = load_tf_targets(args.min_edge_score)
    targets = {t: g for t, g in targets.items() if len(g) >= args.min_targets}
    log(f"TF prior: {len(targets)} TFs with >= {args.min_targets} targets "
        f"at edge score >= {args.min_edge_score}")

    factors = scan.factor_table(args.organism)
    accs = sorted(factors["id.accession"].unique(), key=lambda a: int(a.split("-")[1]))

    rows, meta_rows, skipped = [], [], {}
    null_check = None

    for acc in accs:
        cons = scan.contrasts_for(acc, factors)
        if not cons:
            skipped[acc] = "no factor with an identifiable control level"
            continue
        path = scan.fetch_counts(acc)
        if path is None:
            skipped[acc] = "no expression matrix"
            continue
        try:
            counts = scan.read_expression(path)
        except Exception as e:  # noqa: BLE001
            skipped[acc] = f"unreadable ({type(e).__name__})"
            continue
        if counts is None or counts.empty:
            skipped[acc] = "no usable gene index"
            continue

        counts.index = counts.index.astype(str).str.strip().str.upper()
        counts = counts.loc[~counts.index.duplicated(keep="first")]
        counts = counts.apply(pd.to_numeric, errors="coerce").dropna(how="all")

        is_array = float(np.nanmax(counts.to_numpy())) < 100
        platform = "microarray" if is_array else "rna-seq"
        if is_array:
            logx = counts.loc[counts.notna().sum(axis=1) >= 2]
        else:
            lib = counts.sum(axis=0).replace(0, np.nan)
            cpm = counts.divide(lib, axis=1) * 1e6
            keep = (cpm >= 1.0).sum(axis=1) >= max(2, cpm.shape[1] // 4)
            logx = np.log2(cpm.loc[keep] + 1.0)

        mission, program = mission_of(acc)
        for factor, level, treated, control in cons:
            t = [s for s in treated if s in logx.columns]
            c = [s for s in control if s in logx.columns]
            if len(t) < 2 or len(c) < 2:
                continue
            lfc = (logx[t].mean(axis=1) - logx[c].mean(axis=1))
            lfc = lfc.replace([np.inf, -np.inf], np.nan).dropna()
            if len(lfc) < 1000:
                continue

            order = lfc.rank(ascending=False, method="average").to_numpy()
            ranks = 1.0 - (order - 1) / (len(order) - 1)
            pos = {g: i for i, g in enumerate(lfc.index)}

            if args.verify_null and null_check is None:
                null_check = verify_null(ranks, [20, 100, 500, 1500, 3000],
                                         args.permutations, rng)

            key = f"{acc}|{factor}|{level}"
            row = {"contrast": key}
            for tf, genes in targets.items():
                idx = np.fromiter((pos[g] for g in genes if g in pos), dtype=int)
                row[tf] = analytic_z(ranks, idx)
            rows.append(row)
            meta_rows.append({
                "contrast": key, "accession": acc, "factor": factor, "level": level,
                "mission": mission or "(unknown)", "space_program": program,
                "platform": platform, "n_treated": len(t), "n_control": len(c),
                "n_genes": int(len(lfc)),
                "is_flight": factor in ("spaceflight", "altered gravity"),
                "is_radiation": factor in ("ionizing radiation", "absorbed radiation dose"),
            })
        log(f"  {acc}: {platform}, mission {mission or '?'}, "
            f"{sum(1 for m in meta_rows if m['accession'] == acc)} contrasts")

    if args.verify_null:
        log("\nanalytic vs permutation null:")
        for r in (null_check or []):
            log(f"  n={r['set_size']:<5} analytic {r['z_analytic']:>8}  "
                f"permutation {r['z_permutation']:>8}  |diff| {r['abs_difference']}")
        return 0

    if not rows:
        raise SystemExit("no contrasts scored")

    act = pd.DataFrame(rows).set_index("contrast")
    meta = pd.DataFrame(meta_rows).set_index("contrast")
    OUT.mkdir(parents=True, exist_ok=True)
    act.round(4).to_csv(OUT / "tf_activity.tsv", sep="\t")
    meta.to_csv(OUT / "contrast_meta.tsv", sep="\t")

    # Positive control: in genuine gamma irradiation, SOG1 must be among the most active
    # TFs. If it is not, the features are wrong and nothing downstream is worth running.
    gamma = meta[(meta["is_radiation"]) & (meta["level"].str.contains("gamma", case=False))]
    sog1_ranks = {}
    if SOG1 in act.columns:
        for k in gamma.index:
            r = act.loc[k].rank(ascending=False)
            sog1_ranks[k] = {"sog1_z": round(float(act.loc[k, SOG1]), 2),
                             "sog1_rank": int(r[SOG1]), "of": int(act.shape[1])}
    ok = bool(sog1_ranks) and np.median([v["sog1_rank"] for v in sog1_ranks.values()]) <= 20

    write_json(RESULTS / "qc" / "tf_activity_qc.json", {
        "generated": dt.date.today().isoformat(),
        "null": "analytic (finite-population correction), verified against permutation "
                "via --verify-null",
        "n_tfs": int(act.shape[1]), "n_contrasts": int(act.shape[0]),
        "n_studies": int(meta["accession"].nunique()),
        "n_missions": int(meta.loc[meta["is_flight"], "mission"].nunique()),
        "missions": sorted(meta.loc[meta["is_flight"], "mission"].unique().tolist()),
        "platforms": meta["platform"].value_counts().to_dict(),
        "n_flight_contrasts": int(meta["is_flight"].sum()),
        "n_radiation_contrasts": int(meta["is_radiation"].sum()),
        "positive_control_sog1_in_gamma": sog1_ranks,
        "positive_control_pass": ok,
        "skipped": skipped,
    })
    log(f"\nTF activity: {act.shape[0]} contrasts x {act.shape[1]} TFs "
        f"({meta['accession'].nunique()} studies)")
    log(f"  flight contrasts {int(meta['is_flight'].sum())} over "
        f"{meta.loc[meta['is_flight'], 'mission'].nunique()} missions; "
        f"platforms {meta['platform'].value_counts().to_dict()}")
    log(f"  positive control (SOG1 rank in gamma): "
        f"{'PASS' if ok else 'FAIL'} {[v['sog1_rank'] for v in sog1_ranks.values()]}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
