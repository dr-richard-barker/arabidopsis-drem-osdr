#!/usr/bin/env python3
"""Project spaceflight contrasts onto the DREM radiation trajectory.

The interpretive step, and the one that answers "can flight be explained in terms of the
DREM time-series model". The radiation model spans 10 min to 24 h of wild-type response.
Each of its timepoints has a TF-activity profile, computed exactly as the query contrasts
are. Correlating a flight contrast against all eight asks: if this looked like any phase
of the radiation response, which phase would it be?

Three things keep the answer honest.

**A best match always exists, and counting them is not enough.** Correlating against eight
profiles and reporting the largest always returns something. A per-contrast shuffle test
helps, but it is still too weak on its own: it passed 26 of 35 flight contrasts whose mean
correlation profile is flat to three decimal places, because a single noisy contrast can
beat its own shuffle by chance. Two population-level statistics decide the result instead:
the SHAPE of the mean correlation profile across the trajectory, and whether the
best-matching timepoints CONCENTRATE near that peak. Genuine irradiation gives a peaked
profile with its argmaxes at the peak; noise gives a flat profile with argmaxes scattered
across the trajectory.

(A plain chi-square test of uniformity was tried and discarded: it rejects uniformity for
both groups, for opposite reasons -- irradiation because its argmaxes pile up at the peak,
flight because several timepoints happen to draw none at all -- so it does not separate
signal from noise here. What separates them is the height of the peak and whether the
argmaxes sit near it.)

**The irradiation contrasts calibrate the scale.** Real gamma-irradiation contrasts are
projected too. They should match strongly and early; if they do not, the projection is not
measuring what it claims and no flight result from it is worth reading.

**A flat profile is an answer.** If flight contrasts correlate with no timepoint above
null, that is a finding --- spaceflight is not a time-shifted radiation response --- and it
is reported as such rather than as a weak match to whichever timepoint scored highest.

  results/tf_activity/drem_projection.tsv    contrast x timepoint correlations
  results/qc/drem_projection_qc.json         calibration, matches, and the null

  python3 scripts/22_drem_projection.py
"""
from __future__ import annotations

import argparse
import datetime as dt
import importlib
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import (DATA, RESULTS, log,  # noqa: E402
                         quiet_accelerate_blas_warnings, write_json)

quiet_accelerate_blas_warnings()

matrix = importlib.import_module("19_tf_activity_matrix")

ACT = RESULTS / "tf_activity"
SERIES = RESULTS / "pseudotimeseries"


def drem_tf_profiles(targets: dict[str, list[str]], arm: str) -> pd.DataFrame:
    """TF activity at each DREM timepoint, built the same way as the query features."""
    expr = pd.read_csv(SERIES / f"{arm}_expression.tsv", sep="\t", index_col=0)
    expr.index = expr.index.astype(str).str.upper()
    expr = expr[sorted(expr.columns, key=float)]

    out = {}
    for tp in expr.columns:
        lfc = expr[tp].replace([np.inf, -np.inf], np.nan).dropna()
        order = lfc.rank(ascending=False, method="average").to_numpy()
        ranks = 1.0 - (order - 1) / (len(order) - 1)
        pos = {g: i for i, g in enumerate(lfc.index)}
        col = {}
        for tf, genes in targets.items():
            idx = np.fromiter((pos[g] for g in genes if g in pos), dtype=int)
            col[tf] = matrix.analytic_z(ranks, idx)
        out[tp] = col
    return pd.DataFrame(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", default="A_primary_WildType")
    ap.add_argument("--permutations", type=int, default=1000)
    ap.add_argument("--min-edge-score", type=float, default=0.5)
    ap.add_argument("--min-targets", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1260)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    act = pd.read_csv(ACT / "tf_activity.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(ACT / "contrast_meta.tsv", sep="\t", index_col=0)

    targets = matrix.load_tf_targets(args.min_edge_score)
    targets = {t: g for t, g in targets.items()
               if len(g) >= args.min_targets and t in act.columns}
    ref = drem_tf_profiles(targets, args.arm).dropna(how="any")
    tfs = [t for t in ref.index if t in act.columns]
    ref = ref.loc[tfs]
    log(f"DREM reference: {len(tfs)} TFs x {ref.shape[1]} timepoints ({args.arm})")

    rows = []
    for key in act.index:
        v = act.loc[key, tfs].to_numpy(dtype=float)
        if not np.isfinite(v).all():
            continue
        cors = {tp: float(stats.spearmanr(v, ref[tp].to_numpy())[0]) for tp in ref.columns}
        best_tp = max(cors, key=lambda t: cors[t])
        best = cors[best_tp]

        # Null: shuffle the TF labels of the query, keeping the reference intact. This
        # asks whether THIS contrast's particular pattern matches, not whether any two
        # TF vectors of this length correlate.
        null = np.array([max(float(stats.spearmanr(rng.permutation(v),
                                                   ref[tp].to_numpy())[0])
                             for tp in ref.columns)
                         for _ in range(args.permutations // 10)])
        p = float((np.sum(null >= best) + 1) / (len(null) + 1))

        m = meta.loc[key]
        rows.append({
            "contrast": key, "accession": m["accession"], "factor": m["factor"],
            "level": m["level"], "mission": m["mission"], "platform": m["platform"],
            "is_flight": bool(m["is_flight"]), "is_radiation": bool(m["is_radiation"]),
            "best_timepoint_min": float(best_tp), "best_rho": round(best, 4),
            "p_vs_shuffled_tfs": round(p, 5),
            "matched": bool(p < 0.05),
            **{f"rho_{tp}": round(cors[tp], 4) for tp in ref.columns},
        })

    df = pd.DataFrame(rows)
    df.to_csv(ACT / "drem_projection.tsv", sep="\t", index=False)

    radn = df[df["is_radiation"]]
    fl = df[df["is_flight"]]

    rho_cols = [f"rho_{tp}" for tp in ref.columns]
    times = [float(tp) for tp in ref.columns]

    def summarise(sub: pd.DataFrame) -> dict:
        if sub.empty:
            return {}
        profile = [round(float(sub[c].mean()), 4) for c in rho_cols]
        peak_i = int(np.argmax(profile))
        # Uniformity of the argmax. A real phase match concentrates at the peak; noise
        # spreads the argmax evenly, which is what a flat profile produces.
        counts = sub["best_timepoint_min"].value_counts().reindex(times).fillna(0)
        # Do the argmaxes sit at the profile's own peak, or wander? Measured as the
        # fraction landing on the peak timepoint or an immediate neighbour.
        near = {times[i] for i in (peak_i - 1, peak_i, peak_i + 1) if 0 <= i < len(times)}
        frac_near = float(sum(v for k, v in counts.items() if k in near) / max(len(sub), 1))
        return {
            "n": int(len(sub)),
            "n_matched_per_contrast_shuffle": int(sub["matched"].sum()),
            "mean_rho_profile": dict(zip([f"{t:g}" for t in times], profile)),
            "peak_timepoint_min": times[peak_i],
            "peak_mean_rho": profile[peak_i],
            "profile_range": round(max(profile) - min(profile), 4),
            "best_timepoint_counts": {f"{k:g}": int(v) for k, v in counts.items()},
            "frac_argmax_near_peak": round(frac_near, 3),
        }

    report = {
        "generated": dt.date.today().isoformat(),
        "reference_arm": args.arm,
        "n_tfs": len(tfs),
        "timepoints_min": [float(c) for c in ref.columns],
        "calibration_radiation": summarise(radn),
        "flight": summarise(fl),
        "interpretation": None,
    }
    # State the conclusion in the artefact, not only in the prose that cites it, and
    # base it on the population statistics rather than on a count of per-contrast argmaxes.
    cal, flt = report["calibration_radiation"], report["flight"]
    MIN_PEAK, MIN_NEAR = 0.25, 0.5
    if not cal or cal["peak_mean_rho"] < MIN_PEAK:
        report["interpretation"] = (
            "CALIBRATION FAILED: the irradiation contrasts themselves do not trace the "
            "DREM trajectory, so no flight result from this projection is interpretable.")
    elif (flt.get("peak_mean_rho", 0) < MIN_PEAK
          or flt.get("frac_argmax_near_peak", 0) < MIN_NEAR):
        report["interpretation"] = (
            f"Irradiation contrasts trace the trajectory with a clear peak "
            f"(mean rho {cal['peak_mean_rho']} at {cal['peak_timepoint_min']:g} min, "
            f"profile range {cal['profile_range']}, "
            f"{cal['frac_argmax_near_peak']:.0%} of argmaxes at or beside that peak). "
            f"Flight contrasts do not: their mean profile is essentially flat "
            f"(peak {flt.get('peak_mean_rho')}, range {flt.get('profile_range')}) and "
            f"only {flt.get('frac_argmax_near_peak', 0):.0%} of their best-matching "
            f"timepoints fall near it, so the argmax is arbitrary. Spaceflight is NOT a "
            f"time-shifted radiation response. The per-contrast shuffle test passes "
            f"{flt.get('n_matched_per_contrast_shuffle')} of {flt.get('n')} contrasts "
            f"and is too permissive to be read on its own.")
    else:
        report["interpretation"] = (
            f"Flight contrasts trace the trajectory with a peak of "
            f"{flt['peak_mean_rho']} at {flt['peak_timepoint_min']:g} min, with "
            f"{flt['frac_argmax_near_peak']:.0%} of argmaxes at or beside that peak.")

    write_json(RESULTS / "qc" / "drem_projection_qc.json", report)
    for lab, r in (("radiation", cal), ("flight", flt)):
        log(f"  {lab:<10} peak rho {r.get('peak_mean_rho')} at "
            f"{r.get('peak_timepoint_min')} min | profile range {r.get('profile_range')} "
            f"| argmax near peak {r.get('frac_argmax_near_peak')}")
    log(f"  -> {report['interpretation']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
