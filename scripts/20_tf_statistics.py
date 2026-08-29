#!/usr/bin/env python3
"""Which DREM transcription factors move in spaceflight? Per-TF tests with FDR.

For every TF, the question is whether its activity shifts consistently across flight
contrasts. Two things make that harder than a t-test over 41 contrasts.

**Studies are not independent.** SpaceX-5 contributes five studies and SpaceX-2 four;
their samples flew on the same vehicle, through the same launch and the same ground
handling. Counting them as nine independent observations would inflate every statistic.
So contrasts are averaged within a mission first, and the test runs over MISSIONS. The
effective n is 14, not 41, and that is the number the manuscript reports.

**Direction matters and cancels.** A TF that goes up in half the missions and down in the
other half has a mean near zero, which is a real answer, not a null one. Both the mean
mission z and the fraction of missions agreeing in sign are reported, so a consistent
weak shift is distinguishable from an inconsistent strong one.

Tests, all Benjamini-Hochberg corrected across the TFs tested:
  - flight vs zero: one-sample t over mission means
  - flight vs radiation: Welch t between the two contrast classes, which asks whether a
    TF behaves differently in flight than under the irradiation the model was built on
  - gravity dose-response: Spearman rho against g-level within OSD-251, the one study
    with an ordered gravity series and an onboard 1G control

  results/tf_activity/tf_statistics.tsv    per TF: effects, p, q, sign consistency
  results/qc/tf_statistics_qc.json         counts, and the TFs surviving FDR

  python3 scripts/20_tf_statistics.py
"""
from __future__ import annotations

import argparse
import datetime as dt
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import RESULTS, log, quiet_accelerate_blas_warnings, write_json  # noqa: E402

quiet_accelerate_blas_warnings()

ACT = RESULTS / "tf_activity"
KEY_TFS = {"AT1G25580": "SOG1", "AT4G32730": "MYB3R1", "AT3G09370": "MYB3R3",
           "AT5G11510": "MYB3R4", "AT2G46770": "ANAC043", "AT5G13330": "ERF115"}


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    x = p[ok]
    if x.size == 0:
        return q
    order = np.argsort(x)
    ranked = x[order]
    n = x.size
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def g_level(level: str) -> float | None:
    """'0.36G by centrifugation' -> 0.36; 'uG' -> 0.0."""
    s = str(level)
    if re.fullmatch(r"\s*uG\s*", s, re.I):
        return 0.0
    m = re.match(r"\s*([\d.]+)\s*G\b", s)
    return float(m.group(1)) if m else None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alpha", type=float, default=0.05)
    args = ap.parse_args()

    act = pd.read_csv(ACT / "tf_activity.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(ACT / "contrast_meta.tsv", sep="\t", index_col=0)
    tfs = list(act.columns)

    flight = meta[meta["is_flight"]]
    rad = meta[meta["is_radiation"]]
    log(f"{len(flight)} flight contrasts over {flight['mission'].nunique()} missions; "
        f"{len(rad)} radiation contrasts")

    # Collapse to one value per mission before testing.
    fl_mission = act.loc[flight.index].groupby(flight["mission"]).mean()
    rad_vals = act.loc[rad.index]
    log(f"  collapsed to {fl_mission.shape[0]} mission-level observations")

    rows = []
    for tf in tfs:
        x = fl_mission[tf].dropna().to_numpy()
        if len(x) >= 3:
            t, p = stats.ttest_1samp(x, 0.0)
            mean_z, frac = float(x.mean()), float(np.mean(np.sign(x) == np.sign(x.mean())))
        else:
            t = p = mean_z = frac = np.nan

        r = rad_vals[tf].dropna().to_numpy()
        if len(x) >= 3 and len(r) >= 3:
            t2, p2 = stats.ttest_ind(x, r, equal_var=False)
            rad_mean = float(r.mean())
        else:
            t2 = p2 = rad_mean = np.nan

        rows.append({"TF": tf, "name": KEY_TFS.get(tf, ""),
                     "n_missions": int(len(x)),
                     "flight_mean_z": round(mean_z, 4) if mean_z == mean_z else None,
                     "sign_consistency": round(frac, 3) if frac == frac else None,
                     "flight_t": round(float(t), 3) if t == t else None,
                     "flight_p": float(p) if p == p else np.nan,
                     "radiation_mean_z": round(rad_mean, 4) if rad_mean == rad_mean else None,
                     "flight_vs_radiation_t": round(float(t2), 3) if t2 == t2 else None,
                     "flight_vs_radiation_p": float(p2) if p2 == p2 else np.nan})

    df = pd.DataFrame(rows)
    df["flight_q"] = bh(df["flight_p"].to_numpy())
    df["flight_vs_radiation_q"] = bh(df["flight_vs_radiation_p"].to_numpy())

    # Gravity dose-response inside the one study that has an ordered g-series with an
    # onboard 1G control. This distinguishes a graded response from an isolated blip --
    # OSD-346's 0.38G contrast clears a threshold while its neighbours on both sides do
    # not, which is what an isolated blip looks like.
    dose_rows = []
    for acc in meta["accession"].unique():
        sub = meta[(meta["accession"] == acc) & (meta["factor"] == "altered gravity")]
        if len(sub) < 4:
            continue
        g = sub["level"].map(g_level)
        keep = sub.index[g.notna()]
        if len(keep) < 4:
            continue
        gv = g.loc[keep].to_numpy(dtype=float)
        for tf in tfs:
            y = act.loc[keep, tf].to_numpy(dtype=float)
            if np.isfinite(y).sum() < 4:
                continue
            rho, p = stats.spearmanr(gv, y)
            dose_rows.append({"accession": acc, "TF": tf, "name": KEY_TFS.get(tf, ""),
                              "n_levels": int(len(gv)), "spearman_rho": round(float(rho), 3),
                              "p": float(p)})
    dose = pd.DataFrame(dose_rows)
    if not dose.empty:
        dose["q"] = np.concatenate([bh(dose.loc[dose["accession"] == a, "p"].to_numpy())
                                    for a in dose["accession"].unique()])
        dose.sort_values("p").to_csv(ACT / "gravity_dose_response.tsv", sep="\t", index=False)

    df = df.sort_values("flight_p")
    df.to_csv(ACT / "tf_statistics.tsv", sep="\t", index=False)

    sig = df[df["flight_q"] < args.alpha]
    diff = df[df["flight_vs_radiation_q"] < args.alpha]
    write_json(RESULTS / "qc" / "tf_statistics_qc.json", {
        "generated": dt.date.today().isoformat(),
        "design": "contrasts averaged within mission; tests run over missions, so the "
                  "effective n is the number of missions, not contrasts",
        "n_tfs_tested": int(len(df)),
        "n_flight_contrasts": int(len(flight)),
        "n_missions": int(fl_mission.shape[0]),
        "n_radiation_contrasts": int(len(rad)),
        "n_significant_flight": int(len(sig)),
        "n_significant_flight_vs_radiation": int(len(diff)),
        "alpha": args.alpha,
        "top_flight_tfs": sig.head(15)[
            ["TF", "name", "flight_mean_z", "sign_consistency", "flight_q"]
        ].to_dict("records"),
        "key_tfs": df[df["TF"].isin(KEY_TFS)][
            ["TF", "name", "flight_mean_z", "radiation_mean_z", "flight_q",
             "flight_vs_radiation_q"]].to_dict("records"),
        "gravity_dose_response_top": (
            dose.nsmallest(10, "p")[["accession", "TF", "name", "spearman_rho", "q"]]
            .to_dict("records") if not dose.empty else []),
    })

    log(f"  TFs significant in flight (q<{args.alpha}): {len(sig)}/{len(df)}")
    log(f"  TFs differing flight vs radiation (q<{args.alpha}): {len(diff)}/{len(df)}")
    log("\n  key DREM regulators:")
    for r in df[df["TF"].isin(KEY_TFS)].itertuples():
        log(f"    {r.name or r.TF:<9} flight z={r.flight_mean_z:<8} "
            f"radiation z={r.radiation_mean_z:<8} q={r.flight_q:.3g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
