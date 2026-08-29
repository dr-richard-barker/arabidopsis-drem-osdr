#!/usr/bin/env python3
"""Turn the signature scores into a calibrated radiation decoder, and report predictions.

Scoring one gene set is not enough. Several non-irradiation contrasts elevate the
SOG1-dependent set --- a constitutive DNA-damage state raises those genes too --- so a
decoder built on that set alone would call DNA-damage mutants "irradiated".

What separates acute irradiation is the *joint* pattern, and it is the pattern the source
study described: SOG1 activates its targets while the MYB3R repressors shut down the G2/M
programme. Genuine irradiation therefore shows SOG1-dependent genes UP and MYB3R targets
DOWN simultaneously. A constitutive DDR state shows the first without the second.

    radiation_index = min( z(sog1_dependent),  -z(myb3r_repressed) )

The index is a CONJUNCTION, deliberately. The obvious alternative --- the difference
z(sog1_dependent) - z(myb3r_repressed) --- was tried first and rejected: it scored six
spaceflight and altered-gravity contrasts above the threshold on the strength of the
MYB3R arm alone, with SOG1-dependent z of roughly zero. Spaceflight represses G2/M genes
because growth slows, which has nothing to do with DNA damage, and a difference lets
either arm carry the score by itself. Taking the minimum requires both arms to move, and
reports the weaker of the two required moves, so no contrast can qualify on one arm.

The threshold is statistical, not fitted. Both arms are z-scores against a permutation
null, so requiring the weaker arm to clear z >= 1.96 asks that each arm be individually
significant at two-sided p < 0.05. Fitting a threshold to the labelled data was tried and
abandoned: a midpoint between the labelled-positive minimum and the negative maximum is
hostage to a single mislabelled example at each end, and OSDR's positives include an
exposure -- OSD-782 at 10 cGy -- two to three orders of magnitude below the dose the
signature was derived at, which genuinely produces no response. Calibrating on it drags
the threshold into the noise. A fixed statistical cutoff cannot be gamed by one outlier,
and the labelled contrasts are then reported as a *test* of the decoder rather than as
its calibration set.

  results/decoder/predictions.tsv    every contrast, ranked, with its call
  results/decoder/decoder_report.json  calibration, separation and the predictions
  figures/fig5_decoder.pdf/.png      the fingerprint plot

  python3 scripts/17_decoder_report.py
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import FIGURES, RESULTS, ROOT, log, write_json  # noqa: E402

OUT = RESULTS / "decoder"
MANU_FIGS = ROOT / "manuscript" / "figures"
SCORES = OUT / "study_scores.tsv"

plt.rcParams.update({"figure.dpi": 150, "savefig.bbox": "tight", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})

# A contrast is "known irradiation" when its own OSDR factor names ionizing radiation or
# an absorbed dose. Taken from the metadata, never from a curated accession list, so the
# calibration cannot be quietly tuned by editing a set of study IDs.
RADIATION_FACTORS = ("ionizing radiation", "absorbed radiation dose")


def load() -> pd.DataFrame:
    if not SCORES.exists():
        raise SystemExit("run 16_scan_osdr_plants.py first")
    d = pd.read_csv(SCORES, sep="\t")
    p = d.pivot_table(index=["accession", "factor", "level"],
                      columns="gene_set", values="z").reset_index()
    p.columns.name = None
    return p


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-index", type=float, default=1.96,
                    help="index threshold: the weaker arm's minimum z (default 1.96, "
                         "i.e. each arm individually significant at two-sided p<0.05)")
    args = ap.parse_args()

    d = load()
    for col in ("sog1_dependent", "myb3r_repressed"):
        if col not in d:
            raise SystemExit(f"missing gene set {col} in the scores")

    # Conjunction, not difference — see the module docstring for why.
    d["radiation_index"] = np.minimum(d["sog1_dependent"], -d["myb3r_repressed"])
    d["sog1_arm"] = d["sog1_dependent"]
    d["myb3r_arm"] = -d["myb3r_repressed"]
    d["is_radiation_factor"] = d["factor"].isin(RADIATION_FACTORS)
    # A mutant-vs-wildtype contrast is not an exposure. sog1-1 in an irradiation study
    # is the negative control, and calling it "irradiated" would be a category error.
    d["is_genotype_contrast"] = d["factor"] == "genotype"

    pos = d[d["is_radiation_factor"]]["radiation_index"]
    neg = d[~d["is_radiation_factor"] & ~d["is_genotype_contrast"]]["radiation_index"]
    threshold = float(args.min_index)

    # Four calls, not two. Splitting the "no signal" bucket matters because most
    # spaceflight contrasts are not featureless: they repress the G2/M arm without
    # activating the SOG1 arm, which is proliferation slowing rather than DNA damage.
    # Collapsing that into "no signal" would discard a real and consistent pattern, and
    # a decoder that used the difference statistic would instead have mislabelled it as
    # radiation.
    radiation_like = d["radiation_index"] >= threshold
    inverse = np.minimum(-d["sog1_dependent"], d["myb3r_repressed"]) >= threshold
    proliferation_only = (~radiation_like & ~inverse
                          & (d["myb3r_arm"] >= threshold)
                          & (d["sog1_arm"] < threshold))
    d["call"] = np.select(
        [radiation_like, inverse, proliferation_only],
        ["radiation-like", "inverse (DDR suppressed)",
         "G2/M repressed, no DDR (proliferation slowing)"],
        default="no signal")
    d["novel_prediction"] = (d["call"] == "radiation-like") & ~d["is_radiation_factor"]

    d = d.sort_values("radiation_index", ascending=False)
    keep = ["accession", "factor", "level", "sog1_arm", "myb3r_arm",
            "sog1_dependent", "myb3r_repressed", "ddr_core", "phase_early", "phase_mid",
            "radiation_index", "is_radiation_factor", "call", "novel_prediction"]
    d[keep].round(3).to_csv(OUT / "predictions.tsv", sep="\t", index=False)

    # Separation: does the index rank every known-irradiation contrast above the rest?
    labels = d["is_radiation_factor"].to_numpy()
    scores = d["radiation_index"].to_numpy()
    order = np.argsort(-scores)
    y = labels[order]
    n_pos, n_neg = int(y.sum()), int((~y).sum())
    # AUC by rank-sum, over exposure contrasts only (genotype contrasts excluded above
    # from calibration but retained here so the negative controls are visible).
    ranks = np.arange(1, len(y) + 1)
    auc = ((ranks[y].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)) if n_pos and n_neg else float("nan")
    auc = 1.0 - auc  # ranks are descending, so invert

    report = {
        "generated": dt.date.today().isoformat(),
        "index": "min( z(sog1_dependent), -z(myb3r_repressed) )",
        "index_rejected": {
            "form": "z(sog1_dependent) - z(myb3r_repressed)",
            "why_rejected": ("Scored six spaceflight/altered-gravity contrasts above "
                             "threshold on the MYB3R arm alone, with SOG1-dependent z "
                             "near zero. Spaceflight represses G2/M genes because growth "
                             "slows, not because DNA is damaged; a difference lets either "
                             "arm carry the score unaided."),
        },
        "rationale": ("Acute irradiation moves the SOG1-activated and MYB3R-repressed "
                      "arms in opposite directions. A constitutive DNA-damage state "
                      "raises the SOG1 arm without repressing the G2/M arm, so the "
                      "difference separates exposure from chronic stress."),
        "threshold": {
            "value": threshold,
            "basis": "fixed statistical cutoff: both arms individually significant at "
                     "two-sided p<0.05 (z>=1.96). Not fitted to these data.",
        },
        "performance": {
            "positive_definition": f"factor in {RADIATION_FACTORS}",
            "n_labelled_positive": int(len(pos)),
            "n_other_exposure": int(len(neg)),
            "labelled_positive_index_range": [round(float(pos.min()), 2),
                                              round(float(pos.max()), 2)],
            "other_exposure_index_range": [round(float(neg.min()), 2),
                                           round(float(neg.max()), 2)],
            "labelled_positives_detected": int((pos >= threshold).sum()),
            "labelled_positives_missed": int((pos < threshold).sum()),
            "other_exposures_above_threshold": int((neg >= threshold).sum()),
            "auc_vs_all_contrasts": round(float(auc), 4),
        },
        "known_radiation_ranked": [
            {"accession": r.accession, "level": r.level,
             "index": round(r.radiation_index, 2), "call": r.call}
            for r in d[d["is_radiation_factor"]].itertuples()],
        "novel_predictions": [
            {"accession": r.accession, "factor": r.factor, "level": r.level,
             "index": round(r.radiation_index, 2),
             "sog1_arm_z": round(r.sog1_arm, 2),
             "myb3r_arm_z": round(r.myb3r_arm, 2)}
            for r in d[d["novel_prediction"]].itertuples()],
        "call_counts": None,  # filled below
        "proliferation_only": [
            {"accession": r.accession, "factor": r.factor, "level": r.level,
             "sog1_arm_z": round(r.sog1_arm, 2), "myb3r_arm_z": round(r.myb3r_arm, 2)}
            for r in d[d["call"].str.startswith("G2/M")].itertuples()],
        "negative_controls": [
            {"accession": r.accession, "level": r.level,
             "index": round(r.radiation_index, 2), "call": r.call}
            for r in d[d["is_genotype_contrast"]
                       & d["level"].str.contains("sog1|myb3r", case=False, na=False)].itertuples()],
    }
    report["call_counts"] = {k: int(v) for k, v in d["call"].value_counts().items()}
    write_json(OUT / "decoder_report.json", report)

    # ------------------------------------------------------------------ figure
    fig = plt.figure(figsize=(12.2, 4.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15])
    ax = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    grp = np.where(d["is_radiation_factor"], "irradiation (labelled)",
                   np.where(d["is_genotype_contrast"], "genotype contrast", "other exposure"))
    colours = {"irradiation (labelled)": "#c53030", "genotype contrast": "#805ad5",
               "other exposure": "#2b6cb0"}
    lim = max(abs(d["sog1_dependent"]).max(), abs(d["myb3r_repressed"]).max()) * 1.12

    # Decision region: SOG1 arm up AND MYB3R arm down, both past threshold.
    ax.add_patch(plt.Rectangle((threshold, -lim), lim - threshold, lim - threshold,
                               facecolor="#c53030", alpha=0.07, zorder=0))
    ax.plot([threshold, threshold], [-lim, -threshold], ls="--", lw=1, color="#c53030", zorder=1)
    ax.plot([threshold, lim], [-threshold, -threshold], ls="--", lw=1, color="#c53030", zorder=1)
    ax.axhline(0, color="0.88", lw=0.8, zorder=0)
    ax.axvline(0, color="0.88", lw=0.8, zorder=0)
    for g, c in colours.items():
        m = grp == g
        ax.scatter(d.loc[m, "sog1_dependent"], d.loc[m, "myb3r_repressed"],
                   s=36, alpha=0.85, label=g, color=c, edgecolor="white", linewidth=0.6,
                   zorder=3)
    # Label only the extremes, offset apart so they do not collide.
    for r, dy in zip(d.head(4).itertuples(), (10, -4, -14, -24)):
        ax.annotate(r.accession, (r.sog1_dependent, r.myb3r_repressed), fontsize=6.5,
                    xytext=(-34, dy), textcoords="offset points", color="#742a2a")
    for r in d.tail(2).itertuples():
        ax.annotate(f"{r.accession} {r.level}", (r.sog1_dependent, r.myb3r_repressed),
                    fontsize=6.5, xytext=(8, 2), textcoords="offset points", color="#553c9a")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("SOG1-dependent set  (z)")
    ax.set_ylabel("MYB3R-repressed set  (z)")
    ax.set_title("Both arms must move; shaded corner is the decision region", fontsize=9.5)
    ax.legend(fontsize=7, frameon=False, loc="lower left")

    top = d.head(14).iloc[::-1]
    bar_c = ["#c53030" if p_ else "#2b6cb0" for p_ in top["is_radiation_factor"]]
    ax2.barh(range(len(top)), top["radiation_index"], color=bar_c, height=0.7, zorder=2)
    ax2.axvline(threshold, ls="--", lw=1, color="#c53030", zorder=3)
    ax2.set_yticks([])
    # Labels inside the axes, so they cannot push into the neighbouring panel.
    xmax = float(top["radiation_index"].max()) * 1.05
    for y, r in enumerate(top.itertuples()):
        ax2.text(xmax * 0.012, y, f"{r.accession} · {str(r.level)[:30]}",
                 va="center", ha="left", fontsize=6.8,
                 color="white" if r.radiation_index > xmax * 0.42 else "0.25", zorder=4)
    ax2.set_xlim(0, xmax)
    ax2.set_xlabel("radiation index   min(z SOG1-dep, −z MYB3R-rep)")
    ax2.set_title(f"Top 14 contrasts (dashed = threshold {threshold:g})", fontsize=9.5)

    FIGURES.mkdir(parents=True, exist_ok=True)
    MANU_FIGS.mkdir(parents=True, exist_ok=True)
    for t in (FIGURES / "fig5_decoder.pdf", MANU_FIGS / "fig5_decoder.pdf"):
        fig.savefig(t)
    for t in (FIGURES / "fig5_decoder.png", MANU_FIGS / "fig5_decoder.png"):
        fig.savefig(t, dpi=200)
    plt.close(fig)

    perf = report["performance"]
    log(f"threshold {threshold} (fixed, z>=1.96 on both arms)   AUC {auc:.3f}")
    log(f"  labelled irradiation : {perf['labelled_positives_detected']}"
        f"/{perf['n_labelled_positive']} detected, range {perf['labelled_positive_index_range']}")
    log(f"  other exposures      : {perf['other_exposures_above_threshold']}"
        f"/{perf['n_other_exposure']} above threshold, range {perf['other_exposure_index_range']}")
    for k, v in report["call_counts"].items():
        log(f"  {v:>3}  {k}")
    log(f"  novel radiation predictions: {len(report['novel_predictions'])}")
    for p_ in report["novel_predictions"]:
        log(f"    {p_['accession']:<9} {p_['factor']}={p_['level']}  index={p_['index']}")
    log("  wrote predictions.tsv, decoder_report.json, fig5_decoder")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
