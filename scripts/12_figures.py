#!/usr/bin/env python3
"""Render the manuscript figures from results/.

  figures/fig1_cohort_design.pdf      timepoint coverage per study, both cohorts
  figures/fig2_sentinel_kinetics.pdf  DDR sentinel trajectories, WT vs sog1-1
  figures/fig3_latent_trajectory.pdf  auto-decoder latent code over pseudo-time
  figures/fig4_prior_ablation.pdf     TF attribution under the three priors

Every figure is also written to manuscript/figures/ as PDF, which is what main.tex
includes. A figure whose inputs are missing is skipped with a message rather than
rendered empty, so a missing panel is visible as a missing panel.

  python3 scripts/12_figures.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cohorts import COHORTS  # noqa: E402
from lib_sources import FIGURES, RESULTS, ROOT, log  # noqa: E402

MANU_FIGS = ROOT / "manuscript" / "figures"
plt.rcParams.update({"figure.dpi": 150, "savefig.bbox": "tight",
                     "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False})


def save(fig, name: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    MANU_FIGS.mkdir(parents=True, exist_ok=True)
    # PDF for LaTeX (vector), PNG for the Word build: pandoc embeds whatever file it
    # finds, and Word cannot render an embedded PDF — it shows a placeholder.
    for target in (FIGURES / f"{name}.pdf", MANU_FIGS / f"{name}.pdf"):
        fig.savefig(target)
    for target in (FIGURES / f"{name}.png", MANU_FIGS / f"{name}.png"):
        fig.savefig(target, dpi=200)
    plt.close(fig)
    log(f"  wrote {name}")


def fig_cohort_design() -> None:
    path = RESULTS / "qc" / "metadata_qc.json"
    if not path.exists():
        log("  fig1 skipped: metadata_qc.json missing")
        return
    md = pd.read_csv(ROOT / "data" / "metadata_master.csv")
    md = md[md["time_minutes"].notna()]

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))
    for ax, key in zip(axes, ("A", "B")):
        sub = md[md["cohort"] == key]
        studies = sorted(sub["accession"].unique(), key=lambda a: int(a.split("-")[1]))
        for y, acc in enumerate(studies):
            t = sorted(sub[sub["accession"] == acc]["time_minutes"].unique())
            n = [int((sub[(sub["accession"] == acc) & (sub["time_minutes"] == x)]).shape[0])
                 for x in t]
            ax.scatter(t, [y] * len(t), s=[12 * v for v in n],
                       color="#2b6cb0" if key == "A" else "#2f855a",
                       alpha=0.8, edgecolor="white", linewidth=0.5, zorder=3)
        ax.set_yticks(range(len(studies)))
        ax.set_yticklabels(studies)
        ax.set_xscale("log")
        ax.set_xlabel("time after treatment (min, log scale)")
        ax.grid(axis="x", alpha=0.25, zorder=0)
        ax.set_title(f"Cohort {key}: {COHORTS[key]['label']}", fontsize=9)
    fig.suptitle("Pseudo-time-series coverage (marker area = replicate count)",
                 fontsize=10)
    save(fig, "fig1_cohort_design")


def fig_sentinels() -> None:
    path = RESULTS / "pseudotimeseries" / "series_qc.json"
    if not path.exists():
        log("  fig2 skipped: series_qc.json missing")
        return
    q = json.loads(path.read_text())
    arms = q.get("cohorts", {}).get("A", {}).get("arms", {})
    panels = [("A_primary_WildType", "wild type"), ("A_primary_sog1-1", "sog1-1")]
    if not all(a in arms for a, _ in panels):
        log("  fig2 skipped: sentinel arms missing")
        return

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.2), sharey=True)
    for ax, (slug, label) in zip(axes, panels):
        genes = arms[slug]["sentinels"]["genes"]
        for name, v in sorted(genes.items()):
            traj = v.get("trajectory_log2fc")
            if not traj:
                continue
            xs = sorted(float(k) for k in traj)
            ys = [traj[f"{x:g}"] for x in xs]
            ax.plot(xs, ys, marker="o", ms=3, lw=1.2, label=name)
        ax.set_xscale("log")
        ax.axhline(0, color="0.6", lw=0.8, zorder=0)
        ax.set_xlabel("time after irradiation (min)")
        ax.set_title(label, fontsize=9)
    axes[0].set_ylabel(r"$\log_2$ fold-change vs control")
    axes[1].legend(fontsize=7, ncol=2, frameon=False, loc="upper left")
    fig.suptitle("DNA-damage sentinel kinetics: the response is SOG1-dependent",
                 fontsize=10)
    save(fig, "fig2_sentinel_kinetics")


def fig_latent() -> None:
    files = sorted((RESULTS / "celltypes").glob("*_latent.tsv"))
    if not files:
        log("  fig3 skipped: no latent trajectories")
        return
    fig, ax = plt.subplots(figsize=(6, 3.4))
    cmap = plt.get_cmap("viridis")
    for i, path in enumerate(files):
        d = pd.read_csv(path, sep="\t", index_col=0).sort_index(
            key=lambda ix: [float(x) for x in ix])
        # Distance travelled in latent space from the earliest timepoint: one
        # interpretable scalar out of 32 dimensions.
        base = d.iloc[0].to_numpy()
        dist = np.linalg.norm(d.to_numpy() - base, axis=1)
        ax.plot([float(x) for x in d.index], dist, marker="o", ms=3.5, lw=1.3,
                color=cmap(i / max(len(files) - 1, 1)),
                label=path.name.replace("_latent.tsv", ""))
    ax.set_xscale("log")
    ax.set_xlabel("pseudo-time (min, log scale)")
    ax.set_ylabel("latent displacement from t$_0$")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title("Auto-decoder latent trajectory over pseudo-time", fontsize=10)
    save(fig, "fig3_latent_trajectory")


def fig_ablation() -> None:
    path = RESULTS / "comparison" / "prior_ablation.tsv"
    if not path.exists():
        log("  fig4 skipped: prior_ablation.tsv missing (run 10_compare_models.py)")
        return
    d = pd.read_csv(path, sep="\t")
    if d.empty:
        log("  fig4 skipped: ablation table empty")
        return
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    arms = sorted(d["arm"].unique())
    counts = [(int(d[d["arm"] == a]["attribution_changed"].sum()),
               int((d["arm"] == a).sum())) for a in arms]
    pct = [100 * c / t if t else 0 for c, t in counts]
    # Shade by how many bifurcations the percentage rests on. Three of these arms
    # resolved a single split, and a bare "100%" from one split reads as though it
    # carried the same weight as 50% from six.
    colors = ["#2b6cb0" if t >= 4 else "#9db8d2" for _, t in counts]
    ax.barh(arms, pct, color=colors, height=0.55)
    for y, (v, (c, t)) in enumerate(zip(pct, counts)):
        ax.text(v + 1.5, y, f"{c}/{t}", va="center", fontsize=8)
    ax.set_xlabel("% of bifurcations whose top-ranked TF changes with the prior")
    ax.set_xlim(0, 118)
    ax.set_title("Effect of the cell-type-weighted prior on TF attribution", fontsize=10)
    ax.text(0.99, -0.34, "labels are changed/total splits; pale bars rest on <4 splits",
            transform=ax.transAxes, ha="right", fontsize=7, color="0.35")
    save(fig, "fig4_prior_ablation")


def fig_tf_flight() -> None:
    """Figures 7-10: the spaceflight TF analysis."""
    act_dir = RESULTS / "tf_activity"
    if not (act_dir / "tf_activity.tsv").exists():
        log("  figs 7-10 skipped: run 19-22 first")
        return
    act = pd.read_csv(act_dir / "tf_activity.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(act_dir / "contrast_meta.tsv", sep="\t", index_col=0)
    stats_df = pd.read_csv(act_dir / "tf_statistics.tsv", sep="\t")

    # ---- Fig 7: TF activity heatmap, contrasts x most variable TFs
    sub = meta[meta["is_flight"] | meta["is_radiation"]].copy()
    sub = sub.sort_values(["is_radiation", "mission", "accession"])
    top = act.loc[sub.index].std().nlargest(40).index
    M = act.loc[sub.index, top].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(11, 6.4), constrained_layout=True)
    v = np.nanpercentile(np.abs(M), 98)
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-v, vmax=v)
    ax.set_xticks(range(len(top)))
    ax.set_xticklabels(top, rotation=90, fontsize=5)
    ax.set_yticks(range(len(sub)))
    ax.set_yticklabels([f"{r.accession} {str(r.level)[:22]}" for r in sub.itertuples()],
                       fontsize=5)
    # Mark where radiation contrasts begin, so the two blocks are visually separable.
    nrad = int((~sub["is_radiation"]).sum())
    ax.axhline(nrad - 0.5, color="black", lw=1.2)
    ax.text(len(top) * 0.5, nrad - 1.2, "flight  /  radiation", ha="center",
            fontsize=7, color="black")
    fig.colorbar(im, ax=ax, shrink=0.6, label="TF activity (z)")
    ax.set_title("TF activity across flight and radiation contrasts "
                 "(40 most variable TFs)", fontsize=10)
    save(fig, "fig7_tf_activity_heatmap")

    # ---- Fig 8: per-TF flight effect vs flight-vs-radiation difference
    fig, ax = plt.subplots(figsize=(6.6, 4.4), constrained_layout=True)
    d = stats_df.dropna(subset=["flight_mean_z", "radiation_mean_z"])
    sig = d["flight_vs_radiation_q"] < 0.05
    ax.scatter(d.loc[~sig, "radiation_mean_z"], d.loc[~sig, "flight_mean_z"],
               s=14, color="#a0aec0", alpha=0.7, label="not different (q>=0.05)")
    ax.scatter(d.loc[sig, "radiation_mean_z"], d.loc[sig, "flight_mean_z"],
               s=20, color="#c53030", alpha=0.85, label="differs, q<0.05")
    for r in d[d["name"].astype(str).str.len() > 0].itertuples():
        ax.annotate(r.name, (r.radiation_mean_z, r.flight_mean_z), fontsize=6.5,
                    xytext=(4, 3), textcoords="offset points")
    lim = float(np.nanmax(np.abs(d[["radiation_mean_z", "flight_mean_z"]].to_numpy()))) * 1.1
    ax.plot([-lim, lim], [-lim, lim], ls=":", lw=1, color="0.6", zorder=0)
    ax.axhline(0, color="0.85", lw=0.8, zorder=0)
    ax.axvline(0, color="0.85", lw=0.8, zorder=0)
    ax.set_xlabel("TF activity under irradiation (mean z)")
    ax.set_ylabel("TF activity in spaceflight (mean z, per mission)")
    ax.set_title("The same TFs respond to radiation but not to flight", fontsize=10)
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    save(fig, "fig8_tf_flight_vs_radiation")

    # ---- Fig 9: classifier
    rep_path = act_dir / "classifier_report.json"
    if rep_path.exists():
        rep = json.loads(rep_path.read_text())
        fig, ax = plt.subplots(figsize=(6.2, 3.8), constrained_layout=True)
        labels = ["flight vs\nradiation\n(logistic)", "random\nforest",
                  "permutation\nnull", "platform\ncontrol"]
        vals = [rep.get("auc_logistic"), rep.get("auc_random_forest"),
                rep.get("null_mean"), rep.get("platform_control_auc")]
        cols = ["#2b6cb0", "#4a7fb5", "#a0aec0", "#805ad5"]
        bars = ax.bar(labels, [v or 0 for v in vals], color=cols, width=0.6)
        if rep.get("null_sd"):
            ax.errorbar(2, rep["null_mean"], yerr=rep["null_sd"], color="0.3",
                        capsize=4, lw=1.2)
        ax.axhline(0.5, ls="--", lw=1, color="0.5")
        for b, v in zip(bars, vals):
            if v is not None:
                ax.text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.3f}",
                        ha="center", fontsize=8)
        ax.set_ylim(0, 1.12)
        ax.set_ylabel("AUC (leave-one-mission-out)")
        ax.set_title(f"Flight and radiation are separable; platform is not a confound\n"
                     f"(p = {rep.get('p_permutation')})", fontsize=9.5)
        save(fig, "fig9_classifier")

    # ---- Fig 10: DREM trajectory projection
    proj_path = act_dir / "drem_projection.tsv"
    if proj_path.exists():
        pr = pd.read_csv(proj_path, sep="\t")
        rho_cols = [c for c in pr.columns if c.startswith("rho_")]
        times = [float(c[4:]) for c in rho_cols]
        fig, ax = plt.subplots(figsize=(6.6, 4.0), constrained_layout=True)
        for lab, mask, col in (("irradiation", pr["is_radiation"], "#c53030"),
                               ("spaceflight", pr["is_flight"], "#2b6cb0")):
            sub2 = pr[mask]
            if sub2.empty:
                continue
            m = [sub2[c].mean() for c in rho_cols]
            sd = [sub2[c].std() / max(np.sqrt(len(sub2)), 1) for c in rho_cols]
            ax.plot(times, m, marker="o", ms=4, lw=1.6, color=col, label=lab)
            ax.fill_between(times, np.array(m) - np.array(sd), np.array(m) + np.array(sd),
                            color=col, alpha=0.18)
        ax.set_xscale("log")
        ax.axhline(0, color="0.7", lw=0.8)
        ax.set_xlabel("DREM trajectory timepoint (min, log scale)")
        ax.set_ylabel("Spearman rho with the timepoint's TF profile")
        ax.set_title("Irradiation traces the radiation trajectory; spaceflight is flat",
                     fontsize=10)
        ax.legend(fontsize=8, frameon=False)
        save(fig, "fig10_drem_projection")


def fig_dose_response() -> None:
    """Figure 11: the detection floor against the ISS dose."""
    qc = RESULTS / "qc" / "dose_response_qc.json"
    dr_path = RESULTS / "decoder" / "dose_response.tsv"
    if not (qc.exists() and dr_path.exists()):
        log("  fig11 skipped: run 23_dose_response.py first")
        return
    rep = json.loads(qc.read_text())
    dr = pd.read_csv(dr_path, sep="\t").sort_values("dose_cgy")

    fig, ax = plt.subplots(figsize=(7.4, 4.4), constrained_layout=True)
    iss_lo, iss_hi = rep["iss"]["mission_dose_range_cgy"]
    ax.axvspan(iss_lo, iss_hi, color="#2b6cb0", alpha=0.16, zorder=0)
    ax.text(np.sqrt(iss_lo * iss_hi), ax.get_ylim()[1], " ISS missions\n (11-70 d)",
            fontsize=7.5, color="#2b6cb0", ha="center", va="top")

    br = rep["observed_bracket_cgy"]
    ax.axvspan(br["highest_dose_not_detected"], br["lowest_dose_detected"],
               color="0.85", alpha=0.6, zorder=0)
    ax.axhline(rep["z_threshold"], ls="--", lw=1, color="#c53030")
    ax.text(dr["dose_cgy"].max(), rep["z_threshold"] + 0.35, "z = 1.96",
            fontsize=7.5, color="#c53030", ha="right")
    ax.axhline(0, color="0.85", lw=0.8, zorder=0)

    for arm, col, lab in (("sog1_arm", "#c53030", "SOG1 / DNA-damage arm"),
                          ("myb3r_arm", "#2f855a", "MYB3R / G2-M arm")):
        ax.plot(dr["dose_cgy"], dr[arm], marker="o", ms=6, lw=1.8, color=col, label=lab)
    for r in dr.itertuples():
        ax.annotate(r.accession.replace("OSD-", ""), (r.dose_cgy, r.sog1_arm),
                    fontsize=6, xytext=(0, -12), textcoords="offset points",
                    ha="center", color="0.4")

    ax.set_xscale("log")
    ax.set_xlim(iss_lo * 0.5, dr["dose_cgy"].max() * 1.6)
    ax.set_xlabel("absorbed dose (cGy, log scale)")
    ax.set_ylabel("arm activity (z)")
    ax.set_title("The assay is blind at ISS doses: the DNA-damage arm needs "
                 f"{br['lowest_dose_detected']:g} cGy,\nISS delivers "
                 f"{iss_lo:.2f}-{iss_hi:.2f} cGy", fontsize=9.5)
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    save(fig, "fig11_dose_response")


def fig_radiation_quality() -> None:
    """Figure 12: the matched HZE-vs-gamma pair against the gamma-vs-gamma baseline."""
    qc = RESULTS / "qc" / "radiation_quality_qc.json"
    cpath = RESULTS / "radiation_quality" / "quality_correlations.tsv"
    if not (qc.exists() and cpath.exists()):
        log("  fig12 skipped: run 24_radiation_quality.py first")
        return
    rep = json.loads(qc.read_text())
    c = pd.read_csv(cpath, sep="\t")
    r = c[c["scope"] == "responsive_tfs"]
    base = r[(~r["same_study"]) & (r["same_quality"])
             & (r["quality_a"] == "low-LET photon")]["spearman_rho"]
    matched = rep["matched_test"]["rho_responsive_tfs"]

    fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
    ax.hist(base, bins=10, range=(0, 0.75), color="#a0aec0", alpha=0.85,
            label=f"gamma vs gamma, different studies (n={len(base)})")
    ax.axvline(matched, color="#c53030", lw=2.2,
               label=f"matched gamma vs Fe-56 HZE, same study ({matched:.2f})")
    ax.axvline(float(base.mean()), color="0.35", ls="--", lw=1.2,
               label=f"baseline mean ({base.mean():.2f})")
    ax.set_xlabel("Spearman correlation of the full TF-activity profile")
    ax.set_ylabel("number of study pairs")
    pct = rep["matched_pair_vs_baseline"]["percentile"]
    pv = rep["matched_pair_vs_baseline"]["one_sided_p"]
    ax.set_title("Radiation quality does not explain the flight null:\n"
                 f"the matched HZE pair sits at the {pct:g}th percentile of the "
                 f"gamma baseline (p = {pv})", fontsize=9.5)
    ax.legend(fontsize=7, frameon=False)
    save(fig, "fig12_radiation_quality")


def main() -> int:
    log("rendering figures")
    fig_cohort_design()
    fig_sentinels()
    fig_latent()
    fig_ablation()
    fig_tf_flight()
    fig_dose_response()
    fig_radiation_quality()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
