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


def main() -> int:
    log("rendering figures")
    fig_cohort_design()
    fig_sentinels()
    fig_latent()
    fig_ablation()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
