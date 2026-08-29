#!/usr/bin/env python3
"""Where does the radiation response live? Tissue attribution for the signature and the prior.

Three analyses, in increasing order of how much they can be trusted:

1. **Deconvolution stability (a caveat, measured).** The per-timepoint NNLS fractions are
   NOT stable. The atlas has 183 partly collinear clusters, so one bulk profile admits
   several near-equivalent sparse fits and small changes flip which one wins:
   `rosette_21d_10` swings between 0.00 and 0.47 across adjacent timepoints.

   Aggregating to organ level was the obvious remedy and it does NOT work --- measured
   here, the median relative swing is essentially unchanged (organ 0.61 vs cluster 0.60),
   because the flipping happens between clusters of the same organ as readily as across
   organs. That is reported rather than quietly dropped, because the tempting move is to
   present organ-level fractions as though coarser meant steadier.

   What IS stable is the TIME-AVERAGED composition, which is what the weighting actually
   consumes: averaging over eight timepoints cancels the flip-flopping, and a gene-level
   bootstrap measures how much. So per-timepoint fractions are never interpreted as
   compositional dynamics anywhere in this pipeline, and tissue attribution comes from
   analysis 2, which does not use the deconvolution at all.

2. **Signature-gene tissue enrichment (robust).** Where the deconvolution is a fit, this
   is a lookup: for each atlas cluster, is the SOG1-dependent gene set more highly
   expressed there than genes generally? No NNLS, no collinearity, no per-timepoint
   instability --- just the atlas's own expression matrix, tested against a permutation
   null. This is the analysis that answers "which tissues does the response belong to".

3. **What the weighting actually did.** For each TF, the ratio of its weighted to its
   binding-only edge mass, paired with the cluster where that TF is most expressed. This
   turns the lever from an abstract multiplication into a list of which regulators were
   promoted or demoted, and on account of which tissue.

  results/celltypes/tissue_enrichment.tsv     per cluster x gene set: z, p, fold
  results/celltypes/tf_reweighting.tsv        per TF: weight change and its top cluster
  results/qc/tissue_association_qc.json       stability measurements and summary
  figures/fig6_tissue_association.pdf/.png

  python3 scripts/18_tissue_association.py
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
from lib_sources import (DATA, FIGURES, RESULTS, ROOT, log,  # noqa: E402
                         quiet_accelerate_blas_warnings, write_json)

quiet_accelerate_blas_warnings()

CELLTYPES = RESULTS / "celltypes"
DECODER = RESULTS / "decoder"
PRIORS = DATA / "tf_prior" / "weighted"
MANU_FIGS = ROOT / "manuscript" / "figures"
DEFAULT_ATLAS = Path.home() / "Documents" / "tropism_atlas"

plt.rcParams.update({"figure.dpi": 150, "savefig.bbox": "tight", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})


def symbol_map() -> dict[str, str]:
    """AGI locus -> gene symbol, from the TAIR alias table 04 already downloaded.

    A bar chart of AT3G26790 tells a reader nothing; FUS3 tells them something.
    """
    import re
    path = DATA / "tf_prior" / "tair_gene_aliases.txt"
    out: dict[str, str] = {}
    if not path.exists():
        return out
    for line in path.read_text(errors="replace").splitlines()[1:]:
        f = line.split("\t")
        if len(f) >= 2 and re.fullmatch(r"AT[1-5CM]G\d{5}", f[0].strip(), re.I):
            out.setdefault(f[0].strip().upper(), f[1].strip())
    return out


def organ_of(cluster: str) -> str:
    """`rosette_21d_10` -> `rosette_21d`. The resolution the fit actually supports."""
    return cluster.rsplit("_", 1)[0]


def load_signatures(atlas: Path) -> pd.DataFrame:
    for name in ("cell_type_signatures_full.csv", "cell_type_signatures.csv"):
        p = atlas / name
        if p.exists():
            sig = pd.read_csv(p, index_col=0)
            sig.index = sig.index.astype(str).str.strip().str.upper()
            return np.expm1(sig).clip(lower=0.0)
    raise SystemExit(f"no signature matrix under {atlas}")


def stability(fractions: pd.DataFrame) -> dict:
    """How much do fractions move between adjacent timepoints, per resolution?

    A response cannot be attributed to a compartment whose estimate is not stable enough
    to distinguish from its neighbour, so the number that matters is the size of the
    swings relative to the size of the estimate.
    """
    def swing(df: pd.DataFrame, min_mass: float = 0.02) -> float:
        # Only compartments carrying real mass. Including near-empty ones measures the
        # noise floor rather than the estimate's stability: a fraction moving from
        # 0.0001 to 0.0006 is a 500% relative swing and means nothing.
        keep = df.columns[df.mean() >= min_mass]
        if len(keep) == 0:
            return float("nan")
        sub = df[keep]
        d = sub.diff().abs().iloc[1:]
        return float(np.nanmedian((d / sub.mean()).to_numpy()))

    by_organ = fractions.T.groupby(organ_of).sum().T
    cl, org = swing(fractions), swing(by_organ)
    return {
        "what_is_stable": ("Neither per-timepoint resolution is stable; the "
                           "time-averaged composition used by the weighting is. See "
                           "timeaveraged_bootstrap_cv."),
        "cluster_level_median_relative_swing": round(cl, 3),
        "organ_level_median_relative_swing": round(org, 3),
        "organ_aggregation_helps": bool(org < cl),
        "n_clusters_above_2pct": int((fractions.mean() >= 0.02).sum()),
        "n_organs_above_2pct": int((by_organ.mean() >= 0.02).sum()),
        "n_clusters": int(fractions.shape[1]),
        "n_organs": int(by_organ.shape[1]),
        "organ_means": {k: round(float(v), 4)
                        for k, v in by_organ.mean().sort_values(ascending=False).items()
                        if v > 0.001},
    }


def enrich(sig: pd.DataFrame, genes: list[str], n_perm: int,
           rng: np.random.Generator) -> pd.DataFrame:
    """Per cluster: is this gene set expressed above the atlas background?"""
    present = [g for g in genes if g in sig.index]
    if len(present) < 5:
        return pd.DataFrame()
    # Scale each gene to its own maximum so the test is about WHERE a gene is expressed,
    # not how abundant it is; otherwise a few high-expressors dominate every cluster.
    scaled = sig.div(sig.max(axis=1).replace(0, np.nan), axis=0).dropna()
    present = [g for g in present if g in scaled.index]
    obs = scaled.loc[present].mean(axis=0)

    null = np.stack([scaled.iloc[rng.choice(len(scaled), len(present), replace=False)]
                     .mean(axis=0).to_numpy() for _ in range(n_perm)])
    mu, sd = null.mean(axis=0), null.std(axis=0)
    z = (obs.to_numpy() - mu) / np.where(sd > 0, sd, np.nan)
    p = (np.sum(null >= obs.to_numpy(), axis=0) + 1) / (n_perm + 1)
    return pd.DataFrame({
        "cluster": sig.columns, "organ": [organ_of(c) for c in sig.columns],
        "mean_scaled_expression": obs.to_numpy().round(4),
        "z": np.round(z, 3), "p_empirical": np.round(p, 5),
        "n_genes_tested": len(present),
    })


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    ap.add_argument("--permutations", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1260)
    ap.add_argument("--arm", default="A_primary_WildType")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    sig = load_signatures(args.atlas)
    fractions = pd.read_csv(CELLTYPES / f"{args.arm}_fractions.tsv", sep="\t", index_col=0)
    fractions = fractions.loc[sorted(fractions.index, key=float)]

    # ---------------------------------------------------------------- 1. stability
    st = stability(fractions)
    # The weighting consumes the time-averaged fraction, so that is the quantity whose
    # stability matters. Bootstrap the timepoints to measure it rather than assert it.
    fbar = fractions.mean(axis=0)
    boot = np.stack([fractions.iloc[rng.integers(0, len(fractions), len(fractions))]
                     .mean(axis=0).to_numpy() for _ in range(500)])
    mass = fbar.to_numpy() >= 0.02
    cv = np.divide(boot.std(axis=0), np.where(fbar.to_numpy() > 0, fbar.to_numpy(), np.nan))
    st["timeaveraged_bootstrap_cv"] = round(float(np.nanmedian(cv[mass])), 3)
    st["n_compartments_in_cv"] = int(mass.sum())
    log(f"  time-averaged composition bootstrap CV (compartments >2% mass): "
        f"{st['timeaveraged_bootstrap_cv']}")
    log(f"deconvolution stability: cluster-level relative swing "
        f"{st['cluster_level_median_relative_swing']}, "
        f"organ-level {st['organ_level_median_relative_swing']}")

    # ---------------------------------------------------------------- 2. enrichment
    sets = json.loads((DECODER / "radiation_signature.json").read_text())["sets"]
    frames = []
    for name in ("sog1_dependent", "myb3r_repressed", "ddr_core", "sog1_independent"):
        e = enrich(sig, sets.get(name, []), args.permutations, rng)
        if e.empty:
            continue
        e.insert(0, "gene_set", name)
        frames.append(e)
        best = e.nlargest(3, "z")
        log(f"  {name:<18} top clusters: " +
            ", ".join(f"{r.cluster} (z={r.z:.1f})" for r in best.itertuples()))
    enrichment = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if not enrichment.empty:
        enrichment.to_csv(CELLTYPES / "tissue_enrichment.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------- 3. reweighting
    binding = pd.read_csv(PRIORS / "prior_binding.tsv", sep="\t")
    weighted = pd.read_csv(PRIORS / f"prior_weighted_{args.arm}.tsv", sep="\t")
    m = binding.merge(weighted, on=["TF", "Gene"], suffixes=("_binding", "_weighted"))
    per_tf = m.groupby("TF").agg(
        binding_mass=("Input_binding", "sum"),
        weighted_mass=("Input_weighted", "sum"),
        n_targets=("Gene", "size")).reset_index()
    per_tf["log2_weight_change"] = np.log2(
        (per_tf["weighted_mass"] + 1) / (per_tf["binding_mass"] + 1))
    # Cell-type context is a multiplier in [0,1], so every TF's edge mass shrinks and an
    # absolute change would report "all demoted" regardless of what the atlas said. What
    # re-ranks the TFs against each other is shrinking more or less than the typical TF,
    # so the effect is centred on the median atlas-informed change.

    scaled = sig.div(sig.max(axis=1).replace(0, np.nan), axis=0)
    top_cluster, top_organ, tf_expressed = [], [], []
    for tf in per_tf["TF"]:
        if tf in scaled.index and np.isfinite(scaled.loc[tf]).any():
            c = scaled.loc[tf].idxmax()
            top_cluster.append(c); top_organ.append(organ_of(c)); tf_expressed.append(True)
        else:
            # Not in the atlas gene set, so cell-type context could not inform it and
            # its edges were passed through unchanged by design.
            top_cluster.append(""); top_organ.append("(not in atlas set)")
            tf_expressed.append(False)
    per_tf["top_cluster"] = top_cluster
    per_tf["top_organ"] = top_organ
    per_tf["atlas_informed"] = tf_expressed
    syms = symbol_map()
    per_tf["symbol"] = [syms.get(t, "") for t in per_tf["TF"]]
    per_tf.round(4).to_csv(CELLTYPES / "tf_reweighting.tsv", sep="\t", index=False)

    med = per_tf.loc[per_tf["atlas_informed"], "log2_weight_change"].median()
    per_tf["relative_weight_change"] = per_tf["log2_weight_change"] - med
    per_tf = per_tf.sort_values("relative_weight_change")
    informed = per_tf[per_tf["atlas_informed"]]
    log(f"  reweighting: {len(informed)}/{len(per_tf)} TFs atlas-informed; "
        f"relatively promoted {int((informed['relative_weight_change'] > 0.1).sum())}, "
        f"relatively demoted {int((informed['relative_weight_change'] < -0.1).sum())} "
        f"(median absolute shrinkage {med:.2f} log2, removed as a constant)")

    # ---------------------------------------------------------------- figure
    fig = plt.figure(figsize=(12.4, 4.5), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.15, 1.0])

    ax0 = fig.add_subplot(gs[0, 0])
    by_organ = fractions.T.groupby(organ_of).sum().T
    keep = by_organ.mean().sort_values(ascending=False).head(6).index
    for o in keep:
        ax0.plot([float(t) for t in by_organ.index], by_organ[o],
                 marker="o", ms=3.5, lw=1.4, label=o)
    ax0.set_xscale("log")
    ax0.set_xlabel("time after irradiation (min)")
    ax0.set_ylabel("deconvolved fraction")
    ax0.set_title("Per-timepoint fractions are unstable at\nboth resolutions — not dynamics",
                  fontsize=9)
    ax0.legend(fontsize=6.5, frameon=False, ncol=2)

    ax1 = fig.add_subplot(gs[0, 1])
    if not enrichment.empty:
        e = enrichment[enrichment["gene_set"] == "sog1_dependent"].copy()
        e = e.sort_values("z", ascending=False).head(14).iloc[::-1]
        cols = ["#c53030" if p < 0.05 else "#a0aec0" for p in e["p_empirical"]]
        ax1.barh(range(len(e)), e["z"], color=cols, height=0.7)
        ax1.set_yticks(range(len(e)))
        ax1.set_yticklabels(e["cluster"], fontsize=6.5)
        ax1.axvline(0, color="0.7", lw=0.8)
        ax1.set_xlabel("enrichment of SOG1-dependent genes  (z)")
        ax1.set_title("Where the SOG1-dependent response\nis expressed in the atlas",
                      fontsize=9)

    ax2 = fig.add_subplot(gs[0, 2])
    if not informed.empty:
        ends = pd.concat([informed.head(8), informed.tail(8)])
        cols = ["#2b6cb0" if v > 0 else "#c53030" for v in ends["relative_weight_change"]]
        ax2.barh(range(len(ends)), ends["relative_weight_change"], color=cols, height=0.7)
        ax2.set_yticks(range(len(ends)))
        ax2.set_yticklabels([f"{r.symbol or r.TF} · {r.top_organ}"
                             for r in ends.itertuples()], fontsize=6.2)
        ax2.axvline(0, color="0.7", lw=0.8)
        ax2.set_xlabel("edge-mass change vs the median TF (log2)")
        ax2.set_title("Which regulators the atlas\npromoted or demoted", fontsize=9)

    FIGURES.mkdir(parents=True, exist_ok=True)
    MANU_FIGS.mkdir(parents=True, exist_ok=True)
    for t in (FIGURES / "fig6_tissue_association.pdf", MANU_FIGS / "fig6_tissue_association.pdf"):
        fig.savefig(t)
    for t in (FIGURES / "fig6_tissue_association.png", MANU_FIGS / "fig6_tissue_association.png"):
        fig.savefig(t, dpi=200)
    plt.close(fig)

    write_json(RESULTS / "qc" / "tissue_association_qc.json", {
        "generated": dt.date.today().isoformat(),
        "arm": args.arm,
        "deconvolution_stability": st,
        "stability_caveat": (
            "Per-timepoint fractions are not stable enough to read as compositional "
            "dynamics, at cluster OR organ resolution — aggregating to organ does not "
            "help, which was measured rather than assumed. The time-averaged composition "
            "that the weighting consumes is stable (see timeaveraged_bootstrap_cv). "
            "Tissue attribution in this study therefore rests on the signature-gene "
            "enrichment, which uses the atlas expression matrix directly and involves no "
            "deconvolution."),
        "top_enriched": (
            enrichment.sort_values("z", ascending=False)
            .groupby("gene_set").head(5)[["gene_set", "cluster", "organ", "z", "p_empirical"]]
            .to_dict("records") if not enrichment.empty else []),
        "reweighting": {
            "n_tfs": int(len(per_tf)),
            "n_atlas_informed": int(len(informed)),
            "median_absolute_shrinkage_log2": round(float(med), 3),
            "note": ("Context is a multiplier in [0,1], so absolute edge mass can only "
                     "shrink; TFs are re-ranked by shrinking more or less than the "
                     "median, which is what changes a DREM attribution."),
            "n_relatively_promoted": int((informed["relative_weight_change"] > 0.1).sum()),
            "n_relatively_demoted": int((informed["relative_weight_change"] < -0.1).sum()),
            "most_promoted": informed.tail(6)[
                ["TF", "symbol", "top_cluster", "relative_weight_change"]].to_dict("records"),
            "most_demoted": informed.head(6)[
                ["TF", "symbol", "top_cluster", "relative_weight_change"]].to_dict("records"),
        },
    })
    log("  wrote tissue_enrichment.tsv, tf_reweighting.tsv, fig6_tissue_association")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
