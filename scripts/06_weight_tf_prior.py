#!/usr/bin/env python3
"""Re-weight the TF prior by cell-type context — the auto-decoder lever.

This is the methodological contribution. iDREM can annotate a finished model with
single-cell data, but its own documentation states the cell-type data "is not used when
predicting the iDREM model". Here the single-cell atlas instead *enters* the model,
through the one channel DREM already provides: the third column of the TF-gene
interaction file, which DREM reads as a score in [0,1] rather than a binary flag.

For each arm, given cell-type fractions f_c(t) from `05_deconvolve_celltypes.py` and
per-cluster expression e(g,c) from the atlas signature matrix:

    context(TF) = sum_c fbar_c * ehat(TF, c)          fbar = fraction averaged over t
    context(gene) = sum_c fbar_c * ehat(gene, c)
    w(TF, gene) = binding(TF, gene) * context(TF)^alpha * context(gene)^beta

`ehat` is each gene's cluster profile scaled to its own maximum, so context measures
"is this gene expressed in the cell types that are actually present", independent of
how highly expressed the gene is overall. Weights are renormalised to [0,1] per TF so
they stay in DREM's expected range and so a TF is never globally advantaged simply for
being abundant.

Three priors are written, and the difference between them IS the experiment:

  flat      every edge = 1.0                    the binary prior DREM is usually given
  binding   the DAP-seq/ChIP score              what 04/04b produced
  weighted  binding x cell-type context         the lever

Coverage note, reported every run: the atlas signature matrix covers 4,000 highly
variable genes. A gene outside that set is, by construction, one whose expression does
not vary much across the 183 clusters, so its cell-type context is close to uniform and
the correct weight is its unmodified binding score. Such edges are therefore passed
through unchanged rather than dropped, and the fraction affected is measured so the
ablation's effect size can be read against the number of edges that could move at all.
Run `scripts/build_full_signatures.R` to widen the atlas beyond the HVG set.

  python3 scripts/06_weight_tf_prior.py
  python3 scripts/06_weight_tf_prior.py --alpha 1.0 --beta 0.5
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import DATA, RESULTS, log, quiet_accelerate_blas_warnings, write_json  # noqa: E402

quiet_accelerate_blas_warnings()

PRIOR = DATA / "tf_prior" / "tf_gene_edges.tsv"
CELLTYPES = RESULTS / "celltypes"
OUT = DATA / "tf_prior" / "weighted"
QC = RESULTS / "qc" / "weighted_prior_qc.json"
DEFAULT_ATLAS = Path.home() / "Documents" / "tropism_atlas"


def load_signatures(atlas: Path) -> pd.DataFrame:
    """Per-cluster expression, each gene scaled to its own maximum across clusters."""
    for name in ("cell_type_signatures_full.csv", "cell_type_signatures.csv"):
        path = atlas / name
        if path.exists():
            sig = pd.read_csv(path, index_col=0)
            sig.index = sig.index.astype(str).str.strip().str.upper()
            sig = np.expm1(sig).clip(lower=0.0)
            scaled = sig.div(sig.max(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
            log(f"signatures: {path.name} — {scaled.shape[0]} genes x "
                f"{scaled.shape[1]} clusters")
            return scaled
    raise SystemExit(f"no signature matrix in {atlas}")


def load_prior() -> pd.DataFrame:
    if not PRIOR.exists():
        raise SystemExit("run 04_fetch_tf_prior.py first")
    df = pd.read_csv(PRIOR, sep="\t")
    df["TF"] = df["TF"].astype(str).str.upper()
    df["Gene"] = df["Gene"].astype(str).str.upper()
    return df


def write_prior(path: Path, tf: np.ndarray, gene: np.ndarray, score: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["TF", "Gene", "Input"])
        for t, g, s in zip(tf, gene, score):
            w.writerow([t, g, f"{s:g}"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    ap.add_argument("--alpha", type=float, default=1.0,
                    help="exponent on the TF's cell-type context (default 1.0)")
    ap.add_argument("--beta", type=float, default=0.5,
                    help="exponent on the target's cell-type context (default 0.5); the "
                         "TF side dominates because DREM's attribution asks which TF "
                         "explains a split, not which target")
    ap.add_argument("--floor", type=float, default=0.05,
                    help="minimum weight for an edge with real binding evidence, so "
                         "cell-type context can down-weight an edge but never silently "
                         "delete it (default 0.05)")
    args = ap.parse_args()

    sig = load_signatures(args.atlas)
    prior = load_prior()
    log(f"prior: {len(prior)} edges, {prior['TF'].nunique()} TFs, "
        f"{prior['Gene'].nunique()} target genes")

    OUT.mkdir(parents=True, exist_ok=True)
    report = {"generated": dt.date.today().isoformat(),
              "parameters": {"alpha": args.alpha, "beta": args.beta, "floor": args.floor},
              "atlas": str(args.atlas),
              "n_edges": int(len(prior)),
              "arms": {}}

    # The flat and binding priors do not depend on an arm — written once.
    write_prior(OUT / "prior_flat.tsv", prior["TF"].to_numpy(), prior["Gene"].to_numpy(),
                np.ones(len(prior)))
    write_prior(OUT / "prior_binding.tsv", prior["TF"].to_numpy(),
                prior["Gene"].to_numpy(), prior["Input"].to_numpy())
    log("wrote prior_flat.tsv and prior_binding.tsv")

    for frac_path in sorted(CELLTYPES.glob("*_fractions.tsv")):
        slug = frac_path.name.replace("_fractions.tsv", "")
        fracs = pd.read_csv(frac_path, sep="\t", index_col=0)

        usable = [c for c in fracs.columns if c in sig.columns]
        fbar = fracs[usable].mean(axis=0)
        fbar = fbar / (fbar.sum() or 1.0)

        # context per gene: how present is it in the cell types this arm is made of
        ctx = sig[usable].to_numpy(dtype=float) @ fbar.to_numpy(dtype=float)
        ctx = pd.Series(ctx, index=sig.index)
        ctx = ctx / (ctx.max() or 1.0)

        tf_ctx = prior["TF"].map(ctx)
        gene_ctx = prior["Gene"].map(ctx)
        covered = tf_ctx.notna() | gene_ctx.notna()

        # Genes outside the atlas set are cell-type-invariant by construction, so their
        # correct context is neutral. Neutral is 1.0, not 0.
        tf_term = tf_ctx.fillna(1.0).clip(lower=1e-6) ** args.alpha
        gene_term = gene_ctx.fillna(1.0).clip(lower=1e-6) ** args.beta

        weighted = prior["Input"].to_numpy() * tf_term.to_numpy() * gene_term.to_numpy()

        # Renormalise within each TF so weights span [0,1] per regulator.
        out = pd.DataFrame({"TF": prior["TF"], "Gene": prior["Gene"], "w": weighted})
        peak = out.groupby("TF")["w"].transform("max").replace(0, np.nan)
        out["w"] = (out["w"] / peak).fillna(0.0)
        out.loc[(out["w"] < args.floor) & (prior["Input"].to_numpy() > 0), "w"] = args.floor

        path = OUT / f"prior_weighted_{slug}.tsv"
        write_prior(path, out["TF"].to_numpy(), out["Gene"].to_numpy(),
                    out["w"].round(4).to_numpy())

        changed = float(np.mean(np.abs(out["w"].to_numpy() - prior["Input"].to_numpy()) > 0.01))
        report["arms"][slug] = {
            "prior_file": path.name,
            "n_clusters_used": len(usable),
            "n_edges_atlas_informed": int(covered.sum()),
            "frac_edges_atlas_informed": round(float(covered.mean()), 4),
            "frac_edges_changed_by_0.01": round(changed, 4),
            "n_tfs_with_atlas_context": int(tf_ctx.notna().groupby(prior["TF"]).any().sum()),
            "top_context_cell_types": {k: round(float(v), 4)
                                       for k, v in fbar.sort_values(ascending=False)
                                       .head(8).items()},
            "weight_summary": {
                "mean": round(float(out["w"].mean()), 4),
                "median": round(float(out["w"].median()), 4),
                "p90": round(float(out["w"].quantile(0.9)), 4),
            },
        }
        log(f"  {slug}: {covered.sum():,} atlas-informed edges "
            f"({100 * covered.mean():.0f}%), {100 * changed:.0f}% moved by >0.01")

    write_json(QC, report)
    log(f"\nwrote {QC}")
    return 0 if report["arms"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
