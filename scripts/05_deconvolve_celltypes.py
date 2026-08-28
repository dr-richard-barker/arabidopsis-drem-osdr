#!/usr/bin/env python3
"""Deconvolve each pseudo-timepoint into cell-type fractions, and project the latent code.

This is the first half of the auto-decoder lever. It takes the trained Salk-atlas
auto-decoder from the sibling Tropism project and asks, of every timepoint in the
pseudo-time-series: which cell types is this bulk signal coming from, and where does
that put the timepoint in the atlas's 32-dimensional latent space?

Inputs (from ~/Documents/tropism_atlas, override with --atlas):
  cell_type_signatures.csv              4000 HVGs x 183 atlas clusters
  autodecoder/cluster_stimulus_codes.json  183 clusters x 32 latent dimensions

Outputs:
  results/celltypes/<arm>_fractions.tsv     timepoint x 183 cell-type fractions
  results/celltypes/<arm>_latent.tsv        timepoint x 32 latent dimensions
  results/qc/deconvolution_qc.json          fit residuals and coverage

Method: non-negative least squares of each timepoint's expression profile against the
signature matrix, fractions constrained to sum to one. NNLS rather than a regression
with an intercept because a negative cell-type fraction has no meaning and would let
the fit trade one cell type off against another to chase noise.

Deconvolution runs on ABSOLUTE expression (linear CPM averaged per timepoint from the
raw counts), never on the log-ratio series that DREM consumes. A cell-type fraction is
a statement about composition, so it needs a composition to fit: a log fold-change has
already divided the composition out, and mixing signatures in log space is not a
mixture at all. The atlas signatures are Seurat log-normalised cluster means, so they
are returned to linear space with expm1 before the fit for the same reason.

An important caveat, recorded in the QC and repeated in the manuscript: the atlas is a
*developmental* whole-plant atlas, not an irradiated one. Deconvolution therefore
assumes cell-type expression signatures are stable under gamma irradiation. That is an
assumption the design cannot test, not a finding.

  python3 scripts/05_deconvolve_celltypes.py
  python3 scripts/05_deconvolve_celltypes.py --atlas /path/to/tropism_atlas
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import DATA, RESULTS, log, quiet_accelerate_blas_warnings, write_json  # noqa: E402

quiet_accelerate_blas_warnings()

OUT = RESULTS / "celltypes"
QC = RESULTS / "qc" / "deconvolution_qc.json"
SERIES = RESULTS / "pseudotimeseries"
DEFAULT_ATLAS = Path.home() / "Documents" / "tropism_atlas"


def load_atlas(atlas: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    sig_path = atlas / "cell_type_signatures.csv"
    code_path = atlas / "autodecoder" / "cluster_stimulus_codes.json"
    for p in (sig_path, code_path):
        if not p.exists():
            raise SystemExit(
                f"missing {p}\nThe auto-decoder artefacts come from the sibling "
                f"Tropism_autodecoder_2026 project. Pass --atlas to point at them.")

    sig = pd.read_csv(sig_path, index_col=0)
    sig.index = sig.index.astype(str).str.strip().str.upper()
    # Seurat log-normalised cluster means -> linear space. A mixture of cell types is
    # additive in expression, not in log expression.
    sig = np.expm1(sig).clip(lower=0.0)
    codes = pd.DataFrame(json.load(open(code_path))).T
    codes.index = codes.index.astype(str)
    log(f"atlas: {sig.shape[0]} signature genes x {sig.shape[1]} clusters; "
        f"latent {codes.shape[0]} clusters x {codes.shape[1]} dims")
    return sig, codes


def deconvolve(profile: pd.Series, sig: pd.DataFrame) -> tuple[np.ndarray, float]:
    """NNLS fractions for one timepoint, renormalised to sum to 1.

    Both sides are scaled to unit norm first so the fit is about profile *shape* rather
    than about matching the arbitrary absolute magnitudes of CPM against atlas means.
    """
    shared = sig.index.intersection(profile.index)
    A = sig.loc[shared].to_numpy(dtype=float)
    b = profile.loc[shared].to_numpy(dtype=float)

    finite = np.isfinite(b) & np.isfinite(A).all(axis=1)
    A, b = A[finite], b[finite]
    if b.size == 0 or not np.any(b > 0):
        return np.zeros(sig.shape[1]), float("nan")

    nb = np.linalg.norm(b) or 1.0
    b = b / nb
    A = A / (np.linalg.norm(A, axis=0, keepdims=True) + 1e-12)

    w, residual = nnls(A, b)
    total = w.sum()
    frac = w / total if total > 0 else w
    return frac, float(residual)


def timepoint_profiles(slug: str, md: pd.DataFrame,
                       series_cols: list[str]) -> pd.DataFrame | None:
    """Mean linear-CPM profile per timepoint, from the raw counts of this arm's studies.

    Reconstructed from counts rather than reused from the DREM series because the
    series holds log-ratios; see the module docstring.
    """
    cohort_key, group = slug.split("_")[0], slug.split("_")[1]
    genotype = slug.split("_", 2)[2]

    rows = md[(md["cohort"] == cohort_key) & (md["arm"] == "treated")
              & md["time_minutes"].notna()].copy()
    if genotype != "pooled":
        rows = rows[rows["genotype"].str.replace(" ", "") == genotype]
    if rows.empty:
        return None

    per_tp: dict[str, list[pd.Series]] = {}
    for acc, block in rows.groupby("accession"):
        path = DATA / "counts" / f"{acc}_counts.csv"
        if not path.exists():
            continue
        counts = pd.read_csv(path, index_col=0)
        counts.index = counts.index.astype(str).str.strip().str.upper()
        counts = counts.loc[~counts.index.duplicated(keep="first")]
        lib = counts.sum(axis=0).replace(0, np.nan)
        cpm = counts.divide(lib, axis=1) * 1e6
        for _, r in block.iterrows():
            s = r["sample_name"]
            if s in cpm.columns:
                per_tp.setdefault(f"{r['time_minutes']:g}", []).append(cpm[s])

    keep = {tp: pd.concat(v, axis=1).mean(axis=1)
            for tp, v in per_tp.items() if tp in series_cols}
    return pd.DataFrame(keep) if keep else None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS,
                    help=f"directory holding the auto-decoder artefacts "
                         f"(default {DEFAULT_ATLAS})")
    args = ap.parse_args()

    sig, codes = load_atlas(args.atlas)
    OUT.mkdir(parents=True, exist_ok=True)

    files = sorted(SERIES.glob("*_expression.tsv"))
    if not files:
        raise SystemExit("no pseudo-time-series found — run 03_build_pseudotimeseries.py")

    report = {"generated": dt.date.today().isoformat(),
              "atlas": str(args.atlas),
              "atlas_shape": {"signature_genes": int(sig.shape[0]),
                              "clusters": int(sig.shape[1]),
                              "latent_dims": int(codes.shape[1])},
              "caveat": "The Salk atlas (GSE226097) is a developmental whole-plant atlas, "
                        "not an irradiated one. Deconvolution assumes cell-type signatures "
                        "are stable under the treatment. This is an assumption, not a result.",
              "arms": {}}

    md = pd.read_csv(DATA / "metadata_master.csv", dtype=str).fillna("")
    md["time_minutes"] = pd.to_numeric(md["time_minutes"], errors="coerce")

    for path in files:
        slug = path.name.replace("_expression.tsv", "")
        series_cols = list(pd.read_csv(path, sep="\t", index_col=0, nrows=0).columns)
        expr = timepoint_profiles(slug, md, series_cols)
        if expr is None:
            log(f"  {slug}: no count profiles resolved — skipped")
            continue
        expr = expr[[c for c in series_cols if c in expr.columns]]
        shared = sig.index.intersection(expr.index)
        if len(shared) < 100:
            log(f"  {slug}: only {len(shared)} signature genes present — skipped")
            continue

        fracs, resids = {}, {}
        for tp in expr.columns:
            f, r = deconvolve(expr[tp], sig)
            fracs[tp] = f
            resids[tp] = r

        frac_df = pd.DataFrame(fracs, index=sig.columns).T
        frac_df.to_csv(OUT / f"{slug}_fractions.tsv", sep="\t", index_label="timepoint")

        # Latent trajectory: the fraction-weighted mean of the clusters' stimulus codes.
        # Clusters are matched by name; any signature column without a code is dropped
        # from the projection and counted, never silently treated as a zero vector.
        usable = [c for c in frac_df.columns if c in codes.index]
        w = frac_df[usable].to_numpy(dtype=float)
        w = w / np.clip(w.sum(axis=1, keepdims=True), 1e-12, None)
        latent = pd.DataFrame(w @ codes.loc[usable].to_numpy(dtype=float),
                              index=frac_df.index,
                              columns=[f"z{i}" for i in range(codes.shape[1])])
        latent.round(5).to_csv(OUT / f"{slug}_latent.tsv", sep="\t", index_label="timepoint")

        top = frac_df.mean(axis=0).sort_values(ascending=False).head(8)
        report["arms"][slug] = {
            "n_signature_genes_matched": int(len(shared)),
            "n_clusters_with_latent_code": len(usable),
            "n_clusters_without_latent_code": int(frac_df.shape[1] - len(usable)),
            "timepoints": list(expr.columns),
            "relative_residual": {k: round(v, 4) for k, v in resids.items()},
            "median_relative_residual": round(float(np.median(list(resids.values()))), 4),
            "top_cell_types_mean_fraction": {k: round(float(v), 4) for k, v in top.items()},
        }
        log(f"  {slug}: {len(shared)} genes matched, {len(usable)} clusters coded, "
            f"median residual {np.median(list(resids.values())):.3f}")
        log(f"      top cell types: {', '.join(top.index[:4])}")

    write_json(QC, report)
    log(f"\nwrote {QC}")
    return 0 if report["arms"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
