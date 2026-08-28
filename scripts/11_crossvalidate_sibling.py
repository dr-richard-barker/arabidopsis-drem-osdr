#!/usr/bin/env python3
"""Cross-validate DREM paths against the sibling repository's WGCNA modules.

`Plant_response_to_radiation` analyses an overlapping set of OSDR studies with an
entirely different method: a Gaussian-process autoencoder and WGCNA co-expression
modules, with no DREM and no TF attribution. Its `blue` module is classified as
early-response and its `grey` module carries the DSB-repair core.

The two analyses share raw data but almost no methodology, so agreement between a DREM
path and a WGCNA module is genuine external corroboration rather than a restatement.
Disagreement is informative too and is reported as such: DREM partitions on trajectory
*shape* under a regulatory prior, WGCNA on pairwise co-expression, and they are not
obliged to agree.

Enrichment of each module within each DREM node is tested with a hypergeometric tail
probability over the genes both analyses modelled --- the background is the intersection,
not the genome, because a gene absent from one analysis could never have matched.

  results/comparison/sibling_crossvalidation.tsv   node x module overlap and p-value
  results/comparison/sibling_summary.json          best module per node

  python3 scripts/11_crossvalidate_sibling.py
  python3 scripts/11_crossvalidate_sibling.py --sibling /path/to/Plant_response_to_radiation
"""
from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import RESULTS, log, write_json  # noqa: E402

PARSED = RESULTS / "drem" / "parsed"
OUT = RESULTS / "comparison"
DEFAULT_SIBLING = Path.home() / "Documents" / "Plant_response_to_radiation"


def load_modules(sibling: Path) -> tuple[pd.DataFrame, dict]:
    mods = sibling / "results" / "wgcna" / "modules.csv"
    classes = sibling / "results" / "wgcna" / "module_classification.csv"
    if not mods.exists():
        raise SystemExit(
            f"sibling WGCNA modules not found at {mods}\n"
            f"Pass --sibling to point at a checkout of Plant_response_to_radiation.")
    df = pd.read_csv(mods)
    df["Gene"] = df["Gene"].astype(str).str.strip().str.upper()

    labels = {}
    if classes.exists():
        c = pd.read_csv(classes)
        labels = dict(zip(c["Module"], c["Classification"]))
    return df, labels


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sibling", type=Path, default=DEFAULT_SIBLING)
    ap.add_argument("--run", default="A_primary_WildType__weighted",
                    help="which DREM run to cross-validate (default the weighted WT model)")
    args = ap.parse_args()

    mods, labels = load_modules(args.sibling)
    genes_path = PARSED / f"{args.run}_genes.tsv"
    if not genes_path.exists():
        # Fall back to any wild-type run rather than failing: which prior was fitted
        # first should not decide whether this check can run at all.
        alt = sorted(PARSED.glob("A_primary_WildType__*_genes.tsv"))
        if not alt:
            raise SystemExit(f"no parsed DREM genes for {args.run} — run 09 first")
        genes_path = alt[0]
        log(f"  {args.run} not found; using {genes_path.stem}")

    assign = pd.read_csv(genes_path, sep="\t")
    assign["gene"] = assign["gene"].astype(str).str.strip().str.upper()

    background = sorted(set(assign["gene"]) & set(mods["Gene"]))
    if len(background) < 50:
        log(f"  only {len(background)} genes shared with the sibling analysis — "
            f"too few to test")
        write_json(OUT / "sibling_summary.json",
                   {"generated": dt.date.today().isoformat(),
                    "n_shared_genes": len(background),
                    "status": "insufficient overlap"})
        return 0

    bg = set(background)
    mods = mods[mods["Gene"].isin(bg)]
    assign = assign[assign["gene"].isin(bg)]
    N = len(bg)

    rows = []
    for node, block in assign.groupby("node"):
        node_genes = set(block["gene"])
        n = len(node_genes)
        for module, mblock in mods.groupby("Module"):
            module_genes = set(mblock["Gene"])
            K = len(module_genes)
            k = len(node_genes & module_genes)
            # P(X >= k): survival function at k-1.
            p = float(hypergeom.sf(k - 1, N, K, n)) if k else 1.0
            rows.append({
                "run": genes_path.stem.replace("_genes", ""),
                "node": int(node), "node_size": n,
                "module": module, "module_class": labels.get(module, ""),
                "module_size": K, "overlap": k,
                "expected": round(n * K / N, 2),
                "fold_enrichment": round((k / n) / (K / N), 3) if n and K else None,
                "p_hypergeometric": p,
            })

    df = pd.DataFrame(rows).sort_values(["node", "p_hypergeometric"])
    # Bonferroni over every node x module test actually performed.
    df["p_bonferroni"] = (df["p_hypergeometric"] * len(df)).clip(upper=1.0)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "sibling_crossvalidation.tsv", sep="\t", index=False)

    best = df.loc[df.groupby("node")["p_hypergeometric"].idxmin()]
    sig = best[best["p_bonferroni"] < 0.05]
    summary = {
        "generated": dt.date.today().isoformat(),
        "run": genes_path.stem.replace("_genes", ""),
        "sibling": str(args.sibling),
        "n_shared_genes": N,
        "n_nodes": int(assign["node"].nunique()),
        "n_modules": int(mods["Module"].nunique()),
        "n_nodes_with_significant_module": int(len(sig)),
        "best_module_per_node": [
            {"node": int(r.node), "module": r.module, "class": r.module_class,
             "overlap": int(r.overlap), "fold": r.fold_enrichment,
             "p_bonferroni": round(float(r.p_bonferroni), 6)}
            for r in best.itertuples()],
        "note": ("DREM partitions on trajectory shape under a regulatory prior and WGCNA "
                 "on pairwise co-expression. Agreement is external corroboration; "
                 "disagreement is not by itself evidence against either."),
    }
    write_json(OUT / "sibling_summary.json", summary)
    log(f"  {N} genes shared, {summary['n_nodes']} DREM nodes vs "
        f"{summary['n_modules']} WGCNA modules")
    log(f"  {summary['n_nodes_with_significant_module']} nodes enrich a module "
        f"at Bonferroni p < 0.05")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
