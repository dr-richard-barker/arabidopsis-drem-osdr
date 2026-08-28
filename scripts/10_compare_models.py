#!/usr/bin/env python3
"""Compare the DREM models: prior ablation, genotype contrast, and the published baseline.

Three questions, in the order they should be asked:

1. **Reproduction.** Bourbousse et al. (2018) built a DREM model over the same
   accessions and reported 11 coexpressed gene groups. Our flat-prior wild-type model is
   the closest analogue, so its node count is compared with theirs. This is a
   *calibration* check, not a claim of identity: our series pools three accessions and
   adds two early timepoints (10 and 45 min) the published model did not use, so a
   different count is an expected consequence of a denser series, not a failure. The
   number is reported either way.

2. **Prior ablation.** flat vs binding vs cell-type-weighted, same arm, same seed, same
   data. If the weighted prior changes nothing, that is the finding and it is reported
   as such — a negative result here is publishable and must not be buried.

3. **Genotype contrast.** The cohort's built-in positive control: SOG1 should be
   attributed to early splits in wild type and lose that attribution in sog1-1. If the
   method cannot see a knockout of the master regulator, its TF attributions are not
   trustworthy anywhere else.

  results/comparison/prior_ablation.tsv      per split: top TF under each prior
  results/comparison/genotype_contrast.tsv   SOG1 and MYB3R ranks, WT vs sog1-1
  results/comparison/summary.json            all three questions, answered with numbers

  python3 scripts/10_compare_models.py
"""
from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cohorts import BOURBOUSSE_DOI, BOURBOUSSE_GROUPS  # noqa: E402
from lib_sources import RESULTS, log, write_json  # noqa: E402

PARSED = RESULTS / "drem" / "parsed"
OUT = RESULTS / "comparison"

WATCHED = {
    "AT1G25580": "SOG1",
    "AT4G32730": "MYB3R1",
    "AT3G09370": "MYB3R3",
    "AT5G11510": "MYB3R4",
    "AT2G46770": "ANAC043",
    "AT5G13330": "ERF115",
}
PRIORS = ["flat", "binding", "weighted"]


def load(run: str, kind: str | None = None) -> pd.DataFrame:
    path = PARSED / f"{run}_tfscores.tsv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    if kind and "kind" in df:
        return df[df["kind"] == kind]
    return df


def common_kind(runs: list[str]) -> str | None:
    """The most specific TF test every one of these runs produced.

    DREM does not emit the same test types for every model: a run whose splits never
    resolved a one-vs-others comparison has only two-way files. Ranking a flat-prior
    run by `vsOthers` against a binding-prior run by `2way` would compare different
    statistics and attribute the difference to the prior, so the ablation is run on a
    test type all three priors share, and which one that was is recorded.
    """
    available = []
    for run in runs:
        df = load(run)
        available.append(set(df["kind"].unique()) if not df.empty else set())
    if not available or not all(available):
        return None
    shared = set.intersection(*available)
    for kind in ("vsOthers", "2way", "path"):
        if kind in shared:
            return kind
    return None


def arms_present() -> list[str]:
    arms = set()
    for p in PARSED.glob("*_tfscores.tsv"):
        stem = p.name.replace("_tfscores.tsv", "")
        if "__" in stem:
            arms.add(stem.rsplit("__", 1)[0])
    return sorted(arms)


def rank_table(df: pd.DataFrame) -> pd.DataFrame:
    """Rank TFs within each split. Lower DREM Score = more significant, so rank ascending."""
    out = df.copy()
    out["rank"] = out.groupby("split")["score"].rank(method="min", ascending=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top", type=int, default=1,
                    help="how many top TFs define 'the attribution' for a split")
    args = ap.parse_args()

    if not PARSED.exists():
        raise SystemExit("run 09_parse_drem_model.py first")
    OUT.mkdir(parents=True, exist_ok=True)

    report = {"generated": dt.date.today().isoformat(),
              "score_direction": "lower Score = more significant",
              "reproduction": {}, "ablation": {}, "genotype_contrast": {}}

    # ---------------------------------------------------------------- 1. reproduction
    from json import loads
    summary_path = RESULTS / "drem" / "parsed_summary.json"
    summary = loads(summary_path.read_text()) if summary_path.exists() else {"runs": {}}
    wt_flat = summary["runs"].get("A_primary_WildType__flat")
    report["reproduction"] = {
        "published": {"doi": BOURBOUSSE_DOI, "reported_gene_groups": BOURBOUSSE_GROUPS,
                      "timepoints_used": [20, 90, 180, 360, 720, 1440]},
        "ours": ({"run": "A_primary_WildType__flat",
                  "n_nodes": wt_flat["n_nodes"], "n_splits": wt_flat["n_splits"],
                  "timepoints_used": [0, 10, 20, 45, 90, 180, 360, 720, 1440]}
                 if wt_flat else None),
        "interpretation": (
            "Our series pools OSD-498/508/510 and adds the 10 and 45 min timepoints the "
            "published model did not use, so node counts are not expected to match "
            "exactly. Reported for calibration."),
    }

    # ---------------------------------------------------------------- 2. ablation
    ablation_rows = []
    for arm in arms_present():
        kind = common_kind([f"{arm}__{p}" for p in PRIORS])
        if kind is None:
            log(f"  {arm}: no TF test type shared by all three priors — skipped")
            report["ablation"][arm] = {"status": "no shared test type across priors"}
            continue
        tables = {p: rank_table(load(f"{arm}__{p}", kind)) for p in PRIORS}
        tables = {p: t for p, t in tables.items() if not t.empty}
        if len(tables) < 2:
            continue

        splits = sorted(set.intersection(*(set(t["split"]) for t in tables.values())))
        changed = 0
        for split in splits:
            row = {"arm": arm, "test_type": kind, "split": split}
            tops = {}
            for p, t in tables.items():
                best = t[t["split"] == split].nsmallest(args.top, "score")
                tops[p] = tuple(best["TF"])
                row[f"top_{p}"] = ",".join(best["TF"])
                row[f"score_{p}"] = round(float(best["score"].iloc[0]), 4) if len(best) else None
            distinct = len(set(tops.values()))
            row["attribution_changed"] = distinct > 1
            changed += distinct > 1
            ablation_rows.append(row)

        report["ablation"][arm] = {
            "test_type": kind,
            "priors_compared": sorted(tables),
            "n_splits": len(splits),
            "n_splits_attribution_changed": int(changed),
            "frac_splits_changed": round(changed / len(splits), 3) if splits else None,
        }
        log(f"  {arm}: {changed}/{len(splits)} splits changed top-TF attribution "
            f"across priors")

    if ablation_rows:
        pd.DataFrame(ablation_rows).to_csv(OUT / "prior_ablation.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------- 3. genotype
    contrast_rows = []
    for prior in PRIORS:
        kind = common_kind([f"A_primary_WildType__{prior}", f"A_primary_sog1-1__{prior}"])
        if kind is None:
            continue
        wt = rank_table(load(f"A_primary_WildType__{prior}", kind))
        mut = rank_table(load(f"A_primary_sog1-1__{prior}", kind))
        if wt.empty or mut.empty:
            continue
        for locus, name in WATCHED.items():
            def best_rank(df):
                sub = df[df["TF"] == locus]
                return (int(sub["rank"].min()), round(float(sub["score"].min()), 4),
                        int(sub["TF"].count())) if not sub.empty else (None, None, 0)

            wr, ws, wn = best_rank(wt)
            mr, ms, mn = best_rank(mut)
            contrast_rows.append({
                "prior": prior, "test_type": kind, "TF": name, "locus": locus,
                "wt_best_rank": wr, "wt_best_score": ws, "wt_splits_scored": wn,
                "sog1_best_rank": mr, "sog1_best_score": ms, "sog1_splits_scored": mn,
                "rank_worsens_in_mutant": (None if wr is None or mr is None else mr > wr),
            })

    if contrast_rows:
        cdf = pd.DataFrame(contrast_rows)
        cdf.to_csv(OUT / "genotype_contrast.tsv", sep="\t", index=False)
        sog1 = cdf[cdf["TF"] == "SOG1"]
        report["genotype_contrast"] = {
            "rows": len(cdf),
            "sog1_present_in_prior": bool(sog1["wt_splits_scored"].sum() > 0)
            if not sog1.empty else False,
            "sog1_rank_worsens_in_mutant": (
                sog1["rank_worsens_in_mutant"].dropna().tolist() if not sog1.empty else []),
            "note": ("SOG1 edges come from OSD-496/GSE112529 ChIP-seq (04b), not DAP-seq, "
                     "which does not assay SOG1. Without 04b this contrast is undefined."),
        }
        log(f"  genotype contrast: {len(cdf)} TF x prior rows written")

    write_json(OUT / "summary.json", report)
    log(f"\nwrote {OUT / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
