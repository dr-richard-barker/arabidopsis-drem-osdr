#!/usr/bin/env python3
"""Assemble the cross-study pseudo-time-series DREM consumes.

DREM models each gene as a trajectory of log-ratios anchored at the first timepoint, so
this script turns nine per-study count matrices into one gene x timepoint matrix per
cohort arm, on a common log2 scale, with between-study offsets removed.

  results/pseudotimeseries/<cohort>_<arm>_expression.tsv   gene x timepoint means
  results/pseudotimeseries/<cohort>_<arm>_repeats.tsv      gene x replicate (DREM repeats)
  results/pseudotimeseries/<cohort>_series_qc.json         batch + sanity diagnostics

Pipeline per cohort:

  1. CPM + log2 per study (library-size normalisation is only ever within-study; a raw
     count is not comparable across library preps).
  2. Intersect genes across studies. TAIR locus IDs are stable, so this is a real join,
     not a name-matching exercise.
  3. Express every sample as a log2 ratio against its own study's control mean. This is
     the step that does most of the batch removal: a per-study offset cancels when both
     numerator and denominator come from the same study.
  4. Fit and remove a residual per-study offset on the *shared* timepoints only, so a
     study is only ever adjusted using timepoints another study also measured.
  5. Average replicates into the timepoint grid; keep the per-replicate matrix too,
     because DREM takes a repeats file and using it beats pre-medianing.

Held-out studies (OSD-782) are written to their own files and never pooled into the
primary series: they cross a dose regime, so a shared offset would be fitting away real
biology.

  python3 scripts/03_build_pseudotimeseries.py
  python3 scripts/03_build_pseudotimeseries.py --cohort A --min-cpm 1.0
"""
from __future__ import annotations

import argparse
import datetime as dt
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cohorts import COHORTS  # noqa: E402
from lib_sources import DATA, RESULTS, log, write_json  # noqa: E402

OUT = RESULTS / "pseudotimeseries"

# DNA-damage-response genes with well-established gamma induction in Arabidopsis. They
# are the post-correction sanity check: if batch correction flattens these, the
# correction is wrong, not the genes. Loci from TAIR10.
DDR_SENTINELS = {
    "AT4G21070": "BRCA1",
    "AT5G20850": "RAD51",
    "AT4G02390": "PARP2",
    "AT3G27630": "SMR7",
    "AT5G24280": "GMI1",
    "AT5G48720": "XRI1",
}


def load_counts(acc: str) -> pd.DataFrame:
    df = pd.read_csv(DATA / "counts" / f"{acc}_counts.csv", index_col=0)
    df.index = df.index.astype(str).str.strip().str.upper()
    return df.loc[~df.index.duplicated(keep="first")]


def cpm_log2(counts: pd.DataFrame, min_cpm: float, min_samples: int) -> pd.DataFrame:
    """Within-study CPM then log2. Genes too low to be measurable anywhere are dropped."""
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.divide(lib, axis=1) * 1e6
    keep = (cpm >= min_cpm).sum(axis=1) >= min_samples
    return np.log2(cpm.loc[keep] + 1.0)


def load_metadata(cohort_key: str) -> pd.DataFrame:
    md = pd.read_csv(DATA / "metadata_master.csv", dtype=str).fillna("")
    md = md[md["cohort"] == cohort_key].copy()
    md["time_minutes"] = pd.to_numeric(md["time_minutes"], errors="coerce")
    return md


def study_logratios(acc: str, md: pd.DataFrame, min_cpm: float) -> pd.DataFrame | None:
    """Per-sample log2 ratio against this study's own control mean.

    Returns None when a study has no usable control: without an internal reference its
    values cannot be placed on the same scale as the others, and inventing one from
    another study is exactly the batch artefact this pipeline exists to avoid.
    """
    rows = md[md["accession"] == acc]
    counts = load_counts(acc)
    present = [s for s in rows["sample_name"] if s in counts.columns]
    if not present:
        log(f"  {acc}: no metadata sample joins the count matrix — skipped")
        return None

    logx = cpm_log2(counts[present], min_cpm, min_samples=max(2, len(present) // 4))
    ctrl = [s for s in present
            if rows.set_index("sample_name").loc[s, "arm"] == "control"]

    if ctrl:
        ref = logx[ctrl].mean(axis=1)
        ref_desc = f"{len(ctrl)} in-study control samples"
    else:
        # Studies like OSD-498 label mock samples per timepoint; where no sample is
        # flagged 'control' at all, fall back to the earliest timepoint as t0, which is
        # DREM's own anchoring convention. Recorded in the QC so it is never invisible.
        earliest = rows["time_minutes"].min()
        t0 = [s for s in present
              if rows.set_index("sample_name").loc[s, "time_minutes"] == earliest]
        if not t0:
            log(f"  {acc}: no control and no resolvable t0 — skipped")
            return None
        ref = logx[t0].mean(axis=1)
        ref_desc = f"earliest timepoint ({earliest:g} min, {len(t0)} samples) as t0 anchor"

    log(f"  {acc}: {len(present)} samples, {logx.shape[0]} genes, ref = {ref_desc}")
    lr = logx.subtract(ref, axis=0)
    lr.attrs["reference"] = ref_desc
    return lr


def remove_study_offsets(per_study: dict[str, pd.DataFrame],
                         sample_time: dict[str, float]) -> tuple[pd.DataFrame, dict]:
    """Remove a residual per-study offset, fitted only on shared timepoints.

    A study is adjusted by a single scalar: the median gene-wise difference between its
    own profile and the cross-study consensus, computed *only* at timepoints more than
    one study measured. Studies with no shared timepoint are left untouched rather than
    shifted against a reference they never overlap.
    """
    times_by_study = {acc: {sample_time[s] for s in df.columns if s in sample_time}
                      for acc, df in per_study.items()}
    counts = defaultdict(int)
    for ts in times_by_study.values():
        for t in ts:
            counts[t] += 1
    shared = {t for t, n in counts.items() if n > 1}

    genes = sorted(set.intersection(*(set(df.index) for df in per_study.values())))
    log(f"  shared timepoints: {sorted(shared)}  |  genes in common: {len(genes)}")

    def profile_at(df: pd.DataFrame, times: set[float]) -> pd.Series | None:
        cols = [s for s in df.columns if sample_time.get(s) in times]
        return df.loc[genes, cols].mean(axis=1) if cols else None

    consensus = pd.concat(
        [p for df in per_study.values() if (p := profile_at(df, shared)) is not None],
        axis=1).mean(axis=1)

    offsets, adjusted = {}, {}
    for acc, df in per_study.items():
        own = profile_at(df, shared & times_by_study[acc])
        if own is None:
            offsets[acc] = {"offset": 0.0, "reason": "no timepoint shared with another study"}
            adjusted[acc] = df.loc[genes]
            continue
        off = float(np.median(own - consensus))
        offsets[acc] = {"offset": off,
                        "fitted_on_timepoints": sorted(shared & times_by_study[acc])}
        adjusted[acc] = df.loc[genes] - off

    return pd.concat(adjusted.values(), axis=1), {
        "shared_timepoints": sorted(shared), "per_study_offset": offsets,
        "n_genes_common": len(genes),
    }


def sentinel_report(series: pd.DataFrame, times: list[float], enforce: bool) -> dict:
    """Do the canonical DDR genes still respond after correction?

    `enforce` is False for arms where a flat response is the expected biology rather
    than a defect — the low-dose holdout, and any genotype arm whose whole purpose is
    to lack the response. A gate that fails for a correct reason teaches people to
    ignore the gate.
    """
    out = {}
    for locus, name in DDR_SENTINELS.items():
        if locus not in series.index:
            out[name] = {"locus": locus, "status": "not measured in the common gene set"}
            continue
        row = series.loc[locus]
        # Column labels are formatted with %g, so look them up the same way rather than
        # via str(float) — str(10.0) is '10.0' and would match nothing.
        traj = {f"{t:g}": round(float(row[f"{t:g}"]), 3)
                for t in times if f"{t:g}" in row.index}
        if not traj:
            out[name] = {"locus": locus,
                         "status": f"no timepoint column matched (have {list(row.index)[:5]})"}
            continue
        out[name] = {"locus": locus, "trajectory_log2fc": traj,
                     "max_abs_log2fc": round(float(np.nanmax(np.abs(row.values))), 3)}

    responding = [n for n, v in out.items() if v.get("max_abs_log2fc", 0) >= 1.0]
    measured = [n for n, v in out.items() if "max_abs_log2fc" in v]
    if not measured:
        # Every sentinel failing to resolve is a plumbing bug, not biology — always fatal.
        return {"genes": out, "n_responding_at_2fold": 0, "responding": [],
                "enforced": True, "pass": False,
                "note": "no sentinel resolved a trajectory — check the timepoint columns"}
    return {"genes": out, "n_responding_at_2fold": len(responding),
            "responding": responding, "enforced": enforce,
            "pass": (len(responding) >= 3) if enforce else True}


def build_cohort(key: str, min_cpm: float) -> dict:
    cohort = COHORTS[key]
    md = load_metadata(key)
    sample_time = dict(zip(md["sample_name"], md["time_minutes"]))
    report: dict = {"cohort": key, "name": cohort["name"], "arms": {}}

    for group, accs in (("primary", cohort["primary_studies"]),
                        ("holdout", cohort["holdout_studies"])):
        if not accs:
            continue
        log(f"\n[{key}/{group}] {', '.join(accs)}")
        per_study = {}
        for acc in accs:
            lr = study_logratios(acc, md, min_cpm)
            if lr is not None:
                per_study[acc] = lr
        if not per_study:
            continue

        if len(per_study) > 1:
            combined, batch = remove_study_offsets(per_study, sample_time)
        else:
            only = next(iter(per_study.values()))
            combined, batch = only, {"per_study_offset": {}, "n_genes_common": len(only.index),
                                     "shared_timepoints": [],
                                     "note": "single study — no cross-study offset to fit"}

        # How a cohort splits into arms is a scientific choice, not a default.
        #   Cohort A: genotype IS the experiment (WT vs sog1-1), so it must split.
        #   Cohort B: ecotype is a nuisance variable — the contrast (flight vs ground)
        #     already lives in the log-ratio. Splitting on it would shatter the cohort
        #     into single-study arms of two timepoints each, too few for DREM to place
        #     a bifurcation. So pool, and carry ecotype as an annotation.
        cohort_md = md[md["accession"].isin(per_study)].copy()
        if cohort.get("arm_by", "genotype") == "pooled":
            cohort_md["_arm"] = "pooled"
        else:
            cohort_md["_arm"] = cohort_md["genotype"]

        for arm_name, arm_md in cohort_md.groupby("_arm"):
            timed = arm_md[arm_md["time_minutes"].notna() & (arm_md["arm"] == "treated")]
            cols = [s for s in timed["sample_name"] if s in combined.columns]
            if len(cols) < 2:
                continue
            times = sorted(timed[timed["sample_name"].isin(cols)]["time_minutes"].unique())

            means = pd.DataFrame({
                f"{t:g}": combined[[s for s in cols if sample_time[s] == t]].mean(axis=1)
                for t in times})
            reps = combined[cols].copy()
            seen: dict[str, int] = {}
            labels = []
            for s in cols:
                tp = f"{sample_time[s]:g}"
                seen[tp] = seen.get(tp, 0) + 1
                labels.append(f"{tp}|{seen[tp]}")
            reps.columns = labels

            slug = f"{key}_{group}_{arm_name.replace(' ', '')}"
            OUT.mkdir(parents=True, exist_ok=True)
            means.round(4).to_csv(OUT / f"{slug}_expression.tsv", sep="\t",
                                  index_label="Gene")
            reps.round(4).to_csv(OUT / f"{slug}_repeats.tsv", sep="\t", index_label="Gene")

            report["arms"][slug] = {
                "studies": sorted(per_study),
                "genotype": arm_name,
                "genotypes_pooled": sorted(
                    arm_md[arm_md["sample_name"].isin(cols)]["genotype"].unique()),
                "group": group,
                "n_genes": int(means.shape[0]),
                "timepoints": [float(t) for t in times],
                "n_samples": len(cols),
                "batch": batch,
                # Enforced only where a canonical high-dose DDR is genuinely expected:
                # the wild-type primary arm. sog1-1 is *supposed* to lose the response,
                # and the OSD-782 holdout is 10-100 cGy, 2-3 orders of magnitude below
                # the dose these sentinels were chosen for.
                "sentinels": sentinel_report(
                    means, times,
                    enforce=(key == "A" and group == "primary"
                             and arm_name.lower().startswith("wild"))
                ) if key == "A" else None,
                "expression_file": f"{slug}_expression.tsv",
                "repeats_file": f"{slug}_repeats.tsv",
            }
            log(f"  -> {slug}: {means.shape[0]} genes x {len(times)} timepoints "
                f"({len(cols)} samples)")
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", choices=["A", "B"], help="restrict to one cohort")
    ap.add_argument("--min-cpm", type=float, default=1.0,
                    help="drop genes below this CPM in most samples (default 1.0)")
    args = ap.parse_args()

    reports = {}
    for key in COHORTS:
        if args.cohort and key != args.cohort:
            continue
        reports[key] = build_cohort(key, args.min_cpm)

    qc = {"generated": dt.date.today().isoformat(), "min_cpm": args.min_cpm,
          "cohorts": reports}
    write_json(OUT / "series_qc.json", qc)
    log(f"\nwrote {OUT / 'series_qc.json'}")

    # The sentinel check is the one that can fail the run: it is the evidence that batch
    # correction preserved the biology rather than erasing it.
    bad = [slug for r in reports.values() for slug, a in r["arms"].items()
           if a.get("sentinels") and not a["sentinels"]["pass"]]
    if bad:
        log(f"\nSENTINEL FAILURE in {bad}: fewer than 3 canonical DDR genes reach 2-fold "
            f"after correction. Inspect {OUT / 'series_qc.json'} before trusting the series.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
