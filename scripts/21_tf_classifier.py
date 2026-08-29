#!/usr/bin/env python3
"""Can a classifier separate flight from radiation on TF activity alone? With honest CV.

The per-TF tests found no single TF significantly altered in spaceflight. A classifier
asks the weaker but different question: is there a *multivariate* pattern that individual
tests would miss.

Three design choices carry the result.

**Leave-one-mission-out, always.** SpaceX-5 contributes five studies and SpaceX-2 four.
Studies sharing a mission share hardware, launch and ground handling, so a random split
would put near-duplicates on both sides and report a score that measures memorisation.
Every fold holds out an entire mission.

**The null permutes GROUP labels, not contrast labels.** The first implementation shuffled
labels within each mission and produced a null AUC of 0.993 +/- 0.000 --- identical to the
observed value, because the label is constant within a group (a mission is either a flight
mission or it is not), so shuffling inside it changes nothing at all. That null was
vacuous: it would have "validated" any result whatsoever. The exchangeable unit in a
grouped design with group-constant labels is the GROUP, so the null reassigns labels
across groups, preserving how many groups carry each label.

**A platform control decides whether any of it means anything.** The corpus mixes
microarray and RNA-seq. If a classifier predicts *platform* from the same features as well
as it predicts the biological label, the biological result is confounded and is reported
as confounded.

The primary task is flight vs radiation --- two exposure classes with known different
biology, which the per-TF tests already show is detectable. Flight vs ground is NOT a
classification task here: every contrast is already a treated-vs-control comparison, so
"ground" is the reference inside each feature, not a separate class.

  results/tf_activity/classifier_report.json   AUCs, nulls, controls
  results/tf_activity/coefficient_stability.tsv  per-TF coefficient across folds

  python3 scripts/21_tf_classifier.py
  python3 scripts/21_tf_classifier.py --permutations 500
"""
from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import RESULTS, log, quiet_accelerate_blas_warnings, write_json  # noqa: E402

quiet_accelerate_blas_warnings()

ACT = RESULTS / "tf_activity"
KEY_TFS = {"AT1G25580": "SOG1", "AT4G32730": "MYB3R1", "AT3G09370": "MYB3R3",
           "AT5G11510": "MYB3R4", "AT2G46770": "ANAC043"}


def loso_auc(X, y, groups, seed: int, model: str = "logit"):
    """Leave-one-group-out AUC, pooled over held-out predictions."""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import LeaveOneGroupOut
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    logo = LeaveOneGroupOut()
    preds = np.full(len(y), np.nan)
    coefs = []
    for tr, te in logo.split(X, y, groups):
        # A fold whose training set has only one class cannot be fitted; skipping it is
        # correct, and how many were skipped is reported rather than hidden.
        if len(np.unique(y[tr])) < 2:
            continue
        if model == "logit":
            clf = make_pipeline(StandardScaler(),
                                LogisticRegression(penalty="l2", C=0.1, max_iter=5000,
                                                   class_weight="balanced",
                                                   random_state=seed))
            clf.fit(X[tr], y[tr])
            preds[te] = clf.predict_proba(X[te])[:, 1]
            coefs.append(clf[-1].coef_.ravel())
        else:
            clf = RandomForestClassifier(n_estimators=400, random_state=seed,
                                         class_weight="balanced", n_jobs=-1)
            clf.fit(X[tr], y[tr])
            preds[te] = clf.predict_proba(X[te])[:, 1]
            coefs.append(clf.feature_importances_)

    ok = np.isfinite(preds)
    if ok.sum() < 4 or len(np.unique(y[ok])) < 2:
        return float("nan"), np.array(coefs), int(ok.sum())
    return float(roc_auc_score(y[ok], preds[ok])), np.array(coefs), int(ok.sum())


def permutation_null(X, y, groups, n_perm: int, seed: int) -> np.ndarray:
    """Null AUCs from labels reassigned ACROSS groups.

    Each group gets one label (verified below), so the group is the exchangeable unit.
    Permuting the group->label map preserves the number of groups per class and the
    within-group dependence, while destroying the association with the features.
    """
    rng = np.random.default_rng(seed)
    uniq = np.unique(groups)
    # A group carrying both labels would make group-level permutation ill-defined.
    group_label = {}
    for g in uniq:
        lab = np.unique(y[groups == g])
        if len(lab) != 1:
            raise SystemExit(
                f"group {g!r} carries both labels; group-level permutation assumes one "
                f"label per group. Re-check the grouping before trusting any null.")
        group_label[g] = int(lab[0])

    labels = np.array([group_label[g] for g in uniq])
    out = []
    for i in range(n_perm):
        shuffled = rng.permutation(labels)
        mapping = dict(zip(uniq, shuffled))
        yp = np.array([mapping[g] for g in groups])
        if len(np.unique(yp)) < 2:
            continue
        a, _, _ = loso_auc(X, yp, groups, seed + i)
        if a == a:
            out.append(a)
    return np.array(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--permutations", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1260)
    args = ap.parse_args()

    act = pd.read_csv(ACT / "tf_activity.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(ACT / "contrast_meta.tsv", sep="\t", index_col=0)
    act = act.dropna(axis=1, how="any")
    tfs = list(act.columns)

    sub = meta[meta["is_flight"] | meta["is_radiation"]]
    X = act.loc[sub.index].to_numpy(dtype=float)
    y = sub["is_radiation"].to_numpy().astype(int)
    groups = sub["mission"].fillna("(unknown)").to_numpy()
    # Radiation studies carry no mission; give each its own group so they are never
    # split across a fold boundary with a flight study.
    groups = np.array([g if g not in ("(unknown)", "nan") else f"rad::{a}"
                       for g, a in zip(groups, sub["accession"])])

    log(f"flight vs radiation: {len(y)} contrasts "
        f"({int(y.sum())} radiation, {int((1 - y).sum())} flight), "
        f"{len(np.unique(groups))} CV groups, {len(tfs)} features")

    auc, coefs, n_pred = loso_auc(X, y, groups, args.seed)
    auc_rf, _, _ = loso_auc(X, y, groups, args.seed, model="rf")
    log(f"  logistic AUC {auc:.3f}   random forest AUC {auc_rf:.3f}  "
        f"({n_pred}/{len(y)} contrasts predicted)")

    log(f"  permutation null ({args.permutations} draws, group labels permuted) ...")
    null = permutation_null(X, y, groups, args.permutations, args.seed)
    p = float((np.sum(null >= auc) + 1) / (len(null) + 1)) if len(null) else float("nan")
    log(f"  null AUC {null.mean():.3f} +/- {null.std():.3f}   p = {p:.4g}")

    # Platform control: same features, same CV, different label.
    plat = (sub["platform"] == "microarray").to_numpy().astype(int)
    auc_plat = float("nan")
    if len(np.unique(plat)) == 2:
        auc_plat, _, _ = loso_auc(X, plat, groups, args.seed)
    log(f"  PLATFORM control AUC {auc_plat:.3f} "
        f"({'confounded' if auc_plat == auc_plat and auc_plat >= auc else 'ok'})")

    stab = pd.DataFrame()
    if coefs.size:
        m, s = coefs.mean(axis=0), coefs.std(axis=0)
        stab = pd.DataFrame({
            "TF": tfs, "name": [KEY_TFS.get(t, "") for t in tfs],
            "mean_coefficient": np.round(m, 4), "sd_across_folds": np.round(s, 4),
            # A coefficient that changes sign between folds is not a finding.
            "sign_stability": np.round(np.mean(np.sign(coefs) == np.sign(m), axis=0), 3),
        }).sort_values("mean_coefficient", key=abs, ascending=False)
        stab.to_csv(ACT / "coefficient_stability.tsv", sep="\t", index=False)

    write_json(ACT / "classifier_report.json", {
        "generated": dt.date.today().isoformat(),
        "task": "flight vs radiation exposure, on 474-dim TF activity",
        "cv": "leave-one-mission-out",
        "n_contrasts": int(len(y)), "n_radiation": int(y.sum()),
        "n_flight": int((1 - y).sum()), "n_groups": int(len(np.unique(groups))),
        "n_features": len(tfs),
        "auc_logistic": round(auc, 4) if auc == auc else None,
        "auc_random_forest": round(auc_rf, 4) if auc_rf == auc_rf else None,
        "null_mean": round(float(null.mean()), 4) if len(null) else None,
        "null_sd": round(float(null.std()), 4) if len(null) else None,
        "p_permutation": round(p, 5) if p == p else None,
        "platform_control_auc": round(auc_plat, 4) if auc_plat == auc_plat else None,
        "platform_confounded": bool(auc_plat == auc_plat and auc == auc and auc_plat >= auc),
        "caveat": ("Feature values are TF target-set activities and the sets overlap "
                   "heavily (median 2,467 targets), so a coefficient identifies a target "
                   "set, not a transcription factor acting alone."),
        "top_stable_coefficients": (
            stab[stab["sign_stability"] >= 0.8].head(15).to_dict("records")
            if not stab.empty else []),
    })
    log(f"  wrote classifier_report.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
