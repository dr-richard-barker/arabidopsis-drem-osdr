#!/usr/bin/env python3
"""Detection floor: what dose does this assay actually need, and what does the ISS deliver?

A null result is only interpretable next to the sensitivity that produced it. This fits
the dose-response of both signature arms, estimates the dose at which each becomes
detectable, and places the radiation dose an ISS plant experiment actually receives
against that floor.

Three things the analysis has to keep straight.

**The two arms have different thresholds.** The MYB3R/G2-M arm responds at the lowest dose
tested and barely moves thereafter; the SOG1/DNA-damage arm is absent at 10 cGy and rises
steeply above it. Quoting one floor for "the radiation signature" would be wrong: the
diagnostic arm and the sensitive arm are different arms.

**The dose points come from different studies.** OSD-782 irradiated 4-week-old plants with
a 137Cs source at 1.4 cGy/s; OSD-658 irradiated flask-grown seedlings with sequential ion
beams simulating the GCR spectrum at NSRL, harvesting three hours later. The
fit is therefore a cross-study estimate, not a titration, and the confidence interval is
bootstrapped over the points rather than reported as a regression standard error.

**Suborbital is not orbital.** Two "missions" in the corpus are minutes-long suborbital
flights that accumulate essentially no dose. Pooling them with 11-70 day ISS missions
would understate the exposure of the orbital ones and overstate that of the suborbital.

ISS dose rate is taken from published dosimetry, not assumed: 0.355 +/- 0.04 mGy/day
measured by Bio-PADLES over 1,584 days (Yoshida et al., Heliyon 2022;8(8):e10266).

  results/decoder/dose_response.tsv    per dose point, both arms
  results/decoder/iss_dose.tsv         per mission: duration, dose, and the gap
  results/qc/dose_response_qc.json     fitted floors, bootstrap CIs, the verdict

  python3 scripts/23_dose_response.py
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import RESULTS, log, quiet_accelerate_blas_warnings, write_json  # noqa: E402

quiet_accelerate_blas_warnings()

OUT = RESULTS / "decoder"
ACT = RESULTS / "tf_activity"
BIODATA = "https://visualization.osdr.nasa.gov/biodata/api/v2"

# Yoshida et al., Heliyon 2022;8(8):e10266 — Bio-PADLES, MELFI freezer, 1,584 days.
ISS_DOSE_MGY_PER_DAY = 0.355
ISS_DOSE_SD = 0.04
ISS_DOSE_SOURCE = ("Yoshida et al., Heliyon 2022;8(8):e10266 "
                   "(Bio-PADLES, ISS MELFI, 1,584 days): 0.355 +/- 0.04 mGy/day")

# Significance threshold, matching the decoder's.
Z_THRESHOLD = 1.96

# Suborbital vehicles: minutes above the Karman line, negligible accumulated dose. Named
# rather than inferred, because "0 days" in the metadata is also what a missing date looks
# like and the two must not be conflated.
SUBORBITAL = {"VSS Unity", "Unity 22"}


def parse_date(s: str):
    for f in ("%m/%d/%Y", "%m/%d/%y", "%d-%b-%Y", "%m-%d-%Y", "%Y-%m-%d"):
        try:
            return dt.datetime.strptime(str(s).strip(), f)
        except (ValueError, TypeError):
            continue
    return None


def mission_days(acc: str) -> tuple[str, float | None]:
    u = f"{BIODATA}/dataset/{acc}/metadata/"
    try:
        d = json.load(urllib.request.urlopen(
            urllib.request.Request(u, headers={"User-Agent": "drem"}), timeout=60))
        m = (d.get(acc) or {}).get("metadata", {})
    except Exception:  # noqa: BLE001
        return "", None
    mis = m.get("mission")
    if isinstance(mis, list):
        mis = mis[0] if mis else {}
    if not isinstance(mis, dict):
        return "", None
    # Some studies span two resupply flights and record the mission name as a LIST
    # (e.g. ['SpaceX-2', 'SpaceX-4']). Joining them keeps the study identifiable instead
    # of stringifying a Python list into the output table.
    raw = mis.get("name")
    if isinstance(raw, list):
        name = " + ".join(str(x).strip() for x in raw if x)
    else:
        name = str(raw or "").strip()
    a, b = parse_date(mis.get("start date")), parse_date(mis.get("end date"))
    return name, ((b - a).days if a and b else None)


def fit_floor(dose: np.ndarray, z: np.ndarray, n_boot: int,
              rng: np.random.Generator) -> dict:
    """Dose at which a linear fit of z against log10(dose) crosses Z_THRESHOLD.

    Bootstrapped over the dose points, because with four points from two studies a
    regression standard error would understate the real uncertainty.
    """
    x = np.log10(dose)
    if len(x) < 3 or np.ptp(z) == 0:
        return {"floor_cgy": None, "note": "too few points to fit"}

    def cross(xx, zz):
        slope, intercept = np.polyfit(xx, zz, 1)
        if slope <= 0:
            return np.nan
        return float(10 ** ((Z_THRESHOLD - intercept) / slope))

    point = cross(x, z)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(x), len(x))
        if len(np.unique(x[idx])) < 2:
            continue
        c = cross(x[idx], z[idx])
        if np.isfinite(c) and 0 < c < 1e5:
            boots.append(c)
    boots = np.array(boots)
    return {
        "floor_cgy": round(point, 2) if np.isfinite(point) else None,
        "ci95_cgy": [round(float(np.percentile(boots, 2.5)), 2),
                     round(float(np.percentile(boots, 97.5)), 2)] if len(boots) > 20 else None,
        "n_bootstrap": int(len(boots)),
        "slope_per_log10_dose": round(float(np.polyfit(x, z, 1)[0]), 3),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bootstrap", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1260)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    pred = pd.read_csv(OUT / "predictions.tsv", sep="\t")

    # ---------------------------------------------------------------- dose-response
    dr = pred[pred["level"].astype(str).str.contains("centigray", na=False)].copy()
    dr["dose_cgy"] = dr["level"].str.extract(r"([\d.]+)").astype(float)
    dr = dr[dr["dose_cgy"] > 0].sort_values("dose_cgy")
    if dr.empty:
        raise SystemExit("no dose-labelled contrasts — run 16 and 17 first")
    dr[["accession", "level", "dose_cgy", "sog1_arm", "myb3r_arm",
        "radiation_index", "call"]].round(3).to_csv(
        OUT / "dose_response.tsv", sep="\t", index=False)

    log("dose-response:")
    for r in dr.itertuples():
        log(f"  {r.dose_cgy:>6g} cGy  {r.accession:<9} "
            f"SOG1 {r.sog1_arm:>7.2f}   MYB3R {r.myb3r_arm:>7.2f}")

    floors = {arm: fit_floor(dr["dose_cgy"].to_numpy(), dr[arm].to_numpy(),
                             args.bootstrap, rng)
              for arm in ("sog1_arm", "myb3r_arm")}
    # The bracket is more trustworthy than the fit: the highest dose that fails and the
    # lowest that passes, read straight off the data with no model.
    below = dr[dr["sog1_arm"] < Z_THRESHOLD]["dose_cgy"]
    above = dr[dr["sog1_arm"] >= Z_THRESHOLD]["dose_cgy"]
    bracket = [float(below.max()) if len(below) else None,
               float(above.min()) if len(above) else None]

    # ---------------------------------------------------------------- ISS dose
    meta = pd.read_csv(ACT / "contrast_meta.tsv", sep="\t")
    flight = meta[meta["is_flight"]]
    rows = []
    seen = {}
    for acc in sorted(flight["accession"].unique(), key=lambda a: int(a.split("-")[1])):
        name, days = mission_days(acc)
        if not name or name in seen:
            continue
        seen[name] = days
        # Three categories, not two. A mission with no recorded dates is NOT the same
        # as a suborbital hop: the first is unknown, the second is known to be ~zero
        # dose. Collapsing them would have reported five "suborbital" flights when only
        # two are, and quietly dropped three ISS missions from the dose range.
        if name in SUBORBITAL:
            category = "suborbital"
        elif days is None or days <= 0:
            category = "orbital, duration unknown"
        else:
            category = "orbital"
        rows.append({
            "mission": name, "duration_days": days, "category": category,
            "orbital": category == "orbital",
            "iss_dose_cgy": (round(days * ISS_DOSE_MGY_PER_DAY / 10, 3)
                             if category == "orbital" else None),
        })
    iss = pd.DataFrame(rows).sort_values("iss_dose_cgy")
    iss.to_csv(OUT / "iss_dose.tsv", sep="\t", index=False)

    orb = iss[iss["orbital"]]
    lo, hi = float(orb["iss_dose_cgy"].min()), float(orb["iss_dose_cgy"].max())
    cats = iss["category"].value_counts().to_dict()
    log(f"\nmissions by category: {cats}")
    log(f"  mission dose range: {lo:.2f} - {hi:.2f} cGy "
        f"({int(orb['duration_days'].min())}-{int(orb['duration_days'].max())} days)")

    floor = floors["sog1_arm"].get("floor_cgy")
    ci = floors["sog1_arm"].get("ci95_cgy")
    gap_lo = (bracket[0] / hi) if bracket[0] and hi else None
    gap_hi = (bracket[0] / lo) if bracket[0] and lo else None

    # The verdict is computed, not asserted: the claim stands only if the observed
    # non-detection bracket sits entirely above the ISS dose range.
    supported = bool(bracket[0] and lo and bracket[0] > hi)

    report = {
        "generated": dt.date.today().isoformat(),
        "z_threshold": Z_THRESHOLD,
        "dose_points": dr[["accession", "dose_cgy", "sog1_arm", "myb3r_arm"]]
        .round(3).to_dict("records"),
        "fitted_floor": floors,
        "observed_bracket_cgy": {
            "highest_dose_not_detected": bracket[0],
            "lowest_dose_detected": bracket[1],
            "note": "read directly off the data; more trustworthy than the fit, which "
                    "rests on four points from two studies with different material",
        },
        "iss": {
            "dose_rate_mgy_per_day": ISS_DOSE_MGY_PER_DAY,
            "dose_rate_sd": ISS_DOSE_SD,
            "source": ISS_DOSE_SOURCE,
            "n_orbital_missions": int(len(orb)),
            "mission_categories": cats,
            "suborbital": sorted(set(iss.loc[iss["category"] == "suborbital", "mission"])),
            "orbital_duration_unknown": sorted(set(
                iss.loc[iss["category"] == "orbital, duration unknown", "mission"])),
            "mission_dose_range_cgy": [round(lo, 3), round(hi, 3)],
        },
        "gap": {
            "iss_max_dose_cgy": round(hi, 3),
            "lowest_tested_non_detection_cgy": bracket[0],
            "fold_below_non_detection": (round(bracket[0] / hi, 1) if bracket[0] and hi else None),
            "fold_below_range": ([round(gap_lo, 1), round(gap_hi, 1)]
                                 if gap_lo and gap_hi else None),
            # Acute vs chronic delivery, the other half of the story.
            "acute_dose_rate_cgy_per_s": 1.4,
            "iss_dose_rate_cgy_per_s": round(ISS_DOSE_MGY_PER_DAY / 10 / 86400, 12),
            "dose_rate_ratio": round((1.4) / (ISS_DOSE_MGY_PER_DAY / 10 / 86400), -3),
        },
        "claim_supported": supported,
        "verdict": None,
    }
    report["verdict"] = (
        f"The SOG1/DNA-damage arm is undetected at {bracket[0]:g} cGy and detected at "
        f"{bracket[1]:g} cGy. ISS missions in this corpus deliver "
        f"{lo:.2f}-{hi:.2f} cGy, {report['gap']['fold_below_non_detection']}x below the "
        f"highest dose that already fails to register, at a dose rate ~"
        f"{report['gap']['dose_rate_ratio']:.0e} times lower. The absence of a "
        f"radiation-damage signature in flight is therefore what this assay predicts at "
        f"ISS doses, whether or not damage occurred."
    ) if supported else (
        "The ISS dose range is NOT below the observed non-detection bracket; the "
        "detection-limit argument does not hold as stated and must be rewritten.")

    write_json(RESULTS / "qc" / "dose_response_qc.json", report)
    log(f"\n  SOG1 arm floor: fitted {floor} cGy "
        f"(95% CI {ci}), observed bracket {bracket[0]}-{bracket[1]} cGy")
    log(f"  MYB3R arm floor: fitted {floors['myb3r_arm'].get('floor_cgy')} cGy")
    log(f"  claim supported: {supported}")
    log(f"  {report['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
