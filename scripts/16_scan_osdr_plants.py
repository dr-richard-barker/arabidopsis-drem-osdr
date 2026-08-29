#!/usr/bin/env python3
"""Score every OSDR plant study against the radiation signature.

For each study with an expression matrix, every contrast its metadata supports is
scored against the gene sets from `15_build_radiation_signature.py`. The question is
which studies carry a SOG1-dependent DNA-damage signature --- including studies that were
never about radiation.

Contrast definition is deliberately mechanical rather than clever. For each factor
column, a level is treated as the control if it matches a control vocabulary (ground
control, non-irradiated, mock, untreated, 0 dose); every other level of that factor is
scored against it. A study therefore yields one contrast per (factor, non-control level)
pair, each labelled, and nothing is chosen to make a result look better.

Scoring is a rank statistic, not a fold-change threshold, so it survives the platform and
depth differences between studies: genes are ranked by log2 fold change, the mean
normalised rank of a set is compared with 2,000 random gene sets of the same size drawn
from the same ranking, and the result is reported as a z-score with an empirical p-value.
A positive z means the set sits towards the induced end.

  results/decoder/study_scores.tsv     one row per (study, contrast, gene set)
  results/decoder/scan_qc.json         coverage, failures and skipped studies

  python3 scripts/16_scan_osdr_plants.py
  python3 scripts/16_scan_osdr_plants.py --permutations 5000 --organism "Arabidopsis thaliana"
"""
from __future__ import annotations

import argparse
import datetime as dt
import io
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import (DATA, RESULTS, download, get_text, log,  # noqa: E402
                         osdr_file_index, osdr_url, quiet_accelerate_blas_warnings,
                         write_json)

quiet_accelerate_blas_warnings()

OUT = RESULTS / "decoder"
SCAN_COUNTS = DATA / "scan_counts"
SIGNATURE = OUT / "radiation_signature.json"

# A level counts as the control arm of its factor if it matches any of these. Kept
# explicit and auditable: a mis-set control silently inverts every score in that study.
CONTROL_PATTERNS = [
    r"^ground\s*control", r"\bground control\b", r"\bnon[- ]?irradiat", r"\bmock\b",
    r"\buntreated\b", r"\bcontrol\b", r"^0\s*\{", r"^0\s*(gray|gy|cgy)\b",
    r"\bwild\s*type\b", r"\b1G on Earth\b", r"\bunexposed\b", r"\bno treatment\b",
    # Onboard centrifuge. OSD-251 and OSD-346 are all-flight g-gradient studies whose
    # only control is a 1G centrifuge on orbit; without this they contribute nothing,
    # and they are the only studies carrying a within-flight gravity dose-response.
    r"\b1\s*G by centrifugation\b",
]
CONTROL_RE = re.compile("|".join(CONTROL_PATTERNS), re.I)

# Factors that describe what a sample IS rather than what was done to it. Contrasting
# across them measures tissue or genotype, not treatment, and would drown the signal.
SKIP_FACTORS = {
    "study.factor value.organism part", "study.factor value.tissue",
    "study.factor value.tissue segment", "study.factor value.cell cycle phase",
    "study.factor value.sample preservation method", "study.factor value.age",
    "study.factor value.generation", "study.factor value.ecotype",
    "study.factor value.cultivar", "study.factor value.spacecraft",
}
NA = {"NaN", "nan", "", "{Not Applicable}", "Not Applicable", None}

# The cross-factor control fallback applies ONLY to factors describing an exposure
# magnitude. Applied to any factor it does real damage: on OSD-508/510/320 it paired each
# post-irradiation timepoint against the mock controls, generating 26 extra contrasts that
# are genuine irradiation responses (they score 4-21) but duplicate the study's pooled
# radiation contrast and are not labelled as radiation, so they enter the decoder's
# calibration as false positives and dropped its AUC from 0.95 to 0.85. The fallback
# exists to recover a dose axis; restricting it to dose factors is the scope it earns.
DOSE_FACTOR_RE = re.compile(r"absorbed.*dose|radiation dose|\bdose\b", re.I)


def load_signature() -> dict[str, list[str]]:
    if not SIGNATURE.exists():
        raise SystemExit("run 15_build_radiation_signature.py first")
    return json.loads(SIGNATURE.read_text())["sets"]


def factor_table(organism: str) -> pd.DataFrame:
    cache = DATA / "osdr" / f"factors_{organism.replace(' ', '_')}.csv"
    if not cache.exists():
        from urllib.parse import urlencode
        q = urlencode({"study.characteristics.organism": organism})
        text = get_text("https://visualization.osdr.nasa.gov/biodata/api/v2/"
                        f"query/metadata/?{q}&study.factor%20value&format=csv")
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(text)
    return pd.read_csv(cache, dtype=str).fillna("NaN")


def fetch_counts(acc: str) -> Path | None:
    """Unnormalized counts if the study has them, else normalized, else microarray.

    GeneLab's Affymetrix tables are usable without a probeset map: they carry a `TAIR`
    column of AGI loci alongside the probeset ID, so they join the same gene vocabulary
    as the RNA-seq matrices. Mixing platforms is safe here only because the downstream
    statistic is a within-contrast rank, which carries no platform-specific scale.
    """
    dest = SCAN_COUNTS / f"{acc}_counts.csv"
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    try:
        pairs = osdr_file_index(acc)
    except RuntimeError:
        return None
    prefs = [
        lambda n: n.endswith("RSEM_Unnormalized_Counts_GLbulkRNAseq.csv") and "rRNArm" not in n,
        lambda n: n.endswith("STAR_Unnormalized_Counts_GLbulkRNAseq.csv") and "rRNArm" not in n,
        lambda n: bool(re.search(r"Unnormalized_Counts.*\.csv$", n)) and "rRNArm" not in n,
        lambda n: bool(re.search(r"Normalized_Counts.*\.csv$", n)) and "rRNArm" not in n,
        lambda n: bool(re.search(r"(Unnormalized|Normalized)_Counts.*\.csv$", n)),
        lambda n: bool(re.search(r"array_normalized_expression_probeset.*\.csv$", n)),
    ]
    for pred in prefs:
        hits = sorted((n, u) for n, u in pairs if pred(n))
        if hits:
            try:
                return download(osdr_url(hits[0][1]), dest)
            except RuntimeError:
                return None
    return None


# Annotation columns carried by GeneLab microarray tables ahead of the sample columns.
ARRAY_ANNOTATION = {"TAIR", "SYMBOL", "GENENAME", "REFSEQ", "ENTREZID", "STRING_id",
                    "GOSLIM_IDS", "ProbesetID", "count_ENSEMBL_mappings"}


def read_expression(path: Path) -> pd.DataFrame | None:
    """Load an RNA-seq count matrix or a GeneLab microarray table into gene x sample."""
    df = pd.read_csv(path)
    if "TAIR" in df.columns:
        # Microarray: index on the AGI locus and drop every annotation column, or the
        # gene-name strings would be read as samples.
        df = df[df["TAIR"].notna()].copy()
        df = df.set_index("TAIR")
        df = df.drop(columns=[c for c in df.columns if c in ARRAY_ANNOTATION],
                     errors="ignore")
        # Several probesets can map to one locus; keep their mean rather than an
        # arbitrary first, which would depend on file ordering.
        df = df.apply(pd.to_numeric, errors="coerce")
        return df.groupby(level=0).mean()
    return df.set_index(df.columns[0])


def contrasts_for(acc: str, factors: pd.DataFrame) -> list[tuple[str, str, list[str], list[str]]]:
    """(factor, level, treated_samples, control_samples) for every scoreable contrast.

    A dose factor often has no control level of its own: OSD-658 records
    `{Not Applicable}` for its unirradiated samples rather than `0 cGy`, so its 40 and 80
    cGy arms have no within-factor reference and the whole dose axis is lost -- the study
    collapses to a single pooled "mixed radiation" contrast. Where that happens the
    controls are still identifiable, just from a different factor (`non-irradiated` under
    `ionizing radiation`), so fall back to the study's controls found in any other factor.
    That recovers the dose series without inventing a reference.
    """
    rows = factors[factors["id.accession"] == acc]
    cols = [c for c in factors.columns
            if c.startswith("study.factor value.") and c not in SKIP_FACTORS]

    # Controls identifiable anywhere in this study, for the cross-factor fallback.
    study_controls: set[str] = set()
    for col in cols:
        levels = [v for v in rows[col].unique() if v not in NA]
        if len(levels) < 2:
            continue
        for lv in levels:
            if CONTROL_RE.search(lv):
                study_controls.update(rows[rows[col] == lv]["id.sample name"])

    out = []
    for col in cols:
        levels = [v for v in rows[col].unique() if v not in NA]
        if len(levels) < 2:
            continue
        controls = [v for v in levels if CONTROL_RE.search(lv := v)]
        if controls:
            ctrl_samples = list(rows[rows[col].isin(controls)]["id.sample name"])
            cross = False
        elif not DOSE_FACTOR_RE.search(col):
            continue
        else:
            # No control level in this factor: use the study's controls, minus anything
            # this factor itself assigns a non-control level (a sample cannot be both).
            treated_any = set(rows[rows[col].isin(levels)]["id.sample name"])
            ctrl_samples = sorted(study_controls - treated_any)
            cross = True
            if len(ctrl_samples) < 2:
                continue
        for lv in levels:
            if lv in controls:
                continue
            treated = list(rows[rows[col] == lv]["id.sample name"])
            if len(treated) >= 2 and len(ctrl_samples) >= 2:
                out.append((col.replace("study.factor value.", "")
                            + ("*" if cross else ""), lv, treated, ctrl_samples))
    return out


def score_set(ranks: np.ndarray, idx: np.ndarray, n_perm: int,
              rng: np.random.Generator) -> tuple[float, float, float]:
    """Mean normalised rank of a gene set vs a size-matched permutation null."""
    n = len(idx)
    if n < 5:
        return float("nan"), float("nan"), float("nan")
    obs = float(ranks[idx].mean())
    null = np.array([ranks[rng.choice(len(ranks), n, replace=False)].mean()
                     for _ in range(n_perm)])
    sd = null.std()
    z = (obs - null.mean()) / sd if sd > 0 else 0.0
    # Two-sided empirical p with the +1 correction, so p is never reported as zero.
    p = (np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean())) + 1) / (n_perm + 1)
    return obs, float(z), float(p)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organism", default="Arabidopsis thaliana",
                    help="organism to scan (default Arabidopsis thaliana)")
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--min-cpm", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=1260)
    args = ap.parse_args()

    sets = load_signature()
    factors = factor_table(args.organism)
    accs = sorted(factors["id.accession"].unique(), key=lambda a: int(a.split("-")[1]))
    log(f"{len(accs)} {args.organism} studies in the metadata API")

    rng = np.random.default_rng(args.seed)
    rows, skipped = [], {}

    for acc in accs:
        cons = contrasts_for(acc, factors)
        if not cons:
            skipped[acc] = "no factor with an identifiable control level"
            continue
        path = fetch_counts(acc)
        if path is None:
            skipped[acc] = "no expression matrix in the OSDR file listing"
            continue
        try:
            counts = read_expression(path)
        except Exception as e:  # noqa: BLE001 - a malformed table is a skip, not a crash
            skipped[acc] = f"unreadable count matrix ({type(e).__name__})"
            continue
        if counts is None or counts.empty:
            skipped[acc] = "expression table had no usable gene index"
            continue
        counts.index = counts.index.astype(str).str.strip().str.upper()
        counts = counts.loc[~counts.index.duplicated(keep="first")]
        counts = counts.apply(pd.to_numeric, errors="coerce").dropna(how="all")

        # GeneLab microarray tables are already normalised log2 intensities (values in
        # roughly 2-14). Applying CPM and a second log to them would be nonsense, so the
        # scale is detected rather than assumed.
        is_log_scale = float(np.nanmax(counts.to_numpy())) < 100
        if is_log_scale:
            logx = counts.loc[counts.notna().sum(axis=1) >= 2]
        else:
            lib = counts.sum(axis=0).replace(0, np.nan)
            cpm = counts.divide(lib, axis=1) * 1e6
            keep = (cpm >= args.min_cpm).sum(axis=1) >= max(2, cpm.shape[1] // 4)
            logx = np.log2(cpm.loc[keep] + 1.0)

        scored_any = False
        for factor, level, treated, control in cons:
            t = [s for s in treated if s in logx.columns]
            c = [s for s in control if s in logx.columns]
            if len(t) < 2 or len(c) < 2:
                continue
            lfc = logx[t].mean(axis=1) - logx[c].mean(axis=1)
            lfc = lfc.replace([np.inf, -np.inf], np.nan).dropna()
            if len(lfc) < 1000:
                continue

            order = lfc.rank(ascending=False, method="average").to_numpy()
            ranks = 1.0 - (order - 1) / (len(order) - 1)  # 1 = most induced
            gene_pos = {g: i for i, g in enumerate(lfc.index)}

            scored_any = True
            for set_name, genes in sets.items():
                idx = np.array([gene_pos[g] for g in genes if g in gene_pos], dtype=int)
                obs, z, p = score_set(ranks, idx, args.permutations, rng)
                rows.append({
                    "accession": acc, "factor": factor, "level": level,
                    "n_treated": len(t), "n_control": len(c),
                    "gene_set": set_name,
                    "n_set_genes_measured": int(len(idx)),
                    "n_set_genes_total": len(genes),
                    "mean_rank": round(obs, 4) if obs == obs else None,
                    "z": round(z, 3) if z == z else None,
                    "p_empirical": round(p, 5) if p == p else None,
                })
        if not scored_any:
            skipped[acc] = "contrasts found but samples did not join the count matrix"
        else:
            log(f"  {acc}: {len(cons)} contrast(s) scored")

    if not rows:
        raise SystemExit("nothing scored")
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "study_scores.tsv", sep="\t", index=False)

    write_json(OUT / "scan_qc.json", {
        "generated": dt.date.today().isoformat(),
        "organism": args.organism,
        "permutations": args.permutations,
        "n_studies_in_api": len(accs),
        "n_studies_scored": int(df["accession"].nunique()),
        "n_contrasts": int(df.groupby(["accession", "factor", "level"]).ngroups),
        "skipped": skipped,
        "gene_sets": {k: len(v) for k, v in sets.items()},
    })
    log(f"\nscored {df['accession'].nunique()} studies, "
        f"{df.groupby(['accession','factor','level']).ngroups} contrasts, "
        f"{len(skipped)} skipped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
