#!/usr/bin/env python3
"""Derive a radiation-response signature from the DREM model and the kinetic landscape.

The DREM analysis says more than "these genes changed". It says which genes travel
together, when they diverge, which TF is credited with the divergence, and --- because
the cohort carries a knockout of the master regulator --- which of that depends on SOG1.
This script turns that into gene sets a decoder can score an unseen study against.

Six sets are written. The first is the one that matters:

  sog1_dependent   induced >=2-fold in wild type, a SOG1 ChIP target, and losing that
                   induction in sog1-1. Genetic dependency, not just correlation: a
                   gene that responds in both genotypes is responding to irradiation
                   through some other route and is not diagnostic of THIS pathway.
  sog1_independent induced in wild type but equally induced in sog1-1 --- the contrast
                   set. A study scoring high here but low on sog1_dependent is
                   responding to something other than the canonical DDR.
  myb3r_repressed  repressed in wild type and a MYB3R1 target (the G2/M arm).
  ddr_core         canonical DSB-repair genes, the narrowest and most specific set.
  early / mid / late
                   partitioned by when each gene's |log2FC| peaks (<=90 min,
                   180-360 min, >=720 min) --- the kinetic phases.

Sets are disjoint where that matters and every gene carries its evidence, so a hit can
be traced back to why the gene was in the set at all.

  results/decoder/radiation_signature.json   the sets, with per-gene evidence
  results/decoder/signature_qc.json          sizes, overlaps, and the dependency split

  python3 scripts/15_build_radiation_signature.py
  python3 scripts/15_build_radiation_signature.py --min-lfc 1.0 --dependency-ratio 0.4
"""
from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import DATA, RESULTS, log, write_json  # noqa: E402

SERIES = RESULTS / "pseudotimeseries"
OUT = RESULTS / "decoder"

# Canonical DSB-repair / DNA-damage genes. Same loci the pipeline uses as batch-correction
# sentinels, so the decoder and the QC gate agree on what "the DDR" means.
DDR_CORE = {
    "AT4G21070": "BRCA1", "AT5G20850": "RAD51", "AT4G02390": "PARP2",
    "AT3G27630": "SMR7", "AT5G24280": "GMI1", "AT5G48720": "XRI1",
    "AT1G07500": "SMR5", "AT2G31320": "PARP1", "AT5G66130": "RAD17",
    "AT3G02680": "NBS1", "AT5G64520": "XRCC2", "AT1G80420": "XRCC1",
}
MYB3R1 = "AT4G32730"


def load_arm(name: str) -> pd.DataFrame:
    path = SERIES / f"{name}_expression.tsv"
    if not path.exists():
        raise SystemExit(f"missing {path} — run 03_build_pseudotimeseries.py")
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str).str.upper()
    return df[sorted(df.columns, key=float)]


def tf_targets(locus: str) -> set[str]:
    path = DATA / "tf_prior" / "tf_gene_edges.tsv"
    if not path.exists():
        return set()
    out = set()
    with open(path) as fh:
        next(fh)
        for line in fh:
            tf, gene, _score = line.rstrip("\n").split("\t")
            if tf == locus:
                out.add(gene)
    return out


def phase_of(row: pd.Series, times: list[float]) -> str:
    """Which kinetic phase a gene peaks in."""
    peak_t = times[int(np.argmax(np.abs(row.values)))]
    if peak_t <= 90:
        return "early"
    if peak_t <= 360:
        return "mid"
    return "late"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-lfc", type=float, default=1.0,
                    help="|log2FC| a gene must reach in wild type (default 1.0)")
    ap.add_argument("--dependency-ratio", type=float, default=0.4,
                    help="a gene is SOG1-dependent when its mutant response is below "
                         "this fraction of its wild-type response (default 0.4)")
    ap.add_argument("--max-set", type=int, default=400,
                    help="cap each set at the N strongest genes (default 400)")
    args = ap.parse_args()

    wt = load_arm("A_primary_WildType")
    mut = load_arm("A_primary_sog1-1")
    shared = wt.index.intersection(mut.index)
    wt, mut = wt.loc[shared], mut.loc[shared]
    times = [float(c) for c in wt.columns]
    log(f"arms: {len(shared)} shared genes over {len(times)} timepoints")

    wt_amp = wt.abs().max(axis=1)
    mut_amp = mut.abs().max(axis=1)
    wt_signed = wt.apply(lambda r: r.iloc[int(np.argmax(np.abs(r.values)))], axis=1)

    responsive = wt_amp >= args.min_lfc
    induced = responsive & (wt_signed > 0)
    repressed = responsive & (wt_signed < 0)
    # Ratio of mutant to wild-type amplitude: near 0 means the response needs SOG1.
    ratio = (mut_amp / wt_amp.replace(0, np.nan)).fillna(1.0)
    dependent = ratio <= args.dependency_ratio

    sog1_t = tf_targets("AT1G25580")
    myb_t = tf_targets(MYB3R1)
    log(f"SOG1 ChIP targets: {len(sog1_t)}   MYB3R1 targets: {len(myb_t)}")

    is_sog1_target = pd.Series(shared.isin(sog1_t), index=shared)
    is_myb_target = pd.Series(shared.isin(myb_t), index=shared)

    def take(mask: pd.Series, by: pd.Series) -> list[str]:
        sel = by[mask].sort_values(ascending=False)
        return list(sel.head(args.max_set).index)

    sets = {
        "sog1_dependent": take(induced & dependent & is_sog1_target, wt_amp),
        "sog1_independent": take(induced & ~dependent & is_sog1_target, wt_amp),
        "myb3r_repressed": take(repressed & is_myb_target, wt_amp),
        "ddr_core": [g for g in DDR_CORE if g in shared],
    }
    for phase in ("early", "mid", "late"):
        ph = wt[responsive].apply(lambda r: phase_of(r, times), axis=1)
        sets[f"phase_{phase}"] = take(
            responsive & responsive.index.isin(ph[ph == phase].index), wt_amp)

    evidence = {}
    for name, genes in sets.items():
        evidence[name] = {
            g: {"wt_peak_log2fc": round(float(wt_signed.get(g, np.nan)), 3),
                "mutant_amplitude": round(float(mut_amp.get(g, np.nan)), 3),
                "mutant_over_wt": round(float(ratio.get(g, np.nan)), 3),
                "sog1_target": bool(is_sog1_target.get(g, False)),
                "myb3r1_target": bool(is_myb_target.get(g, False)),
                "name": DDR_CORE.get(g, "")}
            for g in genes}

    OUT.mkdir(parents=True, exist_ok=True)
    write_json(OUT / "radiation_signature.json", {
        "generated": dt.date.today().isoformat(),
        "derived_from": {
            "arms": ["A_primary_WildType", "A_primary_sog1-1"],
            "sog1_edges": "data/tf_prior/sog1_edges.tsv (OSD-496 / GSE112529 ChIP-seq)",
            "note": "SOG1 dependency is measured by comparing the wild-type and sog1-1 "
                    "amplitudes of the same gene, so it is a genetic statement, not a "
                    "correlation with irradiation.",
        },
        "parameters": vars(args),
        "sets": sets,
        "evidence": evidence,
    })

    qc = {
        "generated": dt.date.today().isoformat(),
        "n_genes_considered": int(len(shared)),
        "n_responsive_wt": int(responsive.sum()),
        "n_induced": int(induced.sum()),
        "n_repressed": int(repressed.sum()),
        "n_sog1_dependent_any": int((induced & dependent).sum()),
        "set_sizes": {k: len(v) for k, v in sets.items()},
        "overlaps": {
            f"{a}|{b}": len(set(sets[a]) & set(sets[b]))
            for i, a in enumerate(sets) for b in list(sets)[i + 1:]
        },
        "median_mutant_over_wt": {
            k: round(float(np.median([evidence[k][g]["mutant_over_wt"] for g in v])), 3)
            for k, v in sets.items() if v
        },
    }
    write_json(OUT / "signature_qc.json", qc)

    for k, v in sets.items():
        med = qc["median_mutant_over_wt"].get(k)
        log(f"  {k:<18} n={len(v):<4} median mutant/WT amplitude = {med}")
    return 0 if sets["sog1_dependent"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
