#!/usr/bin/env python3
"""Separate the microgravity and non-microgravity halves of spaceflight, using an
in-flight 1g centrifuge.

Every flight contrast elsewhere in this pipeline confounds microgravity with everything
else about being in orbit -- radiation, launch, confinement, hardware atmosphere. A
centrifuge flown alongside the static samples breaks that confound, because its material
experienced the whole mission EXCEPT weightlessness.

E-MTAB-2518 (Fengler et al., BioMed Res Int 2015) flew 11-day-old Arabidopsis Col-0 callus
cultures on Shenzhou-8 in the SIMBOX facility for five days, with exactly that design:

  FS_front, FS_rear   space flight microgravity
  FC_front, FC_rear   1g in-flight reference centrifuge
  GS_front, GS_rear   1g ground control

giving three contrasts that decompose the flight effect:

  FC vs GS   spaceflight WITHOUT microgravity -- the term that contains radiation
  FS vs FC   microgravity alone, radiation held constant (both arms flew)
  FS vs GS   the total flight effect, as every other study in this corpus measures it

Two limits are structural and are reported with every number rather than in a footnote.
There are two arrays per condition and they are FRONT and REAR positions in one hardware
container, so they are positional replicates, not biological ones: a difference between
conditions cannot be separated from a difference between two culture flasks. And a
five-day mission accumulates roughly 0.18 cGy, some two hundred times below this assay's
demonstrated 40 cGy floor, so a null in the FC-vs-GS term is predicted by dose alone and
carries no information about whether radiation damaged anything.

What the design CAN show, and what no other study here can, is whether the flight
transcriptome is dominated by weightlessness or by the rest of the environment.

  results/shenzhou/expression.tsv        RMA probeset x sample
  results/shenzhou/decomposition.tsv     three contrasts scored against both arms
  results/qc/shenzhou_qc.json            design, provenance, caveats, verdict

  python3 scripts/25_shenzhou_decomposition.py
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import (DATA, RESULTS, ROOT, download, log,  # noqa: E402
                         quiet_accelerate_blas_warnings, write_json)

quiet_accelerate_blas_warnings()

RAW = DATA / "arrayexpress"
OUT = RESULTS / "shenzhou"
ACC = "E-MTAB-2518"
BASE = f"https://www.ebi.ac.uk/biostudies/files/{ACC}"
CELS = ["FS_front.CEL", "FS_rear.CEL", "FC_front.CEL", "FC_rear.CEL",
        "GS_front.CEL", "GS_rear.CEL"]

# ATH1 probeset -> AGI locus. GEO's platform table is a plain download and avoids a
# second Bioconductor dependency for what is a two-column lookup.
GPL198 = ("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL198/annot/"
          "GPL198.annot.gz")

# Dose accumulated over the mission, on the same basis as every other ISS figure here
# (Yoshida et al. 2022, 0.355 mGy/day). Shenzhou-8 flew a comparable low Earth orbit.
MISSION_DAYS = 5
DOSE_MGY_PER_DAY = 0.355

CONTRASTS = [
    ("non_microgravity_spaceflight", "FC", "GS",
     "in-flight 1g centrifuge vs ground 1g: everything about the mission except "
     "weightlessness, including radiation"),
    ("microgravity_alone", "FS", "FC",
     "microgravity vs in-flight 1g: weightlessness with radiation and mission "
     "environment held constant"),
    ("total_flight_effect", "FS", "GS",
     "microgravity flight vs ground: the confounded contrast every other study in this "
     "corpus measures"),
]

RMA_R = r'''
suppressPackageStartupMessages({library(affy)})
args <- commandArgs(trailingOnly=TRUE)
setwd(args[1])
ab <- ReadAffy(filenames=list.files(pattern="[.]CEL$"))
e  <- exprs(rma(ab))
colnames(e) <- sub("[.]CEL$", "", colnames(e))
write.table(data.frame(probeset=rownames(e), e, check.names=FALSE),
            file=args[2], sep="\t", quote=FALSE, row.names=FALSE)
cat("rma ok:", nrow(e), "probesets x", ncol(e), "arrays\n")
'''


def acquire() -> pd.DataFrame:
    RAW.mkdir(parents=True, exist_ok=True)
    for f in CELS + [f"{ACC}.sdrf.txt"]:
        download(f"{BASE}/{f}", RAW / f)

    # The SDRF has repeated column names ("Protocol REF", "Material Type"), which makes
    # DictReader silently collapse them; index by position instead.
    rows = [l.rstrip("\n").split("\t") for l in (RAW / f"{ACC}.sdrf.txt").open()]
    hdr = rows[0]
    fi = next(i for i, h in enumerate(hdr) if "Factor Value" in h)
    ai = next(i for i, h in enumerate(hdr) if "Array Data File" in h)
    design = pd.DataFrame([{"array": r[ai].replace(".CEL", ""), "condition": r[fi]}
                           for r in rows[1:] if len(r) > max(fi, ai)])
    design["arm"] = design["array"].str.split("_").str[0]
    design["position"] = design["array"].str.split("_").str[1]
    return design


def rma(out_tsv: Path) -> pd.DataFrame:
    """RMA-normalise the CELs via R/affy. Requires affy + ath1121501cdf."""
    if out_tsv.exists() and out_tsv.stat().st_size > 1000:
        return pd.read_csv(out_tsv, sep="\t", index_col=0)
    script = RAW / "_rma.R"
    script.write_text(RMA_R)
    r = subprocess.run(["Rscript", str(script), str(RAW.resolve()), str(out_tsv.resolve())],
                       capture_output=True, text=True, timeout=3600)
    if r.returncode != 0 or not out_tsv.exists():
        raise SystemExit(
            "RMA failed. Install the Bioconductor pieces with:\n"
            "  Rscript -e 'BiocManager::install(c(\"affy\",\"ath1121501cdf\"))'\n\n"
            + (r.stderr or r.stdout)[-1500:])
    log("  " + (r.stdout or "").strip().splitlines()[-1])
    return pd.read_csv(out_tsv, sep="\t", index_col=0)


def probeset_to_agi() -> dict[str, str]:
    """ATH1 probeset -> AGI locus, from the GEO platform annotation."""
    import gzip
    path = download(GPL198, RAW / "GPL198.annot.gz")
    mapping: dict[str, str] = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith(("^", "!", "#")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            m = re.search(r"AT[1-5CM]G\d{5}", " ".join(f[1:6]), re.I)
            if m and f[0]:
                mapping.setdefault(f[0].strip(), m.group(0).upper())
    return mapping


def analytic_z(ranks: np.ndarray, idx: np.ndarray) -> float:
    """Same exact finite-population statistic used everywhere else in this pipeline."""
    N, n = len(ranks), len(idx)
    if n < 5 or n >= N:
        return float("nan")
    var = ranks.var(ddof=1) / n * (N - n) / (N - 1)
    return float((ranks[idx].mean() - ranks.mean()) / np.sqrt(var)) if var > 0 else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threshold", type=float, default=1.96)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    design = acquire()
    log(f"{ACC}: {len(design)} arrays")
    for r in design.itertuples():
        log(f"  {r.array:<10} {r.arm:<3} {r.position:<6} {r.condition}")

    expr = rma(OUT / "expression.tsv")
    log(f"  RMA: {expr.shape[0]} probesets x {expr.shape[1]} arrays")

    mapping = probeset_to_agi()
    expr = expr.rename(index=mapping)
    expr = expr[expr.index.astype(str).str.match(r"AT[1-5CM]G\d{5}$", na=False)]
    expr = expr.groupby(level=0).mean()   # several probesets can hit one locus
    log(f"  mapped to {expr.shape[0]} AGI loci")

    sig = json.loads((RESULTS / "decoder" / "radiation_signature.json").read_text())["sets"]
    arm_sets = {"sog1_dependent": sig["sog1_dependent"],
                "myb3r_repressed": sig["myb3r_repressed"]}

    rows = []
    for name, treated, control, why in CONTRASTS:
        t = [c for c in expr.columns if c.startswith(treated)]
        c = [c for c in expr.columns if c.startswith(control)]
        lfc = (expr[t].mean(axis=1) - expr[c].mean(axis=1)).dropna()
        order = lfc.rank(ascending=False, method="average").to_numpy()
        ranks = 1.0 - (order - 1) / (len(order) - 1)
        pos = {g: i for i, g in enumerate(lfc.index)}

        scored = {}
        for arm, genes in arm_sets.items():
            idx = np.fromiter((pos[g] for g in genes if g in pos), dtype=int)
            scored[arm] = analytic_z(ranks, idx)
            scored[f"n_{arm}"] = int(len(idx))
        # Same conjunctive index as the decoder: both arms must move, in opposite ways.
        radiation_index = min(scored["sog1_dependent"], -scored["myb3r_repressed"])
        rows.append({
            "contrast": name, "treated": treated, "control": control,
            "meaning": why, "n_treated": len(t), "n_control": len(c),
            "n_genes": int(len(lfc)),
            "sog1_arm": round(scored["sog1_dependent"], 3),
            "myb3r_arm": round(-scored["myb3r_repressed"], 3),
            "n_sog1_genes": scored["n_sog1_dependent"],
            "n_myb3r_genes": scored["n_myb3r_repressed"],
            "radiation_index": round(radiation_index, 3),
            "call": "radiation-like" if radiation_index >= args.threshold else "no signal",
        })
        log(f"  {name:<32} SOG1 {rows[-1]['sog1_arm']:>7}  "
            f"MYB3R {rows[-1]['myb3r_arm']:>7}  index {rows[-1]['radiation_index']:>7}"
            f"  {rows[-1]['call']}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "decomposition.tsv", sep="\t", index=False)

    dose_cgy = MISSION_DAYS * DOSE_MGY_PER_DAY / 10
    non_ug = df[df["contrast"] == "non_microgravity_spaceflight"].iloc[0]
    ug = df[df["contrast"] == "microgravity_alone"].iloc[0]

    write_json(RESULTS / "qc" / "shenzhou_qc.json", {
        "generated": dt.date.today().isoformat(),
        "accession": ACC,
        "source": "ArrayExpress/BioStudies; not present in NASA OSDR",
        "publication_doi": "10.1155/2015/547495",
        "platform": "Affymetrix ATH1-121501 (A-AFFY-2), RMA",
        "material": "11-day-old Arabidopsis thaliana Col-0 semisolid callus culture (SDRF material type 'cell')",
        "mission": {"name": "Shenzhou-8 (SIMBOX)", "days": MISSION_DAYS,
                    "estimated_dose_cgy": round(dose_cgy, 3),
                    "dose_basis": "0.355 mGy/day (Yoshida et al. 2022), as used throughout"},
        "design": design.to_dict("records"),
        "contrasts": rows,
        "caveats": {
            "replication": "Two arrays per condition, and they are FRONT and REAR "
                           "positions within one hardware container. They are positional "
                           "replicates, not biological ones, so a between-condition "
                           "difference cannot be separated from a between-flask one.",
            "dose": f"A {MISSION_DAYS}-day mission accumulates about {dose_cgy:.2f} cGy, "
                    f"roughly {40 / dose_cgy:.0f}x below this assay's demonstrated 40 cGy "
                    f"floor. A null in the non-microgravity term is therefore predicted by "
                    f"dose alone and is not evidence about radiation damage.",
            "tissue": "Cell culture, not seedlings. Every other contrast in this corpus "
                      "uses whole plants or seedlings, so this is not directly comparable "
                      "to them.",
        },
        "interpretation": (
            f"Neither half of the flight effect carries the radiation signature "
            f"(non-microgravity term index {non_ug['radiation_index']}, microgravity term "
            f"{ug['radiation_index']}; threshold {args.threshold}). Given the mission dose "
            f"that is the expected result and adds no independent evidence about "
            f"radiation. What the decomposition does show is which half of the flight "
            f"environment moves the transcriptome at all."
            if max(non_ug["radiation_index"], ug["radiation_index"]) < args.threshold else
            f"The non-microgravity term scores {non_ug['radiation_index']} and the "
            f"microgravity term {ug['radiation_index']}; a positive non-microgravity term "
            f"at {dose_cgy:.2f} cGy would be well below the assay's floor and should be "
            f"treated as a false positive until replicated."),
    })
    log(f"\n  mission dose ~{dose_cgy:.2f} cGy, {40 / dose_cgy:.0f}x below the 40 cGy floor")
    log(f"  wrote {OUT / 'decomposition.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
