"""Cohort definitions — which OSD accessions enter which pseudo-time-series, and why.

Every membership decision carries its reason in the data structure rather than in a
comment, because the reasons are the scientifically contestable part and they end up
verbatim in the manuscript's Methods and in `docs/data/cohorts.json`.

The timepoints listed here are the ones the OSDR biodata API reported on 2026-08-28.
They are *expectations*, not inputs: `02_harmonise_metadata.py` re-derives them from the
live API and fails if they have changed, because OSDR deposits do get revised.
"""
from __future__ import annotations

# --------------------------------------------------------------------------- Cohort A

# Bourbousse C, Vegesna N, Law JA. PNAS 2018;115(52):E12453-E12462.
# doi:10.1073/pnas.1810582115 — three OSD accessions, ONE experiment. The published
# paper built its own DREM model over these, giving us a reproduction target.
BOURBOUSSE_DOI = "10.1073/pnas.1810582115"
BOURBOUSSE_PMID = "30541889"
BOURBOUSSE_GROUPS = 11  # "revealing 11 coexpressed gene groups" — from the abstract.

COHORT_A = {
    "name": "gamma_irradiation",
    "label": "γ-irradiation DNA damage response",
    "time_unit": "minute",
    "studies": {
        "OSD-498": {
            "role": "core",
            "expected_timepoints": [10, 20, 90, 1440],
            "genotypes": ["Wild Type"],
            "why": "Bourbousse core, deposited as [RNA-seq gIR vs mock wt].",
        },
        "OSD-508": {
            "role": "core",
            "expected_timepoints": [10, 20, 45, 90, 180, 360, 720, 1440],
            "genotypes": ["Wild Type", "sog1-1"],
            "why": "Bourbousse core, deposited as [RNA-seq t vs t0]. Densest arm; carries "
                   "both genotypes.",
        },
        "OSD-510": {
            "role": "core",
            "expected_timepoints": [20, 90, 180, 360, 720, 1440],
            "genotypes": ["Wild Type", "sog1-1"],
            "why": "Bourbousse core, deposited as [RNA-seq DREMmodel] — the accession the "
                   "published DREM model was built from.",
        },
        "OSD-782": {
            "role": "extension_holdout",
            "expected_timepoints": [60, 180, 1440, 4320],
            "genotypes": ["Wild Type"],
            "why": "Independent lab, low-dose (10/100 cGy) Cs-137. Extends the series to "
                   "72 h and is the only genuinely cross-study arm, so it is held out for "
                   "validation rather than pooled: it crosses a dose regime.",
        },
        "OSD-502": {
            "role": "anchor",
            "expected_timepoints": [],
            "genotypes": ["Wild Type", "myb3r135"],
            "why": "Single-timepoint myb3r triple mutant. No time axis, so it cannot join "
                   "the series; used to anchor the repressor arm of the TF attribution.",
        },
    },
    # Timepoints the primary (pooled core) model is built on.
    # Genotype IS the experiment here (WT vs sog1-1), so arms split on it.
    "arm_by": "genotype",
    "primary_studies": ["OSD-498", "OSD-508", "OSD-510"],
    "holdout_studies": ["OSD-782"],
}

# --------------------------------------------------------------------------- Cohort B

COHORT_B = {
    "name": "spaceflight_development",
    "label": "Spaceflight seedling development",
    "time_unit": "day",
    "studies": {
        "OSD-193": {
            "role": "core",
            "expected_timepoints": [4, 8],
            "genotypes": ["Col-0", "Sku6"],
            "why": "Sku6 vs Col-0 roots, flight vs ground.",
        },
        "OSD-218": {
            "role": "core",
            "expected_timepoints": [4, 8],
            "genotypes": ["Col-0", "WS"],
            "why": "Col-0 vs WS, flight vs ground. Kept over its duplicate OSD-219.",
        },
        "OSD-281": {
            "role": "core",
            "expected_timepoints": [4, 8],
            "genotypes": ["Sku5", "WS"],
            "why": "Sku5 vs WS roots, flight vs ground.",
        },
        "OSD-437": {
            "role": "core",
            "expected_timepoints": [4, 6],
            "genotypes": ["Wild Type"],
            "why": "Adds day 6 and, uniquely, an onboard 1G-by-centrifugation control that "
                   "separates microgravity from the spaceflight environment.",
        },
    },
    # Ecotype is a nuisance variable here, not the contrast: splitting on it would leave
    # five single-study arms of two timepoints each. Pool, and record which ecotypes went in.
    "arm_by": "pooled",
    "primary_studies": ["OSD-193", "OSD-218", "OSD-281", "OSD-437"],
    "holdout_studies": [],
}

# --------------------------------------------------------------------------- exclusions

EXCLUSIONS = {
    "OSD-219": "Duplicate deposit: all 32 BioSample IDs are identical to OSD-218's "
               "(assay named gse95373). Including both would double-count every sample. "
               "Verified against the biodata sample-name listing, not assumed.",
    "OSD-615": "Glycomic ELISA profiling, not transcriptomic — no gene-level matrix to "
               "feed DREM.",
    "OSD-45": "Affymetrix microarray (4/8/12 d). Would extend Cohort B to a fourth "
              "timepoint but crosses platforms; available via --include-microarray.",
    "OSD-46": "Affymetrix microarray gamma/Fe-56 time course (atm-1). Cross-platform; "
              "reserved for the optional validation arm.",
    "OSD-320": "Affymetrix microarray gamma/Fe-56 time course (atm-1). Cross-platform.",
    "OSD-329": "Microarray gamma time course (atm, atr mutants). Cross-platform.",
    "OSD-496": "SOG1 ChIP-seq from the same Bourbousse system, but OSDR holds raw files "
               "and ISA only — no processed peaks. Used as an optional TF-prior "
               "refinement, not as a cohort member.",
}

OPTIONAL_MICROARRAY = ["OSD-45", "OSD-46", "OSD-320", "OSD-329"]

COHORTS = {"A": COHORT_A, "B": COHORT_B}


def all_accessions(include_microarray: bool = False) -> list[str]:
    accs = [a for c in COHORTS.values() for a in c["studies"]]
    if include_microarray:
        accs += OPTIONAL_MICROARRAY
    return sorted(set(accs), key=lambda a: int(a.split("-")[1]))


def cohort_of(accession: str) -> str | None:
    for key, c in COHORTS.items():
        if accession in c["studies"]:
            return key
    return None
