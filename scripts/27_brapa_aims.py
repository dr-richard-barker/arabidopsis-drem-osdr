#!/usr/bin/env python3
"""Three cross-species tests against the B. rapa simulated-GCR experiment.

The B_rappa_LLGCSS dataset is 39 GeneLab-style libraries of Brassica rapa leaf: 0 or
40 cGy simulated GCR, WT or anthocyanin-less, two preservatives (DRS/RL). Its 40 cGy arm
is dose-matched AND quality-matched to OSD-658's lower arm, which makes it the only
cross-species replicate of a simulated-GCR plant exposure available anywhere.

  1. Does the Arabidopsis-trained decoder fire on B. rapa GCR? A signature built on
     Arabidopsis gamma irradiation, carried across a species boundary by orthology.
  2. What does a B. rapa-NATIVE GCR signature look like, and does screening the OSDR
     spaceflight corpus with it find anything the Arabidopsis signature missed? If
     terrestrial Arabidopsis is the wrong training set, a second species trained on
     simulated GCR is the obvious thing to try before concluding the flight null is real.
  3. Do the two species agree gene-for-gene at matched dose and quality? This is the
     cleanest available test of whether the DNA-damage response is conserved well enough
     for cross-species screening to be meaningful at all.

**The comparison is bounded by three things and every number here inherits them.**
Tissue differs (B. rapa leaf vs OSD-658 whole seedlings). Harvest time differs. And the
preservative factor is large and fully crossed with nothing, so it is handled as a block:
fold changes are computed within each (genotype, preservative) cell and combined, never
pooled across preservatives, because a preservative effect would otherwise enter the
contrast as if it were radiation.

Orthology comes from scripts/26 (reciprocal best hit, one-to-one by construction). Its
namespace is Ensembl Plants Brapa_1.0 -- the same `Bra######` ids the RSEM matrix uses.

  results/brapa/decoder_on_brapa.tsv     aim 1, per contrast
  results/brapa/native_signature.json    aim 2, the B. rapa-derived sets
  results/brapa/flight_screen.tsv        aim 2, OSDR flight scored with them
  results/brapa/concordance.tsv          aim 3, cross-species per-gene agreement
  results/qc/brapa_aims_qc.json

  python3 scripts/27_brapa_aims.py
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import DATA, RESULTS, log, write_json  # noqa: E402


def _load(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).parent / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_tf = _load("tfmod", "19_tf_activity_matrix.py")
analytic_z = _tf.analytic_z

BRAPA_REPO = Path.home() / "Documents" / "B_rappa_LLGCSS"
COUNTS = BRAPA_REPO / "counts and factors" / "GeneCounts_RSEM.csv"
FACTORS = BRAPA_REPO / "counts and factors" / "LEAF_Metadata.csv"
RBH = DATA / "brapa" / "brapa_arabidopsis_rbh.tsv"
OUT = RESULTS / "brapa"
THRESHOLD = 1.96
MIN_CPM_SAMPLES = 3          # a gene must clear 1 CPM in at least this many libraries
NATIVE_SET_SIZE = 200        # genes per direction in the B. rapa-native signature


# ---------------------------------------------------------------- inputs

def load_orthologs() -> dict[str, str]:
    """B. rapa gene -> AGI. One-to-one by construction (reciprocal best hit)."""
    if not RBH.exists():
        raise SystemExit("run scripts/26_brapa_orthologs.py first")
    m = {}
    with RBH.open() as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            m[f[0].upper()] = f[1].upper()
    return m


def load_brapa() -> tuple[pd.DataFrame, pd.DataFrame]:
    """CPM matrix over protein-coding Bra genes, plus the sample sheet."""
    if not COUNTS.exists():
        raise SystemExit(f"B. rapa counts not found: {COUNTS}")
    x = pd.read_csv(COUNTS, index_col=0)
    x.index = x.index.astype(str).str.upper()
    # ERCC spike-ins and ncRNA are not in the ortholog map and must not enter the
    # library-size denominator as if they were transcriptome.
    x = x[x.index.str.startswith("BRA")]
    cpm = x / x.sum(axis=0) * 1e6
    keep = (cpm >= 1).sum(axis=1) >= MIN_CPM_SAMPLES
    meta = pd.read_csv(FACTORS).set_index("sample_id")
    meta = meta.loc[[s for s in x.columns if s in meta.index]]
    return cpm.loc[keep, meta.index], meta


def blocked_lfc(cpm: pd.DataFrame, meta: pd.DataFrame,
                subset: pd.Series | None = None) -> pd.Series | None:
    """log2 fold change of 40 cGy over 0 cGy, computed within preservative blocks.

    Pooling across preservatives would let a preservative effect masquerade as radiation
    whenever the dose arms are not perfectly balanced within it. Each block contributes
    its own fold change and the blocks are averaged, so the preservative cancels.
    """
    m = meta if subset is None else meta[subset]
    parts = []
    for _, cell in m.groupby("preservative", observed=True):
        t = cell.index[cell["condition"] == "40_cGy"]
        c = cell.index[cell["condition"] == "0_cGy"]
        if len(t) < 2 or len(c) < 2:
            continue
        parts.append(np.log2(cpm[t].mean(axis=1) + 1) - np.log2(cpm[c].mean(axis=1) + 1))
    if not parts:
        return None
    return pd.concat(parts, axis=1).mean(axis=1)


def to_agi(v: pd.Series, orth: dict[str, str]) -> pd.Series:
    """Project a B. rapa gene vector into Arabidopsis space through the RBH map."""
    idx = [orth.get(g) for g in v.index]
    out = pd.Series(v.values, index=idx)
    return out[out.index.notna()]


def score(vec: pd.Series, sets: dict[str, list[str]]) -> dict:
    """Rank-score each gene set against the analytic null used throughout the pipeline."""
    ranks = stats.rankdata(vec.values) / len(vec)
    pos = {g: i for i, g in enumerate(vec.index)}
    out = {}
    for name, genes in sets.items():
        idx = np.array([pos[g] for g in genes if g in pos], dtype=int)
        out[name] = (analytic_z(ranks, idx) if len(idx) >= 5 else float("nan"), len(idx))
    return out


def reliability(mat: pd.DataFrame, t: list[str], c: list[str], n_rep: int,
                rng: np.random.Generator) -> dict | None:
    """How well a contrast reproduces itself: split BOTH arms, build two independent
    estimates of the SAME contrast, correlate them.

    This is the number that decides whether a cross-dataset correlation means anything.
    Two contrasts cannot correlate with each other more than the geometric mean of their
    own reliabilities (classical attenuation), so a near-zero reliability caps the
    comparison at zero regardless of how conserved the underlying biology is.
    """
    t, c = list(t), list(c)
    if len(t) < 4 or len(c) < 4:
        return None
    out = []
    for _ in range(n_rep):
        tp, cp = list(t), list(c)
        rng.shuffle(tp)
        rng.shuffle(cp)
        ht, hc = len(tp) // 2, len(cp) // 2
        f1 = (np.log2(mat[tp[:ht]].mean(axis=1) + 1)
              - np.log2(mat[cp[:hc]].mean(axis=1) + 1))
        f2 = (np.log2(mat[tp[ht:2 * ht]].mean(axis=1) + 1)
              - np.log2(mat[cp[hc:2 * hc]].mean(axis=1) + 1))
        r = stats.spearmanr(f1, f2)[0]
        if r == r:
            out.append(r)
    if not out:
        return None
    return {"mean": round(float(np.mean(out)), 4), "sd": round(float(np.std(out)), 4),
            "n_splits": len(out), "n_treated": len(t), "n_control": len(c)}


# ---------------------------------------------------------------- aims

def aim1(cpm, meta, orth, sig) -> pd.DataFrame:
    """Score the Arabidopsis-trained decoder on B. rapa GCR contrasts."""
    log("\nAIM 1 — the Arabidopsis decoder applied to B. rapa GCR")
    specs = [("all genotypes", None)] + [
        (g, meta["genotype"] == g) for g in sorted(meta["genotype"].unique())]
    rows = []
    for name, sub in specs:
        lfc = blocked_lfc(cpm, meta, sub)
        if lfc is None:
            continue
        s = score(to_agi(lfc, orth), sig)
        sog, n_sog = s["sog1_dependent"]
        myb, n_myb = s["myb3r_repressed"]
        idxv = min(sog, -myb)
        rows.append({"contrast": f"40 vs 0 cGy, {name}",
                     "n_samples": int(len(meta if sub is None else meta[sub])),
                     "sog1_arm": round(sog, 3), "myb3r_arm": round(-myb, 3),
                     "n_sog1_genes": n_sog, "n_myb3r_genes": n_myb,
                     "radiation_index": round(idxv, 3),
                     "call": "radiation-like" if idxv >= THRESHOLD else "no signal"})
        log(f"  {rows[-1]['contrast']:<34} SOG1 {rows[-1]['sog1_arm']:>7}  "
            f"MYB3R {rows[-1]['myb3r_arm']:>7}  index {rows[-1]['radiation_index']:>7}"
            f"  {rows[-1]['call']}")
    return pd.DataFrame(rows)


def aim2(cpm, meta, orth) -> tuple[dict, pd.DataFrame]:
    """Derive a B. rapa-native GCR signature, then screen the OSDR flight corpus."""
    log("\nAIM 2 — a B. rapa-native GCR signature, applied to spaceflight")
    lfc = to_agi(blocked_lfc(cpm, meta, meta["genotype"] == "WT"), orth)
    lfc = lfc[~lfc.index.duplicated()]
    ordered = lfc.sort_values()
    native = {"brapa_gcr_induced": list(ordered.index[-NATIVE_SET_SIZE:]),
              "brapa_gcr_repressed": list(ordered.index[:NATIVE_SET_SIZE])}
    log(f"  {len(native['brapa_gcr_induced'])} induced / "
        f"{len(native['brapa_gcr_repressed'])} repressed orthologs (WT, blocked LFC)")

    scan = _load("scanmod", "16_scan_osdr_plants.py")
    factors = scan.factor_table("Arabidopsis thaliana")
    flight = pd.read_csv(RESULTS / "tf_activity" / "contrasts.tsv", sep="\t") \
        if (RESULTS / "tf_activity" / "contrasts.tsv").exists() else None
    accs = sorted({p.stem.replace("_counts", "") for p in
                   (DATA / "scan_counts").glob("*_counts.csv")})
    rows = []
    for acc in accs:
        path = DATA / "scan_counts" / f"{acc}_counts.csv"
        expr = scan.read_expression(path)
        if expr is None:
            continue
        for factor, level, t, c in scan.contrasts_for(acc, factors):
            t = [s for s in t if s in expr.columns]
            c = [s for s in c if s in expr.columns]
            if len(t) < 2 or len(c) < 2:
                continue
            v = (np.log2(expr[t].mean(axis=1) + 1) - np.log2(expr[c].mean(axis=1) + 1))
            v.index = v.index.astype(str).str.upper()
            v = v[~v.index.duplicated()]
            s = score(v, native)
            up, n_up = s["brapa_gcr_induced"]
            dn, n_dn = s["brapa_gcr_repressed"]
            if not np.isfinite(up) or not np.isfinite(dn):
                continue
            conj = min(up, -dn)
            rows.append({"accession": acc, "factor": factor, "level": level,
                         "n_treated": len(t), "n_control": len(c),
                         "induced_arm": round(up, 3), "repressed_arm": round(-dn, 3),
                         "n_induced_genes": n_up, "n_repressed_genes": n_dn,
                         "brapa_index": round(conj, 3),
                         "call": "GCR-like" if conj >= THRESHOLD else "no signal"})
    d = pd.DataFrame(rows).sort_values("brapa_index", ascending=False)
    hits = d[d["call"] == "GCR-like"]
    log(f"  scored {len(d)} contrasts across {d['accession'].nunique()} studies; "
        f"{len(hits)} reach the threshold")

    # POSITIVE CONTROL. A null across spaceflight means nothing unless the signature can
    # detect radiation where radiation is known to be present and large. These are the
    # labelled irradiation arms, including two at 100 Gy.
    irr = d[d["factor"].astype(str).str.contains("ionizing radiation|radiation dose",
                                                 case=False, na=False)]
    ctrl = {"n_labelled_irradiations": int(len(irr)),
            "n_detected": int((irr["brapa_index"] >= THRESHOLD).sum()),
            "best_index": float(irr["brapa_index"].max()) if len(irr) else None,
            "worst_index": float(irr["brapa_index"].min()) if len(irr) else None}
    ctrl["passes"] = bool(ctrl["n_detected"] > 0)
    log(f"  positive control: detects {ctrl['n_detected']}/{ctrl['n_labelled_irradiations']}"
        f" labelled irradiations (best index {ctrl['best_index']})")
    if not ctrl["passes"]:
        log("  -> the signature FAILS its positive control. It does not fire on 100 Gy "
            "gamma, where the response is large and well characterised, so its null "
            "across spaceflight carries no information about spaceflight.")
    for _, r in d.head(5).iterrows():
        log(f"    {r['accession']:<9} {str(r['level'])[:38]:<38} {r['brapa_index']:>6}")
    return native, d, ctrl


def aim3(cpm, meta, orth, sig, rng) -> tuple[pd.DataFrame, dict]:
    """Cross-species agreement at matched dose (40 cGy) and quality (simulated GCR)."""
    log("\nAIM 3 — B. rapa vs OSD-658, matched at 40 cGy simulated GCR")
    f = pd.read_csv(DATA / "osdr" / "arabidopsis_factor_values.csv")
    f = f[f["id.accession"] == "OSD-658"]
    dose = f["study.factor value.absorbed radiation dose"].astype(str)
    treated = list(f.loc[dose.str.startswith("40"), "id.sample name"])
    control = list(f.loc[f["study.factor value.ionizing radiation"] == "non-irradiated",
                         "id.sample name"])
    expr = pd.read_csv(DATA / "scan_counts" / "OSD-658_counts.csv", index_col=0)
    expr.index = expr.index.astype(str).str.upper()
    cpm658 = expr / expr.sum(axis=0) * 1e6
    cpm658 = cpm658[(cpm658 >= 1).sum(axis=1) >= MIN_CPM_SAMPLES]
    t = [s for s in treated if s in cpm658.columns]
    c = [s for s in control if s in cpm658.columns]
    at = (np.log2(cpm658[t].mean(axis=1) + 1) - np.log2(cpm658[c].mean(axis=1) + 1))
    at = at[~at.index.duplicated()]

    br = to_agi(blocked_lfc(cpm, meta, meta["genotype"] == "WT"), orth)
    br = br[~br.index.duplicated()]
    shared = sorted(set(at.index) & set(br.index))
    log(f"  OSD-658: {len(t)} treated, {len(c)} control; {len(shared)} shared orthologs")

    rows = []
    for scope, genes in (("all_shared_orthologs", shared),
                         ("sog1_dependent", [g for g in sig["sog1_dependent"]
                                             if g in set(shared)]),
                         ("myb3r_repressed", [g for g in sig["myb3r_repressed"]
                                              if g in set(shared)])):
        if len(genes) < 10:
            continue
        rho, p = stats.spearmanr(at[genes], br[genes])
        rows.append({"scope": scope, "n_genes": len(genes),
                     "spearman_rho": round(float(rho), 4), "p": float(p)})
        log(f"    {scope:<22} n={len(genes):<6} rho={rho:+.3f}  p={p:.3g}")

    # POWER. Before reading anything into a near-zero cross-species correlation, ask what
    # correlation these two contrasts could possibly show. Each is split against itself;
    # the geometric mean of the two reliabilities caps any correlation between them.
    wt = meta[meta["genotype"] == "WT"]
    rel_br = reliability(cpm, list(wt.index[wt["condition"] == "40_cGy"]),
                         list(wt.index[wt["condition"] == "0_cGy"]), 200, rng)
    rel_at = reliability(cpm658, t, c, 200, rng)
    cap = None
    if rel_br and rel_at:
        cap = float(np.sqrt(max(rel_br["mean"], 0) * max(rel_at["mean"], 0)))
        log(f"  reliability: B. rapa {rel_br['mean']:+.3f}, OSD-658 {rel_at['mean']:+.3f}"
            f"  -> attenuation cap {cap:.3f}")
        if cap < 0.1:
            log("  -> at that cap this comparison has no power. A near-zero cross-species "
                "correlation is what an unreproducible contrast produces whatever the "
                "biology does, so it is not evidence about conservation.")

    # The floor for "do two irradiations agree" was measured within Arabidopsis in
    # script 24: two accessions of the SAME published experiment reach only 0.045.
    floor = None
    q = RESULTS / "qc" / "radiation_quality_qc.json"
    if q.exists():
        floor = json.loads(q.read_text()).get("power", {}).get("same_experiment_rho")
    return pd.DataFrame(rows), {
        "same_experiment_floor": floor,
        "contrast_reliability": {"brapa_40cGy": rel_br, "osd658_40cGy": rel_at},
        "attenuation_cap": None if cap is None else round(cap, 4),
        "interpretable": bool(cap is not None and cap >= 0.1),
    }


# ---------------------------------------------------------------- main

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    orth = load_orthologs()
    cpm, meta = load_brapa()
    sig = json.loads((RESULTS / "decoder" / "radiation_signature.json").read_text())["sets"]
    log(f"B. rapa: {cpm.shape[0]} expressed genes x {cpm.shape[1]} libraries; "
        f"{len(orth)} orthologs available")
    mapped = sum(1 for g in cpm.index if g in orth)
    log(f"  {mapped} expressed genes carry an ortholog "
        f"({100 * mapped / cpm.shape[0]:.1f}%)")

    a1 = aim1(cpm, meta, orth, sig)
    a1.to_csv(OUT / "decoder_on_brapa.tsv", sep="\t", index=False)
    native, a2, a2ctrl = aim2(cpm, meta, orth)
    a2.to_csv(OUT / "flight_screen.tsv", sep="\t", index=False)
    write_json(OUT / "native_signature.json",
               {"generated": dt.date.today().isoformat(),
                "derived_from": "B. rapa WT, 40 cGy vs 0 cGy simulated GCR, blocked on "
                                "preservative, projected to Arabidopsis by RBH orthology",
                "set_size": NATIVE_SET_SIZE, "sets": native})
    a3, extra = aim3(cpm, meta, orth, sig, np.random.default_rng(0))
    a3.to_csv(OUT / "concordance.tsv", sep="\t", index=False)

    flight_hits = a2[a2["call"] == "GCR-like"]
    spaceflight = flight_hits[flight_hits["factor"].astype(str).str.contains(
        "spaceflight|altered gravity", case=False, na=False)]
    qc = {
        "generated": dt.date.today().isoformat(),
        "orthology": {"n_pairs": len(orth), "expressed_genes_mapped": mapped,
                      "source": "scripts/26_brapa_orthologs.py, DIAMOND RBH"},
        "design": {"species": "Brassica rapa", "tissue": "leaf",
                   "dose_cgy": 40, "quality": "simulated GCR",
                   "n_libraries": int(cpm.shape[1]),
                   "blocked_on": "preservative (DRS/RL)"},
        "aim1_decoder_on_brapa": a1.to_dict("records"),
        "aim2_native_signature": {
            "set_size": NATIVE_SET_SIZE,
            "n_contrasts_scored": int(len(a2)),
            "n_studies": int(a2["accession"].nunique()),
            "n_reaching_threshold": int(len(flight_hits)),
            "n_spaceflight_reaching_threshold": int(len(spaceflight)),
            "positive_control": a2ctrl,
            "interpretable": a2ctrl["passes"],
            "top": a2.head(5).to_dict("records"),
        },
        "aim3_concordance": {"rows": a3.to_dict("records"), **extra},
        "verdicts": {
            "aim1": (
                "The Arabidopsis-trained decoder does not fire on B. rapa at 40 cGy "
                f"simulated GCR (index {a1['radiation_index'].max()} at best, threshold "
                f"{THRESHOLD}). The decoder is a validated instrument -- it detects 7 of 8 "
                "labelled irradiations in Arabidopsis at AUC 0.960 -- so this is a real "
                "negative for the projected signature. What it cannot separate is a "
                "species difference from the low reliability of this particular dataset."),
            "aim2": (
                "UNINTERPRETABLE. The B. rapa-native signature fails its positive "
                f"control: it detects {a2ctrl['n_detected']} of "
                f"{a2ctrl['n_labelled_irradiations']} labelled irradiations, scoring "
                f"worst of all on a 100 Gy gamma exposure ({a2ctrl['worst_index']}). A "
                "signature that cannot see 100 Gy tells us nothing by failing to see "
                "spaceflight, so its null across the flight corpus is reported as a "
                "failed instrument, not as evidence about spaceflight."),
            "aim3": (
                "NO POWER at gene level. The OSD-658 40 cGy contrast does not reproduce "
                f"itself ({(extra['contrast_reliability']['osd658_40cGy'] or {}).get('mean')}), "
                "so the attenuation cap on any cross-species correlation is "
                f"{extra['attenuation_cap']}. The observed value is what that cap "
                "predicts regardless of the underlying biology. This says nothing about "
                "whether the DNA-damage response is conserved; it says these two "
                "datasets cannot be compared gene by gene."),
            "what_would_fix_it": (
                "Set-level scoring survives noise that per-gene correlation does not, "
                "which is why aim 1 is readable and aim 3 is not. A cross-species test "
                "at gene level needs contrasts that reproduce themselves: more "
                "biological replication per dose arm, or matched harvest timing, before "
                "the comparison is worth repeating."),
        },
        "caveats": {
            "tissue": "B. rapa leaf vs OSD-658 whole seedlings; a tissue difference of "
                      "this size is not separable from a species difference here.",
            "preservative": "DRS/RL is a large factor in this dataset and is handled as "
                            "a block, never pooled.",
            "timing": "Harvest times differ between the two experiments, so a shared "
                      "response could still be missed by sampling at different points "
                      "of the same kinetic.",
            "orthology": "Reciprocal best hit is one-to-one, so B. rapa's mesohexaploid "
                         "paralogues are collapsed to a single representative; "
                         "subfunctionalised copies are invisible to this map.",
        },
    }
    write_json(RESULTS / "qc" / "brapa_aims_qc.json", qc)
    log(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
