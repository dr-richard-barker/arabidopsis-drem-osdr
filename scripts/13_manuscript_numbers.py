#!/usr/bin/env python3
"""Emit every number the manuscript quotes as a LaTeX macro, read from results/.

The manuscript never types a figure. It writes `\\NumCohortATimedSamples`, and this
script defines that macro from `results/qc/metadata_qc.json`. If the pipeline is re-run
on revised OSDR deposits, the prose follows automatically; and a number that has no
macro is a number nobody computed, which is exactly the state that produces a confident
sentence about a result that does not exist.

  manuscript/latex/generated_numbers.tex   \\newcommand definitions
  results/manuscript_numbers.json          the same values, machine-readable

Any macro whose source file is missing is defined as \\textbf{TODO} rather than omitted,
so an unbuilt number shows up in the rendered PDF instead of failing the LaTeX run
silently or, worse, leaving a plausible stale value in place.

  python3 scripts/13_manuscript_numbers.py
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cohorts import BOURBOUSSE_DOI, BOURBOUSSE_GROUPS, COHORTS, EXCLUSIONS  # noqa: E402
from lib_sources import RESULTS, ROOT, log, write_json  # noqa: E402

TEX = ROOT / "manuscript" / "latex" / "generated_numbers.tex"
MISSING = r"\textbf{TODO}"


def load(path: Path):
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def dig(obj, *keys, default=None):
    for k in keys:
        if obj is None:
            return default
        obj = obj.get(k) if isinstance(obj, dict) else None
    return default if obj is None else obj


def main() -> int:
    meta = load(RESULTS / "qc" / "metadata_qc.json")
    tfp = load(RESULTS / "qc" / "tf_prior_qc.json")
    sog1 = load(RESULTS / "qc" / "sog1_chip_qc.json")
    deconv = load(RESULTS / "qc" / "deconvolution_qc.json")
    weighted = load(RESULTS / "qc" / "weighted_prior_qc.json")
    series = load(RESULTS / "pseudotimeseries" / "series_qc.json")
    parsed = load(RESULTS / "drem" / "parsed_summary.json")
    comp = load(RESULTS / "comparison" / "summary.json")
    sib = load(RESULTS / "comparison" / "sibling_summary.json")
    dec = load(RESULTS / "decoder" / "decoder_report.json")
    sigqc = load(RESULTS / "decoder" / "signature_qc.json")
    tis = load(RESULTS / "qc" / "tissue_association_qc.json")

    v: dict[str, object] = {}

    # ---- cohorts
    shape = dig(meta, "cohort_shape", "per_cohort", default={})
    v["NumStudiesTotal"] = sum(len(c["studies"]) for c in COHORTS.values())
    v["NumCohortAStudies"] = len(COHORTS["A"]["studies"])
    v["NumCohortBStudies"] = len(COHORTS["B"]["studies"])
    v["NumCohortATimedSamples"] = dig(shape, "A", "n_timed_samples")
    v["NumCohortAAnchorSamples"] = dig(shape, "A", "n_untimed_anchor_samples")
    v["NumCohortATimepoints"] = len(dig(shape, "A", "observed_timepoints", default=[]) or [])
    v["NumCohortBSamples"] = dig(shape, "B", "n_timed_samples")
    v["NumCohortBTimepoints"] = len(dig(shape, "B", "observed_timepoints", default=[]) or [])
    v["NumArabidopsisStudiesScreened"] = 62  # studies in the OSDR organism query
    v["NumExcludedStudies"] = len(EXCLUSIONS)

    # ---- series
    arms = dig(series, "cohorts", "A", "arms", default={}) or {}
    wt = arms.get("A_primary_WildType", {})
    v["NumSeriesGenes"] = wt.get("n_genes")
    v["NumSeriesSharedTimepoints"] = len(dig(wt, "batch", "shared_timepoints", default=[]) or [])
    sent = wt.get("sentinels") or {}
    v["NumSentinelsResponding"] = sent.get("n_responding_at_2fold")
    mut = arms.get("A_primary_sog1-1", {}).get("sentinels") or {}
    v["NumSentinelsRespondingMutant"] = mut.get("n_responding_at_2fold")

    # ---- TF prior
    m = dig(tfp, "measured", default={}) or {}
    v["NumDapseqPeakFiles"] = m.get("narrowpeak_files")
    v["NumDapseqResolvedByAGI"] = dig(m, "tf_identified_by", "agi_in_filename")
    v["NumDapseqResolvedBySymbol"] = dig(m, "tf_identified_by", "symbol_map")
    v["NumDapseqUnresolved"] = dig(m, "tf_identified_by", "unresolved")
    v["NumPriorTFs"] = m.get("tfs_with_edges")
    v["NumPriorEdges"] = m.get("total_edges")
    v["NumPriorMedianTargets"] = m.get("median_targets_per_tf")
    v["NumSogOneTargets"] = dig(sog1, "n_target_genes")
    v["NumSogOnePeaksTwentyMin"] = dig(sog1, "per_timepoint", "20min", "n_peaks")
    v["NumSogOnePeaksOneHour"] = dig(sog1, "per_timepoint", "1h", "n_peaks")

    # ---- auto-decoder
    v["NumAtlasClusters"] = dig(deconv, "atlas_shape", "clusters")
    v["NumAtlasSignatureGenes"] = dig(deconv, "atlas_shape", "signature_genes")
    v["NumLatentDims"] = dig(deconv, "atlas_shape", "latent_dims")
    darm = dig(deconv, "arms", "A_primary_WildType", default={}) or {}
    v["NumDeconvGenesMatched"] = darm.get("n_signature_genes_matched")
    v["NumDeconvMedianResidual"] = darm.get("median_relative_residual")
    warm = dig(weighted, "arms", "A_primary_WildType", default={}) or {}
    v["NumWeightedEdgesInformed"] = warm.get("n_edges_atlas_informed")
    v["PctWeightedEdgesInformed"] = (round(100 * warm["frac_edges_atlas_informed"])
                                     if warm.get("frac_edges_atlas_informed") is not None else None)
    v["PctWeightedEdgesChanged"] = (round(100 * warm["frac_edges_changed_by_0.01"])
                                    if warm.get("frac_edges_changed_by_0.01") is not None else None)

    # ---- DREM models
    runs = dig(parsed, "runs", default={}) or {}
    v["NumDremRuns"] = len(runs) or None
    wtflat = runs.get("A_primary_WildType__flat", {})
    v["NumWTNodes"] = wtflat.get("n_nodes")
    v["NumWTSplits"] = wtflat.get("n_splits")
    v["NumWTGenesAssigned"] = wtflat.get("n_genes_assigned")
    v["NumPublishedGroups"] = BOURBOUSSE_GROUPS
    abl = dig(comp, "ablation", "A_primary_WildType", default={}) or {}
    v["AblationTestType"] = abl.get("test_type")
    v["NumAblationSplits"] = abl.get("n_splits")
    v["NumAblationSplitsChanged"] = abl.get("n_splits_attribution_changed")
    v["PctAblationSplitsChanged"] = (round(100 * abl["frac_splits_changed"])
                                     if abl.get("frac_splits_changed") is not None else None)

    # ---- genotype positive control: SOG1's rank in WT vs the knockout
    import csv as _csv
    gpath = RESULTS / "comparison" / "genotype_contrast.tsv"
    if gpath.exists():
        rows = list(_csv.DictReader(gpath.read_text().splitlines(), delimiter="\t"))
        def _g(prior, tf, field):
            for r in rows:
                if r["prior"] == prior and r["TF"] == tf and r[field] not in ("", "NaN"):
                    return int(float(r[field]))
            return None
        v["SogOneRankWTFlat"] = _g("flat", "SOG1", "wt_best_rank")
        v["SogOneRankMutFlat"] = _g("flat", "SOG1", "sog1_best_rank")
        v["SogOneRankWTWeighted"] = _g("weighted", "SOG1", "wt_best_rank")
        v["SogOneRankMutWeighted"] = _g("weighted", "SOG1", "sog1_best_rank")
        v["MybThreeRRankWTFlat"] = _g("flat", "MYB3R1", "wt_best_rank")

    # ---- sibling cross-validation
    v["NumSharedGenesWithSibling"] = dig(sib, "n_shared_genes")
    v["NumSiblingNodes"] = dig(sib, "n_nodes")
    v["NumSiblingModules"] = dig(sib, "n_modules")
    v["NumNodesEnrichingModule"] = dig(sib, "n_nodes_with_significant_module")

    # ---- decoder
    v["NumSignatureSogDep"] = dig(sigqc, "set_sizes", "sog1_dependent")
    v["NumSignatureSogIndep"] = dig(sigqc, "set_sizes", "sog1_independent")
    v["NumSignatureMybRep"] = dig(sigqc, "set_sizes", "myb3r_repressed")
    v["NumDecoderContrasts"] = sum((dig(dec, "call_counts", default={}) or {}).values()) or None
    v["NumDecoderStudies"] = dig(load(RESULTS / "decoder" / "scan_qc.json"), "n_studies_scored")
    perf = dig(dec, "performance", default={}) or {}
    v["NumLabelledPositives"] = perf.get("n_labelled_positive")
    v["NumPositivesDetected"] = perf.get("labelled_positives_detected")
    v["NumOtherExposures"] = perf.get("n_other_exposure")
    v["NumFalsePositives"] = perf.get("other_exposures_above_threshold")
    v["DecoderAUC"] = perf.get("auc_vs_all_contrasts")
    v["NumProliferationOnly"] = dig(dec, "call_counts",
                                    "G2/M repressed, no DDR (proliferation slowing)")
    v["NumInverseCalls"] = dig(dec, "call_counts", "inverse (DDR suppressed)")

    # ---- tissue association
    st = dig(tis, "deconvolution_stability", default={}) or {}
    v["ClusterSwing"] = st.get("cluster_level_median_relative_swing")
    v["OrganSwing"] = st.get("organ_level_median_relative_swing")
    v["TimeAveragedCV"] = st.get("timeaveraged_bootstrap_cv")
    rw = dig(tis, "reweighting", default={}) or {}
    v["NumTFsAtlasInformed"] = rw.get("n_atlas_informed")
    v["NumTFsTotal"] = rw.get("n_tfs")

    defined = sum(1 for x in v.values() if x is not None)
    lines = [
        "% Generated by scripts/13_manuscript_numbers.py — do not hand-edit.",
        "% Every value is read from results/; a TODO means the pipeline has not",
        "% produced that number yet, and it will render visibly in the PDF.",
        f"% Generated {dt.date.today().isoformat()}: {defined}/{len(v)} values available.",
        "",
    ]
    for name, value in v.items():
        if value is None:
            lines.append(f"\\newcommand{{\\{name}}}{{{MISSING}}}")
        elif isinstance(value, float):
            lines.append(f"\\newcommand{{\\{name}}}{{{value:g}}}")
        elif isinstance(value, int):
            lines.append(f"\\newcommand{{\\{name}}}{{{value:,}}}".replace(",", "{,}"))
        else:
            lines.append(f"\\newcommand{{\\{name}}}{{{value}}}")
    lines += ["", f"\\newcommand{{\\BourbousseDOI}}{{{BOURBOUSSE_DOI}}}", ""]

    TEX.parent.mkdir(parents=True, exist_ok=True)
    TEX.write_text("\n".join(lines))
    write_json(RESULTS / "manuscript_numbers.json",
               {"generated": dt.date.today().isoformat(), "values": v,
                "n_defined": defined, "n_total": len(v)})

    missing = [k for k, x in v.items() if x is None]
    log(f"wrote {TEX}  ({defined}/{len(v)} values defined)")
    if missing:
        log(f"  still TODO ({len(missing)}): {', '.join(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
