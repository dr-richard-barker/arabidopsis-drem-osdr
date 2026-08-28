#!/usr/bin/env python3
"""Write the expression, repeats and settings files DREM's batch mode consumes.

One run is produced per (arm x prior) combination, because the comparison between
priors on an otherwise identical run is the experiment:

  results/drem/inputs/<arm>__<prior>/
      expression.txt   Gene + timepoint columns, log2 ratios, t=0 column prepended
      repeats.txt      the per-replicate matrix in the same shape
      settings.txt     DREM 2.0 batch settings

Two choices here matter and are easy to get silently wrong:

*t=0 column.* `03_build_pseudotimeseries.py` expresses every timepoint as a log2 ratio
against the untreated control, so the untreated state is by definition 0. DREM anchors
trajectories at the first column, so a literal zero column is prepended: without it the
earliest measured timepoint (10 min, already strongly induced) would be treated as the
origin and the entire early response would vanish from the model.

*Normalisation mode.* For the same reason the mode is `No normalization/add 0` rather
than DREM's usual `Normalize data`. The other modes re-anchor on the first column's
values, which would double-subtract a baseline that has already been removed. This mode
requires log-space input, which is what the series is.

The random seed is pinned so that a difference between two runs can only come from the
prior, never from the search.

  python3 scripts/07_write_drem_inputs.py
  python3 scripts/07_write_drem_inputs.py --min-abs-log-ratio 1.0 --node-penalty 40
"""
from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import DATA, RESULTS, log, write_json  # noqa: E402

SERIES = RESULTS / "pseudotimeseries"
PRIORS = DATA / "tf_prior" / "weighted"
OUT = RESULTS / "drem" / "inputs"

# Settings DREM 2.0 reads in batch mode. Only the keys we deliberately set are listed;
# the template is written in full so a run is reproducible from the file alone.
SETTINGS_TEMPLATE = """#Main Input:
TF-gene_Interaction_Source\tUser Provided
TF-gene_Interactions_File\t{tf_file}
Expression_Data_File\t{expression_file}
Saved_Model_File
Gene_Annotation_Source\tNo annotations
Gene_Annotation_File\t
Cross_Reference_Source\tNo cross references
Cross_Reference_File\t
Normalize_Data[Log normalize data,Normalize data,No normalization/add 0]\tNo normalization/add 0
Spot_IDs_in_the_data_file\tfalse

#Repeat Data:
Repeat_Data_Files(comma delimited list)\t{repeat_file}
Repeat_Data_is_from[Different time periods,The same time period]\tThe same time period

#Filtering:
Filter_Gene_If_It_Has_No_Static_Input_Data\tfalse
Maximum_Number_of_Missing_Values\t0
Minimum_Correlation_between_Repeats\t0
Minimum_Absolute_Log_Ratio_Expression\t{min_abs_log_ratio}
Change_should_be_based_on[Maximum-Minimum,Difference From 0]\tDifference From 0
Pre-filtered_Gene_File

#Search Options
Allow_Path_Merges\tfalse
Maximum_number_of_paths_out_of_split\t{max_paths}
Use_transcription_factor-gene_interaction_data_to_build\ttrue
Saved_Model[Use As Is/Start Search From/Do Not Use]\tDo Not Use
Convergence_Likelihood_%\t0.01
Minimum_Standard_Deviation\t0.0

#Model Selection Options
Model_Selection_Framework[Penalized Likelihood,Train-Test]\tPenalized Likelihood
Penalized_Likelihood_Node_Penalty\t{node_penalty}
Random_Seed\t{seed}
Main_search_score_%\t0
Main_search_difference_threshold\t0
Delete_path_score_%\t0.15
Delete_path_difference_threshold\t0
Delay_split_score_%\t0.15
Delay_split_difference_threshold\t0
Merge_path_score_%\t0.15
Merge_path_difference_threshold\t0

#Gene Annotations:
Include_Biological_Process\tfalse
Include_Molecular_Function\tfalse
Include_Cellular_Process\tfalse
Only_include_annotations_with_these_evidence_codes\t
Only_include_annotations_with_these_taxon_IDs
Category_ID_file

#GO Analysis
Minimum_GO_level\t3
Minimum_number_of_genes\t5
Number_of_samples_for_randomized_multiple_hypothesis_correction\t500
Multiple_hypothesis_correction_method_enrichment[Bonferroni,Randomization]\tRandomization

#Expression Scaling Options
Regulator_Types_Used_For_Activity_Scoring\tBoth
Expression_Scaling_Weight\t1.0
Minimum_TF_Expression_After_Scaling\t0.5

#Interface
X-axis_Scale_Factor\t1
Y-axis_Scale_Factor\t1.2
X-axis_scale[Uniform,Based on Real Time]\tBased on Real Time
Key_Input_X_p-val_10^-X\t3
Minimum_Split_Percent\t0
Scale_Node_Areas_By_The_Factor\t1
Key_Input_Significance_Based_On[Path Significance Conditional on Split,Path Significance Overall,Split Significance]\tPath Significance Conditional on Split
"""


def with_zero_anchor(df: pd.DataFrame) -> pd.DataFrame:
    """Prepend an explicit t=0 column of zeros — the untreated control state."""
    out = df.copy()
    ordered = sorted(out.columns, key=lambda c: float(c))
    out = out[ordered]
    out.insert(0, "0", 0.0)
    return out


def choose_filter(expr: pd.DataFrame, target: int, floor: float, ceiling: float) -> tuple[float, int]:
    """Pick the |log2| threshold that admits about `target` genes.

    A single fixed threshold cannot serve both cohorts. The arms differ by orders of
    magnitude in response size: at |log2| >= 1 the wild-type irradiation arm admits 3,073
    genes and exhausts an 8 GB JVM heap during DREM's model search, while the same
    threshold on the spaceflight cohort admits 870 and at |log2| >= 2 only 75 — too few
    for a bifurcation to mean anything.

    Targeting a gene count instead keeps model complexity comparable across arms, which
    is also what makes cross-arm comparison fair. The chosen threshold is recorded per
    run. Within an arm the three priors see identical expression data, so they get an
    identical threshold and the ablation stays controlled.
    """
    peak = expr.abs().max(axis=1)
    lo, hi = floor, ceiling
    best = floor
    for _ in range(40):
        mid = (lo + hi) / 2
        n = int((peak >= mid).sum())
        if n > target:
            lo = mid
        else:
            hi = mid
        best = mid
    n = int((peak >= best).sum())
    # Never go below the caller's floor, even if that admits more than the target.
    if best < floor:
        best, n = floor, int((peak >= floor).sum())
    return round(best, 3), n


def split_replicates(reps: pd.DataFrame) -> list[pd.DataFrame]:
    """Split the wide per-replicate matrix into one full matrix per replicate.

    DREM does not take a single file holding every replicate column: it requires each
    repeat to be its own gene x timepoint matrix with *exactly* the same columns as the
    main expression file, supplied as a comma-delimited list. Feeding it the wide form
    fails with "Repeat data set must have same number of columns as original".

    Only as many replicate files are emitted as the least-replicated timepoint supports,
    so every emitted file is complete; emitting a ragged one would either drop a
    timepoint or silently recycle another replicate's values.
    """
    by_tp: dict[str, list[pd.Series]] = {}
    for i, col in enumerate(reps.columns):
        # 03 labels replicate columns "<timepoint>|<n>"; fall back to the bare name for
        # a single-replicate timepoint.
        tp = str(col).split("|")[0]
        by_tp.setdefault(tp, []).append(reps.iloc[:, i])
    if not by_tp:
        return []
    k = min(len(v) for v in by_tp.values())
    ordered = sorted(by_tp, key=lambda c: float(c))
    out = []
    for i in range(k):
        frame = pd.DataFrame({tp: by_tp[tp][i] for tp in ordered})
        out.append(frame)
    return out


def priors_for(slug: str) -> dict[str, Path]:
    """The three priors this arm is run against."""
    found = {"flat": PRIORS / "prior_flat.tsv",
             "binding": PRIORS / "prior_binding.tsv",
             "weighted": PRIORS / f"prior_weighted_{slug}.tsv"}
    return {k: v for k, v in found.items() if v.exists()}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target-genes", type=int, default=1200,
                    help="aim for about this many genes past DREM's expression filter; "
                         "the threshold is solved per arm (default 1200)")
    ap.add_argument("--min-abs-log-ratio", type=float, default=1.0,
                    help="floor for the solved threshold (default 1.0 = 2-fold)")
    ap.add_argument("--max-abs-log-ratio", type=float, default=6.0,
                    help="ceiling for the solved threshold (default 6.0)")
    ap.add_argument("--node-penalty", type=int, default=40,
                    help="penalised-likelihood node penalty (default 40, DREM's own)")
    ap.add_argument("--max-paths", type=int, default=3,
                    help="maximum paths out of a split (default 3)")
    ap.add_argument("--seed", type=int, default=1260,
                    help="random seed, pinned so prior is the only difference between runs")
    args = ap.parse_args()

    series = sorted(SERIES.glob("*_expression.tsv"))
    if not series:
        raise SystemExit("run 03_build_pseudotimeseries.py first")

    report = {"generated": dt.date.today().isoformat(),
              "parameters": vars(args) | {"seed": args.seed},
              "runs": {}}

    for expr_path in series:
        slug = expr_path.name.replace("_expression.tsv", "")
        expr = with_zero_anchor(pd.read_csv(expr_path, sep="\t", index_col=0))
        threshold, n_pass = choose_filter(expr, args.target_genes,
                                          args.min_abs_log_ratio, args.max_abs_log_ratio)
        rep_path = SERIES / f"{slug}_repeats.tsv"
        rep_frames = (split_replicates(pd.read_csv(rep_path, sep="\t", index_col=0))
                      if rep_path.exists() else [])
        rep_frames = [with_zero_anchor(r) for r in rep_frames]

        for prior_name, prior_path in priors_for(slug).items():
            run_dir = OUT / f"{slug}__{prior_name}"
            run_dir.mkdir(parents=True, exist_ok=True)

            expr.round(4).to_csv(run_dir / "expression.txt", sep="\t", index_label="Gene")
            rep_paths = []
            for i, frame in enumerate(rep_frames, start=1):
                rp = run_dir / f"repeat{i}.txt"
                frame.round(4).to_csv(rp, sep="\t", index_label="Gene")
                rep_paths.append(str(rp.resolve()))
            repeat_arg = ",".join(rep_paths)

            # DREM resolves paths relative to its working directory; absolute paths in
            # the settings file remove any doubt about where it looked.
            (run_dir / "settings.txt").write_text(SETTINGS_TEMPLATE.format(
                tf_file=prior_path.resolve(),
                expression_file=(run_dir / "expression.txt").resolve(),
                repeat_file=repeat_arg,
                min_abs_log_ratio=threshold,
                max_paths=args.max_paths,
                node_penalty=args.node_penalty,
                seed=args.seed))

            report["runs"][run_dir.name] = {
                "arm": slug,
                "prior": prior_name,
                "prior_file": str(prior_path.relative_to(DATA.parent)),
                "n_genes": int(expr.shape[0]),
                "expression_filter_log2": threshold,
                "n_genes_past_filter": n_pass,
                "timepoints": list(expr.columns),
                "n_replicate_files": len(rep_paths),
            }
            log(f"  {run_dir.name}: {expr.shape[0]} genes x {expr.shape[1]} columns, "
                f"filter |log2|>={threshold} admits {n_pass}")

    write_json(RESULTS / "drem" / "inputs_manifest.json", report)
    log(f"\n{len(report['runs'])} DREM runs prepared in {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
