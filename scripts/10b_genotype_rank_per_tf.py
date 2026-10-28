#!/usr/bin/env python3
"""Rank SOG1 and MYB3R1 among transcription factors, not among TF-branch scores.

`10_compare_models.py` ranks every row of a DREM split table. For a `vsOthers` split
DREM writes one row per TF *per branch*, so a three-way split has 3 x 479 = 1,437 rows,
and a rank taken over those rows can exceed the number of TFs. The genotype contrast
reported SOG1 at rank 574 in sog1-1 (weighted prior) while the manuscript described it
as a rank "among 478 TFs", which is not a number that can exist.

This script reads the same parsed DREM tables (no DREM run) and the same test type
10_compare_models chose for each prior, collapses each split to one score per TF (the
TF's best, i.e. lowest, score across that split's branches), ranks TFs within the split
(method "min", lower = more significant), and reports each watched TF's best rank
across splits. For `2way` splits with one row per TF this equals 10's rank.

  results/comparison/genotype_rank_per_tf.tsv

  python3 scripts/10b_genotype_rank_per_tf.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import RESULTS, log  # noqa: E402

PARSED = RESULTS / "drem" / "parsed"
OUT = RESULTS / "comparison"
WATCHED = {"AT1G25580": "SOG1", "AT4G32730": "MYB3R1"}
ARMS = {"wt": "A_primary_WildType", "sog1": "A_primary_sog1-1"}


def per_tf_rank(run: str, kind: str, locus: str) -> tuple[int | None, int, int]:
    """Best per-TF rank of `locus`, the TF count in that split, and its row count."""
    df = pd.read_csv(PARSED / f"{run}_tfscores.tsv", sep="\t")
    df = df[df["kind"] == kind]
    best = df.groupby(["split", "TF"])["score"].min().reset_index()
    best["rank"] = best.groupby("split")["score"].rank(method="min", ascending=True)
    hit = best[best["TF"] == locus].sort_values("rank")
    if hit.empty:
        return None, 0, 0
    split = hit.iloc[0]["split"]
    n_tfs = int((best["split"] == split).sum())
    n_rows = int((df["split"] == split).sum())
    return int(hit.iloc[0]["rank"]), n_tfs, n_rows


def main() -> int:
    contrast = pd.read_csv(OUT / "genotype_contrast.tsv", sep="\t")
    kinds = contrast.drop_duplicates("prior").set_index("prior")["test_type"].to_dict()
    rows = []
    for prior, kind in kinds.items():
        for locus, name in WATCHED.items():
            row = {"prior": prior, "test_type": kind, "TF": name, "locus": locus}
            for tag, arm in ARMS.items():
                rank, n_tfs, n_rows = per_tf_rank(f"{arm}__{prior}", kind, locus)
                row[f"{tag}_rank_among_tfs"] = rank
                row[f"{tag}_n_tfs_in_split"] = n_tfs
                row[f"{tag}_rows_in_split"] = n_rows
            rows.append(row)
            log(f"{prior:8s} {kind:8s} {name:6s} WT {row['wt_rank_among_tfs']}/"
                f"{row['wt_n_tfs_in_split']}  sog1-1 {row['sog1_rank_among_tfs']}/"
                f"{row['sog1_n_tfs_in_split']}")
    out = OUT / "genotype_rank_per_tf.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    log(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
