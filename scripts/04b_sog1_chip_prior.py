#!/usr/bin/env python3
"""Add SOG1 target edges to the TF prior, called from the study's own ChIP-seq.

Why this script exists: `04_fetch_tf_prior.py` reports that **SOG1 is absent from
DAP-seq**. SOG1 is the NAC master activator of the plant DNA-damage response and the
subject of the very study this cohort comes from — a DREM model built on the DAP-seq
prior alone physically cannot attribute any bifurcation to it. Without this step the
headline biological question is unanswerable.

The gap is filled from OSD-496 / GEO GSE112529, the SOG1-3xFLAG ChIP-seq that Bourbousse
et al. ran on the same seedlings as the RNA-seq. OSDR holds only raw reads and ISA for
OSD-496, so the processed coverage comes from GEO:

  GSM3072266/67  IP    SOG1-3xFLAG, 20 min / 1 h post-gamma
  GSM3072264/65  input SOG1-3xFLAG, 20 min / 1 h
  GSM3072268/69  IP    wild type (no FLAG epitope) — the specificity control

Peaks are called by binned enrichment rather than by shelling out to MACS: a bin is
retained when IP exceeds *both* its matched input and the no-FLAG WT IP by
`--min-fold`, which removes the two distinct artefact classes (chromatin accessibility,
and antibody background) that each control alone would miss. Retained bins are merged
and assigned to genes by TSS proximity, exactly as the DAP-seq edges are, so the two
sources produce commensurable scores.

  data/tf_prior/sog1_edges.tsv        TF, Gene, Input  (SOG1 only)
  data/tf_prior/tf_gene_edges.tsv     rewritten to include SOG1 unless --no-merge
  results/qc/sog1_chip_qc.json        measured bin, peak and target counts

  python3 scripts/04b_sog1_chip_prior.py
  python3 scripts/04b_sog1_chip_prior.py --bin 200 --min-fold 2.0 --window 1000
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import DATA, RESULTS, download, log, write_json  # noqa: E402

OUT = DATA / "tf_prior"
CHIP = DATA / "chip"
QC = RESULTS / "qc" / "sog1_chip_qc.json"

SOG1_LOCUS = "AT1G25580"
# Per-SAMPLE supplementary paths. The series-level suppl/ directory only offers the
# 241 MB GSE112529_RAW.tar; the individual bedGraphs live under geo/samples/.
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3072nnn/{gsm}/suppl/{fname}"
TRACKS = {
    "ip_20min": ("GSM3072266", "GSM3072266_IP_SOG1-3xFLAG_20min.bedGraph.gz"),
    "ip_1h": ("GSM3072267", "GSM3072267_IP_SOG1-3xFLAG_1h.bedGraph.gz"),
    "input_20min": ("GSM3072264", "GSM3072264_input_SOG1-3xFLAG_20min.bedGraph.gz"),
    "input_1h": ("GSM3072265", "GSM3072265_input_SOG1-3xFLAG_1h.bedGraph.gz"),
    "wt_20min": ("GSM3072268", "GSM3072268_IP_wt_20min.bedGraph.gz"),
    "wt_1h": ("GSM3072269", "GSM3072269_IP_wt_1h.bedGraph.gz"),
}

# TAIR10 chromosome lengths, used only to size the bin arrays.
CHROM_LEN = {"1": 30427671, "2": 19698289, "3": 23459830,
             "4": 18585056, "5": 26975502, "Mt": 366924, "Pt": 154478}


def norm_chrom(c: str) -> str:
    c = str(c).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return {"c": "Pt", "m": "Mt", "chrc": "Pt", "chrm": "Mt"}.get(c.lower(), c)


def binned_coverage(path: Path, bin_size: int) -> dict[str, np.ndarray]:
    """Sum bedGraph signal into fixed-width bins per chromosome.

    Streamed line by line: these files are 30-56 MB gzipped and tens of millions of
    intervals, far too large to hold as records.
    """
    arrays = {c: np.zeros(L // bin_size + 2, dtype=np.float64)
              for c, L in CHROM_LEN.items()}
    n_lines = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(("track", "#", "browser")):
                continue
            f = line.split()
            if len(f) < 4:
                continue
            chrom = norm_chrom(f[0])
            arr = arrays.get(chrom)
            if arr is None:
                continue
            try:
                start, end, val = int(f[1]), int(f[2]), float(f[3])
            except ValueError:
                continue
            n_lines += 1
            b0, b1 = start // bin_size, min(end // bin_size, len(arr) - 1)
            if b0 == b1:
                arr[b0] += val * (end - start)
            else:
                # Spread an interval spanning several bins across them by length.
                arr[b0] += val * ((b0 + 1) * bin_size - start)
                if b1 > b0 + 1:
                    arr[b0 + 1:b1] += val * bin_size
                arr[b1] += val * (end - b1 * bin_size)
    log(f"    {path.name}: {n_lines:,} intervals binned")
    return arrays


def normalise(arrays: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Scale to a common total so IP and input are comparable depth-for-depth."""
    total = sum(a.sum() for a in arrays.values())
    if total <= 0:
        return arrays
    factor = 1e9 / total
    return {c: a * factor for c, a in arrays.items()}


def call_peaks(ip, inp, wt, bin_size: int, min_fold: float, pseudo: float):
    """Bins enriched over BOTH the matched input and the no-FLAG WT IP."""
    peaks: dict[str, list[tuple[int, int, float]]] = {}
    stats = {"bins_tested": 0, "bins_over_input": 0, "bins_over_both": 0}
    for chrom in CHROM_LEN:
        a, b, c = ip.get(chrom), inp.get(chrom), wt.get(chrom)
        if a is None or b is None or c is None:
            continue
        n = min(len(a), len(b), len(c))
        a, b, c = a[:n], b[:n], c[:n]
        stats["bins_tested"] += n

        over_input = (a + pseudo) / (b + pseudo)
        over_wt = (a + pseudo) / (c + pseudo)
        stats["bins_over_input"] += int((over_input >= min_fold).sum())
        hit = (over_input >= min_fold) & (over_wt >= min_fold) & (a > 0)
        stats["bins_over_both"] += int(hit.sum())

        # Merge runs of adjacent enriched bins into single peaks.
        idx = np.flatnonzero(hit)
        if idx.size == 0:
            continue
        breaks = np.flatnonzero(np.diff(idx) > 1)
        starts = np.r_[idx[0], idx[breaks + 1]]
        ends = np.r_[idx[breaks], idx[-1]]
        peaks[chrom] = [
            (int(s) * bin_size, (int(e) + 1) * bin_size,
             float(np.max(over_input[int(s):int(e) + 1])))
            for s, e in zip(starts, ends)]
    return peaks, stats


def load_tss() -> dict[str, list[tuple[int, str]]]:
    path = OUT / "gene_tss.tsv"
    if not path.exists():
        raise SystemExit("run 04_fetch_tf_prior.py first — gene_tss.tsv is missing")
    by_chrom: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for line in path.read_text().splitlines()[1:]:
        gene, chrom, tss, _strand = line.split("\t")
        by_chrom[norm_chrom(chrom)].append((int(tss), gene))
    return {c: sorted(v) for c, v in by_chrom.items()}


def assign(peaks, tss_by_chrom, window: int) -> dict[str, float]:
    best: dict[str, float] = {}
    for chrom, plist in peaks.items():
        entries = tss_by_chrom.get(chrom)
        if not entries:
            continue
        positions = np.array([t for t, _ in entries])
        names = [g for _, g in entries]
        for start, end, score in plist:
            mid = (start + end) // 2
            i = int(np.searchsorted(positions, mid))
            for j in (i - 1, i):
                if 0 <= j < len(names) and abs(int(positions[j]) - mid) <= window:
                    if score > best.get(names[j], -np.inf):
                        best[names[j]] = score
    return best


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", type=int, default=200, help="bin width in bp (default 200)")
    ap.add_argument("--min-fold", type=float, default=2.0,
                    help="required enrichment over BOTH controls (default 2.0)")
    ap.add_argument("--window", type=int, default=1000,
                    help="max bp from peak centre to TSS (default 1000, matches DAP-seq)")
    ap.add_argument("--pseudo", type=float, default=1.0,
                    help="pseudocount added before ratios (default 1.0)")
    ap.add_argument("--no-merge", action="store_true",
                    help="write sog1_edges.tsv but leave tf_gene_edges.tsv untouched")
    args = ap.parse_args()

    CHIP.mkdir(parents=True, exist_ok=True)
    log("downloading SOG1 ChIP-seq coverage from GEO GSE112529 (~265 MB) ...")
    cov = {}
    for key, (gsm, fname) in TRACKS.items():
        path = download(GEO.format(gsm=gsm, fname=fname), CHIP / fname)
        log(f"  binning {key}")
        cov[key] = normalise(binned_coverage(path, args.bin))

    all_targets: dict[str, float] = {}
    per_time = {}
    for tp in ("20min", "1h"):
        peaks, stats = call_peaks(cov[f"ip_{tp}"], cov[f"input_{tp}"], cov[f"wt_{tp}"],
                                  args.bin, args.min_fold, args.pseudo)
        n_peaks = sum(len(v) for v in peaks.values())
        targets = assign(peaks, load_tss(), args.window)
        per_time[tp] = {**stats, "n_peaks": n_peaks, "n_target_genes": len(targets)}
        log(f"  {tp}: {n_peaks:,} peaks -> {len(targets):,} target genes")
        for g, s in targets.items():
            all_targets[g] = max(all_targets.get(g, 0.0), s)

    if not all_targets:
        log("no SOG1 targets called — the prior keeps its DAP-seq-only gap")
        write_json(QC, {"generated": dt.date.today().isoformat(), "per_timepoint": per_time,
                        "n_target_genes": 0, "merged": False})
        return 1

    # Rank-normalise to [0,1] so SOG1 edges sit on the same scale as the DAP-seq ones.
    genes = list(all_targets)
    ranks = np.argsort(np.argsort([all_targets[g] for g in genes]))
    denom = max(len(genes) - 1, 1)
    scored = {g: round(float(r) / denom, 4) for g, r in zip(genes, ranks)}

    sog1_path = OUT / "sog1_edges.tsv"
    with open(sog1_path, "w") as fh:
        fh.write("TF\tGene\tInput\n")
        for g, s in sorted(scored.items()):
            fh.write(f"{SOG1_LOCUS}\t{g}\t{s:g}\n")
    log(f"wrote {sog1_path}  ({len(scored)} SOG1 target genes)")

    merged = False
    if not args.no_merge:
        main_path = OUT / "tf_gene_edges.tsv"
        existing = [l for l in main_path.read_text().splitlines()[1:]
                    if l and not l.startswith(SOG1_LOCUS)]
        new = [f"{SOG1_LOCUS}\t{g}\t{s:g}" for g, s in sorted(scored.items())]
        main_path.write_text("TF\tGene\tInput\n" + "\n".join(existing + new) + "\n")
        merged = True
        log(f"merged into {main_path}  ({len(existing) + len(new)} total edges)")

    write_json(QC, {
        "generated": dt.date.today().isoformat(),
        "source": {"osdr": "OSD-496", "geo": "GSE112529",
                   "tracks": {k: v[1] for k, v in TRACKS.items()}},
        "parameters": {"bin_bp": args.bin, "min_fold": args.min_fold,
                       "window_bp": args.window, "pseudocount": args.pseudo},
        "per_timepoint": per_time,
        "n_target_genes": len(scored),
        "merged_into_main_prior": merged,
        "note": "Peaks require enrichment over both the matched input and the no-FLAG "
                "wild-type IP, so accessibility and antibody-background artefacts are "
                "excluded by construction.",
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
