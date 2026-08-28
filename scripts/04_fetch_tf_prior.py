#!/usr/bin/env python3
"""Build the TF -> target-gene prior DREM needs, from public Arabidopsis resources.

DREM's second input is a static TF-gene interaction table. This builds one from
DAP-seq, which is the genome-wide in-vitro binding atlas for Arabidopsis:

  1. PlantTFDB `Ath_TF_list.txt.gz` -> the TF universe (locus -> family, gene name).
  2. Ensembl Plants release-54 `Arabidopsis_thaliana.TAIR10.54.gtf.gz` -> gene models.
     This is deliberately the *same* annotation GeneLab used to produce the count
     matrices (stated in the OSDR processing protocol), so peak-to-gene assignment and
     expression share one coordinate system and one locus vocabulary.
  3. GEO `GSE60143_RAW.tar` (77 MB) -> per-TF narrowPeak files from O'Malley et al.
     2016. The assayed TF's AGI locus is encoded in each filename.

A peak is assigned to a gene when it falls within `--window` bp of that gene's TSS
(default 3000, upstream-weighted). The edge score is the peak's signalValue, rank-
normalised to [0,1] within each TF so a TF with deep sequencing does not dominate one
with shallow — DREM reads column 3 as a score, and an unnormalised score would encode
library depth rather than binding confidence.

  data/tf_prior/tf_list.tsv          TF locus, family, name
  data/tf_prior/gene_tss.tsv         locus, chrom, TSS, strand
  data/tf_prior/tf_gene_edges.tsv    TF, Gene, Input   <- DREM's format, flat prior
  results/qc/tf_prior_qc.json        measured counts at every stage

Everything reported here is measured from the downloads. No edge or TF count is quoted
from a publication.

  python3 scripts/04_fetch_tf_prior.py
  python3 scripts/04_fetch_tf_prior.py --window 5000 --max-peaks-per-tf 5000
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import re
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import DATA, RESULTS, download, get_json, log, write_json  # noqa: E402

OUT = DATA / "tf_prior"
QC = RESULTS / "qc" / "tf_prior_qc.json"

PLANTTFDB_TF_LIST = "https://planttfdb.gao-lab.org/download/TF_list/Ath_TF_list.txt.gz"
ENSEMBL_GTF = ("https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-54/gtf/"
               "arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.54.gtf.gz")
DAPSEQ_TAR = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60143/suppl/"
              "GSE60143_RAW.tar")
TAIR_ALIASES = ("https://www.arabidopsis.org/api/download-files/download?filePath="
                "Public_Data_Releases/TAIR_Data_20140331/gene_aliases_20140331.txt")
CONNECTF_TFS = "https://connectf.org/api/tfs/"

AGI = re.compile(r"AT[1-5CM]G\d{5}", re.I)

# Only 194 of the 814 DAP-seq peak files name the assayed TF by AGI locus. The other
# 620 use a gene symbol (WRKY18, ANAC057, MYB33) — including most of the NAC, MYB, WRKY
# and bZIP families, which is exactly where the DNA-damage regulators live. Resolving
# symbols is therefore not a nicety: skipping them would silently discard three
# quarters of the atlas and the entire repressor arm of this study's biology.
# Filenames look like  GSM1925556_DAPSeq-NAC_tnt-ANAC057_col_a-chr1-5_GEM_events...
DAP_TF_TOKEN = re.compile(r"_DAPSeq-[^_]+_[A-Za-z]+-([^_]+)_", re.I)

# Regulators this study's biology stands or falls on. Their presence in the prior is
# reported every run: DAP-seq assayed a large but incomplete TF panel, and a DREM model
# simply cannot attribute a split to a TF that has no edges. SOG1 in particular is
# absent from DAP-seq, which is why `04b_sog1_chip_prior.py` exists.
KEY_REGULATORS = {
    "AT1G25580": "SOG1",      # NAC master activator of the DDR
    "AT4G32730": "MYB3R1",
    "AT3G09370": "MYB3R3",
    "AT5G11510": "MYB3R4",
    "AT3G01600": "ANAC044",
    "AT5G14490": "ANAC085",
    "AT5G13330": "ERF115",
}


def fetch_tf_list() -> dict[str, dict]:
    path = download(PLANTTFDB_TF_LIST, OUT / "Ath_TF_list.txt.gz")
    tfs: dict[str, dict] = {}
    with gzip.open(path, "rt") as fh:
        header = next(fh).rstrip("\n").split("\t")
        cols = {c.strip().lower(): i for i, c in enumerate(header)}
        gene_i = cols.get("gene_id", 1)
        fam_i = cols.get("family", 2)
        name_i = cols.get("gene_name", gene_i)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(gene_i, fam_i):
                continue
            locus = f[gene_i].split(".")[0].upper()
            if not AGI.fullmatch(locus):
                continue
            tfs[locus] = {"family": f[fam_i],
                          "name": f[name_i] if name_i < len(f) else ""}
    (OUT / "tf_list.tsv").write_text(
        "TF\tfamily\tname\n"
        + "".join(f"{k}\t{v['family']}\t{v['name']}\n" for k, v in sorted(tfs.items())))
    log(f"PlantTFDB: {len(tfs)} Arabidopsis TF loci")
    return tfs


def build_symbol_map(gtf_path: Path) -> tuple[dict[str, str], dict]:
    """symbol -> AGI locus, from three complementary sources.

    None is sufficient alone: the Ensembl GTF carries one current symbol per gene and
    misses older aliases; TAIR's 2014 alias table has the historical names (SOG1,
    ANAC057, CBF3) but predates newer ones; ConnecTF adds a few TF-specific entries.
    Earlier sources win, so the annotation matching the count matrices takes precedence.
    """
    sym: dict[str, str] = {}
    provenance = {"ensembl_gtf": 0, "tair_aliases": 0, "connectf": 0}

    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gid = re.search(r'gene_id "([^"]+)"', f[8])
            gname = re.search(r'gene_name "([^"]+)"', f[8])
            if gid and gname and gname.group(1).upper() not in sym:
                sym[gname.group(1).upper()] = gid.group(1).split(".")[0].upper()
                provenance["ensembl_gtf"] += 1

    try:
        alias_path = download(TAIR_ALIASES, OUT / "tair_gene_aliases.txt")
        for line in alias_path.read_text(errors="replace").splitlines()[1:]:
            f = line.split("\t")
            if len(f) >= 2 and AGI.fullmatch(f[0].strip()):
                key = f[1].upper().strip()
                if key and key not in sym:
                    sym[key] = f[0].upper().strip()
                    provenance["tair_aliases"] += 1
    except RuntimeError as e:
        log(f"  TAIR aliases unavailable ({e}) — continuing without them")

    try:
        for entry in get_json(CONNECTF_TFS):
            name = (entry.get("name") or "").upper().strip()
            locus = (entry.get("value") or "").upper().strip()
            if name and AGI.fullmatch(locus) and name not in sym:
                sym[name] = locus
                provenance["connectf"] += 1
    except RuntimeError as e:
        log(f"  ConnecTF unavailable ({e}) — continuing without it")

    log(f"symbol map: {len(sym)} aliases "
        f"(Ensembl {provenance['ensembl_gtf']}, TAIR {provenance['tair_aliases']}, "
        f"ConnecTF {provenance['connectf']})")
    return sym, provenance


def fetch_gene_tss() -> tuple[dict[str, tuple[str, int, str]], Path]:
    """locus -> (chrom, TSS, strand), from the gene features of the Ensembl GTF."""
    path = download(ENSEMBL_GTF, OUT / "Arabidopsis_thaliana.TAIR10.54.gtf.gz")
    genes: dict[str, tuple[str, int, str]] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = re.search(r'gene_id "([^"]+)"', f[8])
            if not m:
                continue
            locus = m.group(1).split(".")[0].upper()
            start, end, strand = int(f[3]), int(f[4]), f[6]
            genes[locus] = (f[0], start if strand == "+" else end, strand)
    (OUT / "gene_tss.tsv").write_text(
        "gene\tchrom\ttss\tstrand\n"
        + "".join(f"{g}\t{c}\t{t}\t{s}\n" for g, (c, t, s) in sorted(genes.items())))
    log(f"Ensembl TAIR10.54: {len(genes)} gene models")
    return genes, path


def tss_index(genes: dict[str, tuple[str, int, str]]):
    """Per-chromosome sorted TSS arrays for vectorised nearest-gene lookup."""
    by_chrom: dict[str, list[tuple[int, str, str]]] = defaultdict(list)
    for locus, (chrom, tss, strand) in genes.items():
        by_chrom[norm_chrom(chrom)].append((tss, locus, strand))
    return {c: (np.array([t for t, _, _ in sorted(v)]),
                [l for _, l, _ in sorted(v)],
                [s for _, _, s in sorted(v)])
            for c, v in by_chrom.items()}


def norm_chrom(c: str) -> str:
    """'Chr1', 'chr1', '1' -> '1'. DAP-seq and Ensembl disagree on chromosome naming."""
    c = str(c).strip()
    c = re.sub(r"^chr", "", c, flags=re.I)
    return {"c": "Pt", "m": "Mt"}.get(c.lower(), c)


def peaks_to_edges(tf_locus: str, peaks: list[tuple[str, int, float]], idx,
                   window: int, max_peaks: int) -> dict[str, float]:
    """Assign each peak to the nearest TSS within `window`; keep the best score per gene."""
    if max_peaks and len(peaks) > max_peaks:
        peaks = sorted(peaks, key=lambda p: -p[2])[:max_peaks]
    best: dict[str, float] = {}
    for chrom, pos, score in peaks:
        entry = idx.get(norm_chrom(chrom))
        if entry is None:
            continue
        tss_arr, loci, _ = entry
        i = int(np.searchsorted(tss_arr, pos))
        for j in (i - 1, i):
            if 0 <= j < len(loci) and abs(int(tss_arr[j]) - pos) <= window:
                g = loci[j]
                if score > best.get(g, -np.inf):
                    best[g] = score
    if not best:
        return {}
    # Rank-normalise within the TF: an edge score must mean "how confident is this
    # binding site relative to this TF's others", not "how deeply was this TF sequenced".
    genes = list(best)
    order = np.argsort(np.argsort([best[g] for g in genes]))
    denom = max(len(genes) - 1, 1)
    return {g: round(float(r) / denom, 4) for g, r in zip(genes, order)}


def resolve_tf(name: str, sym: dict[str, str]) -> tuple[str | None, str]:
    """The assayed TF's AGI locus, and how it was identified."""
    if (hit := AGI.search(name)):
        return hit.group(0).upper(), "agi_in_filename"
    if (m := DAP_TF_TOKEN.search(name)):
        token = m.group(1).upper()
        if (locus := sym.get(token)):
            return locus, "symbol_map"
        # Constructs are suffixed (VRN1a, ABR1_col); retry on the leading word.
        base = re.split(r"[^A-Z0-9]", token)[0]
        if base != token and (locus := sym.get(base)):
            return locus, "symbol_map_basename"
    return None, "unresolved"


def read_dapseq(window: int, max_peaks: int, idx, tfs: dict,
                sym: dict[str, str]) -> tuple[dict, dict]:
    path = download(DAPSEQ_TAR, OUT / "GSE60143_RAW.tar")
    edges: dict[str, dict[str, float]] = {}
    stats = {"members": 0, "narrowpeak": 0, "tf_in_planttfdb": 0,
             "resolved_by": {"agi_in_filename": 0, "symbol_map": 0,
                             "symbol_map_basename": 0, "unresolved": 0},
             "unresolved_examples": []}

    with tarfile.open(path) as tar:
        for member in tar:
            stats["members"] += 1
            name = member.name
            if "narrowPeak" not in name:
                continue
            stats["narrowpeak"] += 1
            tf, how = resolve_tf(name, sym)
            stats["resolved_by"][how] += 1
            if tf is None:
                # ChIP-seq controls, non-Arabidopsis orthologue assays (ZmARF5) and a
                # residue of post-2014 symbols. Listed in the QC rather than dropped
                # silently, so the coverage gap stays visible.
                if len(stats["unresolved_examples"]) < 25:
                    stats["unresolved_examples"].append(name)
                continue

            fh = tar.extractfile(member)
            if fh is None:
                continue
            raw = fh.read()
            try:
                text = gzip.decompress(raw).decode("utf-8", "replace")
            except (OSError, EOFError):
                text = raw.decode("utf-8", "replace")

            peaks = []
            for line in text.splitlines():
                f = line.split("\t")
                if len(f) < 7:
                    continue
                try:
                    start, end, signal = int(f[1]), int(f[2]), float(f[6])
                except ValueError:
                    continue
                peaks.append((f[0], (start + end) // 2, signal))
            if not peaks:
                continue

            new = peaks_to_edges(tf, peaks, idx, window, max_peaks)
            # A TF can be assayed in several constructs; keep the strongest evidence.
            cur = edges.setdefault(tf, {})
            for g, s in new.items():
                if s > cur.get(g, -1.0):
                    cur[g] = s

    stats["tf_in_planttfdb"] = sum(1 for t in edges if t in tfs)
    return edges, stats


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    # 1 kb, not 3 kb. DAP-seq is promiscuous: at +/-3 kb the median TF acquires ~5,200
    # targets, roughly a fifth of the genome, and a prior that broad makes DREM's TF
    # attribution at a bifurcation nearly uninformative. 1 kb keeps promoter-proximal
    # binding, which is what a regulatory prior should encode.
    ap.add_argument("--window", type=int, default=1000,
                    help="max bp from a peak summit to a TSS (default 1000)")
    ap.add_argument("--max-peaks-per-tf", type=int, default=5000,
                    help="keep only the strongest N peaks per assay (0 = all)")
    ap.add_argument("--min-score", type=float, default=0.0,
                    help="drop edges scoring below this after rank normalisation")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    tfs = fetch_tf_list()
    genes, gtf_path = fetch_gene_tss()
    sym, sym_provenance = build_symbol_map(gtf_path)
    idx = tss_index(genes)

    log("reading DAP-seq peaks (77 MB tar, first run downloads it) ...")
    edges, stats = read_dapseq(args.window, args.max_peaks_per_tf, idx, tfs, sym)

    path = OUT / "tf_gene_edges.tsv"
    n_edges = 0
    with open(path, "w") as fh:
        fh.write("TF\tGene\tInput\n")  # DREM's required header
        for tf in sorted(edges):
            for gene, score in sorted(edges[tf].items()):
                if score >= args.min_score:
                    fh.write(f"{tf}\t{gene}\t{score:g}\n")
                    n_edges += 1

    per_tf = {tf: len(v) for tf, v in edges.items()}
    qc = {
        "generated": dt.date.today().isoformat(),
        "sources": {"tf_list": PLANTTFDB_TF_LIST, "annotation": ENSEMBL_GTF,
                    "binding": DAPSEQ_TAR, "aliases": TAIR_ALIASES,
                    "tf_symbols": CONNECTF_TFS},
        "symbol_map": {"n_aliases": len(sym), "provenance": sym_provenance},
        "parameters": {"window_bp": args.window,
                       "max_peaks_per_tf": args.max_peaks_per_tf,
                       "min_score": args.min_score},
        "measured": {
            "planttfdb_tf_loci": len(tfs),
            "annotated_genes": len(genes),
            "tar_members": stats["members"],
            "narrowpeak_files": stats["narrowpeak"],
            "tf_identified_by": stats["resolved_by"],
            "unresolved_examples": stats["unresolved_examples"],
            "tfs_with_edges": len(edges),
            "tfs_also_in_planttfdb": stats["tf_in_planttfdb"],
            "total_edges": n_edges,
            "median_targets_per_tf": int(np.median(list(per_tf.values()))) if per_tf else 0,
            "max_targets_per_tf": max(per_tf.values()) if per_tf else 0,
            "genes_with_any_edge": len({g for v in edges.values() for g in v}),
        },
        # Surfaced every run so an absent master regulator can never be discovered late.
        "key_regulators": {
            name: ({"present": True, "n_targets": len(edges[locus])} if locus in edges
                   else {"present": False,
                         "note": "not assayed in DAP-seq; no edges available from this source"})
            for locus, name in KEY_REGULATORS.items()
        },
        "output": str(path.relative_to(DATA.parent)),
    }
    write_json(QC, qc)
    log(f"\nwrote {path}")

    # This script REPLACES tf_gene_edges.tsv, so a previous 04b merge is gone. Silently
    # losing SOG1 would silently remove the one regulator the cohort's conclusions depend
    # on, and no downstream step would notice.
    if (OUT / "sog1_edges.tsv").exists() and "AT1G25580" not in {t for t in edges}:
        log("\nNOTE: sog1_edges.tsv exists but this rewrite dropped SOG1 from the prior.\n"
            "      Re-run:  python3 scripts/04b_sog1_chip_prior.py")
    for k, v in qc["measured"].items():
        if k != "unresolved_examples":
            log(f"  {k}: {v}")
    return 0 if n_edges else 1


if __name__ == "__main__":
    raise SystemExit(main())
