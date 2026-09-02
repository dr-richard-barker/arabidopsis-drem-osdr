#!/usr/bin/env python3
"""Write MANIFEST.tsv: provenance, licence and SHA-256 for every input and output.

Large inputs are catalogued but not redistributed — the manifest records where each came
from and what its checksum was, so `run_all.sh` can reconstruct them and a reader can
verify they got the same bytes. `redistributed_here` says plainly which is which.

  MANIFEST.tsv          one row per file
  CHECKSUMS.sha256      sha256sum-compatible, for the files this repo does ship

  python3 scripts/14_manifest.py
"""
from __future__ import annotations

import datetime as dt
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import ROOT, log  # noqa: E402

MANIFEST = ROOT / "MANIFEST.tsv"
CHECKSUMS = ROOT / "CHECKSUMS.sha256"

# path glob -> (description, source, licence, redistributed_here)
SPEC = [
    ("data/counts/*.csv",
     "Unnormalized RSEM count matrix",
     "https://osdr.nasa.gov/bio/repo/data/studies/{acc}",
     "NASA Open Data (cite the original investigators)", "no"),
    ("data/runsheets/*.csv",
     "GeneLab runsheet: sample IDs and factor values",
     "https://osdr.nasa.gov/bio/repo/data/studies/{acc}",
     "NASA Open Data", "no"),
    ("data/osdr/study_catalog.json",
     "Measured per-study file counts, titles and publications",
     "derived from the OSDR file and biodata APIs", "CC-BY-4.0", "yes"),
    ("data/osdr/arabidopsis_factor_values.csv",
     "Sample-level factor values for every OSDR Arabidopsis study",
     "https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/",
     "NASA Open Data", "no"),
    ("data/metadata_master.csv",
     "Harmonised cohort metadata with a numeric time axis",
     "derived from the OSDR biodata API", "CC-BY-4.0", "yes"),
    ("data/tf_prior/GSE60143_RAW.tar",
     "Arabidopsis DAP-seq peaks (O'Malley et al. 2016)",
     "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60143/suppl/GSE60143_RAW.tar",
     "NCBI GEO terms", "no"),
    ("data/tf_prior/Arabidopsis_thaliana.TAIR10.54.gtf.gz",
     "Ensembl Plants release-54 TAIR10 annotation",
     "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-54/gtf/arabidopsis_thaliana/",
     "EMBL-EBI (Apache-2.0 style, no restriction)", "no"),
    ("data/tf_prior/Ath_TF_list.txt.gz",
     "PlantTFDB Arabidopsis TF list",
     "https://planttfdb.gao-lab.org/download/TF_list/Ath_TF_list.txt.gz",
     "PlantTFDB terms", "no"),
    ("data/tf_prior/tair_gene_aliases.txt",
     "TAIR gene alias table (symbol -> AGI locus)",
     "https://www.arabidopsis.org/", "TAIR terms", "no"),
    ("data/chip/*.bedGraph.gz",
     "SOG1-3xFLAG ChIP-seq coverage (Bourbousse et al. 2018)",
     "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112529",
     "NCBI GEO terms", "no"),
    ("data/arrayexpress/*.CEL",
     "Shenzhou-8 SIMBOX raw Affymetrix ATH1 arrays (Fengler et al. 2015, E-MTAB-2518)",
     "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2518",
     "EBI BioStudies terms", "no"),
    ("data/arrayexpress/*.sdrf.txt",
     "Shenzhou-8 sample and data relationship file; carries the growth-condition design",
     "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2518",
     "EBI BioStudies terms", "yes"),
    ("data/arrayexpress/*.annot.gz",
     "GPL198 (ATH1-121501) probeset annotation, used for probeset to AGI mapping",
     "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL198/annot/GPL198.annot.gz",
     "NCBI GEO terms", "no"),
    ("data/brapa/TAIR10.pep.all.fa.gz",
     "TAIR10 reference proteome, the Arabidopsis side of the orthology search",
     "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-54/fasta/"
     "arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz",
     "Ensembl Plants terms", "no"),
    ("data/brapa/brapa_arabidopsis_rbh.tsv",
     "B. rapa to Arabidopsis reciprocal-best-hit orthologs; replaces a BioMart route that "
     "returns no Arabidopsis homolog attributes for brapa_eg_gene",
     "derived by scripts/26_brapa_orthologs.py", "CC-BY-4.0", "yes"),
    ("data/tf_prior/sog1_edges.tsv",
     "SOG1 target edges called from GSE112529 coverage",
     "derived by scripts/04b_sog1_chip_prior.py", "CC-BY-4.0", "yes"),
    ("data/tf_prior/gene_tss.tsv",
     "Gene TSS coordinates from the Ensembl annotation",
     "derived by scripts/04_fetch_tf_prior.py", "CC-BY-4.0", "yes"),
    ("data/tf_prior/tf_list.tsv",
     "TF locus and family table", "derived from PlantTFDB", "CC-BY-4.0", "yes"),
    ("results/pseudotimeseries/*.tsv",
     "Pseudo-time-series expression and replicate matrices",
     "derived by scripts/03_build_pseudotimeseries.py", "CC-BY-4.0", "yes"),
    ("results/celltypes/*.tsv",
     "Cell-type fractions and latent trajectories",
     "derived by scripts/05_deconvolve_celltypes.py", "CC-BY-4.0", "yes"),
    ("results/qc/*.json",
     "Quality-control reports for every pipeline stage",
     "derived by the pipeline", "CC-BY-4.0", "yes"),
    ("results/drem/parsed/*.tsv",
     "Parsed DREM TF scores and gene-to-node assignments",
     "derived by scripts/09_parse_drem_model.py", "CC-BY-4.0", "yes"),
    ("results/drem/runs/*/model.txt",
     "Fitted DREM models",
     "produced by DREM 2.0.7 (GPL-3.0, Ernst lab)", "CC-BY-4.0", "yes"),
    ("results/comparison/*.tsv",
     "Prior ablation and genotype contrast tables",
     "derived by scripts/10_compare_models.py", "CC-BY-4.0", "yes"),
    ("vendor/drem.jar",
     "DREM 2.0.7 engine",
     "https://github.com/jernst98/STEM_DREM",
     "GPL-3.0 (downloaded at run time, not redistributed)", "no"),
    ("manuscript/manuscript.pdf",
     "Manuscript, PDF", "built from manuscript/latex", "CC-BY-4.0", "yes"),
    ("manuscript/manuscript.docx",
     "Manuscript, Word", "built from manuscript/latex", "CC-BY-4.0", "yes"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    today = dt.date.today().isoformat()
    rows, shipped = [], []

    for pattern, desc, source, licence, redistributed in SPEC:
        for path in sorted(ROOT.glob(pattern)):
            if not path.is_file():
                continue
            rel = path.relative_to(ROOT)
            src = source
            if "{acc}" in src:
                acc = rel.name.split("_")[0]
                src = src.format(acc=acc)
            digest = sha256(path)
            rows.append([str(rel), desc, src, licence, redistributed,
                         str(path.stat().st_size), digest, today])
            if redistributed == "yes":
                shipped.append((digest, rel))

    MANIFEST.write_text(
        "path\tdescription\tsource\tlicence\tredistributed_here\tbytes\tsha256\tretrieved\n"
        + "".join("\t".join(r) + "\n" for r in rows))
    CHECKSUMS.write_text("".join(f"{d}  {p}\n" for d, p in shipped))

    n_ship = sum(1 for r in rows if r[4] == "yes")
    log(f"wrote {MANIFEST}  ({len(rows)} files: {n_ship} redistributed, "
        f"{len(rows) - n_ship} catalogued but re-fetchable)")
    log(f"wrote {CHECKSUMS}  ({len(shipped)} entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
