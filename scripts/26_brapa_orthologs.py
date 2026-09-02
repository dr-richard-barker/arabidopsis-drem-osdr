#!/usr/bin/env python3
"""Build a Brassica rapa -> Arabidopsis ortholog map by reciprocal best hit.

Why this exists rather than a database lookup: Ensembl Plants has no Arabidopsis-homolog
attributes for `brapa_eg_gene` in BioMart, and its plants FTP carries no homology TSV for
the species pair. The B_rappa_LLGCSS repository already contains the evidence of that
dead end -- `annotation/brapa_to_arabidopsis_orthologs.tsv` is an empty file holding the
error message

    Query ERROR: caught BioMart::Exception::Usage:
    Attribute athaliana_eg_homolog_ensembl_gene NOT FOUND

and the same query fails the same way today.

The BioMart export that repository does carry (`NewTest/B.rapa_to_Ara_mart_export.txt`)
maps only 3,506 B. rapa genes, which reaches 21% of this pipeline's transcription factors
and neither SOG1 nor MYB3R1 -- so both arms of the radiation decoder are unmappable
through it, and no cross-species scoring is possible.

Reciprocal best hit against the TAIR10 proteome fixes that from sequence. RBH rather than
one-way best hit because B. rapa is a mesohexaploid: a one-way search returns three or
more B. rapa paralogues for a single Arabidopsis gene, and taking the top hit in one
direction silently picks whichever paralogue happened to score highest.

  data/brapa/brapa_arabidopsis_rbh.tsv   brapa_gene, agi, identity, evalue, bitscore
  results/qc/brapa_ortholog_qc.json      coverage, and how far it beats the existing map

  python3 scripts/26_brapa_orthologs.py
  python3 scripts/26_brapa_orthologs.py --brapa-cds /path/to/cds.fa.gz
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import DATA, RESULTS, download, log, write_json  # noqa: E402

OUT = DATA / "brapa"
TAIR_PROT = ("https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-54/fasta/"
             "arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz")
DEFAULT_CDS = (Path.home() / "Documents" / "B_rappa_LLGCSS" /
               "BrapaFPsc_277_v1.3.cds_primaryTranscriptOnly.fa.gz")
EXISTING_MAP = (Path.home() / "Documents" / "B_rappa_LLGCSS" / "NewTest" /
                "B.rapa_to_Ara_mart_export.txt")

AGI = re.compile(r"AT[1-5CM]G\d{5}", re.I)


def need(tool: str) -> str:
    p = shutil.which(tool) or shutil.which(f"/opt/homebrew/bin/{tool}")
    if not p:
        raise SystemExit(f"{tool} not found. Install with:  brew install {tool}")
    return p


CODON = {}
for i, aa in enumerate(
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"):
    CODON["AAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATGATT"
          "CAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTT"
          "GAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTT"
          "TAATACTAGTATTCATCCTCGTCTTGATGCTGGTGTTTATTCTTGTTT"[i * 3:i * 3 + 3]] = aa


def translate(seq: str) -> str:
    """CDS -> protein, frame 0, truncated at the first stop. DIAMOND cannot build a
    database from nucleotides (makedb aborts), and it has no tblastn mode, so the
    reverse half of a reciprocal search is impossible unless both sides are protein."""
    seq = seq.upper().replace("U", "T")
    aas = [CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3)]
    if "*" in aas:
        aas = aas[:aas.index("*")]
    return "".join(aas)


def clean_fasta(src: Path, dest: Path, id_from_header, as_protein: bool = False) -> int:
    """Rewrite a FASTA with short, stable ids, keeping the longest sequence per id.

    Deduplicating matters: TAIR10's proteome carries every splice isoform, so 48,321
    sequences collapse to 27,628 loci. Without this, "best hit" is a race between
    isoforms of the same gene rather than a choice between genes.
    """
    if dest.exists() and dest.stat().st_size > 1000:
        return sum(1 for l in dest.open() if l.startswith(">"))
    opener = gzip.open if src.suffix == ".gz" else open
    seqs: dict[str, str] = {}
    ident, buf = None, []

    def flush() -> None:
        if ident is None:
            return
        s = "".join(buf)
        if as_protein:
            s = translate(s)
        if len(s) > len(seqs.get(ident, "")):
            seqs[ident] = s

    with opener(src, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                ident, buf = id_from_header(line[1:].strip()), []
            else:
                buf.append(line.strip())
    flush()
    with dest.open("w") as out:
        for k, v in seqs.items():
            out.write(f">{k}\n")
            for i in range(0, len(v), 60):
                out.write(v[i:i + 60] + "\n")
    return len(seqs)


def run_diamond(dmd: str, query: Path, db_fa: Path, out: Path, threads: int) -> Path:
    """One DIAMOND search, cached on the output path."""
    if out.exists() and out.stat().st_size > 0:
        return out
    db = db_fa.with_suffix(".dmnd")
    if not db.exists():
        subprocess.run([dmd, "makedb", "--in", str(db_fa), "-d", str(db.with_suffix("")),
                        "--quiet"], check=True)
    subprocess.run(
        [dmd, "blastp", "-q", str(query), "-d", str(db), "-o", str(out),
         "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore",
         "--max-target-seqs", "1", "--evalue", "1e-5", "--sensitive",
         "--threads", str(threads), "--quiet"], check=True)
    return out


def best_hits(path: Path) -> dict[str, tuple[str, float, float, float]]:
    """query -> (subject, identity, evalue, bitscore), keeping the highest bitscore."""
    best: dict[str, tuple[str, float, float, float]] = {}
    with path.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            q, s, pid, ev, bs = f[0], f[1], float(f[2]), float(f[3]), float(f[4])
            if q not in best or bs > best[q][3]:
                best[q] = (s, pid, ev, bs)
    return best


def existing_map_coverage() -> int:
    if not EXISTING_MAP.exists():
        return 0
    seen = set()
    with EXISTING_MAP.open(errors="replace") as fh:
        import csv as _csv
        for r in _csv.DictReader(fh):
            if AGI.search(r.get("Gene description") or ""):
                seen.add(r["Gene stable ID"].strip().upper())
    return len(seen)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--brapa-cds", type=Path, default=DEFAULT_CDS)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-identity", type=float, default=50.0)
    args = ap.parse_args()

    dmd = need("diamond")
    if not args.brapa_cds.exists():
        raise SystemExit(f"B. rapa CDS not found: {args.brapa_cds}")
    OUT.mkdir(parents=True, exist_ok=True)

    tair_gz = download(TAIR_PROT, OUT / "TAIR10.pep.all.fa.gz")
    # Ensembl protein headers carry the locus in `gene:AT1G01010`; fall back to any AGI.
    def tair_id(h: str) -> str:
        m = re.search(r"gene:(AT[1-5CM]G\d{5})", h, re.I) or AGI.search(h)
        return m.group(1).upper() if m and m.lastindex else (m.group(0).upper() if m else h.split()[0])
    n_tair = clean_fasta(tair_gz, OUT / "tair10_pep.fa", tair_id)

    # Phytozome headers are `Brara.K01534.1 pacid=... locus=Brara.K01534 ...`. Splitting
    # the first token on "." yields "Brara" for every gene in the file -- one id for
    # 40,492 sequences, which silently collapses the whole search to a single hit.
    def brapa_id(h: str) -> str:
        m = re.search(r"locus=(\S+)", h)
        return (m.group(1) if m else re.sub(r"\.\d+$", "", h.split()[0])).upper()
    n_brapa = clean_fasta(args.brapa_cds, OUT / "brapa_pep.fa", brapa_id, as_protein=True)
    log(f"sequences: {n_brapa} B. rapa loci (CDS translated), {n_tair} TAIR10 loci")

    log("  forward search (B. rapa -> Arabidopsis) ...")
    fwd = best_hits(run_diamond(dmd, OUT / "brapa_pep.fa", OUT / "tair10_pep.fa",
                                OUT / "fwd.tsv", args.threads))
    log(f"    {len(fwd)} B. rapa genes with a hit")

    log("  reverse search (Arabidopsis -> B. rapa) ...")
    rev = best_hits(run_diamond(dmd, OUT / "tair10_pep.fa", OUT / "brapa_pep.fa",
                                OUT / "rev.tsv", args.threads))
    log(f"    {len(rev)} Arabidopsis genes with a hit")

    # Reciprocal best hit: A's best hit is B, and B's best hit is A.
    rbh = []
    for bra, (agi, pid, ev, bs) in fwd.items():
        back = rev.get(agi)
        if back and back[0] == bra and pid >= args.min_identity:
            rbh.append((bra, agi, pid, ev, bs))
    rbh.sort()

    path = OUT / "brapa_arabidopsis_rbh.tsv"
    with path.open("w") as fh:
        fh.write("brapa_gene\tagi\tidentity\tevalue\tbitscore\n")
        for bra, agi, pid, ev, bs in rbh:
            fh.write(f"{bra}\t{agi}\t{pid:g}\t{ev:g}\t{bs:g}\n")

    agis = {r[1] for r in rbh}
    # What this unlocks, measured against the sets that actually matter.
    tfs = set()
    act = RESULTS / "tf_activity" / "tf_activity.tsv"
    if act.exists():
        tfs = set(act.open().readline().rstrip("\n").split("\t")[1:])
    sig_path = RESULTS / "decoder" / "radiation_signature.json"
    sig = json.loads(sig_path.read_text())["sets"] if sig_path.exists() else {}

    qc = {
        "generated": dt.date.today().isoformat(),
        "method": "reciprocal best hit, DIAMOND blastx/blastp, identity >= "
                  f"{args.min_identity}%, e-value <= 1e-5",
        "why_not_a_database": "Ensembl Plants BioMart has no athaliana homolog attributes "
                              "for brapa_eg_gene and the plants FTP has no homology TSV "
                              "for this pair; the B_rappa_LLGCSS repo contains an empty "
                              "brapa_to_arabidopsis_orthologs.tsv holding that same error",
        "n_brapa_loci": n_brapa, "n_tair_proteins": n_tair,
        "n_forward_hits": len(fwd), "n_reverse_hits": len(rev),
        "n_rbh_pairs": len(rbh), "n_arabidopsis_loci": len(agis),
        "existing_biomart_map_genes": existing_map_coverage(),
        "coverage": {},
    }
    if tfs:
        qc["coverage"]["drem_tfs"] = [len(tfs & agis), len(tfs)]
        for key, name in (("AT1G25580", "SOG1"), ("AT4G32730", "MYB3R1")):
            qc["coverage"][name] = key in agis
    for s in ("sog1_dependent", "myb3r_repressed"):
        if s in sig:
            qc["coverage"][s] = [len(set(sig[s]) & agis), len(sig[s])]

    write_json(RESULTS / "qc" / "brapa_ortholog_qc.json", qc)
    log(f"\n  RBH pairs: {len(rbh)}  ({len(agis)} Arabidopsis loci)")
    log(f"  existing BioMart map covered {qc['existing_biomart_map_genes']} B. rapa genes")
    for k, v in qc["coverage"].items():
        log(f"    {k}: {v}")
    return 0 if rbh else 1


if __name__ == "__main__":
    raise SystemExit(main())
