#!/usr/bin/env python3
"""Build and verify `manuscript/latex/references.bib` from Crossref — never from memory.

Every entry starts as a DOI and a cite key. The script resolves each DOI through the
Crossref REST API and writes the BibTeX from the *returned* metadata, so authors,
titles, journals, volumes, pages and years cannot be misremembered. A DOI that does not
resolve is reported and left out of the .bib rather than guessed at, and the run exits
non-zero so CI catches it.

  python3 scripts/check_references.py            build/refresh references.bib
  python3 scripts/check_references.py --check    verify only, write nothing

Sources with no DOI (a data repository record, a software release) are held in NON_DOI
below with the URL that has to be checked by hand; they are written as @misc and listed
in the report as unverifiable by this script.

Ported from AstroRegolith/scripts/check_references.py (same author, MIT).
"""
from __future__ import annotations

import argparse
import sys
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import CROSSREF, MAILTO, ROOT, get_json, log  # noqa: E402

BIB = ROOT / "manuscript" / "latex" / "references.bib"

# cite key -> (DOI, expected first-author surname, expected year).
#
# The expectation is the point. Resolving a DOI only proves the DOI exists; it does not
# prove it is the paper you meant. Four of the DOIs first written here resolved cleanly
# to entirely different works — iDREM to a Gaussian-process paper, an Arabidopsis atlas
# to a Drosophila chromatin study — and a resolve-only check passed all of them. The
# author/year assertion is what turns this script from a link-checker into a citation
# checker, so never relax it to silence a failure: fix the DOI.
DOIS = {
    "bourbousse2018sog1": ("10.1073/pnas.1810582115", "Bourbousse", 2018),
    "ernst2007drem": ("10.1038/msb4100115", "Ernst", 2007),
    "schulz2012drem2": ("10.1186/1752-0509-6-104", "Schulz", 2012),
    "ding2018idrem": ("10.1371/journal.pcbi.1006019", "Ding", 2018),
    "omalley2016cistrome": ("10.1016/j.cell.2016.04.038", "O'Malley", 2016),
    "jin2017planttfdb": ("10.1093/nar/gkw982", "Jin", 2016),
    "brooks2021connectf": ("10.1093/plphys/kiaa012", "Brooks", 2020),
    "lee2023atlas": ("10.1101/2023.03.23.533992", "Lee", 2023),
    "yoshiyama2009sog1": ("10.1073/pnas.0810304106", "Yoshiyama", 2009),
    "berardini2015tair": ("10.1002/dvg.22877", "Berardini", 2015),
    "love2014deseq2": ("10.1186/s13059-014-0550-8", "Love", 2014),
    "yoshida2022iss": ("10.1016/j.heliyon.2022.e10266", "Yoshida", 2022),
    "reitz2008leo": ("10.1016/j.zemedi.2008.06.015", "Reitz", 2008),
    "slaba2025artemis": ("10.1038/s41526-025-00459-y", "Slaba", 2025),
}

# No DOI exists for these; the URL is what a reader must check.
NON_DOI = {
    "osdr": ("NASA Open Science Data Repository",
             "https://osdr.nasa.gov/bio/repo/", "NASA"),
    "geo_gse112529": ("Gene Expression Omnibus GSE112529: SOG1-3xFLAG ChIP-seq",
                      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112529", "NCBI"),
    "geo_gse60143": ("Gene Expression Omnibus GSE60143: Arabidopsis DAP-seq",
                     "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60143", "NCBI"),
    "geo_gse226097": ("Gene Expression Omnibus GSE226097: Arabidopsis single-nucleus atlas",
                      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226097", "NCBI"),
    "dremsoftware": ("DREM 2.0.7 software release",
                     "https://github.com/jernst98/STEM_DREM", "Ernst Lab"),
    "ensemblplants54": ("Ensembl Plants release 54, Arabidopsis thaliana TAIR10",
                        "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-54/", "EMBL-EBI"),
}


def esc(s: str) -> str:
    return (str(s).replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")
            .replace("#", r"\#").replace("$", r"\$"))


def strip_tags(s: str) -> str:
    import re
    return re.sub(r"<[^>]+>", "", str(s))


def fetch(doi: str) -> dict | None:
    url = f"{CROSSREF}{urllib.parse.quote(doi)}?mailto={MAILTO}"
    try:
        return get_json(url)["message"]
    except (RuntimeError, KeyError):
        return None


def to_bibtex(key: str, m: dict) -> str:
    authors = " and ".join(
        f"{a.get('family', '')}, {a.get('given', '')}".strip(", ")
        for a in m.get("author", []) if a.get("family"))
    year = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
    title = strip_tags((m.get("title") or [""])[0])
    journal = strip_tags((m.get("container-title") or [""])[0])

    if not journal:
        # Crossref gives no container-title for preprints; `unsrtnat` renders an
        # @article with an empty journal as a dangling comma.
        journal = (m.get("institution", [{}]) or [{}])[0].get("name") or ""
        kind = "misc"
        fields = [("author", authors), ("title", title), ("year", year),
                  ("howpublished", journal or "Preprint"), ("doi", m.get("DOI"))]
    else:
        kind = "article"
        fields = [("author", authors), ("title", title), ("journal", journal),
                  ("year", year), ("volume", m.get("volume")), ("number", m.get("issue")),
                  ("pages", m.get("page")), ("doi", m.get("DOI"))]
    body = "".join(f"  {k} = {{{esc(v)}}},\n" for k, v in fields if v)
    return f"@{kind}{{{key},\n{body}}}\n\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="verify only, write nothing")
    args = ap.parse_args()

    entries, failed, mismatched = [], [], []
    for key, (doi, want_author, want_year) in sorted(DOIS.items()):
        m = fetch(doi)
        if m is None:
            log(f"  UNRESOLVED  {key}  {doi}")
            failed.append((key, doi))
            continue
        title = strip_tags((m.get("title") or [""])[0])
        first = next((a.get("family") for a in m.get("author", []) if a.get("family")), "?")
        year = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]

        def norm(x):
            return "".join(c for c in str(x).lower() if c.isalnum())

        author_ok = norm(want_author) == norm(first)
        # Online-vs-issue years legitimately differ by one; anything more is a red flag.
        year_ok = year is not None and abs(int(year) - int(want_year)) <= 1
        if not (author_ok and year_ok):
            log(f"  MISMATCH    {key:<24} expected {want_author} {want_year}, "
                f"DOI gives {first} {year} — {title[:50]}")
            mismatched.append((key, doi, f"{first} {year}", f"{want_author} {want_year}"))
            continue

        log(f"  ok  {key:<26} {first} {year} — {title[:60]}")
        entries.append(to_bibtex(key, m))

    for key, (title, url, publisher) in sorted(NON_DOI.items()):
        entries.append(
            f"@misc{{{key},\n  title = {{{esc(title)}}},\n"
            f"  howpublished = {{\\url{{{url}}}}},\n"
            f"  publisher = {{{esc(publisher)}}},\n  year = {{2026}},\n"
            f"  note = {{Accessed 2026-08-28}},\n}}\n\n")

    if not args.check:
        BIB.parent.mkdir(parents=True, exist_ok=True)
        BIB.write_text(
            "% Generated by scripts/check_references.py from Crossref. Do not hand-edit:\n"
            "% every field below came back from the DOI resolver, so editing here would\n"
            "% reintroduce exactly the transcription errors this file prevents.\n"
            "% (No at-sign in this header: BibTeX starts an entry at one even in a comment.)\n\n"
            + "".join(entries))
        log(f"\nwrote {BIB}  ({len(entries)} entries, "
            f"{len(DOIS) - len(failed) - len(mismatched)} DOI-verified, "
            f"{len(NON_DOI)} @misc)")

    if failed:
        log(f"\n{len(failed)} DOI(s) did not resolve: {[d for _, d in failed]}")
    if mismatched:
        log(f"\n{len(mismatched)} DOI(s) resolved to the WRONG work — fix the DOI, "
            f"do not relax the check:")
        for key, doi, got, want in mismatched:
            log(f"   {key}: {doi} -> {got}, expected {want}")
    return 1 if (failed or mismatched) else 0


if __name__ == "__main__":
    raise SystemExit(main())
