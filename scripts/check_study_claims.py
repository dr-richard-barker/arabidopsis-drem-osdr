#!/usr/bin/env python3
r"""Fail if the manuscript describes a study in a way the OSDR record does not support.

`check_references.py` guards citations and `13_manuscript_numbers.py` guards numbers.
Nothing guarded the third class of factual claim: prose statements about what a study
actually did. Two errors got into the manuscript through that gap and neither was
catchable by any existing check.

  * OSD-658 was described as irradiating "dry dormant seed", and was EXCLUDED from an
    entire analysis on that basis. Its protocol says "The Arabidopsis seedlings were
    sequentially exposed..." and "RNA extracts from whole seedlings". The dry-seed
    sentence belongs to a separate phenotyping sub-experiment in the same deposit; one
    sentence was read and attributed to the wrong assay.
  * OSD-320 was described as irradiating "6-day-old seedlings". Its protocol says "Eight
    days prior to irradiation... seeds were... plated". The figure was carried across
    from OSD-508/510, where 6 days is correct.

The check is TWO-SIDED, and it has to be. A positive check alone --- "does the protocol
support this claim?" --- does not catch the error that motivated the script. The sentence
"Dry seeds were used for GCR experiments" really is in OSD-658's protocol; it just
describes a different sub-experiment. A regex asking whether that string exists returns
true, and the false description passes. Verified by deliberately reverting the pattern:
the positive check still returned 20/20.

So there are two tables. POSITIVE claims must be supported by the protocol text, which
catches invented facts. FORBIDDEN patterns must be absent from the manuscript and README,
which catches descriptions known to be wrong --- seeded from errors actually made, so a
correction cannot silently regress.

**Patterns must be tight.** An earlier audit of exactly these studies used `seedling` as
the pattern for "6-day-old seedlings"; it matched an 8-day-old study and the error
survived its own audit. A pattern that cannot fail is not a check. Every pattern here
therefore encodes the specific quantity claimed, not the topic of the claim, and
`--self-test` proves each one rejects a plausible wrong value.

  python3 scripts/check_study_claims.py
  python3 scripts/check_study_claims.py --self-test
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import log, osdr_meta  # noqa: E402

# accession -> list of (claim as written in the manuscript, tight regex, a wrong value
# the pattern must reject). The third element is what makes the check falsifiable.
CLAIMS: dict[str, list[tuple[str, str, str]]] = {
    "OSD-658": [
        ("seedlings, not seeds, were irradiated for the RNA-seq arm",
         r"Arabidopsis seedlings were sequentially exposed",
         "Arabidopsis dry seeds were sequentially exposed"),
        ("RNA came from whole seedlings",
         r"RNA extracts from whole seedlings",
         "RNA extracts from dry seeds"),
        ("harvested 3 h after irradiation",
         r"Three hours after irradiation",
         "Two days after irradiation"),
        ("40 and 80 cGy simulated GCR",
         r"40 cGy[^.]{0,40}80 cGy[^.]{0,40}GCR irradiation",
         "4 cGy and 8 cGy GCR irradiation"),
        ("sequential ion beams at NSRL",
         r"NSRL GCR and SPE simulations",
         "gamma irradiation at NSRL"),
    ],
    "OSD-320": [
        ("seedlings plated eight days before irradiation",
         r"Eight days prior to irradiation",
         "Six days prior to irradiation"),
        ("Ws and atm-1 background",
         r"wild-?type Ws and atm-1",
         "wild-type Col-0 and atm-1"),
        ("gamma at 100 Gy",
         r"gamma radiation \(100 Gy", "gamma radiation (10 Gy"),
        ("Fe-56 at 30 Gy, 1 GeV/n",
         r"1 GeV/n 56Fe particles \(30 Gy", "1 GeV/n 56Fe particles (3 Gy"),
        ("both beams at 7 Gy/min", r"7 Gy/min", "70 Gy/min"),
    ],
    "OSD-508": [
        ("6-day-old seedlings", r"6-day-old seedlings", "8-day-old seedlings"),
        ("100 Gy Co-60", r"dose of 100 Gy using a Co60", "dose of 10 Gy using a Co60"),
        ("10 Gy/min", r"10 Gy/min", "1 Gy/min"),
    ],
    "OSD-510": [
        ("6-day-old seedlings", r"6-day-old seedlings", "8-day-old seedlings"),
        ("100 Gy Co-60", r"dose of 100 Gy using a Co60", "dose of 10 Gy using a Co60"),
    ],
    "OSD-782": [
        ("4-week-old plants", r"At 4 weeks of growth", "At 4 days of growth"),
        ("137Cs irradiator", r"137Cs irradiator", "Co60 irradiator"),
        ("1.4 cGy per second", r"1\.4 cGy per second", "14 cGy per second"),
        ("7 s = 10 cGy, 71 s = 100 cGy",
         r"7 sec \(10 cGy equivalent\)[^.]{0,40}71 sec \(100 cGy equivalent\)",
         "7 sec (100 cGy equivalent)"),
    ],
    "OSD-496": [
        ("SOG1-3xFLAG ChIP", r"SOG1[- ]3.?FLAG|pSOG1::SOG1", "SOG1-GFP"),
    ],
}


# Descriptions that must NOT appear in the manuscript or README. Each entry records an
# error that was actually made and corrected, so the fix cannot regress unnoticed. The
# positive table cannot express these: the wrong text is often present in the deposit,
# attached to a different sub-experiment.
FORBIDDEN: list[tuple[str, str]] = [
    # Unconditional. No study in this analysis irradiated dry seed for its RNA-seq arm,
    # so the phrase has no legitimate use in the manuscript at all. An earlier version
    # required proximity to "OSD-658" and missed the offending sentence, which said "the
    # GCR arm" instead; it also used [^.] to bound the gap, which breaks on the period
    # inside $^{137}$Cs. Both mistakes made the guard inert -- a proximity window is a
    # liability when the thing you are guarding against can be phrased without the anchor.
    (r"dry\s+(dormant\s+)?seed",
     "No study's RNA-seq arm used dry seed. OSD-658 irradiated flask-grown seedlings; "
     "the dry-seed sentence in that deposit belongs to its phenotyping sub-experiment."),
    # Proximity is unavoidable here because "6-day" is correct for OSD-508/510. Dots are
    # allowed inside the window for the reason above.
    (r"OSD-320.{0,200}?6-day|6-day.{0,200}?OSD-320",
     "OSD-320 plated its seedlings eight days before irradiation, not six. The 6-day "
     "figure belongs to OSD-508/510."),
]

MANUSCRIPT_FILES = ["manuscript/latex/main.tex", "README.md"]

# The scan covers scripts/ too. Both errors it guards against were written into script
# docstrings *before* they reached the manuscript -- 24_radiation_quality.py carried the
# OSD-320 "6-day-old" figure for several revisions after the manuscript was corrected --
# so a manuscript-only scan lets the wrong description survive where the next author will
# read it. This file is the sole exemption: it must hold the patterns to test them.
SCAN_DIRS = ["manuscript/latex/sections", "scripts"]
SCAN_EXEMPT = {"scripts/check_study_claims.py"}


def scan_forbidden() -> list[tuple[str, str, str]]:
    """Manuscript-side check: descriptions known to be wrong must not reappear."""
    from lib_sources import ROOT
    texts = {}
    for rel in MANUSCRIPT_FILES:
        f = ROOT / rel
        if f.exists():
            texts[rel] = f.read_text(errors="replace")
    for d in SCAN_DIRS:
        base = ROOT / d
        if not base.is_dir():
            continue
        for f in sorted(base.iterdir()):
            rel = str(f.relative_to(ROOT))
            if f.suffix not in (".tex", ".py") or rel in SCAN_EXEMPT:
                continue
            texts[rel] = f.read_text(errors="replace")

    hits = []
    for rel, raw in texts.items():
        # Flattened: a line break must not let a forbidden phrase hide, which is how one
        # instance survived an earlier manual sweep.
        flat = re.sub(r"\s+", " ", raw)
        for pat, why in FORBIDDEN:
            for m in re.finditer(pat, flat, re.I):
                hits.append((rel, flat[max(0, m.start() - 70):m.start() + 90], why))
    return hits


def protocol(acc: str) -> str:
    d = osdr_meta(acc).get("study protocol description") or ""
    if isinstance(d, list):
        d = " ".join(map(str, d))
    return re.sub(r"\s+", " ", str(d))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--self-test", action="store_true",
                    help="prove each pattern rejects a plausible wrong value")
    args = ap.parse_args()

    if args.self_test:
        bad = []
        for acc, claims in CLAIMS.items():
            for label, pat, wrong in claims:
                # The pattern must NOT match the counterfactual text. A pattern that
                # matches both the right and the wrong value cannot catch an error.
                if re.search(pat, wrong, re.I):
                    bad.append(f"{acc}: {label!r} also matches the wrong value {wrong!r}")
        for b in bad:
            log(f"  TOO LOOSE  {b}")
        log(f"\nself-test: {sum(len(v) for v in CLAIMS.values()) - len(bad)}"
            f"/{sum(len(v) for v in CLAIMS.values())} patterns reject their counterfactual")
        return 1 if bad else 0

    failures = []
    for acc in sorted(CLAIMS, key=lambda a: int(a.split("-")[1])):
        try:
            text = protocol(acc)
        except RuntimeError as e:
            log(f"  {acc}: OSDR unreachable ({e}) — cannot verify")
            failures.append((acc, "protocol unavailable"))
            continue
        if not text:
            log(f"  {acc}: no protocol description in the OSDR record")
            failures.append((acc, "no protocol text"))
            continue
        for label, pat, _wrong in CLAIMS[acc]:
            if re.search(pat, text, re.I):
                log(f"  ok    {acc}  {label}")
            else:
                log(f"  FAIL  {acc}  {label}")
                failures.append((acc, label))

    forbidden = scan_forbidden()
    for rel, seg, why in forbidden:
        log(f"  FORBIDDEN  {rel}: ...{seg.strip()[:110]}...")
        log(f"             {why}")
    if not forbidden:
        log(f"  ok    no corrected description has regressed "
            f"({len(FORBIDDEN)} patterns checked)")

    total = sum(len(v) for v in CLAIMS.values())
    log(f"\n{total - len(failures)}/{total} manuscript claims verified against OSDR; "
        f"{len(forbidden)} forbidden descriptions found")
    failures += [(rel, why) for rel, _seg, why in forbidden]
    if failures:
        log("\nUnverified claims — correct the manuscript, or the pattern if the claim "
            "is right and the pattern is wrong:")
        for acc, label in failures:
            log(f"  {acc}: {label}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
