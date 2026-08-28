#!/usr/bin/env python3
r"""Fail on undefined citations, undefined cross-references, or missing figure files.

Three failure modes that all survive proofreading because LaTeX only warns:

  * a \cite key references.bib does not define renders as a silent "[?]";
  * a \ref to a label nobody declared renders as "??";
  * an \includegraphics pointing at a file that is not there stops the PDF build but
    is easy to miss in a long log.

Unused bibliography keys are reported but not fatal — a .bib may legitimately carry a
work the current draft does not cite.

  python3 scripts/check_citations.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import ROOT, log  # noqa: E402

LATEX = ROOT / "manuscript" / "latex"
BIB = LATEX / "references.bib"

CITE = re.compile(r"\\cite[tp]?\*?(?:\[[^\]]*\])*\{([^}]*)\}")
ENTRY = re.compile(r"^@\w+\{([^,]+),", re.M)
REF = re.compile(r"\\(?:page)?ref\{([^}]*)\}")
LABEL = re.compile(r"\\label\{([^}]*)\}")
GRAPHIC = re.compile(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}")


def main() -> int:
    if not BIB.exists():
        log("references.bib missing — run scripts/check_references.py")
        return 1
    defined = set(ENTRY.findall(BIB.read_text()))

    used: set[str] = set()
    for tex in sorted(LATEX.rglob("*.tex")):
        for group in CITE.findall(tex.read_text()):
            used.update(k.strip() for k in group.split(",") if k.strip())

    undefined = sorted(used - defined)
    unused = sorted(defined - used)
    log(f"{len(used)} keys cited, {len(defined)} defined")
    if unused:
        log(f"  unused (not fatal): {', '.join(unused)}")

    # Cross-references
    refs: set[str] = set()
    labels: set[str] = set()
    graphics: list[tuple[str, Path]] = []
    for tex in sorted(LATEX.rglob("*.tex")):
        text = tex.read_text()
        for group in REF.findall(text):
            refs.update(k.strip() for k in group.split(",") if k.strip())
        labels.update(LABEL.findall(text))
        for g in GRAPHIC.findall(text):
            graphics.append((g, tex))

    dangling = sorted(refs - labels)
    log(f"{len(refs)} cross-references, {len(labels)} labels")

    # Figure files. Two subtleties:
    #  * LaTeX resolves \includegraphics against the MAIN document's directory, not the
    #    directory of the \input-ed file that contains the directive, so a path written
    #    in sections/results.tex is relative to latex/, not to sections/.
    #  * An extensionless path is matched against a set of known extensions.
    missing_figs = []
    for g, tex in graphics:
        candidates = [LATEX / g, tex.parent / g]
        found = False
        for base in candidates:
            if base.suffix:
                found = found or base.exists()
            else:
                found = found or any(base.with_suffix(e).exists()
                                     for e in (".pdf", ".png", ".jpg", ".eps"))
        if not found:
            missing_figs.append(g)
    log(f"{len(graphics)} figures included, {len(missing_figs)} missing")

    ok = True
    if undefined:
        log(f"  UNDEFINED CITATIONS: {', '.join(undefined)}")
        ok = False
    if dangling:
        log(f"  UNDEFINED REFERENCES: {', '.join(dangling)}")
        ok = False
    if missing_figs:
        log(f"  MISSING FIGURE FILES: {', '.join(missing_figs)}")
        ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
