#!/usr/bin/env python3
"""Check the toolchain this pipeline needs and print the exact command for anything missing.

Run this before anything else. It exits non-zero if a *required* tool is absent so
`run_all.sh` and CI stop early rather than failing halfway through a download.

  python3 scripts/00_env_check.py              # required tools only
  python3 scripts/00_env_check.py --manuscript # also the PDF/DOCX build chain
  python3 scripts/00_env_check.py --json       # machine-readable report
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import log  # noqa: E402

# Homebrew keeps openjdk off the default PATH ("keg-only") so it cannot shadow the
# system Java. Look there explicitly rather than telling the user to edit their shell.
BREW_JDK_BINS = [
    "/opt/homebrew/opt/openjdk@21/bin",
    "/opt/homebrew/opt/openjdk/bin",
    "/usr/local/opt/openjdk@21/bin",
    "/usr/local/opt/openjdk/bin",
]


def _runs(path: str, args: list[str]) -> str | None:
    """Return the tool's first version line, or None if invoking it does not succeed.

    Existence is not enough. macOS ships a /usr/bin/java shim that is present on every
    Mac and exits non-zero with "Unable to locate a Java Runtime" when no JDK is
    installed — trusting shutil.which alone reports Java as available on a machine that
    cannot run a single jar.
    """
    try:
        out = subprocess.run([path, *args], capture_output=True, text=True, timeout=60)
    except Exception:  # noqa: BLE001 - an unrunnable tool is a missing tool here
        return None
    if out.returncode != 0:
        return None
    lines = (out.stdout or out.stderr).strip().splitlines()
    return lines[0] if lines else "(no version output)"


def probe(tool: str, vargs: list[str]) -> tuple[str | None, str | None]:
    """Find a working `tool`, returning (path, version-line) or (None, None)."""
    candidates: list[str] = []
    if tool == "java":
        # Homebrew's openjdk is keg-only, so it is deliberately absent from PATH.
        # Prefer it over the system shim rather than asking the user to edit their shell.
        candidates += [str(Path(d) / "java") for d in BREW_JDK_BINS if (Path(d) / "java").exists()]
    if (hit := shutil.which(tool)):
        candidates.append(hit)

    for path in candidates:
        if (ver := _runs(path, vargs)) is not None:
            return path, ver
    return None, None


REQUIRED = [
    ("python3", ["--version"], None, "Pipeline is Python throughout."),
    ("java", ["--version"], "brew install openjdk@21",
     "Runs the DREM 2.0.7 and iDREM jars in batch mode."),
]

MANUSCRIPT = [
    ("pandoc", ["--version"], "brew install pandoc",
     "Builds manuscript.docx from the same LaTeX source as the PDF."),
    ("latexmk", ["-version"], "brew install texlive",
     "Drives the pdflatex/bibtex passes for manuscript.pdf."),
    ("pdflatex", ["--version"], "brew install texlive",
     "LaTeX engine."),
]

PY_MODULES = [
    ("numpy", "Deconvolution and model comparison."),
    ("scipy", "Non-negative least squares for cell-type deconvolution."),
    ("pandas", "Count-matrix and metadata handling."),
]


def check(rows, optional: bool):
    results, missing = [], []
    for tool, vargs, install, why in rows:
        path, ver = probe(tool, vargs)
        ok = path is not None
        results.append({
            "tool": tool, "found": ok, "path": path, "version": ver,
            "install": install, "why": why, "optional": optional,
        })
        if not ok:
            missing.append((tool, install, why))
    return results, missing


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manuscript", action="store_true",
                    help="also check the PDF/DOCX build chain")
    ap.add_argument("--json", action="store_true", help="machine-readable report")
    args = ap.parse_args()

    results, missing = check(REQUIRED, optional=False)
    soft_missing: list = []
    if args.manuscript:
        m_results, soft_missing = check(MANUSCRIPT, optional=True)
        results += m_results

    for mod, why in PY_MODULES:
        try:
            __import__(mod)
            found, ver = True, sys.modules[mod].__dict__.get("__version__", "?")
        except ImportError:
            found, ver = False, None
        results.append({
            "tool": f"python:{mod}", "found": found, "path": None, "version": ver,
            "install": "pip install -r requirements.txt", "why": why, "optional": False,
        })
        if not found:
            missing.append((f"python:{mod}", "pip install -r requirements.txt", why))

    if args.json:
        print(json.dumps(results, indent=1))
    else:
        for r in results:
            mark = "ok  " if r["found"] else ("MISS" if not r["optional"] else "miss")
            ver = r["version"] or r["install"] or ""
            print(f"[{mark}] {r['tool']:<16} {ver}")

    if missing:
        log("\nMissing required tools:")
        for tool, install, why in missing:
            log(f"  {tool}: {why}\n    -> {install}")
    if soft_missing:
        log("\nMissing manuscript-build tools (pipeline still runs; `make pdf`/`make docx` will not):")
        for tool, install, why in soft_missing:
            log(f"  {tool}: {why}\n    -> {install}")

    if not missing and (jp := probe("java", ["--version"])[0]):
        # Hand run_all.sh / 08_run_drem.py the resolved interpreter so nobody has to
        # export PATH by hand.
        (Path(__file__).resolve().parent.parent / ".java_path").write_text(jp + "\n")
        if not os.environ.get("JAVA_BIN"):
            log(f"\njava resolved to {jp} (written to .java_path)")

    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
