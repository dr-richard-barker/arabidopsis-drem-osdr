#!/usr/bin/env python3
"""Run DREM 2.0 in batch mode over every prepared input directory.

The engine is the real thing — Jason Ernst's `drem.jar`, downloaded and checksummed
rather than vendored, so results stay citable as DREM rather than as a reimplementation.
Its batch interface is:

    java -jar drem.jar -b settings.txt outmodel.txt geneassign.txt TFSCOREDIR

  results/drem/runs/<arm>__<prior>/model.txt        the fitted model
  results/drem/runs/<arm>__<prior>/geneassign.txt   gene -> path assignments
  results/drem/runs/<arm>__<prior>/tfscores/        per-split TF activity scores
  results/drem/runs/<arm>__<prior>/drem.log         stdout+stderr, kept for the record
  results/drem/run_report.json                      exit status and timing per run

Runs are independent, so a failure is recorded and the rest continue; the exit code is
non-zero if any run failed.

  python3 scripts/08_run_drem.py
  python3 scripts/08_run_drem.py --only A_primary_WildType__weighted --timeout 3600
  python3 scripts/08_run_drem.py --smoke-test     # iDREM's bundled example, no real data
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib_sources import RESULTS, ROOT, VENDOR, download, log, write_json  # noqa: E402

INPUTS = RESULTS / "drem" / "inputs"
RUNS = RESULTS / "drem" / "runs"
REPORT = RESULTS / "drem" / "run_report.json"

DREM_JAR_URL = "https://raw.githubusercontent.com/jernst98/STEM_DREM/master/drem.jar"
# Verified 2026-08-28 against the file this pipeline was developed with. A changed
# digest is reported, not fatal: upstream may legitimately publish a new build, and
# silently running a different engine would be worse than a warning.
DREM_JAR_SHA256 = "16312a6a713238a829abb0e1bf3929f8378216ddb17b57fe30e038e892739530"


def java_bin() -> str:
    """The JDK `00_env_check.py` resolved, else whatever is on PATH."""
    marker = ROOT / ".java_path"
    if marker.exists():
        p = marker.read_text().strip()
        if p and Path(p).exists():
            return p
    return shutil.which("java") or "java"


def ensure_jar() -> tuple[Path, dict]:
    jar = VENDOR / "drem.jar"
    download(DREM_JAR_URL, jar)
    digest = hashlib.sha256(jar.read_bytes()).hexdigest()
    info = {"path": str(jar.relative_to(ROOT)), "url": DREM_JAR_URL,
            "sha256": digest, "expected_sha256": DREM_JAR_SHA256,
            "matches_expected": digest == DREM_JAR_SHA256,
            "licence": "GPL-3.0 (Ernst lab, jernst98/STEM_DREM) — downloaded, not redistributed"}
    if not info["matches_expected"]:
        log(f"WARNING: drem.jar sha256 {digest} != expected {DREM_JAR_SHA256}. "
            f"Upstream may have republished; results are from the downloaded build.")
    return jar, info


def run_one(jar: Path, run_dir: Path, out_dir: Path, heap: str, timeout: int) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    tfscores = out_dir / "tfscores"
    tfscores.mkdir(exist_ok=True)
    cmd = [java_bin(), f"-Xmx{heap}", "-jar", str(jar), "-b",
           str(run_dir / "settings.txt"), str(out_dir / "model.txt"),
           str(out_dir / "geneassign.txt"), str(tfscores)]

    started = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                              cwd=str(run_dir),
                              # DREM is a Swing application even in batch mode; without
                              # headless mode it tries to open a display and dies on a
                              # server or in CI.
                              env={**os.environ, "JAVA_TOOL_OPTIONS": "-Djava.awt.headless=true"})
        rc, out, err = proc.returncode, proc.stdout, proc.stderr
        timed_out = False
    except subprocess.TimeoutExpired as e:
        rc, out, err, timed_out = -1, (e.stdout or ""), (e.stderr or ""), True

    elapsed = round(time.time() - started, 1)
    (out_dir / "drem.log").write_text(
        f"$ {' '.join(cmd)}\n\n--- stdout ---\n{out}\n--- stderr ---\n{err}\n")

    model = out_dir / "model.txt"
    assign = out_dir / "geneassign.txt"
    result = {
        "exit_code": rc, "timed_out": timed_out, "seconds": elapsed,
        "model_written": model.exists() and model.stat().st_size > 0,
        "model_bytes": model.stat().st_size if model.exists() else 0,
        "geneassign_written": assign.exists() and assign.stat().st_size > 0,
        "n_tfscore_files": len(list(tfscores.glob("*"))),
    }
    result["ok"] = bool(result["model_written"])
    if not result["ok"]:
        result["stderr_tail"] = (err or out or "")[-800:]
    return result


def smoke_test(jar: Path, heap: str, timeout: int) -> int:
    """Prove the wrapper against DREM's own inputs before trusting it on real data."""
    import urllib.request
    sandbox = RESULTS / "drem" / "smoke"
    sandbox.mkdir(parents=True, exist_ok=True)
    base = "https://raw.githubusercontent.com/phoenixding/idrem/master/example/inputs/"
    files = {"expression.txt": "example_expression_data_file.txt"}
    for local, remote in files.items():
        dest = sandbox / local
        if not dest.exists():
            urllib.request.urlretrieve(base + remote, dest)  # noqa: S310 - fixed URL
            log(f"  fetched {remote}")

    # A minimal all-ones TF file over the example's own genes, so the smoke test
    # exercises the batch path without depending on a species-specific prior.
    genes = [l.split("\t")[0] for l in
             (sandbox / "expression.txt").read_text().splitlines()[1:] if l.strip()]
    (sandbox / "tf.txt").write_text(
        "TF\tGene\tInput\n" + "".join(f"TF1\t{g}\t1\n" for g in genes[:200]))

    from importlib import import_module
    tmpl = import_module("07_write_drem_inputs").SETTINGS_TEMPLATE
    (sandbox / "settings.txt").write_text(tmpl.format(
        tf_file=(sandbox / "tf.txt").resolve(),
        expression_file=(sandbox / "expression.txt").resolve(),
        repeat_file="", min_abs_log_ratio=0.5, max_paths=3, node_penalty=40, seed=1260))

    res = run_one(jar, sandbox, sandbox / "out", heap, timeout)
    log(f"smoke test: {'PASS' if res['ok'] else 'FAIL'} "
        f"(exit {res['exit_code']}, {res['seconds']}s, model {res['model_bytes']} bytes)")
    if not res["ok"]:
        log(res.get("stderr_tail", ""))
    return 0 if res["ok"] else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", help="run just this input directory name")
    ap.add_argument("--heap", default="8g", help="JVM max heap (default 8g)")
    ap.add_argument("--timeout", type=int, default=7200, help="per-run seconds (default 7200)")
    ap.add_argument("--smoke-test", action="store_true",
                    help="run DREM's own example through the wrapper and stop")
    args = ap.parse_args()

    jar, jar_info = ensure_jar()
    log(f"engine: {jar_info['path']} sha256={jar_info['sha256'][:16]}... "
        f"(expected match: {jar_info['matches_expected']})")
    log(f"java: {java_bin()}")

    if args.smoke_test:
        return smoke_test(jar, args.heap, args.timeout)

    dirs = sorted(d for d in INPUTS.iterdir() if d.is_dir())
    if args.only:
        dirs = [d for d in dirs if d.name == args.only]
    if not dirs:
        raise SystemExit("no prepared inputs — run 07_write_drem_inputs.py first")

    report = {"generated": dt.date.today().isoformat(), "engine": jar_info,
              "java": java_bin(), "runs": {}}
    failed = []
    for d in dirs:
        log(f"\n[{d.name}]")
        res = run_one(jar, d, RUNS / d.name, args.heap, args.timeout)
        report["runs"][d.name] = res
        log(f"  {'ok' if res['ok'] else 'FAILED'} in {res['seconds']}s "
            f"(model {res['model_bytes']} bytes, {res['n_tfscore_files']} TF score files)")
        if not res["ok"]:
            failed.append(d.name)
            log(f"  {res.get('stderr_tail', '')[:400]}")

    write_json(REPORT, report)
    log(f"\nwrote {REPORT}  ({len(dirs) - len(failed)}/{len(dirs)} succeeded)")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
