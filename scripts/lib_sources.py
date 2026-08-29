"""Shared helpers for the arabidopsis-drem-osdr acquisition scripts.

Everything here is stdlib-only so `scripts/00_*`–`02_*` can run in a bare CI container.
The analysis scripts (03 onward) use numpy/scipy/pandas.

Why the scripts fetch instead of the browser: both NASA APIs restrict CORS —
`osdr.nasa.gov/osdr/data/*` answers `Access-Control-Allow-Origin: osdr.nasa.gov` and
`visualization.osdr.nasa.gov/biodata/api/v2/` sends no CORS header at all. A page served
from github.io therefore cannot call them. These scripts run offline (or in CI) and bake
JSON into `docs/data/`, which the static site reads.

Ported from AstroRegolith/scripts/lib_sources.py (same author, MIT).
"""
from __future__ import annotations

import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
VENDOR = ROOT / "vendor"
SITE_DATA = ROOT / "docs" / "data"

OSDR_API = "https://osdr.nasa.gov/osdr/data"
OSDR_BIODATA = "https://visualization.osdr.nasa.gov/biodata/api/v2"
OSDR_STUDY_URL = "https://osdr.nasa.gov/bio/repo/data/studies/{acc}"
CROSSREF = "https://api.crossref.org/works/"
MAILTO = "dr.richard.barker@gmail.com"  # Crossref's polite pool

UA = "arabidopsis-drem-osdr/0.1 (https://github.com/dr-richard-barker/arabidopsis-drem-osdr)"


def _open(url: str, timeout: int = 120):
    req = urllib.request.Request(url, headers={"Accept": "*/*", "User-Agent": UA})
    return urllib.request.urlopen(req, timeout=timeout)


def get_json(url: str, retries: int = 3, timeout: int = 120):
    """GET a URL and parse JSON, with a short linear backoff."""
    last = None
    for attempt in range(retries):
        try:
            with _open(url, timeout) as r:
                return json.load(r)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            last = e
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"GET {url} failed after {retries} tries: {last}")


def get_text(url: str, retries: int = 3, timeout: int = 180) -> str:
    last = None
    for attempt in range(retries):
        try:
            with _open(url, timeout) as r:
                return r.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, TimeoutError) as e:
            last = e
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"GET {url} failed after {retries} tries: {last}")


def download(url: str, dest: Path, retries: int = 3, timeout: int = 900) -> Path:
    """Download to `dest` unless it already exists (idempotent re-runs)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    last = None
    for attempt in range(retries):
        try:
            with _open(url, timeout) as r, open(tmp, "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            os.replace(tmp, dest)
            return dest
        except (urllib.error.URLError, TimeoutError) as e:
            last = e
            time.sleep(1.5 * (attempt + 1))
    tmp.unlink(missing_ok=True)
    raise RuntimeError(f"download {url} failed after {retries} tries: {last}")


# --------------------------------------------------------------------------- OSDR

def osdr_files(accession: str) -> dict:
    """The OSDR file listing for e.g. 'OSD-508'."""
    return get_json(f"{OSDR_API}/osd/files/{accession.split('-')[-1]}")


def osdr_isa(accession: str) -> dict:
    """Raw ISA-JSON for a study — carries the publication block (DOI, PMID, authors)."""
    return get_json(f"{OSDR_API}/osd/meta/{accession.split('-')[-1]}")


def osdr_meta(accession: str) -> dict:
    """Flat study metadata (title, organism, factors) from the biodata v2 endpoint.

    The files API carries no title, and `/osdr/data/osd/meta/` returns raw ISA-JSON;
    this endpoint is the one that answers with a stable flat block.
    """
    d = get_json(f"{OSDR_BIODATA}/dataset/{accession}/metadata/")
    return (d.get(accession) or {}).get("metadata", {})


def osdr_factor_table(organism: str = "Arabidopsis thaliana") -> str:
    """Sample-level factor-value table (CSV text) for every study of one organism.

    The `study.factor value` bare parameter asks for *all* factor-value columns; the
    endpoint rejects `format=json`, so CSV is the only machine format it offers.
    """
    q = urllib.parse.urlencode({"study.characteristics.organism": organism})
    return get_text(f"{OSDR_BIODATA}/query/metadata/?{q}&study.factor%20value&format=csv")


def osdr_file_index(accession: str) -> list[tuple[str, str]]:
    """Every (file_name, remote_url) pair in a study's file listing.

    OSDR nests these at varying depths depending on the deposit, so walk rather than
    index into a fixed path.
    """
    pairs: list[tuple[str, str]] = []

    def walk(o):
        if isinstance(o, dict):
            if "file_name" in o and "remote_url" in o:
                pairs.append((o["file_name"], o["remote_url"]))
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk(osdr_files(accession))
    return sorted(set(pairs))


def osdr_url(remote_url: str) -> str:
    """Absolutise a `remote_url` from the file listing.

    The listing gives host-relative paths (`/geode-py/ws/studies/...`); a few entries
    are already absolute, so pass those through untouched.
    """
    if remote_url.startswith("http://") or remote_url.startswith("https://"):
        return remote_url
    return "https://osdr.nasa.gov" + remote_url


def osdr_publication(accession: str) -> dict | None:
    """The first publication block (title/doi/pubMedID/authorList) in a study's ISA-JSON."""
    found: list[dict] = []

    def walk(o):
        if isinstance(o, dict):
            if "title" in o and ("doi" in o or "pubMedID" in o):
                found.append({k: o.get(k) for k in ("title", "doi", "pubMedID", "authorList")})
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk(osdr_isa(accession))
    return found[0] if found else None


# --------------------------------------------------------------------------- misc

def _json_default(o):
    """Coerce numpy scalars the JSON encoder refuses.

    numpy types reach here constantly and invisibly: `bool(x) and np.median(...) <= 20`
    evaluates to np.bool_, not bool, because Python's `and` returns the operand rather
    than a coerced truth value, and `.value_counts().to_dict()` yields np.int64 keys and
    values. Handling it once here beats casting at every call site and forgetting one.
    """
    if hasattr(o, "item"):          # numpy scalar
        return o.item()
    if hasattr(o, "tolist"):        # numpy array
        return o.tolist()
    raise TypeError(f"Object of type {type(o).__name__} is not JSON serializable")


def write_json(path: Path, obj) -> Path:
    """Write pretty, key-sorted JSON so re-runs produce byte-identical files."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, sort_keys=True, ensure_ascii=False,
                               default=_json_default) + "\n")
    return path


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def die(msg: str, code: int = 1):
    log(f"ERROR: {msg}")
    raise SystemExit(code)


def quiet_accelerate_blas_warnings() -> None:
    """Silence the spurious FP warnings NumPy 2.x raises via Apple's Accelerate BLAS.

    On macOS/arm64, `numpy 2.0.x` routes matmul through vecLib, which raises divide-by-
    zero / overflow / invalid flags on the padded lanes of a tiled multiply even when
    every input and output is finite. The scripts that hit this verify their outputs are
    finite instead (see results/qc/*.json), so the warnings are noise that would
    otherwise train a reader to ignore real ones.
    """
    import warnings
    for msg in ("divide by zero encountered in matmul",
                "overflow encountered in matmul",
                "invalid value encountered in matmul"):
        warnings.filterwarnings("ignore", message=msg, category=RuntimeWarning)
