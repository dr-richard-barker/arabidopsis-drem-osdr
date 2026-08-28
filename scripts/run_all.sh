#!/usr/bin/env bash
# Reproduce the whole analysis from public sources.
#
#   bash scripts/run_all.sh              full run (downloads ~400 MB, DREM takes hours)
#   bash scripts/run_all.sh --dry-run    resolve every URL, download and fit nothing
#   bash scripts/run_all.sh --no-drem    everything except the DREM fits
#
# Each step is idempotent: downloads are skipped when the file is already present, so a
# re-run after an interruption resumes rather than restarts.
set -euo pipefail
cd "$(dirname "$0")/.."

DRY=""; NO_DREM=""
for a in "$@"; do
  case "$a" in
    --dry-run) DRY="--dry-run" ;;
    --no-drem) NO_DREM="1" ;;
    *) echo "unknown option: $a" >&2; exit 2 ;;
  esac
done

step() { echo; echo "=== $* ==="; }

step "00 environment"
python3 scripts/00_env_check.py --manuscript

step "01 acquire OSDR"
python3 scripts/01_fetch_osdr.py $DRY

if [ -n "$DRY" ]; then echo; echo "dry run complete — nothing downloaded or fitted."; exit 0; fi

step "02 harmonise metadata"
python3 scripts/02_harmonise_metadata.py

step "03 build pseudo-time-series"
python3 scripts/03_build_pseudotimeseries.py

step "04 TF prior from DAP-seq"
python3 scripts/04_fetch_tf_prior.py

step "04b SOG1 edges from ChIP-seq"
python3 scripts/04b_sog1_chip_prior.py

step "05 deconvolve cell types"
python3 scripts/05_deconvolve_celltypes.py

step "06 cell-type-weighted prior"
python3 scripts/06_weight_tf_prior.py

step "07 write DREM inputs"
python3 scripts/07_write_drem_inputs.py

if [ -z "$NO_DREM" ]; then
  step "08 run DREM (hours)"
  python3 scripts/08_run_drem.py
  step "09 parse models"
  python3 scripts/09_parse_drem_model.py
  step "10 compare models"
  python3 scripts/10_compare_models.py
else
  echo; echo "skipping DREM fits (--no-drem)"
fi

step "12 figures"
python3 scripts/12_figures.py || echo "  (figures step skipped)"

step "13 manuscript numbers"
python3 scripts/13_manuscript_numbers.py

step "14 manifest"
python3 scripts/14_manifest.py

step "references"
python3 scripts/check_references.py

echo; echo "done. Build the manuscript with:  cd manuscript && make pdf docx"
