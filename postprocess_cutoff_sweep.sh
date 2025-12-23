#!/bin/bash
# rc掃引結果から推奨rcを決め、指定inputの COULOMB_CUTOFF を更新する補助スクリプト
#
# 例:
#   bash postprocess_cutoff_sweep.sh results/aor_cutoff_sweep inputs/input_coulomb.dat 0.05 0.5
#
# NOTE:
#   既定では dry-run（更新しない）です。更新したい場合は APPLY=1 を付けて実行してください。
#   APPLY=1 bash postprocess_cutoff_sweep.sh ...

set -euo pipefail

ROOT_DIR="${1:-results/aor_cutoff_sweep}"
INPUT_FILE="${2:-inputs/input_coulomb.dat}"
REF_RC="${3:-0.05}"
TOL_DEG="${4:-0.5}"

SUMMARY_CSV="${ROOT_DIR%/}/sweep_summary.csv"

python3 src/summarize_cutoff_sweep.py \
  --root "${ROOT_DIR}" \
  --out "${SUMMARY_CSV}" \
  --ref-rc "${REF_RC}" \
  --tol-deg "${TOL_DEG}"

DRY="--dry-run"
if [ "${APPLY:-0}" = "1" ]; then
  DRY=""
fi

python3 src/apply_recommended_cutoff.py \
  --summary "${SUMMARY_CSV}" \
  --input "${INPUT_FILE}" \
  --ref-rc "${REF_RC}" \
  --tol-deg "${TOL_DEG}" \
  ${DRY}









