#!/usr/bin/env bash
set -euo pipefail

# End-to-end extended preset with higher spatial/temporal complexity.
# Produces:
# 1) dataset runs + manifest
# 2) model fit run
# 3) closed-loop baseline-vs-MPC run

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

OUTPUTS_ROOT="${REPO_ROOT}/outputs/runs"
mkdir -p "${OUTPUTS_ROOT}"

PLANT_BASE_CFG="$(mktemp /tmp/erd_extended_base_XXXX.yaml)"
MODEL_CFG="$(mktemp /tmp/erd_extended_model_XXXX.yaml)"
PLANT_CTRL_CFG="$(mktemp /tmp/erd_extended_ctrlplant_XXXX.yaml)"
CONTROL_CFG="$(mktemp /tmp/erd_extended_control_XXXX.yaml)"
DATA_LOG="$(mktemp /tmp/erd_extended_data_XXXX.log)"
MODEL_LOG="$(mktemp /tmp/erd_extended_model_XXXX.log)"
CONTROL_LOG="$(mktemp /tmp/erd_extended_control_XXXX.log)"

cleanup() {
  rm -f "${PLANT_BASE_CFG}" "${MODEL_CFG}" "${PLANT_CTRL_CFG}" "${CONTROL_CFG}" \
        "${DATA_LOG}" "${MODEL_LOG}" "${CONTROL_LOG}"
}
trap cleanup EXIT

printf "[extended] Building base plant config...\n"
uv run python - "${PLANT_BASE_CFG}" "${OUTPUTS_ROOT}" <<'PY'
from dataclasses import replace
from pathlib import Path
import sys
from erd_fipy import load_run_config, save_run_config

out_path = Path(sys.argv[1])
outputs_root = sys.argv[2]
cfg = load_run_config(None, preset="report")
cfg = replace(
    cfg,
    domain=replace(cfg.domain, N_r=48, N_phi=96),
    time=replace(cfg.time, dt=0.0015, T_final=1.5, snapshot_every=4),
    output=replace(cfg.output, outputs_root=outputs_root, tag="preset_extended_media"),
)
save_run_config(out_path, cfg)
PY

printf "[extended] Generating dataset...\n"
uv run python erd_fipy/scripts/generate_training_runs.py \
  --config "${PLANT_BASE_CFG}" \
  --preset report \
  --tag preset_extended_media \
  --seed 41 \
  --n-train 12 \
  --n-val 0 \
  --n-test 4 \
  --block-steps 16 \
  --train-time 1.5 \
  --test-time 1.5 | tee "${DATA_LOG}"

MANIFEST_PATH="$(awk -F': ' '/^manifest: /{print $2}' "${DATA_LOG}" | tail -1)"
DATASET_ID="$(awk -F': ' '/^dataset_id: /{print $2}' "${DATA_LOG}" | tail -1)"
if [[ -z "${MANIFEST_PATH}" || ! -f "${MANIFEST_PATH}" ]]; then
  echo "Failed to resolve manifest path from dataset generation output" >&2
  exit 1
fi

printf "[extended] Writing model config...\n"
cat > "${MODEL_CFG}" <<YAML
name: preset_extended_model
seed: 41
data:
  manifest_path: ${MANIFEST_PATH}
  train_split: train
  val_split: val
  test_split: test
  max_train: 12
  max_val: null
  max_test: 4
  time_stride: 1
  time_limit: null
pod:
  rank: 7
  center: true
deriv:
  order: 2
sindy:
  optimizer: trappingsr3
  threshold: 0.05
  alpha: 0.0001
  max_iter: 100
  backoff_factor: 0.1
  backoff_retries: 3
  fallback_optimizer: stlsq
  validation_spike_threshold: 0.6
  trapping:
    method: local
    reg_weight_lam: 0.01
    relax_coeff_nu: 1.0
    tol: 1.0e-6
    max_iter: 240
    eta: 10.0
    stability_alpha: 1.0e-3
    stability_beta: 1.0e-3
    control_threshold: 0.015
    control_alpha: 1.0e-4
    control_max_iter: 6
acceptance:
  min_pod_energy: 0.8
  max_rank: 99
  max_mean_field_rel_l2: 0.20
rollout:
  horizon_steps: 180
plots:
  dpi: 120
  movie_fps: 12
output:
  outputs_root: ${OUTPUTS_ROOT}
  tag: extended_model_${DATASET_ID}
YAML

printf "[extended] Fitting model...\n"
uv run python model/scripts/run_pipeline.py --config "${MODEL_CFG}" | tee "${MODEL_LOG}"
MODEL_RUN_DIR="$(awk -F': ' '/^run_dir: /{print $2}' "${MODEL_LOG}" | tail -1)"
if [[ -z "${MODEL_RUN_DIR}" || ! -d "${MODEL_RUN_DIR}" ]]; then
  echo "Failed to resolve model run_dir from model output" >&2
  exit 1
fi

printf "[extended] Building control plant config...\n"
uv run python - "${MANIFEST_PATH}" "${PLANT_CTRL_CFG}" "${OUTPUTS_ROOT}" <<'PY'
import sys
import yaml
from pathlib import Path

manifest_path = Path(sys.argv[1])
out_path = Path(sys.argv[2])
outputs_root = sys.argv[3]
with manifest_path.open("r", encoding="utf-8") as f:
    manifest = yaml.safe_load(f) or {}
base_cfg = manifest["base_config"]
base_cfg.setdefault("control", {})
base_cfg["control"]["H"] = 12
base_cfg["control"]["N_shoot"] = 128
base_cfg["control"]["rate_penalty"] = 1.0
base_cfg["control"]["weights"] = {
    "w_wob_abs": 5.0e10,
    "w_wob_growth": 2.0e11,
    "w_wall": 0.4,
    "w_coh": 0.3,
    "w_coh_growth": 0.6,
    "w_wall_slope": 0.15,
    "w_pow": 0.02,
    "w_rate": 0.15,
}
base_cfg.setdefault("output", {})
base_cfg["output"]["outputs_root"] = outputs_root
with out_path.open("w", encoding="utf-8") as f:
    yaml.safe_dump(base_cfg, f, sort_keys=False)
PY

printf "[extended] Writing control config...\n"
cat > "${CONTROL_CFG}" <<YAML
model_run_dir: ${MODEL_RUN_DIR}
plant_config_path: ${PLANT_CTRL_CFG}
preset: report
seed: 23
outputs_root: ${OUTPUTS_ROOT}
tag: extended_control_${DATASET_ID}
YAML

printf "[extended] Running closed-loop baseline vs MPC...\n"
uv run python control/scripts/run_closed_loop.py --config "${CONTROL_CFG}" | tee "${CONTROL_LOG}"
CONTROL_RUN_DIR="$(awk -F': ' '/^run_dir: /{print $2}' "${CONTROL_LOG}" | tail -1)"
if [[ -z "${CONTROL_RUN_DIR}" || ! -d "${CONTROL_RUN_DIR}" ]]; then
  echo "Failed to resolve control run_dir from control output" >&2
  exit 1
fi

printf "\n[extended] Complete\n"
printf "  dataset_manifest: %s\n" "${MANIFEST_PATH}"
printf "  model_run_dir:    %s\n" "${MODEL_RUN_DIR}"
printf "  control_run_dir:  %s\n" "${CONTROL_RUN_DIR}"
printf "  summary:          %s\n" "${CONTROL_RUN_DIR}/summary.json"
printf "  plot:             %s\n" "${CONTROL_RUN_DIR}/plots/open_vs_closed_metrics.png"
