# ERD Capstone: Toy Twin + ROM + MPC Guide

This repository implements a complete toy digital-twin workflow for a ring-discharge-like system:

1. `erd_fipy` generates high-fidelity field data from a simple PDE plant.
2. `model/erd_model` fits POD + controlled SINDy reduced dynamics.
3. `control/erd_control` runs baseline-vs-MPC closed-loop comparisons using the learned ROM.

The design target is a **low-dimensional, control-ready toy** where the control path is intentionally indirect:

`u(t) -> f_omega -> omega -> psi -> velocity -> n -> metrics`.

## Repository layout

```text
erd-capstone/
├── erd_fipy/
│   ├── erd_fipy/          # plant, metrics, IO, stepping
│   └── scripts/           # run_demo.py, generate_training_runs.py
├── model/
│   ├── erd_model/         # POD + SINDy pipeline
│   └── scripts/           # run_pipeline.py, default_config.yaml
├── control/
│   ├── erd_control/       # controller + closed-loop orchestration
│   └── scripts/           # run_closed_loop.py, run_rom_case_study.py, default_config.yaml
└── outputs/runs/          # canonical run artifacts
```

## Environment

From repo root:

```bash
uv sync --python 3.10 --group dev
```

Run commands with `uv run ...`.

## Canonical workflow

### 1) Generate training/test data

The dataset generator is the canonical source of manifests and trajectory runs.

```bash
uv run python erd_fipy/scripts/generate_training_runs.py \
  --preset smoke \
  --tag meeting_demo \
  --seed 42
```

Defaults are intentionally simple:

- `smoke`: `1` train run, `0` val runs, `1` test run.
- reproducible seed-incremented trajectories,
- persistently exciting piecewise-constant controls,
- deterministic disturbances with per-run phase jitter (still reproducible by seed).

Useful overrides:

```bash
--n-train --n-val --n-test --block-steps --train-time --test-time --skip-media
```

### 2) Fit POD + controlled SINDy

Copy and edit `model/scripts/default_config.yaml`.

Important: `data.manifest_path` is required and must point to your dataset manifest.

```bash
uv run python model/scripts/run_pipeline.py --config model/scripts/default_config.yaml
```

Model run outputs include:

- `model/basis_U.npy`, `model/mean_state.npy`, `model/sindy_model.pkl`, `model/Xi.npy`
- `model/training_manifest.yaml`
- `model/lineage.json` (dataset linkage)
- `metrics/train_curves.json` and `metrics/test_curves.json` (when available)

### 3) Run closed-loop baseline vs MPC

Copy and edit `control/scripts/default_config.yaml`.

Set:

- `model_run_dir` (required)
- optionally `plant_config_path` (otherwise inferred from model manifest lineage)

Then run:

```bash
uv run python control/scripts/run_closed_loop.py --config control/scripts/default_config.yaml
```

Closed-loop outputs include:

- `stages/open_loop/*`
- `stages/closed_loop/*`
- `plots/open_vs_closed_metrics.png`
- `metrics/relative_deltas.json`
- `model/lineage.json`

## Run-folder and lineage conventions

All pipelines write timestamped folders:

```text
outputs/runs/<YYYYmmdd_HHMMSS>_<tag>/
```

Lineage is explicit via metadata files:

- dataset manifest: `dataset_id` + split entries with `run_dir`
- model run: `model/lineage.json` links to dataset manifest and run dirs
- control run: `model/lineage.json` links to model run and stage run dirs

Recommended tags:

- dataset generation: `--tag <project_or_meeting_name>`
- model/control configs: `tag: auto` (derives lineage-aware tags)

## Config ownership rules

To avoid conflicting inputs:

- **Plant physics/time/grid** are owned by `erd_fipy` config.
- **Model fitting options** are owned by `model` config.
- **Closed-loop orchestration options** are owned by `control` config.
- Cross-stage linkage happens by paths (`manifest_path`, `model_run_dir`), not duplicated knobs.

## Notes on dynamics and demonstration quality

Current defaults are tuned for short, visually nontrivial runs:

- lower damping/diffusion than early versions,
- stronger disturbances and initial perturbations,
- stronger admissible control amplitudes,
- short smoke horizon suitable for meeting demos.

If dynamics are too quiet or too unstable, tune only a few knobs first:

1. `pde.gamma`, `pde.nu`, `pde.D_r`, `pde.D_phi`
2. disturbance amplitudes in `disturbance.mode1/mode2`
3. `forcing.drive_u0_base` and `forcing.u_bounds`

## Minimal smoke commands (copy/paste)

```bash
# 1) dataset
uv run python erd_fipy/scripts/generate_training_runs.py --preset smoke --tag smoke_demo --seed 7

# 2) edit manifest path in model/scripts/default_config.yaml, then fit
uv run python model/scripts/run_pipeline.py --config model/scripts/default_config.yaml

# 3) edit model_run_dir in control/scripts/default_config.yaml, then control
uv run python control/scripts/run_closed_loop.py --config control/scripts/default_config.yaml
```

## Troubleshooting

- `FiPy is required but not installed.`
  - Run `uv sync --python 3.10 --group dev`.
- `data.manifest_path must be set ...`
  - Set `model` config `data.manifest_path` explicitly.
- Non-finite ROM rollouts in MPC:
  - lower model rank (`pod.rank`), increase SINDy threshold, and/or tighten control bounds.

