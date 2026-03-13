"""Generate reproducible ERD identification datasets (train/test by default).

This script is the canonical data factory for the toy twin. It runs the plant
with persistently exciting piecewise-constant controls, writes standard run
artifacts for each trajectory, and produces a manifest consumed by
``model/scripts/run_pipeline.py``.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime
from pathlib import Path
import sys
from typing import Dict, List

import numpy as np
import yaml

if __package__ in (None, ""):
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)

from erd_fipy import ERDPlant, load_run_config, preset_config


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for dataset generation.

    Returns:
        Configured argument parser.
    """

    p = argparse.ArgumentParser(description="Generate ERD train/test trajectories for ROM identification")
    p.add_argument("--config", type=str, default=None, help="Optional YAML override for the plant config")
    p.add_argument("--preset", type=str, default="smoke", choices=["smoke", "report"])
    p.add_argument("--seed", type=int, default=0, help="Base RNG seed; each trajectory increments it")
    p.add_argument("--tag", type=str, default="dataset", help="Human-friendly dataset tag")
    p.add_argument("--manifest", type=str, default=None, help="Optional explicit manifest output path")
    p.add_argument("--n-train", type=int, default=None, help="Override train trajectory count")
    p.add_argument("--n-val", type=int, default=None, help="Override validation trajectory count")
    p.add_argument("--n-test", type=int, default=None, help="Override test trajectory count")
    p.add_argument("--block-steps", type=int, default=None, help="Piecewise-constant control block length")
    p.add_argument("--train-time", type=float, default=None, help="Optional train trajectory duration")
    p.add_argument("--val-time", type=float, default=None, help="Optional validation trajectory duration")
    p.add_argument("--test-time", type=float, default=None, help="Optional test trajectory duration")
    p.add_argument(
        "--skip-media",
        action="store_true",
        help="Skip plot/GIF rendering for faster dataset generation",
    )
    return p


def random_piecewise_schedule(
    rng: np.random.Generator,
    bounds: np.ndarray,
    base_u0: float,
    block_steps: int,
    t_final: float,
    warmup_fraction: float,
    excite_fraction: float,
):
    """Create a piecewise-constant random control schedule.

    Args:
        rng: Random-number generator for reproducible sampling.
        bounds: Per-channel absolute bounds ``[u0, ..., u4]``.
        base_u0: Nominal drive value for ``u0``.
        block_steps: Number of plant steps per control block.
        t_final: Final simulation time used for staged excitation.
        warmup_fraction: Early-time low-excitation fraction.
        excite_fraction: Mid-run high-excitation fraction.

    Returns:
        Schedule callback ``u = f(k, t, n, omega)`` compatible with
        :meth:`erd_fipy.stepping.ERDPlant.run`.
    """

    current = np.array([base_u0, 0.0, 0.0, 0.0, 0.0], dtype=float)

    def schedule(k: int, _t: float, _n: np.ndarray, _omega: np.ndarray) -> np.ndarray:
        """Return the current control block value, resampling on block boundaries."""

        nonlocal current
        if (k % block_steps) == 0:
            sample = rng.uniform(-bounds, bounds)
            frac = float(np.clip(_t / max(t_final, 1e-12), 0.0, 1.0))
            warm = float(np.clip(warmup_fraction, 0.0, 1.0))
            excite = float(np.clip(excite_fraction, 0.0, 1.0 - warm))

            if frac < warm:
                asym_scale = 0.30
                u0_scale = 0.20
            elif frac < warm + excite:
                asym_scale = 1.00
                u0_scale = 0.55
            else:
                asym_scale = 1.20
                u0_scale = 0.50
                sample[1:] += rng.normal(loc=0.0, scale=0.15 * bounds[1:], size=4)

            sample[1:] = np.clip(sample[1:] * asym_scale, -bounds[1:], bounds[1:])
            # Keep axisymmetric drive near operating point while still
            # modulating it enough to avoid near-constant controls.
            sample[0] = np.clip(base_u0 + u0_scale * sample[0], -bounds[0], bounds[0])
            current = sample
        return current

    return schedule


def _tfinal_for_split(base_tfinal: float, split_name: str, args: argparse.Namespace) -> float:
    """Resolve trajectory duration for a split from CLI overrides.

    Args:
        base_tfinal: Default ``T_final`` from the selected run config.
        split_name: Dataset split name (``train``, ``val``, or ``test``).
        args: Parsed CLI arguments.

    Returns:
        Duration to use for trajectories in the given split.
    """

    if split_name == "train" and args.train_time is not None:
        return float(args.train_time)
    if split_name == "val" and args.val_time is not None:
        return float(args.val_time)
    if split_name == "test" and args.test_time is not None:
        return float(args.test_time)
    return float(base_tfinal)


def main() -> None:
    """Generate train/val/test trajectories and write a dataset manifest."""

    args = build_parser().parse_args()
    base_cfg = load_run_config(args.config, preset=args.preset)
    preset = preset_config(args.preset)

    counts = {
        "train": args.n_train if args.n_train is not None else preset.split.n_train,
        "val": args.n_val if args.n_val is not None else preset.split.n_val,
        "test": args.n_test if args.n_test is not None else preset.split.n_test,
    }
    block_steps = int(args.block_steps if args.block_steps is not None else preset.block_steps)
    bounds = np.asarray(base_cfg.forcing.u_bounds, dtype=float)

    outputs_root = Path(base_cfg.output.outputs_root).resolve()
    outputs_root.mkdir(parents=True, exist_ok=True)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    dataset_id = f"{stamp}_{args.tag}_{args.preset}"

    manifest_entries: Dict[str, List[Dict[str, object]]] = {"train": [], "val": [], "test": []}
    seed = int(args.seed)

    if args.skip_media:
        # RunWriter reads this env var to skip plot/GIF generation.
        import os

        os.environ["ERD_SKIP_MEDIA"] = "1"

    for split_name in ("train", "val", "test"):
        n_split = int(counts[split_name])
        split_tfinal = _tfinal_for_split(base_cfg.time.T_final, split_name, args)
        for i in range(n_split):
            local_seed = seed
            seed += 1
            rng = np.random.default_rng(local_seed)

            run_tag = f"{dataset_id}_{split_name}_{i:03d}"

            cfg_i = replace(
                base_cfg,
                time=replace(base_cfg.time, T_final=split_tfinal),
                init=replace(base_cfg.init, seed=local_seed),
                output=replace(base_cfg.output, tag=run_tag, preset=args.preset),
            )

            plant = ERDPlant(cfg_i)
            schedule = random_piecewise_schedule(
                rng=rng,
                bounds=bounds,
                base_u0=cfg_i.forcing.drive_u0_base,
                block_steps=block_steps,
                t_final=cfg_i.time.T_final,
                warmup_fraction=cfg_i.disturbance.warmup_fraction,
                excite_fraction=cfg_i.disturbance.excite_fraction,
            )
            result = plant.run(control_schedule=schedule, write_artifacts=True)

            manifest_entries[split_name].append(
                {
                    "run_dir": result.run_dir,
                    "dataset_id": dataset_id,
                    "seed": local_seed,
                    "tag": run_tag,
                    "split": split_name,
                    "index": i,
                    "n_steps": cfg_i.time.n_steps,
                    "dt": cfg_i.time.dt,
                    "T_final": cfg_i.time.T_final,
                }
            )
            print(f"[{split_name}] {i + 1}/{n_split} -> {result.run_dir}")

    if args.manifest is None:
        manifest_path = outputs_root / f"{dataset_id}_manifest.yaml"
    else:
        manifest_path = Path(args.manifest).resolve()

    payload = {
        "dataset_id": dataset_id,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "preset": args.preset,
        "tag": args.tag,
        "seed": int(args.seed),
        "block_steps": block_steps,
        "counts": counts,
        "base_config": base_cfg.to_dict(),
        "splits": manifest_entries,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False)

    print(f"dataset_id: {dataset_id}")
    print(f"manifest: {manifest_path}")


if __name__ == "__main__":
    main()
