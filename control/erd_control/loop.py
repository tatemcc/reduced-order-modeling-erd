"""Closed-loop orchestration for ERD baseline-vs-MPC comparisons."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime
import json
import logging
from pathlib import Path
import shutil
from typing import Dict

import matplotlib.pyplot as plt
import yaml

from erd_fipy import ERDPlant, default_open_loop_control, load_run_config, run_config_from_dict
from erd_fipy.config import RunConfig

from .config import ControlRunConfig, to_dict
from .controller import RandomShootingMPCController



def _make_top_run_dir(outputs_root: str, tag: str) -> Path:
    """Create canonical closed-loop run directory with required subfolders.

    Args:
        outputs_root: Parent directory for all runs.
        tag: Suffix tag used in the timestamped run folder name.

    Returns:
        Newly created top-level run directory path.
    """

    root = Path(outputs_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base = root / f"{stamp}_{tag}"
    out = base
    idx = 1
    while out.exists():
        out = Path(f"{base}_{idx:02d}")
        idx += 1

    (out / "stages" / "open_loop").mkdir(parents=True, exist_ok=True)
    (out / "stages" / "closed_loop").mkdir(parents=True, exist_ok=True)
    (out / "plots").mkdir(parents=True, exist_ok=True)
    (out / "movies").mkdir(parents=True, exist_ok=True)
    (out / "model").mkdir(parents=True, exist_ok=True)
    (out / "metrics").mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(parents=True, exist_ok=True)
    return out



def _setup_logger(run_dir: Path) -> logging.Logger:
    """Create file + stream logger for closed-loop orchestration.

    Args:
        run_dir: Top-level closed-loop run directory.

    Returns:
        Configured logger instance.
    """

    logger = logging.getLogger(f"erd_control_{run_dir.name}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    fh = logging.FileHandler(run_dir / "logs" / "run.log")
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(fh)

    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(sh)
    return logger



def _load_plant_cfg(cfg: ControlRunConfig) -> RunConfig:
    """Resolve plant configuration either directly or from training manifest.

    Args:
        cfg: Closed-loop run configuration.

    Returns:
        Plant :class:`RunConfig` object.
    """

    if cfg.plant_config_path is not None:
        return load_run_config(cfg.plant_config_path, preset=cfg.preset)

    manifest_path = Path(cfg.model_run_dir) / "model" / "training_manifest.yaml"
    with manifest_path.open("r", encoding="utf-8") as f:
        manifest = yaml.safe_load(f) or {}

    base_cfg = manifest.get("base_config")
    if base_cfg is None:
        raise ValueError("Could not infer plant config from model manifest; provide plant_config_path")
    return run_config_from_dict(base_cfg)



def _copy_model_artifacts(model_run_dir: Path, dst: Path) -> None:
    """Copy serialized model artifacts into the canonical closed-loop folder.

    Args:
        model_run_dir: Source model run directory.
        dst: Destination ``model/`` folder inside closed-loop run.

    Returns:
        None.
    """

    src = model_run_dir / "model"
    for item in src.iterdir():
        if item.is_file():
            shutil.copy2(item, dst / item.name)



def _plot_comparison(run_dir: Path, open_curves: Dict, closed_curves: Dict) -> None:
    """Create summary comparison plots for open-loop vs closed-loop runs.

    Args:
        run_dir: Top-level closed-loop run directory.
        open_curves: Metric curves from the baseline run.
        closed_curves: Metric curves from the MPC run.

    Returns:
        None.
    """

    t_o = open_curves["t"]
    t_c = closed_curves["t"]

    fig, axes = plt.subplots(2, 3, figsize=(14, 7))
    axes = axes.reshape(-1)

    axes[0].plot(t_o, open_curves["E_wob"], label="open")
    axes[0].plot(t_c, closed_curves["E_wob"], label="closed")
    axes[0].set_title("E_wob")

    axes[1].plot(t_o, open_curves["L_w"], label="open")
    axes[1].plot(t_c, closed_curves["L_w"], label="closed")
    axes[1].set_title("L_w")

    axes[2].plot(t_o, open_curves["sigma_r"], label="open")
    axes[2].plot(t_c, closed_curves["sigma_r"], label="closed")
    axes[2].set_title("sigma_r")

    axes[3].plot(t_o, open_curves["P_ctrl"], label="open")
    axes[3].plot(t_c, closed_curves["P_ctrl"], label="closed")
    axes[3].set_title("P_ctrl")

    axes[4].plot(t_o, open_curves.get("L_w_cum", []), label="open")
    axes[4].plot(t_c, closed_curves.get("L_w_cum", []), label="closed")
    axes[4].set_title("L_w_cum")

    axes[5].plot(t_o, open_curves.get("badness_score", []), label="open")
    axes[5].plot(t_c, closed_curves.get("badness_score", []), label="closed")
    axes[5].set_title("badness_score")

    for ax in axes:
        ax.grid(alpha=0.3)
        ax.set_xlabel("t")
        ax.legend()

    fig.tight_layout()
    fig.savefig(run_dir / "plots" / "open_vs_closed_metrics.png", dpi=150)
    plt.close(fig)



def run_closed_loop(cfg: ControlRunConfig) -> str:
    """Run baseline + MPC simulations and write a canonical comparison folder.

    Execution sequence:
    1) Create top-level run directory and logger.
    2) Load plant config and ROM artifacts from the model run.
    3) Run an open-loop baseline stage from identical initial conditions.
    4) Run a closed-loop stage using random-shooting MPC at each plant step.
    5) Copy model artifacts, generate comparison plots, and write summary JSON.

    Args:
        cfg: Closed-loop run configuration.

    Returns:
        Path to the top-level closed-loop run directory.
    """

    model_run_name = Path(cfg.model_run_dir).resolve().name
    run_tag = cfg.tag if cfg.tag.strip().lower() not in {"", "auto"} else f"control_{model_run_name}"
    run_dir = _make_top_run_dir(cfg.outputs_root, run_tag)
    logger = _setup_logger(run_dir)

    logger.info("Loading model and plant configuration")
    plant_cfg = _load_plant_cfg(cfg)

    model_bundle = RandomShootingMPCController.load_model_bundle(cfg.model_run_dir)
    controller = RandomShootingMPCController(model_bundle=model_bundle, plant_cfg=plant_cfg, seed=cfg.seed)

    logger.info("Running open-loop baseline")
    open_cfg = replace(
        plant_cfg,
        output=replace(plant_cfg.output, run_dir=str(run_dir / "stages" / "open_loop"), tag="open_loop"),
    )
    open_plant = ERDPlant(open_cfg)
    open_result = open_plant.run(write_artifacts=True)

    logger.info("Running closed-loop MPC")
    closed_cfg = replace(
        plant_cfg,
        output=replace(plant_cfg.output, run_dir=str(run_dir / "stages" / "closed_loop"), tag="closed_loop"),
    )
    closed_plant = ERDPlant(closed_cfg)

    u0 = default_open_loop_control(closed_cfg)

    def callback(x_k, u_prev, t_k, k):
        """Bridge ERDPlant callback signature to controller API."""

        _ = k
        return controller.select_u_from_state(x_k, u_prev if u_prev is not None else u0, t_k)

    closed_result = closed_plant.run(controller_callback=callback, write_artifacts=True)

    logger.info("Writing canonical comparison artifacts")
    _copy_model_artifacts(Path(cfg.model_run_dir), run_dir / "model")
    _plot_comparison(run_dir, open_result.curves, closed_result.curves)

    for stage in ("open_loop", "closed_loop"):
        src = run_dir / "stages" / stage / "movies" / "n.gif"
        if src.is_file():
            shutil.copy2(src, run_dir / "movies" / f"{stage}_n.gif")

    def _mean(curves: Dict, key: str) -> float:
        vals = curves.get(key, [])
        if not vals:
            return 0.0
        return float(sum(vals) / len(vals))

    def _final(curves: Dict, key: str) -> float:
        vals = curves.get(key, [])
        if not vals:
            return 0.0
        return float(vals[-1])

    metric_keys = [
        "E_wob",
        "L_w",
        "sigma_r",
        "P_ctrl",
        "E_wob_excess",
        "sigma_r_excess",
        "L_w_cum",
        "badness_score",
    ]

    summary = {
        "run_dir": str(run_dir),
        "open_loop_run": open_result.run_dir,
        "closed_loop_run": closed_result.run_dir,
        "dataset_id": model_bundle.dataset_id,
        "model_run_dir": str(Path(cfg.model_run_dir).resolve()),
    }
    for key in metric_keys:
        mean_open = _mean(open_result.curves, key)
        mean_closed = _mean(closed_result.curves, key)
        final_open = _final(open_result.curves, key)
        final_closed = _final(closed_result.curves, key)

        mean_denom = max(abs(mean_open), 1e-12)
        final_denom = max(abs(final_open), 1e-12)

        summary[f"mean_{key}_open"] = mean_open
        summary[f"mean_{key}_closed"] = mean_closed
        summary[f"relative_delta_{key}"] = (mean_closed - mean_open) / mean_denom
        summary[f"final_{key}_open"] = final_open
        summary[f"final_{key}_closed"] = final_closed
        summary[f"relative_final_delta_{key}"] = (final_closed - final_open) / final_denom

    with (run_dir / "metrics" / "relative_deltas.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    with (run_dir / "model" / "lineage.json").open("w", encoding="utf-8") as f:
        json.dump(
            {
                "model_run_dir": str(Path(cfg.model_run_dir).resolve()),
                "open_loop_run_dir": open_result.run_dir,
                "closed_loop_run_dir": closed_result.run_dir,
                "resolved_plant_config_source": cfg.plant_config_path or "model/training_manifest.yaml",
            },
            f,
            indent=2,
        )

    with (run_dir / "config.yaml").open("w", encoding="utf-8") as f:
        yaml.safe_dump({"control": to_dict(cfg), "plant": plant_cfg.to_dict()}, f, sort_keys=False)

    with (run_dir / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Closed-loop run complete: {run_dir}")
    return str(run_dir)
