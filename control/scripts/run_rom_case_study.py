"""Run a compact ROM-space baseline-vs-MPC case study and save artifacts."""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np

if __package__ in (None, ""):
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)

from erd_control.controller import RandomShootingMPCController
from erd_fipy import default_open_loop_control, load_run_config


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for the ROM case-study runner."""
    p = argparse.ArgumentParser(description="Run compact ROM-space MPC case study")
    p.add_argument("--model-run-dir", required=True)
    p.add_argument("--plant-config", required=True)
    p.add_argument("--source-run-dir", required=True)
    p.add_argument("--steps", type=int, default=220)
    p.add_argument(
        "--outputs-root",
        default="outputs/runs",
    )
    p.add_argument("--seed", type=int, default=9)
    return p.parse_args()


def make_out_dir(root: str) -> Path:
    """Create timestamped ROM case-study output directory."""
    out = Path(root).resolve() / f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_rom_case_study"
    for sub in ("raw", "controlled", "metrics", "plots"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    return out


def load_source_initial_state(source_run_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load initial stacked state and source timestamps from a generated run."""
    with h5py.File(source_run_dir / "fields" / "snapshots.h5", "r") as hf:
        n0 = np.asarray(hf["n"][0], dtype=float)
        w0 = np.asarray(hf["omega"][0], dtype=float)
        t = np.asarray(hf["t_snap"], dtype=float)
    return np.concatenate([n0.reshape(-1), w0.reshape(-1)]), n0, t


def summarize(base: dict[str, list[float]], ctrl: dict[str, list[float]]) -> dict[str, float]:
    """Compute aggregate before/after metrics and relative deltas."""
    out: dict[str, float] = {}
    for key in ("E_wob", "L_w", "sigma_r2", "P_ctrl"):
        b = float(np.mean(base[key]))
        c = float(np.mean(ctrl[key]))
        out[f"mean_{key}_baseline"] = b
        out[f"mean_{key}_controlled"] = c
        out[f"relative_delta_{key}"] = (c - b) / max(b, 1e-12)
    return out


def main() -> None:
    """Run baseline/control ROM rollouts, then write data and plots."""
    args = parse_args()
    out_dir = make_out_dir(args.outputs_root)

    model_run = Path(args.model_run_dir).resolve()
    source_run = Path(args.source_run_dir).resolve()
    plant_cfg = load_run_config(args.plant_config, preset="smoke")
    controller = RandomShootingMPCController(
        model_bundle=RandomShootingMPCController.load_model_bundle(model_run),
        plant_cfg=plant_cfg,
        seed=args.seed,
    )

    x0, n0, source_t = load_source_initial_state(source_run)
    a0 = controller.project(x0)
    n_steps = int(args.steps)
    dt = float(controller.model_bundle.dt)
    t = np.asarray([k * dt for k in range(n_steps + 1)], dtype=float)
    u_open = np.asarray(default_open_loop_control(plant_cfg), dtype=float)

    a_b = np.empty((n_steps + 1, a0.size), dtype=float)
    a_c = np.empty((n_steps + 1, a0.size), dtype=float)
    x_b = np.empty((n_steps + 1, x0.size), dtype=float)
    x_c = np.empty((n_steps + 1, x0.size), dtype=float)
    u_b = np.empty((n_steps, 5), dtype=float)
    u_c = np.empty((n_steps, 5), dtype=float)
    m_b = {"E_wob": [], "L_w": [], "sigma_r2": [], "P_ctrl": []}
    m_c = {"E_wob": [], "L_w": [], "sigma_r2": [], "P_ctrl": []}

    a_b[0] = a0
    a_c[0] = a0
    x_b[0] = controller.lift(a0)
    x_c[0] = controller.lift(a0)

    u_prev = u_open.copy()
    for k in range(n_steps):
        # Baseline
        u_b[k] = u_open
        a_b[k + 1] = controller._predict_step(a_b[k], u_open)
        x_b[k + 1] = controller.lift(a_b[k + 1])
        ew, lw, sr2, pc = controller._metrics_from_n(controller._decode_n(x_b[k + 1]), u_open)
        m_b["E_wob"].append(ew)
        m_b["L_w"].append(lw)
        m_b["sigma_r2"].append(sr2)
        m_b["P_ctrl"].append(pc)

        # Controlled
        u = controller.select_u(a_c[k], u_prev, float(t[k]))
        u_c[k] = u
        a_c[k + 1] = controller._predict_step(a_c[k], u)
        x_c[k + 1] = controller.lift(a_c[k + 1])
        ew, lw, sr2, pc = controller._metrics_from_n(controller._decode_n(x_c[k + 1]), u)
        m_c["E_wob"].append(ew)
        m_c["L_w"].append(lw)
        m_c["sigma_r2"].append(sr2)
        m_c["P_ctrl"].append(pc)
        u_prev = u

    np.savez_compressed(out_dir / "raw" / "source_slice.npz", t_source=source_t[: n_steps + 1], n0=n0)
    np.savez_compressed(out_dir / "raw" / "baseline_rom.npz", t=t, a=a_b, x=x_b, u=u_b)
    np.savez_compressed(out_dir / "controlled" / "mpc_rom.npz", t=t, a=a_c, x=x_c, u=u_c)

    summary = summarize(m_b, m_c)
    summary.update(
        {
            "model_run_dir": str(model_run),
            "source_run_dir": str(source_run),
            "output_dir": str(out_dir),
            "steps": n_steps,
            "dt": dt,
        }
    )
    (out_dir / "metrics" / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (out_dir / "metrics" / "baseline_curves.json").write_text(json.dumps({"t": t[1:].tolist(), **m_b}, indent=2))
    (out_dir / "metrics" / "controlled_curves.json").write_text(json.dumps({"t": t[1:].tolist(), **m_c}, indent=2))

    fig, ax = plt.subplots(2, 2, figsize=(10, 6))
    ax = ax.reshape(-1)
    for i, key in enumerate(("E_wob", "L_w", "sigma_r2", "P_ctrl")):
        ax[i].plot(t[1:], m_b[key], label="baseline")
        ax[i].plot(t[1:], m_c[key], label="controlled")
        ax[i].set_title(key)
        ax[i].grid(alpha=0.3)
        ax[i].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "plots" / "baseline_vs_controlled_metrics.png", dpi=150)
    plt.close(fig)

    print(f"output_dir: {out_dir}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
