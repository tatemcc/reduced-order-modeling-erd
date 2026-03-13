"""Analyze variation diagnostics for an existing ERD run folder."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
import yaml

if __package__ in (None, ""):
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)

from erd_fipy.config import run_config_from_dict
from erd_fipy.metrics import compute_variation_diagnostics


def build_parser() -> argparse.ArgumentParser:
    """Build command-line parser for variation analysis.

    Returns:
        Configured argument parser.
    """

    parser = argparse.ArgumentParser(description="Compute ERD variation diagnostics from saved snapshots")
    parser.add_argument("--run-dir", type=str, required=True, help="Run directory containing fields/snapshots.h5")
    parser.add_argument(
        "--baseline-run",
        type=str,
        default=None,
        help="Optional baseline run for delta-median comparison",
    )
    return parser


def _load_run(run_dir: Path) -> tuple[object, np.ndarray, np.ndarray, np.ndarray]:
    """Load run config and density snapshots from one run folder.

    Args:
        run_dir: Run directory path.

    Returns:
        Tuple ``(cfg, n_snapshots, r, t_snap)``.
    """

    with (run_dir / "config.yaml").open("r", encoding="utf-8") as f:
        cfg_raw = yaml.safe_load(f) or {}
    cfg = run_config_from_dict(cfg_raw)

    with h5py.File(run_dir / "fields" / "snapshots.h5", "r") as hf:
        n = np.asarray(hf["n"])
        r = np.asarray(hf["r"])
        t = np.asarray(hf["t_snap"])
    return cfg, n, r, t


def _plot(run_dir: Path, payload: dict) -> None:
    """Render the variation diagnostics figure.

    Args:
        run_dir: Run directory path.
        payload: Diagnostics payload dictionary.

    Returns:
        None.
    """

    diag = payload["diagnostics"]
    acceptance = payload["acceptance"]
    baseline = payload.get("baseline", {})
    summary = diag.get("summary", {})

    t = np.asarray(diag.get("t", []), dtype=float)
    ratio = np.asarray(diag.get("ratio_nonax_over_axisym", []), dtype=float)
    mid = np.asarray(diag.get("band_ratio_mid_to_low", []), dtype=float)
    high = np.asarray(diag.get("band_ratio_high_to_low", []), dtype=float)
    td = np.asarray(diag.get("delta_frame_l2_rel_t", []), dtype=float)
    delta = np.asarray(diag.get("delta_frame_l2_rel", []), dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    axes = axes.reshape(-1)

    axes[0].plot(t, ratio, label="E_nonax / E0")
    axes[0].axhline(0.2, color="r", linestyle="--", label="target 0.2")
    axes[0].set_title("Non-axisymmetry Ratio")
    axes[0].legend()

    axes[1].plot(t, mid, label="mid / low")
    axes[1].plot(t, high, label="high / low")
    axes[1].axhline(0.05, color="r", linestyle="--", label="target 0.05")
    axes[1].set_title("Spectral Band Ratios")
    axes[1].legend()

    axes[2].plot(td, delta, label="frame Δ")
    if "baseline_delta_median" in baseline:
        axes[2].axhline(5.0 * float(baseline["baseline_delta_median"]), color="r", linestyle="--", label="5x baseline")
    axes[2].set_title("Temporal Activity")
    axes[2].legend()

    axes[3].axis("off")
    lines = [
        f"criterion1_nonaxis_growth: {acceptance['criterion1_nonaxis_growth']}",
        f"criterion2_mid_band: {acceptance['criterion2_mid_band']}",
        f"criterion2_high_band: {acceptance['criterion2_high_band']}",
        f"criterion3_temporal_activity: {acceptance['criterion3_temporal_activity']}",
        f"all_passed: {acceptance['all_passed']}",
    ]
    if "delta_median_vs_baseline_factor" in summary:
        lines.append(f"delta_vs_baseline_factor: {summary['delta_median_vs_baseline_factor']:.3f}")
    axes[3].text(0.0, 1.0, "\n".join(lines), va="top", family="monospace")

    for ax in axes[:3]:
        ax.grid(alpha=0.3)
        ax.set_xlabel("t")

    fig.tight_layout()
    fig.savefig(run_dir / "plots" / "variation_diagnostics.png", dpi=150)
    plt.close(fig)


def main() -> None:
    """Compute and save variation diagnostics for one run.

    Returns:
        None.
    """

    args = build_parser().parse_args()
    run_dir = Path(args.run_dir).resolve()
    (run_dir / "metrics").mkdir(parents=True, exist_ok=True)
    (run_dir / "plots").mkdir(parents=True, exist_ok=True)
    (run_dir / "notes").mkdir(parents=True, exist_ok=True)

    cfg, n, r, t = _load_run(run_dir)
    diag = compute_variation_diagnostics(cfg=cfg, n_snapshots=n, r=r, t_snap=t)

    baseline_info: dict[str, object] = {}
    if args.baseline_run:
        _, n_base, r_base, t_base = _load_run(Path(args.baseline_run).resolve())
        base_diag = compute_variation_diagnostics(cfg=cfg, n_snapshots=n_base, r=r_base, t_snap=t_base)
        base_med = float(base_diag["summary"].get("delta_frame_l2_rel_median", 0.0))
        cur_med = float(diag["summary"].get("delta_frame_l2_rel_median", 0.0))
        diag["summary"]["delta_median_vs_baseline_factor"] = cur_med / max(base_med, 1e-14)
        baseline_info = {
            "baseline_run_dir": str(Path(args.baseline_run).resolve()),
            "baseline_delta_median": base_med,
            "delta_median_vs_baseline_factor": diag["summary"]["delta_median_vs_baseline_factor"],
        }

    summary = diag.get("summary", {})
    acceptance = {
        "criterion1_nonaxis_growth": bool(
            summary.get("ratio_nonax_over_axisym_first20_max", 0.0) >= 0.2
            and summary.get("ratio_nonax_over_axisym_first20_std", 0.0) > 1e-4
        ),
        "criterion2_mid_band": bool(summary.get("band_ratio_mid_to_low_max", 0.0) >= 0.05),
        "criterion2_high_band": bool(summary.get("band_ratio_high_to_low_max", 0.0) >= 0.05),
        "criterion3_temporal_activity": bool(summary.get("delta_median_vs_baseline_factor", 0.0) >= 5.0),
    }
    acceptance["all_passed"] = bool(all(acceptance.values()))

    payload = {"diagnostics": diag, "baseline": baseline_info, "acceptance": acceptance}
    with (run_dir / "metrics" / "variation_diagnostics.json").open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    _plot(run_dir, payload)

    note = [
        "# Variation Analysis",
        "",
        f"- run: `{run_dir}`",
        f"- nonaxis first20 max: {summary.get('ratio_nonax_over_axisym_first20_max', 0.0):.4f}",
        f"- mid/low max: {summary.get('band_ratio_mid_to_low_max', 0.0):.4f}",
        f"- high/low max: {summary.get('band_ratio_high_to_low_max', 0.0):.4f}",
        f"- delta median: {summary.get('delta_frame_l2_rel_median', 0.0):.6f}",
    ]
    if "delta_median_vs_baseline_factor" in summary:
        note.append(f"- delta vs baseline factor: {summary['delta_median_vs_baseline_factor']:.3f}")
    note.append(f"- acceptance all_passed: {acceptance['all_passed']}")
    (run_dir / "notes" / "variation_analysis.md").write_text("\n".join(note) + "\n", encoding="utf-8")

    print(f"run_dir: {run_dir}")
    print(f"variation_json: {run_dir / 'metrics' / 'variation_diagnostics.json'}")
    print(f"variation_plot: {run_dir / 'plots' / 'variation_diagnostics.png'}")
    print(f"all_passed: {acceptance['all_passed']}")


if __name__ == "__main__":
    main()
