"""Plot and GIF generation for ERD model-fit runs."""

from __future__ import annotations

from pathlib import Path

import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np

from .config import PlotConfig
from .metrics import EnergyCurves, ErrorCurves
from .pod import PODResult
from .rollout import RolloutResult
from .sindy_model import SINDyControlFitResult


def _savefig(fig, path: Path, dpi: int) -> None:
    """Save and close a Matplotlib figure."""

    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def generate_all_plots_and_movies(
    cfg: PlotConfig,
    run_dir: Path,
    pod: PODResult,
    sindy: SINDyControlFitResult,
    rollout: RolloutResult,
    err: ErrorCurves,
    energy: EnergyCurves,
) -> None:
    """Generate all standard visual artifacts for a model run.

    Args:
        cfg: Plot and movie settings.
        run_dir: Model run directory containing ``plots/`` and ``movies/``.
        pod: POD fit result.
        sindy: Controlled SINDy fit result.
        rollout: Representative rollout result.
        err: Error curves for representative rollout.
        energy: Energy curves for representative rollout.

    Returns:
        None.
    """

    plots = run_dir / "plots"
    movies = run_dir / "movies"

    s = pod.s
    cum = np.cumsum(s**2) / np.sum(s**2)

    fig = plt.figure()
    plt.semilogy(s, marker="o")
    plt.title("POD Singular Values")
    plt.xlabel("mode")
    _savefig(fig, plots / "pod_singular_values.png", cfg.dpi)

    fig = plt.figure()
    plt.plot(cum, marker="o")
    plt.title("POD Cumulative Energy")
    plt.xlabel("mode")
    _savefig(fig, plots / "pod_cumulative_energy.png", cfg.dpi)

    fig = plt.figure(figsize=(10, 4))
    plt.imshow(np.abs(sindy.coefficient_matrix), aspect="auto")
    plt.colorbar()
    plt.title("abs(Xi)")
    plt.xlabel("feature")
    plt.ylabel("equation")
    _savefig(fig, plots / "sindy_Xi_abs_heatmap.png", cfg.dpi)

    fig = plt.figure(figsize=(10, 4))
    for i in range(rollout.A_true.shape[1]):
        plt.plot(rollout.A_true[:, i], color=f"C{i}", label=f"a{i} true")
        plt.plot(rollout.A_pred[:, i], "--", color=f"C{i}", label=f"a{i} pred")
    plt.title("Coefficient Rollout")
    plt.xlabel("step")
    if rollout.A_true.shape[1] <= 8:
        plt.legend(ncol=2, fontsize=8)
    _savefig(fig, plots / "coeff_time_series_true_vs_pred.png", cfg.dpi)

    fig = plt.figure(figsize=(10, 4))
    plt.plot(err.field_rel_l2, label="field_rel_l2")
    plt.plot(err.coeff_mse, label="coeff_mse")
    plt.title("Validation Curves")
    plt.xlabel("step")
    plt.legend()
    _savefig(fig, plots / "metrics_validation_curves.png", cfg.dpi)

    fig = plt.figure(figsize=(10, 4))
    plt.plot(energy.energy_true, label="energy_true")
    plt.plot(energy.energy_pred, label="energy_pred")
    plt.title("Energy True vs Pred")
    plt.xlabel("step")
    plt.legend()
    _savefig(fig, plots / "metrics_energy_true_vs_pred.png", cfg.dpi)

    frames = []
    n_true = rollout.fields_true[:, 0]
    n_pred = rollout.fields_pred[:, 0]
    for k in range(n_true.shape[0]):
        fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
        vmin = float(min(np.min(n_true[k]), np.min(n_pred[k])))
        vmax = float(max(np.max(n_true[k]), np.max(n_pred[k])))
        axes[0].imshow(n_true[k], origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
        im1 = axes[1].imshow(n_pred[k], origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
        axes[0].set_title("n true")
        axes[1].set_title("n pred")
        fig.colorbar(im1, ax=axes.ravel().tolist(), fraction=0.03)
        fig.tight_layout()

        fig.canvas.draw()
        frame = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
        frames.append(frame)
        plt.close(fig)

    imageio.mimsave(movies / "rollout_n_true_vs_pred.gif", frames, fps=cfg.movie_fps)
