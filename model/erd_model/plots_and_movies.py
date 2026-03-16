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


def _write_contact_sheet(
    data: np.ndarray,
    path: Path,
    title_prefix: str,
    dpi: int,
    cmap: str = "viridis",
) -> None:
    """Write a small contact sheet for a rollout field movie.

    Args:
        data: Field movie array with shape ``(T, N_r, N_phi)``.
        path: Output image path.
        title_prefix: Title stem used for each sampled frame.
        dpi: Figure export DPI.
        cmap: Matplotlib colormap name.

    Returns:
        None.
    """

    if data.shape[0] == 0:
        return
    ncols = min(6, data.shape[0])
    idx = np.linspace(0, data.shape[0] - 1, ncols, dtype=int)
    vmin = float(np.percentile(data, 2.0))
    vmax = float(np.percentile(data, 98.0))
    fig, axes = plt.subplots(1, ncols, figsize=(2.5 * ncols, 3.0), squeeze=False)
    for ax, k in zip(axes[0], idx, strict=True):
        im = ax.imshow(data[k], origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(f"{title_prefix} k={k}")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(im, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
    fig.tight_layout()
    _savefig(fig, path, dpi)


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

    _write_contact_sheet(rollout.fields_true[:, 0], plots / "rollout_n_true_contact.png", "n true", cfg.dpi)
    _write_contact_sheet(rollout.fields_pred[:, 0], plots / "rollout_n_pred_contact.png", "n pred", cfg.dpi)

    frames = []
    n_true = rollout.fields_true[:, 0]
    n_pred = rollout.fields_pred[:, 0]
    vmin = float(np.percentile(np.concatenate([n_true.reshape(-1), n_pred.reshape(-1)]), 2.0))
    vmax = float(np.percentile(np.concatenate([n_true.reshape(-1), n_pred.reshape(-1)]), 98.0))
    diff_lim = float(np.percentile(np.abs((n_pred - n_true).reshape(-1)), 98.0))
    for k in range(n_true.shape[0]):
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))
        axes[0].imshow(n_true[k], origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
        im1 = axes[1].imshow(n_pred[k], origin="lower", cmap="viridis", vmin=vmin, vmax=vmax)
        im2 = axes[2].imshow(
            n_pred[k] - n_true[k],
            origin="lower",
            cmap="coolwarm",
            vmin=-diff_lim,
            vmax=diff_lim,
        )
        axes[0].set_title("n true")
        axes[1].set_title("n pred")
        axes[2].set_title("pred - true")
        fig.colorbar(im1, ax=axes[:2].ravel().tolist(), fraction=0.03, pad=0.02)
        fig.colorbar(im2, ax=[axes[2]], fraction=0.046, pad=0.04)
        fig.tight_layout()

        fig.canvas.draw()
        frame = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
        frames.append(frame)
        plt.close(fig)

    imageio.mimsave(movies / "rollout_n_true_vs_pred.gif", frames, fps=cfg.movie_fps)
