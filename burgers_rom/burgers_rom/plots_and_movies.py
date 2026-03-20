"""
Plot and movie generation for Burgers ROM + SINDy runs.

This module is called from run_pipeline after a run completes.
It writes:
- POD singular value and cumulative energy plots
- POD basis modes in field space for each component
- SINDy coefficient matrix heatmap
- Coefficient time series plots
- Coefficient phase portraits (a_i(t), a_j(t))
- Rollout true vs predicted figures and GIF movies per component
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio

from .config import PlotConfig
from .pod import PODResult
from .rollout import RolloutResult
from .snapshot import SnapshotLayout, state_vec_to_fields
from .sindy_model import SINDyFitResult


def _safe_sympy_symbols(r: int, style: str) -> Optional[List[str]]:
    """
    Build sympy-style symbol strings for coefficient labels.

    Parameters
    ----------
    r : int
        Number of coefficients.
    style : str
        'a_i' or 'x_i'.

    Returns
    -------
    list[str] or None
        Label strings, or None if sympy is unavailable.
    """
    try:
        import sympy as sp
    except Exception:
        return None

    prefix = "a" if style == "a_i" else "x"
    syms = sp.symbols(f"{prefix}0:{r}")
    return [str(s) for s in syms]


def _savefig(fig, path: Path, dpi: int) -> None:
    """
    Save a matplotlib figure.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure handle.
    path : Path
        Output path.
    dpi : int
        DPI for saving.
    """
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_pod_singular_values(
    pod: PODResult,
    out_dir: Path,
    dpi: int,
) -> None:
    """
    Plot singular value decay and cumulative energy.

    Parameters
    ----------
    pod : PODResult
        POD results.
    out_dir : Path
        Output directory for figures.
    dpi : int
        DPI for saving.
    """
    s = np.asarray(pod.s)
    energy = s**2
    cum = np.cumsum(energy) / np.sum(energy)

    fig = plt.figure()
    plt.semilogy(s, marker="o")
    plt.xlabel("mode index")
    plt.ylabel("singular value")
    plt.title("POD singular values")
    _savefig(fig, out_dir / "pod_singular_values.png", dpi=dpi)

    fig = plt.figure()
    plt.plot(cum, marker="o", markersize=4)
    plt.axhline(0.99, linestyle="--")
    plt.xlabel("mode index")
    plt.ylabel("cumulative energy fraction")
    plt.title("POD cumulative energy")
    _savefig(fig, out_dir / "pod_cumulative_energy.png", dpi=dpi)


def plot_pod_basis_fields(
    U: np.ndarray,
    layout: SnapshotLayout,
    out_dir: Path,
    n_modes: int,
    cmap: str,
    dpi: int,
) -> None:
    """
    Plot POD basis modes in field space for each component.

    Parameters
    ----------
    U : np.ndarray
        POD basis, shape (n_state, r).
    layout : SnapshotLayout
        Layout for reshaping basis vectors into (C, ny, nx).
    out_dir : Path
        Output directory for figures.
    n_modes : int
        Number of modes to plot.
    cmap : str
        Matplotlib colormap name.
    dpi : int
        DPI for saving.
    """
    r = U.shape[1]
    # plot all modes if n_modes is greater than available modes
    n_plot = min(n_modes, r)

    for k in range(n_plot):
        fields = state_vec_to_fields(U[:, k], layout)

        vmin = float(np.min(fields))
        vmax = float(np.max(fields))
        vmax_abs = max(abs(vmin), abs(vmax))
        vmin, vmax = -vmax_abs, vmax_abs

        fig, axes = plt.subplots(1, layout.n_components, figsize=(5.0 * layout.n_components, 4.0), dpi=dpi)
        if layout.n_components == 1:
            axes = [axes]

        for c in range(layout.n_components):
            ax = axes[c]
            im = ax.imshow(fields[c], origin="lower", aspect="equal", cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_title(f"mode {k} component {c}")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        fig.suptitle(f"POD basis mode {k}")
        fig.tight_layout()
        _savefig(fig, out_dir / f"pod_basis_mode_{k:03d}.png", dpi=dpi)


def plot_sindy_coeff_matrix(
    sindy: SINDyFitResult,
    out_dir: Path,
    dpi: int,
) -> None:
    """
    Plot absolute value heatmap of the SINDy coefficient matrix.

    Parameters
    ----------
    sindy : SINDyFitResult
        SINDy fit results.
    out_dir : Path
        Output directory.
    dpi : int
        DPI for saving.
    """
    Xi = np.asarray(sindy.coefficient_matrix)

    fig = plt.figure(figsize=(10, 4))
    plt.imshow(np.abs(Xi), aspect="auto")
    plt.xlabel("feature index")
    plt.ylabel("equation index")
    plt.title("abs(Xi) heatmap")
    plt.colorbar()
    _savefig(fig, out_dir / "sindy_Xi_abs_heatmap.png", dpi=dpi)


def plot_coeff_time_series(
    A_true: np.ndarray,
    A_pred: np.ndarray,
    out_dir: Path,
    dpi: int,
    sympy_labels: bool,
    sympy_style: str,
) -> None:
    """
    Plot coefficient time series for true and predicted trajectories.

    Parameters
    ----------
    A_true : np.ndarray
        True coefficients, shape (T, r).
    A_pred : np.ndarray
        Predicted coefficients, shape (T, r).
    out_dir : Path
        Output directory.
    dpi : int
        DPI for saving.
    sympy_labels : bool
        If True, attempt sympy-style labels.
    sympy_style : str
        Label style: 'a_i' or 'x_i'.
    """
    A_true = np.asarray(A_true)
    A_pred = np.asarray(A_pred)
    T, r = A_true.shape

    labels = None
    if sympy_labels:
        labels = _safe_sympy_symbols(r, sympy_style)

    fig = plt.figure(figsize=(10, 5))
    for i in range(r):
        name = labels[i] if labels is not None else f"a{i}"
        # sync the color across true and pred for the same coefficient
        plt.plot(A_true[:, i], linewidth=1.0, label=f"{name} true", color=f"C{i}")
        plt.plot(A_pred[:, i], linewidth=1.0, linestyle="--", label=f"{name} pred", color=f"C{i}")
    plt.xlabel("t index")
    plt.ylabel("coefficient value")
    plt.title("POD coefficients: true vs predicted")
    if r <= 10:
        plt.legend(ncol=2, fontsize=8)
    _savefig(fig, out_dir / "coeff_time_series_true_vs_pred.png", dpi=dpi)


def plot_coeff_phase_portraits(
    A: np.ndarray,
    out_dir: Path,
    dpi: int,
    max_pairs: int,
    sympy_labels: bool,
    sympy_style: str,
) -> None:
    """
    Plot phase portraits (a_i(t), a_j(t)) for coefficient trajectories.

    Parameters
    ----------
    A : np.ndarray
        Coefficients, shape (T, r).
    out_dir : Path
        Output directory.
    dpi : int
        DPI for saving.
    max_pairs : int
        Maximum number of pairs to plot.
    sympy_labels : bool
        If True, attempt sympy-style labels.
    sympy_style : str
        Label style: 'a_i' or 'x_i'.
    """
    A = np.asarray(A)
    T, r = A.shape

    labels = None
    if sympy_labels:
        labels = _safe_sympy_symbols(r, sympy_style)

    pairs: List[Tuple[int, int]] = []
    for i in range(r):
        for j in range(i + 1, r):
            pairs.append((i, j))
            if len(pairs) >= max_pairs:
                break
        if len(pairs) >= max_pairs:
            break

    for i, j in pairs:
        xi = A[:, i]
        xj = A[:, j]
        li = labels[i] if labels is not None else f"a{i}"
        lj = labels[j] if labels is not None else f"a{j}"

        fig = plt.figure(figsize=(4.5, 4.5))
        plt.plot(xi, xj, linewidth=1.0)
        plt.xlabel(li)
        plt.ylabel(lj)
        plt.title(f"phase portrait: {li} vs {lj}")
        _savefig(fig, out_dir / f"phase_{i:02d}_{j:02d}.png", dpi=dpi)


def _frame_true_pred(
    true_field: np.ndarray,
    pred_field: np.ndarray,
    title_left: str,
    title_right: str,
    cmap: str,
    dpi: int,
) -> np.ndarray:
    """
    Render a side-by-side true vs predicted heatmap frame.

    Parameters
    ----------
    true_field : np.ndarray
        True field, shape (ny, nx).
    pred_field : np.ndarray
        Predicted field, shape (ny, nx).
    title_left : str
        Left title.
    title_right : str
        Right title.
    cmap : str
        Colormap.
    dpi : int
        DPI.

    Returns
    -------
    np.ndarray
        RGB image array.
    """
    vmin = float(min(true_field.min(), pred_field.min()))
    vmax = float(max(true_field.max(), pred_field.max()))

    fig, axes = plt.subplots(1, 2, figsize=(8, 3), dpi=dpi)
    axes[0].imshow(true_field, origin="lower", aspect="equal", vmin=vmin, vmax=vmax, cmap=cmap)
    axes[0].set_title(title_left)
    axes[0].axis("off")

    axes[1].imshow(pred_field, origin="lower", aspect="equal", vmin=vmin, vmax=vmax, cmap=cmap)
    axes[1].set_title(title_right)
    axes[1].axis("off")

    fig.tight_layout()
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    frame = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(h, w, 4)[:, :, :3]
    plt.close(fig)
    return frame


def save_rollout_movies(
    rollout: RolloutResult,
    out_dir: Path,
    fps: int,
    every: int,
    cmap: str,
    dpi: int,
    components: Optional[List[int]],
) -> None:
    """
    Save GIF movies for rollout true vs predicted fields.

    Parameters
    ----------
    rollout : RolloutResult
        Rollout outputs including fields_true and fields_pred.
    out_dir : Path
        Output directory.
    fps : int
        Frames per second for GIF.
    every : int
        Take every k-th frame.
    cmap : str
        Colormap.
    dpi : int
        DPI for rendering frames.
    components : list[int], optional
        Components to render. If None, render all.
    """
    ft = rollout.fields_true
    fp = rollout.fields_pred

    T, C, _, _ = ft.shape
    comp_list = list(range(C)) if components is None else components

    duration = 1.0 / max(1, fps)

    for c in comp_list:
        frames = []
        for t_idx in range(0, T, max(1, every)):
            frames.append(
                _frame_true_pred(
                    ft[t_idx, c],
                    fp[t_idx, c],
                    title_left=f"true comp {c} t={t_idx}",
                    title_right=f"pred comp {c} t={t_idx}",
                    cmap=cmap,
                    dpi=dpi,
                )
            )
        gif_path = out_dir / f"rollout_true_vs_pred_comp{c}.gif"
        imageio.mimsave(gif_path, frames, duration=duration, loop=0)


def save_rollout_snapshot_figures(
    rollout: RolloutResult,
    out_dir: Path,
    cmap: str,
    dpi: int,
    components: Optional[List[int]],
) -> None:
    """
    Save a small set of true vs predicted snapshots.

    Parameters
    ----------
    rollout : RolloutResult
        Rollout outputs.
    out_dir : Path
        Output directory.
    cmap : str
        Colormap.
    dpi : int
        DPI for saving.
    components : list[int], optional
        Components to render. If None, render all.
    """
    ft = rollout.fields_true
    fp = rollout.fields_pred
    T, C, _, _ = ft.shape

    comp_list = list(range(C)) if components is None else components
    times = [0, T // 3, (2 * T) // 3, T - 1]

    for c in comp_list:
        for t_idx in times:
            true_field = ft[t_idx, c]
            pred_field = fp[t_idx, c]

            vmin = float(min(true_field.min(), pred_field.min()))
            vmax = float(max(true_field.max(), pred_field.max()))

            fig, axes = plt.subplots(1, 2, figsize=(8, 3), dpi=dpi)
            im0 = axes[0].imshow(true_field, origin="lower", aspect="equal", vmin=vmin, vmax=vmax, cmap=cmap)
            axes[0].set_title(f"true comp {c} t={t_idx}")
            axes[0].axis("off")

            im1 = axes[1].imshow(pred_field, origin="lower", aspect="equal", vmin=vmin, vmax=vmax, cmap=cmap)
            axes[1].set_title(f"pred comp {c} t={t_idx}")
            axes[1].axis("off")

            fig.colorbar(im1, ax=axes, fraction=0.046, pad=0.04)
            # fig.tight_layout() # seems to break for this one due to colourbar
            _savefig(fig, out_dir / f"rollout_snapshot_comp{c}_t{t_idx:04d}.png", dpi=dpi)


def plot_metrics_curves_from_artifacts(
    rundir: Path,
    out_dir: Path,
    dpi: int,
) -> None:
    """
    Plot key rollout curves from metrics/curves.json and save to figures/.

    Curves:
    - Coefficient MSE: MSE(t) = (1/r) sum_i (a_i(t) - ahat_i(t))^2
    - Relative field L2 error: ||u_hat - u||_2 / ||u||_2
    - Energy: E(t) = 0.5 * integral (u^2 + v^2) dA
    - Energy drift: E_hat(t) - E(t)

    Parameters
    ----------
    rundir : Path
        Run directory containing metrics/curves.json.
    out_dir : Path
        Figures output directory.
    dpi : int
        DPI for saving.
    """
    import json

    curves_path = rundir / "metrics" / "curves.json"
    if not curves_path.exists():
        return

    with curves_path.open("r", encoding="utf-8") as f:
        curves = json.load(f)

    coeff_mse = np.asarray(curves.get("coeff_mse", []), dtype=float)
    field_rel = np.asarray(curves.get("field_rel_l2", []), dtype=float)
    energy_true = np.asarray(curves.get("energy_true", []), dtype=float)
    energy_pred = np.asarray(curves.get("energy_pred", []), dtype=float)

    t = np.arange(len(coeff_mse), dtype=int)

    if coeff_mse.size:
        fig = plt.figure(figsize=(7.5, 4.0), dpi=dpi)
        plt.plot(t, coeff_mse, linewidth=2.0, label="MSE")
        plt.grid(True, alpha=0.3)
        plt.xlabel("time index")
        plt.ylabel("MSE")
        plt.title(r"Coefficient error: $\mathrm{MSE}(t)=\frac{1}{r}\sum_{i=0}^{r-1}(a_i-\hat a_i)^2$")
        plt.legend()
        _savefig(fig, out_dir / "metrics_coeff_mse_vs_time.png", dpi=dpi)

    t2 = np.arange(len(field_rel), dtype=int)
    if field_rel.size:
        fig = plt.figure(figsize=(7.5, 4.0), dpi=dpi)
        plt.plot(t2, field_rel, linewidth=2.0, label="relative L2")
        plt.grid(True, alpha=0.3)
        plt.xlabel("time index")
        plt.ylabel("relative L2 error")
        plt.title(r"Field error: $\frac{\|\hat u-u\|_2}{\|u\|_2}$ (stacked components)")
        plt.legend()
        _savefig(fig, out_dir / "metrics_field_rel_l2_vs_time.png", dpi=dpi)

    t3 = np.arange(len(energy_true), dtype=int)
    if energy_true.size and energy_pred.size and energy_true.shape == energy_pred.shape:
        fig = plt.figure(figsize=(7.5, 4.0), dpi=dpi)
        plt.plot(t3, energy_true, linewidth=2.0, label="true")
        plt.plot(t3, energy_pred, linewidth=2.0, linestyle="--", label="pred")
        plt.grid(True, alpha=0.3)
        plt.xlabel("time index")
        plt.ylabel("energy")
        plt.title(r"Energy: $E(t)=\frac{1}{2}\int_{\Omega}(u^2+v^2)\,dA$")
        plt.legend()
        _savefig(fig, out_dir / "metrics_energy_true_vs_pred.png", dpi=dpi)

        fig = plt.figure(figsize=(7.5, 4.0), dpi=dpi)
        plt.plot(t3, energy_pred - energy_true, linewidth=2.0, label=r"$\hat E - E$")
        plt.grid(True, alpha=0.3)
        plt.xlabel("time index")
        plt.ylabel("energy drift")
        plt.title(r"Energy drift: $\Delta E(t)=\hat E(t)-E(t)$")
        plt.legend()
        _savefig(fig, out_dir / "metrics_energy_drift.png", dpi=dpi)

def generate_all_plots_and_movies(
    cfg: PlotConfig,
    rundir: Path,
    layout: SnapshotLayout,
    pod: PODResult,
    sindy: SINDyFitResult,
    rollout: RolloutResult,
) -> None:
    """
    Generate all configured plots and movies for a completed run.

    Parameters
    ----------
    cfg : PlotConfig
        Plot configuration.
    rundir : Path
        Run output directory.
    layout : SnapshotLayout
        Snapshot layout for reshaping basis modes.
    pod : PODResult
        POD results.
    sindy : SINDyFitResult
        SINDy results.
    rollout : RolloutResult
        Rollout results.
    """
    if not cfg.enabled:
        return

    fig_dir = rundir / cfg.figures_subdir
    mov_dir = rundir / cfg.movies_subdir
    fig_dir.mkdir(parents=True, exist_ok=True)
    mov_dir.mkdir(parents=True, exist_ok=True)

    (fig_dir / "pod_basis").mkdir(parents=True, exist_ok=True)
    plot_pod_singular_values(pod=pod, out_dir=fig_dir, dpi=cfg.dpi)
    plot_pod_basis_fields(
        U=pod.U,
        layout=layout,
        out_dir=fig_dir / "pod_basis",
        n_modes=cfg.basis_n_modes,
        cmap=cfg.basis_cmap,
        dpi=cfg.dpi,
    )

    plot_sindy_coeff_matrix(sindy=sindy, out_dir=fig_dir, dpi=cfg.dpi)

    if cfg.coeff_time_series:
        plot_coeff_time_series(
            A_true=rollout.A_true,
            A_pred=rollout.A_pred,
            out_dir=fig_dir,
            dpi=cfg.dpi,
            sympy_labels=cfg.sympy_labels,
            sympy_style=cfg.sympy_label_style,
        )

    if cfg.coeff_pair_phase:
        (fig_dir / "phase_portraits").mkdir(parents=True, exist_ok=True)
        plot_coeff_phase_portraits(
            A=rollout.A_true,
            out_dir=fig_dir / "phase_portraits",
            dpi=cfg.dpi,
            max_pairs=cfg.coeff_pair_max_pairs,
            sympy_labels=cfg.sympy_labels,
            sympy_style=cfg.sympy_label_style,
        )

    comps = cfg.movie_components
    (fig_dir / "rollout_snapshots").mkdir(parents=True, exist_ok=True)
    save_rollout_snapshot_figures(
        rollout=rollout,
        out_dir=fig_dir / "rollout_snapshots",
        cmap=cfg.rollout_cmap,
        dpi=cfg.dpi,
        components=comps,
    )

    save_rollout_movies(
        rollout=rollout,
        out_dir=mov_dir,
        fps=cfg.movie_fps,
        every=cfg.movie_every,
        cmap=cfg.rollout_cmap,
        dpi=cfg.dpi,
        components=comps,
    )

    if getattr(cfg, "metrics_curves", True):
        plot_metrics_curves_from_artifacts(
            rundir=rundir,
            out_dir=fig_dir,
            dpi=cfg.dpi,
        )