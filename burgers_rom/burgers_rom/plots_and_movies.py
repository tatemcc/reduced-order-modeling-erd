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
import matplotlib.animation as animation
import imageio.v2 as imageio
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import RectBivariateSpline
from matplotlib.colors import LightSource, hsv_to_rgb, Normalize

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

def compute_vorticity(q: np.ndarray, L: float = 1.0) -> np.ndarray:
    """
    Compute scalar vorticity field from velocity field.

    omega = dv/dx - du/dy

    Parameters
    ----------
    q : np.ndarray
        Velocity field of shape (T, 2, ny, nx).
        q[:, 0] is u (x-velocity), q[:, 1] is v (y-velocity).
    L : float
        Domain size length (assumed square [0, L]x[0, L]).

    Returns
    -------
    omega : np.ndarray
        Vorticity field of shape (T, ny, nx).
    """
    if q.ndim != 4 or q.shape[1] != 2:
        raise ValueError("Vorticity requires velocity field of shape (T, 2, ny, nx)")

    T, C, ny, nx = q.shape
    dx = L / (nx - 1) if nx > 1 else 1.0
    dy = L / (ny - 1) if ny > 1 else 1.0

    u = q[:, 0]
    v = q[:, 1]

    # du/dy: gradient of u along axis 1 (y-axis)
    dudy = np.gradient(u, axis=1) / dy

    # dv/dx: gradient of v along axis 2 (x-axis)
    dvdx = np.gradient(v, axis=2) / dx

    omega = dvdx - dudy
    return omega

def compute_gasdynamics_vars(q: np.ndarray, gamma: float = 1.4) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute pressure and velocity magnitude from conserved gas dynamics variables.

    Assumes q = [rho, rho*u_x, rho*u_y, E]

    Parameters
    ----------
    q : np.ndarray
        Conserved variables of shape (T, 4, ny, nx).
    gamma : float
        Adiabatic index.

    Returns
    -------
    pressure : np.ndarray
        Pressure field of shape (T, ny, nx).
    velocity_mag : np.ndarray
        Velocity magnitude field of shape (T, ny, nx).
    """
    if q.ndim != 4 or q.shape[1] != 4:
        raise ValueError("Gas dynamics requires field of shape (T, 4, ny, nx)")

    rho = q[:, 0]
    rho_ux = q[:, 1]
    rho_uy = q[:, 2]
    E = q[:, 3]

    # Add a small epsilon to rho to avoid division by zero in low-density regions
    rho_stable = rho + 1e-9

    # Kinetic energy: 0.5 * rho * (ux^2 + uy^2) = 0.5 * ( (rho*ux)^2 + (rho*uy)^2 ) / rho
    kinetic_energy = 0.5 * (rho_ux**2 + rho_uy**2) / rho_stable
    
    # Pressure: p = (gamma - 1) * (E - KE)
    pressure = (gamma - 1) * (E - kinetic_energy)

    # Velocity magnitude: |u| = sqrt(ux^2 + uy^2) = sqrt( (rho*ux/rho)^2 + (rho*uy/rho)^2 )
    ux = rho_ux / rho_stable
    uy = rho_uy / rho_stable
    velocity_mag = np.sqrt(ux**2 + uy**2)

    return pressure, velocity_mag

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

def plot_mse_comparison_with_paper(
    rundir: Path,
    out_dir: Path,
    dpi: int,
    equation: str,
) -> None:
    """
    Plot coefficient MSE vs. time and compare with data from Arora et al. 2023.

    This function reproduces Figure 3 from https://arxiv.org/abs/2306.05805
    for the relevant equation and overlays the current run's results.

    Parameters
    ----------
    rundir : Path
        Run directory containing metrics/curves.json.
    out_dir : Path
        Figures output directory.
    dpi : int
        DPI for saving.
    equation : str
        Name of the equation, e.g., 'burgers'.
    """
    import json

    # Data digitized from Figure 3 of Arora et al., "Trapping SINDy...",
    # arXiv:2306.05805.
    paper_data = {
        "kuramotosivashinsky": {
            "SINDy (paper)": ([1, 10, 20, 40, 60, 80, 100], [2e-6, 1e-3, 1e-1, 1e1, 1e2, 1e3, 1e4]),
            "Trapping SINDy (paper)": ([1, 10, 20, 40, 60, 80, 100], [2e-6, 8e-5, 1e-4, 1.2e-4, 1.5e-4, 1.8e-4, 2e-4]),
            "Galerkin (paper)": ([1, 10, 20, 40, 60, 80, 100], [2e-6, 2e-3, 2e-1, 2e1, 2e2, 2e3, 2e4]),
            "ML-Galerkin (paper)": ([1, 10, 20, 40, 60, 80, 100], [2e-6, 8e-6, 1e-5, 1.2e-5, 1.5e-5, 1.8e-5, 2e-5]),
        },
        "burgers": {
            "SINDy (paper)": ([1, 10, 20, 30, 40, 50, 60], [1e-5, 1e-3, 1e-2, 1e-1, 1, 10, 100]),
            "Trapping SINDy (paper)": ([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], [1e-5, 2e-4, 4e-4, 6e-4, 1e-3, 2e-3, 4e-3, 8e-3, 2e-2, 4e-2, 8e-2]),
            "Galerkin (paper)": ([1, 10, 20, 30, 40, 50, 60], [2e-5, 2e-3, 2e-2, 0.2, 2, 20, 200]),
            "ML-Galerkin (paper)": ([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], [1e-5, 1e-4, 2e-4, 3e-4, 4e-4, 6e-4, 1e-3, 2e-3, 4e-3, 8e-3, 2e-2]),
        },
        "reactiondiffusion": {
            "SINDy (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 1e-4, 1e-3, 1e-1, 1, 5, 10]),
            "Trapping SINDy (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 8e-5, 1e-4, 1.2e-4, 1.5e-4, 1.8e-4, 2e-4]),
            "Galerkin (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 1e-3, 1e-2, 1, 10, 50, 100]),
            "ML-Galerkin (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 8e-6, 1e-5, 1.2e-5, 1.5e-5, 1.8e-5, 2e-5]),
        },
        "gasdynamics": {
            "SINDy (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 1e-3, 1e-1, 1e1, 1e2, 5e2, 1e3]),
            "Trapping SINDy (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 8e-5, 1e-4, 1.2e-4, 1.5e-4, 1.8e-4, 2e-4]),
            "Galerkin (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 1e-2, 1, 1e2, 1e3, 5e3, 1e4]),
            "ML-Galerkin (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 8e-6, 1e-5, 1.2e-5, 1.5e-5, 1.8e-5, 2e-5]),
        },
        "advection": {
            "SINDy (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 8e-6, 1e-5, 1.2e-5, 1.5e-5, 1.8e-5, 2e-5]),
            "Trapping SINDy (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 1.1e-6, 1.2e-6, 1.3e-6, 1.4e-6, 1.5e-6, 1.6e-6]),
            "Galerkin (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 9e-6, 1.1e-5, 1.3e-5, 1.6e-5, 1.9e-5, 2.1e-5]),
            "ML-Galerkin (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 1.2e-6, 1.3e-6, 1.4e-6, 1.5e-6, 1.6e-6, 1.7e-6]),
        },
        "wave": {
            "SINDy (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 8e-6, 1e-5, 1.2e-5, 1.5e-5, 1.8e-5, 2e-5]),
            "Trapping SINDy (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 1.1e-6, 1.2e-6, 1.3e-6, 1.4e-6, 1.5e-6, 1.6e-6]),
            "Galerkin (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 9e-6, 1.1e-5, 1.3e-5, 1.6e-5, 1.9e-5, 2.1e-5]),
            "ML-Galerkin (paper)": ([1, 10, 20, 40, 60, 80, 100], [1e-6, 1.2e-6, 1.3e-6, 1.4e-6, 1.5e-6, 1.6e-6, 1.7e-6]),
        },
    }

    equation_titles = {
        "kuramotosivashinsky": "Kuramoto-Sivashinsky",
        "burgers": "2D Burgers",
        "reactiondiffusion": "2D Reaction-Diffusion",
        "gasdynamics": "1D Gas Dynamics",
        "advection": "2D Advection",
        "wave": "2D Wave",
    }

    # Axis limits for each equation, matching the paper's Figure 3.
    axis_limits = {
        "kuramotosivashinsky": ([0.8, 120], [1e-7, 1e5]),
        "burgers": ([0.8, 120], [1e-6, 1e3]),
        "reactiondiffusion": ([0.8, 120], [1e-7, 1e2]),
        "gasdynamics": ([0.8, 120], [1e-7, 1e4]),
        "advection": ([0.8, 120], [5e-8, 5e-5]),
        "wave": ([0.8, 120], [5e-8, 5e-5]),
    }

    if equation not in paper_data:
        print(f"MSE comparison plot not available for equation '{equation}'. Skipping.")
        return

    curves_path = rundir / "metrics" / "curves.json"
    if not curves_path.exists():
        return

    with curves_path.open("r", encoding="utf-8") as f:
        curves = json.load(f)

    coeff_mse = np.asarray(curves.get("coeff_mse", []), dtype=float)
    dt = curves.get("dt")

    if not coeff_mse.size or dt is None:
        print("MSE comparison plot requires 'coeff_mse' and 'dt' in curves.json. Skipping.")
        return

    # Add markers to match the style of the paper's plots.
    styles = {
        "SINDy (paper)": {"color": "blue", "linestyle": "solid", "marker": "o", "markersize": 5},
        "Trapping SINDy (paper)": {"color": "orange", "linestyle": "dashed", "marker": "s", "markersize": 5},
        "Galerkin (paper)": {"color": "green", "linestyle": "dotted", "marker": "^", "markersize": 5},
        "ML-Galerkin (paper)": {"color": "red", "linestyle": "dashdot", "marker": "x", "markersize": 6},
    }

    fig, ax = plt.subplots(figsize=(8, 6), dpi=dpi)

    for name, (t_paper, mse_paper) in paper_data[equation].items():
        ax.plot(t_paper, mse_paper, label=name, **styles[name], markevery=1)

    t_run = np.arange(len(coeff_mse)) * float(dt)
    ax.plot(t_run, coeff_mse, color="black", linewidth=2.5, label="SINDy")

    ax.set_xscale("log"); ax.set_yscale("log")

    # Apply axis limits to match the paper's formatting
    if equation in axis_limits:
        xlim, ylim = axis_limits[equation]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    ax.set_xlabel("Time"); ax.set_ylabel("Test MSE")
    title = equation_titles.get(equation, equation.capitalize())
    ax.set_title(f"Coefficient MSE Comparison for {title}")
    ax.legend(loc="upper left"); ax.grid(True, which="both", ls="--", alpha=0.5)
    
    _savefig(fig, out_dir / "metrics_coeff_mse_comparison.png", dpi=dpi)

def animate_3d_surface(
    q_true: np.ndarray,
    q_pred: np.ndarray,
    equation: str,
    output_path: Path,
    fps: int = 15,
    dpi: int = 100,
    interp_factor: int = 3,
) -> None:
    """
    Animate a 3D surface plot of fields, with behavior depending on the equation.

    For 'burgers', this plots vorticity as height and velocity direction as color.
    Uses cubic spline interpolation for a smooth appearance.

    Parameters
    ----------
    q_true : np.ndarray
        Ground truth trajectory of shape (T, C, ny, nx).
    q_pred : np.ndarray
        Predicted trajectory of shape (T, C, ny, nx).
    equation : str
        Name of the equation, e.g., 'burgers'.
    output_path : Path
        Destination path for the movie file.
    fps : int
        Frames per second for the animation.
    dpi : int
        DPI for rendering frames.
    interp_factor : int
        Factor by which to increase grid resolution for interpolation.
    """
    T, C, ny, nx = q_true.shape

    # --- Data preparation based on equation ---
    if equation == "burgers":
        if C != 2:
            print(f"Skipping 3D surface plot for 'burgers': requires 2 components, found {C}.")
            return
        
        # Height data: vorticity
        z_true = compute_vorticity(q_true)
        z_pred = compute_vorticity(q_pred)

        # Color data (velocity angle) is handled inside the update loop
        # by interpolating u and v components.
        z_label = "Vorticity"
        color_mode = "hsv"
        cmap = "hsv"  # For burgers, this is handled by hsv_to_rgb

    elif equation == "kuramotosivashinsky":
        if C != 1:
            print(f"Skipping 3D surface plot for 'kuramotosivashinsky': requires 1 component, found {C}.")
            return
        
        # Height is the field value itself
        z_true = q_true[:, 0]
        z_pred = q_pred[:, 0]

        # Color is the magnitude of the spatial gradient
        grad_true_y, grad_true_x = np.gradient(z_true, axis=(1, 2))
        color_data_true = np.sqrt(grad_true_x**2 + grad_true_y**2)
        grad_pred_y, grad_pred_x = np.gradient(z_pred, axis=(1, 2))
        color_data_pred = np.sqrt(grad_pred_x**2 + grad_pred_y**2)

        z_label = "Field Value u(x,y)"
        color_mode = "colormap"
        cmap = "viridis"

    elif equation == "gasdynamics":
        if C != 4:
            print(f"Skipping 3D surface plot for 'gasdynamics': requires 4 components, found {C}.")
            return
        
        # Height is pressure, color is velocity magnitude
        z_true, color_data_true = compute_gasdynamics_vars(q_true)
        z_pred, color_data_pred = compute_gasdynamics_vars(q_pred)

        z_label = "Pressure"
        color_mode = "colormap"
        cmap = "plasma"

    else:
        print(f"3D surface plot not implemented for equation '{equation}'. Skipping.")
        return

    # --- Interpolation setup ---
    x = np.arange(nx)
    y = np.arange(ny)
    x_fine = np.linspace(0, nx - 1, nx * interp_factor)
    y_fine = np.linspace(0, ny - 1, ny * interp_factor)
    xx_fine, yy_fine = np.meshgrid(x_fine, y_fine)

    # --- Plotting setup ---
    fig = plt.figure(figsize=(16, 7), dpi=dpi)
    ax_true = fig.add_subplot(1, 2, 1, projection='3d')
    ax_pred = fig.add_subplot(1, 2, 2, projection='3d')
    title = fig.suptitle("Time: 0")

    # Create a light source for manual shading
    light = LightSource(azdeg=225, altdeg=45)

    # Determine shared z-axis limits for consistency
    # Determine shared z-axis limits for consistency, but constrain the bounds
    # to be no more than double the true data's range to prevent outliers in
    # the prediction from making the plot unreadable.
    z_true_min, z_true_max = z_true.min(), z_true.max()
    z_pred_min, z_pred_max = z_pred.min(), z_pred.max()

    z_min_unbounded = min(z_true_min, z_pred_min)
    z_max_unbounded = max(z_true_max, z_pred_max)

    # Define the capped range based on the true data's range.
    z_true_range = z_true_max - z_true_min
    # If the true range is zero, use a small default range to avoid a collapsed axis.
    # The range is set relative to the magnitude of the data.
    if z_true_range < 1e-9:
        z_true_range = max(1.0, abs(z_true_max) * 0.2)

    z_true_center = (z_true_max + z_true_min) / 2

    # The allowed range is centered on the true data, with twice the span.
    z_min_limit = z_true_center - z_true_range
    z_max_limit = z_true_center + z_true_range

    z_min = max(z_min_unbounded, z_min_limit)
    z_max = min(z_max_unbounded, z_max_limit)

    # Pre-calculate color normalization for colormap-based plots
    if color_mode == "colormap":
        color_min = min(color_data_true.min(), color_data_pred.min())
        color_max = max(color_data_true.max(), color_data_pred.max())
        # Handle case where max == min to avoid division by zero
        if color_max - color_min < 1e-9:
            color_max = color_min + 1.0
        color_norm = Normalize(vmin=color_min, vmax=color_max)
        cmap_obj = plt.get_cmap(cmap)

    def update(frame: int):
        ax_true.clear()
        ax_pred.clear()
        
        title.set_text(f"Time: {frame}")

        # --- Process and plot TRUE data ---
        # Interpolate height (vorticity)
        spline_z_true = RectBivariateSpline(y, x, z_true[frame], kx=3, ky=3)
        z_fine_true = spline_z_true(y_fine, x_fine)

        # --- Manually combine color (from angle) and shading (from height) ---
        illumination_true = light.hillshade(z_fine_true, vert_exag=0.8)
        rescaled_illumination_true = 0.5 + 0.5 * illumination_true

        if color_mode == "hsv":
            # Interpolate color (velocity components) and compute angle
            u_true_frame, v_true_frame = q_true[frame, 0], q_true[frame, 1]
            spline_u_true = RectBivariateSpline(y, x, u_true_frame, kx=3, ky=3)
            spline_v_true = RectBivariateSpline(y, x, v_true_frame, kx=3, ky=3)
            u_fine_true = spline_u_true(y_fine, x_fine)
            v_fine_true = spline_v_true(y_fine, x_fine)
            color_fine_true = np.arctan2(v_fine_true, u_fine_true)

            # 1. Calculate per-vertex hue from velocity angle.
            hue_true = (color_fine_true + np.pi) / (2 * np.pi)
            # 2. Calculate per-vertex saturation (full saturation).
            saturation_true = np.ones_like(hue_true)
            # 3. Combine into an HSV array and convert to RGB. This gives per-vertex colors.
            hsv_vertex_true = np.stack([hue_true, saturation_true, rescaled_illumination_true], axis=-1)
            rgb_vertex_true = hsv_to_rgb(hsv_vertex_true)

        elif color_mode == "colormap":
            # Interpolate scalar color data
            spline_color_true = RectBivariateSpline(y, x, color_data_true[frame], kx=3, ky=3)
            color_fine_true = spline_color_true(y_fine, x_fine)

            # 1. Get RGB values from colormap using pre-calculated normalization.
            colors_from_map = cmap_obj(color_norm(color_fine_true))[:, :, :3]  # Drop alpha
            # 2. Modulate brightness with illumination to get per-vertex colors.
            rgb_vertex_true = colors_from_map * rescaled_illumination_true[..., np.newaxis]

        else:
            # Fallback to simple shading if color mode is unknown
            rgb_vertex_true = light.shade(z_fine_true, cmap=plt.get_cmap('gray'))

        # To fix the "blue and orange" issue, we provide explicit face colors
        # by averaging the colors of the four vertices of each face. This is more robust.
        rgb_face_true = (rgb_vertex_true[:-1, :-1, :] + rgb_vertex_true[1:, :-1, :] + rgb_vertex_true[:-1, 1:, :] + rgb_vertex_true[1:, 1:, :]) / 4.0

        ax_true.set_title("True Field")
        ax_true.set_zlim(z_min, z_max)
        ax_true.set_xlabel("x"); ax_true.set_ylabel("y"); ax_true.set_zlabel(z_label)
        ax_true.plot_surface(xx_fine, yy_fine, z_fine_true, facecolors=rgb_face_true, rstride=1, cstride=1, antialiased=True, shade=False)

        # --- Process and plot PRED data ---
        # Interpolate height (vorticity)
        spline_z_pred = RectBivariateSpline(y, x, z_pred[frame], kx=3, ky=3)
        z_fine_pred = spline_z_pred(y_fine, x_fine)

        # --- Manually combine color and shading for predicted data ---
        illumination_pred = light.hillshade(z_fine_pred, vert_exag=0.8)
        rescaled_illumination_pred = 0.5 + 0.5 * illumination_pred

        if color_mode == "hsv":
            u_pred_frame, v_pred_frame = q_pred[frame, 0], q_pred[frame, 1]
            spline_u_pred = RectBivariateSpline(y, x, u_pred_frame, kx=3, ky=3)
            spline_v_pred = RectBivariateSpline(y, x, v_pred_frame, kx=3, ky=3)
            u_fine_pred = spline_u_pred(y_fine, x_fine)
            v_fine_pred = spline_v_pred(y_fine, x_fine)
            color_fine_pred = np.arctan2(v_fine_pred, u_fine_pred)

            hue_pred = (color_fine_pred + np.pi) / (2 * np.pi)
            saturation_pred = np.ones_like(hue_pred)
            hsv_vertex_pred = np.stack([hue_pred, saturation_pred, rescaled_illumination_pred], axis=-1)
            rgb_vertex_pred = hsv_to_rgb(hsv_vertex_pred)

        elif color_mode == "colormap":
            spline_color_pred = RectBivariateSpline(y, x, color_data_pred[frame], kx=3, ky=3)
            color_fine_pred = spline_color_pred(y_fine, x_fine)

            colors_from_map = cmap_obj(color_norm(color_fine_pred))[:, :, :3]
            rgb_vertex_pred = colors_from_map * rescaled_illumination_pred[..., np.newaxis]

        else:
            rgb_vertex_pred = light.shade(z_fine_pred, cmap=plt.get_cmap('gray'))

        rgb_face_pred = (rgb_vertex_pred[:-1, :-1, :] + rgb_vertex_pred[1:, :-1, :] + rgb_vertex_pred[:-1, 1:, :] + rgb_vertex_pred[1:, 1:, :]) / 4.0

        ax_pred.set_title("Predicted Field")
        ax_pred.set_zlim(z_min, z_max)
        ax_pred.set_xlabel("x"); ax_pred.set_ylabel("y"); ax_pred.set_zlabel(z_label)
        ax_pred.plot_surface(xx_fine, yy_fine, z_fine_pred, facecolors=rgb_face_pred, rstride=1, cstride=1, antialiased=True, shade=False)

        return [ax_true, ax_pred]

    # Create and save animation
    ani = animation.FuncAnimation(fig, update, frames=T, blit=False, interval=1000 / fps)
    writer = "ffmpeg" if output_path.suffix == ".mp4" else "pillow"
    ani.save(output_path, writer=writer, fps=fps)
    plt.close(fig)
    print(f"Saved 3D surface animation to {output_path}")

def generate_all_plots_and_movies(
    cfg: PlotConfig,
    rundir: Path,
    layout: SnapshotLayout,
    pod: PODResult,
    sindy: SINDyFitResult,
    rollout: RolloutResult,
    equation: str,
) -> None:
    """
    Generate all configured plots and movies for a completed run.

    Parameters
    ----------
    cfg : PlotConfig
        Plot configuration.
    rundir : Path
        Run output directory.
    equation : str
        Name of the equation.
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
        plot_mse_comparison_with_paper(
            rundir=rundir,
            out_dir=fig_dir,
            dpi=cfg.dpi,
            equation=equation,
        )

    if getattr(cfg, "movie_3d_surface", False):
        interp_factor = getattr(cfg, "movie_3d_interp_factor", 3)
        animate_3d_surface(
            q_true=rollout.fields_true,
            q_pred=rollout.fields_pred,
            equation=equation,
            output_path=mov_dir / "rollout_3d_surface.mp4",
            fps=cfg.movie_fps,
            dpi=cfg.dpi,
            interp_factor=interp_factor,
        )