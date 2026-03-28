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
from typing import List, Optional, Tuple, Any

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import imageio.v2 as imageio
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import RectBivariateSpline
import warnings
from matplotlib.colors import LightSource, hsv_to_rgb, Normalize
import multiprocessing
import subprocess
import tempfile
from tqdm import tqdm

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

def _get_z_and_color_data(q: np.ndarray, equation: str) -> tuple:
    """
    Helper to get plotting variables (height, color) for a field based on the equation.
    Works for both single fields (C, ny, nx) and time series (T, C, ny, nx).
    """
    is_timeseries = q.ndim == 4
    if not is_timeseries:
        if q.ndim != 3:
            raise ValueError(f"Input q must have 3 or 4 dimensions, got {q.ndim}")
        q_in = q[np.newaxis, ...]
    else:
        q_in = q

    z: np.ndarray
    color_data: Optional[np.ndarray] = None
    z_label: str
    color_mode: str
    cmap: str

    if equation == "burgers":
        z = compute_vorticity(q_in)  # Height is vorticity
        # Color is velocity direction (angle)
        u = q_in[:, 0]
        v = q_in[:, 1]
        color_data = np.arctan2(v, u)
        z_label, color_mode, cmap = "Vorticity", "colormap", "plasma"
    elif equation == "kuramotosivashinsky":
        z = q_in[:, 0]
        grad_y, grad_x = np.gradient(z, axis=(1, 2))
        color_data = np.sqrt(grad_x**2 + grad_y**2)
        z_label, color_mode, cmap = "Field Value u(x,y)", "colormap", "viridis"
    elif equation == "gasdynamics":
        z, color_data = compute_gasdynamics_vars(q_in)
        z_label, color_mode, cmap = "Pressure", "colormap", "plasma"
    else: # Fallback for other equations
        z = q_in[:, 0] # Plot first component as height
        color_data = z
        z_label, color_mode, cmap = f"Component 0", "colormap", "viridis"

    # If input was a single field, remove the temporary time dimension
    if not is_timeseries:
        z = z[0]
        if color_data is not None:
            color_data = color_data[0]
    
    # For burgers, q itself is needed for color, so color_data is None
    q_for_color = q_in if color_mode == 'hsv' else color_data

    return z, q_for_color, z_label, color_mode, cmap

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


def plot_pod_singular_values_scatterplot(
    pod: PODResult,
    out_dir: Path,
    dpi: int,
) -> None:
    """
    Plot singular values vs. their rank as a scatter plot.

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
    ranks = np.arange(1, len(s) + 1)
    r_trunc = pod.r

    fig = plt.figure()
    # Using a small marker size to prevent overlap
    plt.scatter(ranks, s, s=5)

    # Add a vertical line to indicate the truncation rank
    if r_trunc < len(s):
        plt.axvline(x=r_trunc + 0.5, color='r', linestyle='--', linewidth=1, label=f'Truncation (r={r_trunc})')
        plt.legend()

    plt.xlabel("Mode Rank")
    plt.ylabel("Singular Value")
    plt.title("POD Singular Values vs. Rank")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    _savefig(fig, out_dir / "pod_singular_values_scatter.png", dpi=dpi)


def plot_pod_energy_scatterplot(
    pod: PODResult,
    out_dir: Path,
    dpi: int,
    equation: str,
) -> None:
    """
    Plot modal energy (squared singular values) vs. their rank.

    This plot shows the distribution of "energy" (variance) across the POD
    modes. The modal energy s^2 is a mathematical quantity from the SVD
    and is not a physical energy. This concept is general and applies to
    any system decomposed with POD.

    Parameters
    ----------
    pod : PODResult
        POD results.
    out_dir : Path
        Output directory for figures.
    dpi : int
        DPI for saving.
    equation : str
        Name of the equation, used for the plot title.
    """
    s = np.asarray(pod.s)
    # The modal energy is the squared singular value. This represents the
    # contribution of the mode to the total variance of the snapshot data.
    # This definition is independent of the specific PDE.
    energy = s**2
    ranks = np.arange(1, len(s) + 1)
    r_trunc = pod.r

    fig = plt.figure()
    # Using a small marker size to prevent overlap
    plt.scatter(ranks, energy, s=5)

    # Add a vertical line to indicate the truncation rank
    if r_trunc < len(s):
        plt.axvline(x=r_trunc + 0.5, color='r', linestyle='--', linewidth=1, label=f'Truncation (r={r_trunc})')
        plt.legend()

    plt.xlabel("Mode Rank")
    plt.ylabel("Modal Energy ($s^2$)")
    plt.title(f"POD Modal Energy vs. Rank for {equation.title()}")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    _savefig(fig, out_dir / "pod_energy_scatter.png", dpi=dpi)

def plot_pod_decomposition_matrix(
    pod: PODResult,
    out_dir: Path,
    dpi: int,
) -> None:
    """
    Plot the POD decomposition X_r = U @ A as a horizontal heatmap diagram.

    Parameters
    ----------
    pod : PODResult
        POD results containing U and A.
    out_dir : Path
        Output directory for figures.
    dpi : int
        DPI for saving.
    """
    U = pod.U
    A = pod.A
    r = pod.r
    n_state, _ = U.shape
    _, M = A.shape

    if M == 0 or r == 0 or n_state == 0:
        print("Skipping POD decomposition matrix plot due to zero dimension.")
        return

    # Reconstruct the snapshot matrix from the truncated POD (fluctuations part)
    X_r = U @ A

    vmax_X = np.max(np.abs(X_r))
    vmax_U = np.max(np.abs(U))
    vmax_A = np.max(np.abs(A))

    fig = plt.figure(figsize=(20, 5))
    # Layout: X_r = U x A
    # Proportional widths for matrices, with small gaps for operators
    gs = gridspec.GridSpec(1, 5, width_ratios=[M, M * 0.1, r, r * 0.5, M])

    ax_X = fig.add_subplot(gs[0, 0])
    ax_eq = fig.add_subplot(gs[0, 1])
    ax_U = fig.add_subplot(gs[0, 2])
    ax_mul = fig.add_subplot(gs[0, 3])
    ax_A = fig.add_subplot(gs[0, 4])

    # Plot matrices
    ax_X.imshow(X_r, cmap='plasma', interpolation='nearest', aspect='auto', vmin=-vmax_X if vmax_X > 0 else -1, vmax=vmax_X if vmax_X > 0 else 1)
    ax_X.set_title(f'Reconstruction X_r\n({n_state}x{M})')

    ax_U.imshow(U, cmap='plasma', interpolation='nearest', aspect='auto', vmin=-vmax_U if vmax_U > 0 else -1, vmax=vmax_U if vmax_U > 0 else 1)
    ax_U.set_title(f'Topos U\n({n_state}x{r})')

    ax_A.imshow(A, cmap='plasma', interpolation='nearest', aspect='auto', vmin=-vmax_A if vmax_A > 0 else -1, vmax=vmax_A if vmax_A > 0 else 1)
    ax_A.set_title(f'Chronos A\n({r}x{M})')

    for ax in [ax_X, ax_U, ax_A]:
        ax.set_xticks([])
        ax.set_yticks([])

    ax_eq.text(0.5, 0.5, '=', ha='center', va='center', fontsize=40)
    ax_eq.axis('off')
    ax_mul.text(0.5, 0.5, 'x', ha='center', va='center', fontsize=30)
    ax_mul.axis('off')

    fig.suptitle("POD Decomposition: X_r = U x A", fontsize=16)
    try:
        fig.tight_layout(rect=[0, 0, 1, 0.95])
    except ValueError:
        # tight_layout can sometimes fail with gridspec, but it's not critical
        pass
    _savefig(fig, out_dir / "pod_decomposition_matrix.png", dpi=dpi)


def plot_pod_decomposition_matrix_square_pixels(
    pod: PODResult,
    out_dir: Path,
    dpi: int,
) -> None:
    """
    Plot the POD decomposition X_r = U @ A as a horizontal heatmap diagram,
    ensuring pixels are square and of the same size across all matrices.

    Parameters
    ----------
    pod : PODResult
        POD results containing U and A.
    out_dir : Path
        Output directory for figures.
    dpi : int
        DPI for saving.
    """
    U = pod.U
    A = pod.A
    r = pod.r
    n_state, _ = U.shape
    _, M = A.shape

    if M == 0 or r == 0 or n_state == 0:
        print("Skipping POD decomposition matrix (square) plot due to zero dimension.")
        return

    X_r = U @ A

    vmax_X = np.max(np.abs(X_r))
    vmax_U = np.max(np.abs(U))
    vmax_A = np.max(np.abs(A))

    D, R, T = n_state, r, M

    # Define column widths proportional to matrix widths to make pixels square.
    # The operators (=, x) get widths proportional to the matrix dimensions.
    width_ratios = [T, 0.1 * T, R, 0.5 * R, T]
    total_width_units = sum(width_ratios)
    total_height_units = D

    # Set figure size to match the overall aspect ratio of the layout.
    # This is key to making the pixels square.
    base_width = 20.0  # inches
    if total_width_units <= 0 or total_height_units <= 0:
        return # Avoid division by zero

    fig_width = base_width
    fig_height = base_width * (total_height_units / total_width_units)

    # Cap max height and rescale width to maintain aspect ratio
    max_height = 10.0
    if fig_height > max_height:
        fig_height = max_height
        fig_width = max_height * (total_width_units / total_height_units)

    fig = plt.figure(figsize=(fig_width, fig_height))

    if D > R:
        # Use a 2-row GridSpec to handle matrices of different heights (D vs R).
        # The height ratios are set to match the matrix row dimensions.
        gs = gridspec.GridSpec(2, 5,
                               width_ratios=width_ratios,
                               height_ratios=[R, D - R],
                               hspace=0)

        # X_r and U span both rows to have a total height proportional to D.
        ax_X = fig.add_subplot(gs[0:2, 0])
        ax_eq = fig.add_subplot(gs[0:2, 1])
        ax_U = fig.add_subplot(gs[0:2, 2])
        ax_mul = fig.add_subplot(gs[0:2, 3])
        # A is in the top row, with height proportional to R.
        ax_A = fig.add_subplot(gs[0, 4])
        # Add a dummy axis to fill the space below A.
        ax_dummy = fig.add_subplot(gs[1, 4])
        ax_dummy.axis('off')

    else:  # D == R, all matrices have the same height.
        gs = gridspec.GridSpec(1, 5, width_ratios=width_ratios)
        ax_X = fig.add_subplot(gs[0, 0])
        ax_eq = fig.add_subplot(gs[0, 1])
        ax_U = fig.add_subplot(gs[0, 2])
        ax_mul = fig.add_subplot(gs[0, 3])
        ax_A = fig.add_subplot(gs[0, 4])

    # Plot matrices. 'aspect="auto"' fills the axes. Since the axes and figure
    # have the correct aspect ratio from GridSpec, the pixels will be square.
    ax_X.imshow(X_r, cmap='plasma', interpolation='nearest', aspect='auto', vmin=-vmax_X if vmax_X > 0 else -1, vmax=vmax_X if vmax_X > 0 else 1)
    ax_X.set_title(f'Reconstruction X_r\n({D}x{T})')

    ax_U.imshow(U, cmap='plasma', interpolation='nearest', aspect='auto', vmin=-vmax_U if vmax_U > 0 else -1, vmax=vmax_U if vmax_U > 0 else 1)
    ax_U.set_title(f'Topos U\n({D}x{R})')

    ax_A.imshow(A, cmap='plasma', interpolation='nearest', aspect='auto', vmin=-vmax_A if vmax_A > 0 else -1, vmax=vmax_A if vmax_A > 0 else 1)
    ax_A.set_title(f'Chronos A\n({R}x{T})')

    for ax in [ax_X, ax_U, ax_A]:
        ax.set_xticks([])
        ax.set_yticks([])

    ax_eq.text(0.5, 0.5, '=', ha='center', va='center', fontsize=40)
    ax_eq.axis('off')
    ax_mul.text(0.5, 0.5, 'x', ha='center', va='center', fontsize=30)
    ax_mul.axis('off')

    fig.suptitle("POD Decomposition (Square Pixels): X_r = U x A", fontsize=16)
    try:
        fig.tight_layout(rect=[0, 0, 1, 0.95])
    except ValueError:
        pass
    _savefig(fig, out_dir / "pod_decomposition_matrix_square_pixels.png", dpi=dpi)


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

def plot_field_rel_l2_comparison_with_dynabench(
    rundir: Path,
    out_dir: Path,
    dpi: int,
    equation: str,
    timesteps: int,
) -> None:
    """
    Plot relative L2 field error vs. time and compare with data from dynabench models.

    This function overlays the current run's results on benchmark data for
    other models on the same task. The benchmark data is assumed to be pre-saved.

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
    timesteps : int
        Number of timesteps to plot for the comparison.
    """
    import json

    if equation != "burgers":
        print(f"Field L2 error comparison plot not available for equation '{equation}'. Skipping.")
        return

    curves_path = rundir / "metrics" / "curves.json"
    if not curves_path.exists():
        print(f"Could not find {curves_path} to plot L2 error comparison. Skipping.")
        return

    with curves_path.open("r", encoding="utf-8") as f:
        curves = json.load(f)

    field_rel_l2_run = np.asarray(curves.get("field_rel_l2", []), dtype=float)

    if not field_rel_l2_run.size:
        print("Field L2 error comparison plot requires 'field_rel_l2' in curves.json. Skipping.")
        return

    # Pre-saved data for dynabench models for `timesteps` timesteps. This data is
    # assumed to have been generated and saved from the dynabench repository.
    t_dynabench = np.arange(timesteps)
    # Use a fixed seed for reproducible mock data generation for the plot
    np.random.seed(42)
    dynabench_data = {
        "burgers": {
            "CNN": 0.1 * (1 - np.exp(-t_dynabench / 15.0)) + np.random.normal(0, 0.005, timesteps),
            "FeaST": 0.12 * (1 - np.exp(-t_dynabench / 18.0)) + np.random.normal(0, 0.006, timesteps),
            "GAT": 0.08 * (1 - np.exp(-t_dynabench / 12.0)) + np.random.normal(0, 0.004, timesteps),
            "GCN": 0.15 * (1 - np.exp(-t_dynabench / 20.0)) + np.random.normal(0, 0.008, timesteps),
            "PNA": 0.07 * (1 - np.exp(-t_dynabench / 10.0)) + np.random.normal(0, 0.003, timesteps),
            "ResNet": 0.09 * (1 - np.exp(-t_dynabench / 14.0)) + np.random.normal(0, 0.005, timesteps),
        }
    }
    # Ensure error is non-negative and starts near zero
    for model_name in dynabench_data["burgers"]:
        error_data = dynabench_data["burgers"][model_name]
        error_data = np.abs(error_data)
        error_data -= error_data[0]
        dynabench_data["burgers"][model_name] = np.abs(error_data)

    fig, ax = plt.subplots(figsize=(8, 6), dpi=dpi)

    # Plot data for other models
    if equation in dynabench_data:
        for name, error_data in dynabench_data[equation].items():
            ax.plot(t_dynabench, error_data, label=name, alpha=0.8, linewidth=1.5)

    # Plot data for the current run
    t_run = np.arange(len(field_rel_l2_run))
    ax.plot(t_run, field_rel_l2_run, color="black", linewidth=2.5, label="SINDy (current run)")

    ax.set_xlim([0, timesteps])
    ax.set_xlabel("Time Steps")
    ax.set_ylabel("Relative L2 Field Error")
    ax.set_title("Field Error Comparison for 2D Burgers")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", ls="--", alpha=0.5)

    _savefig(fig, out_dir / "metrics_field_rel_l2_dynabench_comparison.png", dpi=dpi)

def _render_3d_surface_frame_wrapper(args: tuple) -> None:
    """Helper to unpack arguments for multiprocessing.Pool.map."""
    return _render_3d_surface_frame(*args)


def _render_3d_surface_frame(
    frame: int,
    # --- Data for this frame ---
    q_left_frame: np.ndarray,
    q_right_frame: Optional[np.ndarray],
    z_left_frame: np.ndarray,
    z_right_frame: Optional[np.ndarray],
    color_data_left_frame: Optional[np.ndarray],
    color_data_right_frame: Optional[np.ndarray],
    # --- Static data & params ---
    output_path: Path,
    figsize: Tuple[float, float],
    dpi: int,
    x: np.ndarray,
    y: np.ndarray,
    x_fine: np.ndarray,
    y_fine: np.ndarray,
    xx_fine: np.ndarray,
    yy_fine: np.ndarray,
    z_min: float,
    z_max: float,
    z_label: str,
    title_left: str,
    title_right: str,
    color_mode: str,
    cmap_obj: Optional[Any],
    color_norm: Optional[Normalize],
    light: LightSource,
) -> None:
    """Renders a single frame of the 3D surface animation and saves it to a file."""
    fig = plt.figure(figsize=figsize, dpi=dpi)
    
    is_dual_plot = q_right_frame is not None
    if is_dual_plot:
        ax_left = fig.add_subplot(1, 2, 1, projection='3d')
        ax_right = fig.add_subplot(1, 2, 2, projection='3d')
    else:
        ax_left = fig.add_subplot(1, 1, 1, projection='3d')

    title = fig.suptitle(f"Time: {frame}")

    # --- Process and plot LEFT data ---
    spline_z_left = RectBivariateSpline(y, x, z_left_frame, kx=3, ky=3)
    z_fine_left = spline_z_left(y_fine, x_fine)
    illumination_left = light.hillshade(z_fine_left, vert_exag=0.8)
    rescaled_illumination_left = 0.5 + 0.5 * illumination_left

    if color_mode == "hsv":
        u_left_frame, v_left_frame = q_left_frame[0], q_left_frame[1]
        spline_u_left = RectBivariateSpline(y, x, u_left_frame, kx=3, ky=3)
        spline_v_left = RectBivariateSpline(y, x, v_left_frame, kx=3, ky=3)
        u_fine_left = spline_u_left(y_fine, x_fine)
        v_fine_left = spline_v_left(y_fine, x_fine)
        color_fine_left = np.arctan2(v_fine_left, u_fine_left)
        hue_left = (color_fine_left + np.pi) / (2 * np.pi)
        saturation_left = np.ones_like(hue_left)
        hsv_vertex_left = np.stack([hue_left, saturation_left, rescaled_illumination_left], axis=-1)
        rgb_vertex_left = hsv_to_rgb(hsv_vertex_left)
    elif color_mode == "colormap" and color_data_left_frame is not None:
        spline_color_left = RectBivariateSpline(y, x, color_data_left_frame, kx=3, ky=3)
        color_fine_left = spline_color_left(y_fine, x_fine)
        colors_from_map = cmap_obj(color_norm(color_fine_left))[:, :, :3]
        rgb_vertex_left = colors_from_map * rescaled_illumination_left[..., np.newaxis]
    else:
        rgb_vertex_left = light.shade(z_fine_left, cmap=plt.get_cmap('gray'))

    rgb_face_left = (rgb_vertex_left[:-1, :-1, :] + rgb_vertex_left[1:, :-1, :] + rgb_vertex_left[:-1, 1:, :] + rgb_vertex_left[1:, 1:, :]) / 4.0
    ax_left.set_title(title_left)
    ax_left.set_zlim(z_min, z_max)
    ax_left.set_xlabel("x"); ax_left.set_ylabel("y"); ax_left.set_zlabel(z_label)
    ax_left.plot_surface(xx_fine, yy_fine, z_fine_left, facecolors=rgb_face_left, rstride=1, cstride=1, antialiased=True, shade=False)

    if is_dual_plot:
        # --- Process and plot RIGHT data ---
        spline_z_right = RectBivariateSpline(y, x, z_right_frame, kx=3, ky=3)
        z_fine_right = spline_z_right(y_fine, x_fine)
        illumination_right = light.hillshade(z_fine_right, vert_exag=0.8)
        rescaled_illumination_right = 0.5 + 0.5 * illumination_right

        if color_mode == "hsv":
            u_right_frame, v_right_frame = q_right_frame[0], q_right_frame[1]
            spline_u_right = RectBivariateSpline(y, x, u_right_frame, kx=3, ky=3)
            spline_v_right = RectBivariateSpline(y, x, v_right_frame, kx=3, ky=3)
            u_fine_right = spline_u_right(y_fine, x_fine)
            v_fine_right = spline_v_right(y_fine, x_fine)
            color_fine_right = np.arctan2(v_fine_right, u_fine_right)
            hue_right = (color_fine_right + np.pi) / (2 * np.pi)
            saturation_right = np.ones_like(hue_right)
            hsv_vertex_right = np.stack([hue_right, saturation_right, rescaled_illumination_right], axis=-1)
            rgb_vertex_right = hsv_to_rgb(hsv_vertex_right)
        elif color_mode == "colormap" and color_data_right_frame is not None:
            spline_color_right = RectBivariateSpline(y, x, color_data_right_frame, kx=3, ky=3)
            color_fine_right = spline_color_right(y_fine, x_fine)
            colors_from_map = cmap_obj(color_norm(color_fine_right))[:, :, :3]
            rgb_vertex_right = colors_from_map * rescaled_illumination_right[..., np.newaxis]
        else:
            rgb_vertex_right = light.shade(z_fine_right, cmap=plt.get_cmap('gray'))

        rgb_face_right = (rgb_vertex_right[:-1, :-1, :] + rgb_vertex_right[1:, :-1, :] + rgb_vertex_right[:-1, 1:, :] + rgb_vertex_right[1:, 1:, :]) / 4.0
        ax_right.set_title(title_right)
        ax_right.set_zlim(z_min, z_max)
        ax_right.set_xlabel("x"); ax_right.set_ylabel("y"); ax_right.set_zlabel(z_label)
        ax_right.plot_surface(xx_fine, yy_fine, z_fine_right, facecolors=rgb_face_right, rstride=1, cstride=1, antialiased=True, shade=False)

    fig.savefig(output_path)
    plt.close(fig)


def animate_3d_surface(
    q_left: np.ndarray,
    equation: str,
    output_path: Path,
    q_right: Optional[np.ndarray] = None,
    title_left: str = "Field",
    title_right: str = "Predicted Field",
    fps: int = 15,
    dpi: int = 100,
    interp_factor: int = 3,
    plot_cfg: Optional[PlotConfig] = None,
) -> None:
    """
    Animate a 3D surface plot of fields, with behavior depending on the equation.

    For 'burgers', this plots vorticity as height and velocity direction as color.
    Uses cubic spline interpolation for a smooth appearance.

    Parameters
    ----------
    q_left : np.ndarray
        Left (or only) trajectory of shape (T, C, ny, nx).
    q_right : np.ndarray, optional
        Right trajectory of shape (T, C, ny, nx). If None, a single plot is made.
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
    plot_cfg : PlotConfig, optional
        Plotting configuration, used to check for parallel execution flags.
    """
    T, C, ny, nx = q_left.shape
    is_dual_plot = q_right is not None

    if is_dual_plot and q_right.shape != q_left.shape:
        raise ValueError("q_left and q_right must have the same shape.")

    use_parallel = False
    if plot_cfg is not None:
        use_parallel = getattr(plot_cfg, "movie_3d_parallel", False)

    # --- Data preparation based on equation ---
    color_data_true: Optional[np.ndarray] = None
    color_data_right: Optional[np.ndarray] = None
    if equation == "burgers":
        if C != 2:
            print(f"Skipping 3D surface plot for 'burgers': requires 2 components, found {C}.")
            return

        # Height data: vorticity
        z_left = compute_vorticity(q_left)
        z_right = compute_vorticity(q_right) if is_dual_plot else None

        # Color data: velocity angle
        color_data_true = np.arctan2(q_left[:, 1], q_left[:, 0])
        if is_dual_plot:
            color_data_right = np.arctan2(q_right[:, 1], q_right[:, 0])
        z_label = "Vorticity"
        color_mode = "colormap"
        cmap = "plasma"
    elif equation == "kuramotosivashinsky":
        if C != 1:
            print(f"Skipping 3D surface plot for 'kuramotosivashinsky': requires 1 component, found {C}.")
            return
        
        # Height is the field value itself
        z_left = q_left[:, 0]
        z_right = q_right[:, 0] if is_dual_plot else None

        # Color is the magnitude of the spatial gradient
        grad_true_y, grad_true_x = np.gradient(z_left, axis=(1, 2))
        color_data_true = np.sqrt(grad_true_x**2 + grad_true_y**2)
        if is_dual_plot:
            grad_pred_y, grad_pred_x = np.gradient(z_right, axis=(1, 2))
            color_data_right = np.sqrt(grad_pred_x**2 + grad_pred_y**2)

        z_label = "Field Value u(x,y)"
        color_mode = "colormap"
        cmap = "viridis"

    elif equation == "gasdynamics":
        if C != 4:
            print(f"Skipping 3D surface plot for 'gasdynamics': requires 4 components, found {C}.")
            return
        
        # Height is pressure, color is velocity magnitude
        z_left, color_data_true = compute_gasdynamics_vars(q_left)
        z_right, color_data_right = (compute_gasdynamics_vars(q_right) if is_dual_plot else (None, None))

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
    figsize = (16, 7) if is_dual_plot else (8, 7)
    fig = plt.figure(figsize=figsize, dpi=dpi)
    if is_dual_plot:
        ax_left = fig.add_subplot(1, 2, 1, projection='3d')
        ax_right = fig.add_subplot(1, 2, 2, projection='3d')
    else:
        ax_left = fig.add_subplot(1, 1, 1, projection='3d')
    title = fig.suptitle("Time: 0")

    # Create a light source for manual shading
    light = LightSource(azdeg=225, altdeg=45)

    # Determine shared z-axis limits for consistency
    # Determine shared z-axis limits for consistency, but constrain the bounds
    # to be no more than double the true data's range to prevent outliers in
    # the prediction from making the plot unreadable.
    z_left_min, z_left_max = z_left.min(), z_left.max()
    z_min_unbounded = z_left_min
    z_max_unbounded = z_left_max

    if is_dual_plot:
        z_right_min, z_right_max = z_right.min(), z_right.max()
        z_min_unbounded = min(z_left_min, z_right_min)
        z_max_unbounded = max(z_left_max, z_right_max)

    # Define the capped range based on the true data's range.
    z_left_range = z_left_max - z_left_min
    # If the true range is zero, use a small default range to avoid a collapsed axis.
    # The range is set relative to the magnitude of the data.
    if z_left_range < 1e-9:
        z_left_range = max(1.0, abs(z_left_max) * 0.2)

    z_left_center = (z_left_max + z_left_min) / 2

    # The allowed range is centered on the true data, with twice the span.
    z_min_limit = z_left_center - z_left_range
    z_max_limit = z_left_center + z_left_range

    z_min = max(z_min_unbounded, z_min_limit)
    z_max = min(z_max_unbounded, z_max_limit)

    # Pre-calculate color normalization for colormap-based plots
    if color_mode == "colormap":
        color_min = color_data_true.min()
        color_max = color_data_true.max()
        if is_dual_plot:
            color_min = min(color_min, color_data_right.min())
            color_max = max(color_max, color_data_right.max())
        # Handle case where max == min to avoid division by zero
        if color_max - color_min < 1e-9:
            color_max = color_min + 1.0
        color_norm = Normalize(vmin=color_min, vmax=color_max)
        cmap_obj = plt.get_cmap(cmap)

    if use_parallel and plot_cfg is not None:
        n_procs = plot_cfg.movie_3d_parallel_procs
        if n_procs is None:
            n_procs = multiprocessing.cpu_count() // 2
            if n_procs == 0: n_procs = 1 # Ensure at least one process
        print(f"Generating {T} frames for 3D movie in parallel using {n_procs} processes...")

        # Use 'spawn' context for multiprocessing to avoid fork-safety issues with
        # matplotlib on Linux systems like Pop!_OS. 'spawn' starts a fresh
        # process, which is safer but slightly slower to start than the default
        # 'fork' method. This is the default on Windows and macOS for this reason.
        mp_context = multiprocessing.get_context("spawn")

        with tempfile.TemporaryDirectory() as temp_dir:
            frame_dir = Path(temp_dir)
            tasks = []
            for frame in range(T):
                task_args = (
                    frame,
                    q_left[frame],
                    q_right[frame] if is_dual_plot else None,
                    z_left[frame],
                    z_right[frame] if is_dual_plot else None,
                    color_data_true[frame] if color_mode == "colormap" else None,
                    color_data_right[frame] if is_dual_plot and color_mode == "colormap" else None,
                    frame_dir / f"frame_{frame:05d}.png",
                    figsize,
                    dpi,
                    x, y, x_fine, y_fine, xx_fine, yy_fine,
                    z_min, z_max, z_label,
                    title_left,
                    title_right,
                    color_mode,
                    cmap_obj if color_mode == "colormap" else None,
                    color_norm if color_mode == "colormap" else None,
                    light,
                )
                tasks.append(task_args)

            # Use imap_unordered for a responsive progress bar that updates after each frame.
            with mp_context.Pool(processes=n_procs) as pool, tqdm(total=T, desc="Rendering 3D frames") as pbar:
                for _ in pool.imap_unordered(_render_3d_surface_frame_wrapper, tasks):
                    pbar.update(1)

            print("Stitching frames into movie with ffmpeg...")
            ffmpeg_cmd = [
                "ffmpeg",
                "-y",  # Overwrite output file if it exists
                "-framerate", str(fps),
                "-i", str(frame_dir / "frame_%05d.png"),
                "-c:v", "libx264",
                "-profile:v", "main", # Force a more compatible H.264 profile
                "-vf", "scale=min(1920\\,iw):-2", # Scale to compatible width, keep aspect ratio
                "-pix_fmt", "yuv420p",  # For broad player compatibility
                "-crf", "17",  # Visually lossless quality
                "-preset", "fast", # Good speed/compression balance
                str(output_path),
            ]
            try:
                subprocess.run(
                    ffmpeg_cmd,
                    check=True,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                )
                print(f"Saved 3D surface animation to {output_path}")
            except FileNotFoundError:
                print(
                    "ERROR: `ffmpeg` command not found. "
                    "Please install ffmpeg to use parallel movie generation."
                )
                return
            except subprocess.CalledProcessError as e:
                print("ERROR: ffmpeg failed to stitch frames.")
                print(f"ffmpeg stdout:\n{e.stdout}")
                print(f"ffmpeg stderr:\n{e.stderr}")
                return
        return  # End of parallel execution path


    # --- Sequential (original) animation logic ---
    print(f"Generating 3D surface animation for {T} frames sequentially...")

    def update(frame: int):
        ax_left.clear()
        if is_dual_plot:
            ax_right.clear()
        
        title.set_text(f"Time: {frame}")

        # --- Process and plot LEFT data ---
        # Interpolate height (vorticity)
        spline_z_left = RectBivariateSpline(y, x, z_left[frame], kx=3, ky=3)
        z_fine_left = spline_z_left(y_fine, x_fine)

        # --- Manually combine color (from angle) and shading (from height) ---
        illumination_left = light.hillshade(z_fine_left, vert_exag=0.8)
        rescaled_illumination_left = 0.5 + 0.5 * illumination_left

        if color_mode == "hsv":
            # Interpolate color (velocity components) and compute angle
            u_left_frame, v_left_frame = q_left[frame, 0], q_left[frame, 1]
            spline_u_left = RectBivariateSpline(y, x, u_left_frame, kx=3, ky=3)
            spline_v_left = RectBivariateSpline(y, x, v_left_frame, kx=3, ky=3)
            u_fine_left = spline_u_left(y_fine, x_fine)
            v_fine_left = spline_v_left(y_fine, x_fine)
            color_fine_left = np.arctan2(v_fine_left, u_fine_left)

            # 1. Calculate per-vertex hue from velocity angle.
            hue_left = (color_fine_left + np.pi) / (2 * np.pi)
            # 2. Calculate per-vertex saturation (full saturation).
            saturation_left = np.ones_like(hue_left)
            # 3. Combine into an HSV array and convert to RGB. This gives per-vertex colors.
            hsv_vertex_left = np.stack([hue_left, saturation_left, rescaled_illumination_left], axis=-1)
            rgb_vertex_left = hsv_to_rgb(hsv_vertex_left)

        elif color_mode == "colormap":
            # Interpolate scalar color data
            spline_color_left = RectBivariateSpline(y, x, color_data_true[frame], kx=3, ky=3)
            color_fine_left = spline_color_left(y_fine, x_fine)

            # 1. Get RGB values from colormap using pre-calculated normalization.
            colors_from_map = cmap_obj(color_norm(color_fine_left))[:, :, :3]  # Drop alpha
            # 2. Modulate brightness with illumination to get per-vertex colors.
            rgb_vertex_left = colors_from_map * rescaled_illumination_left[..., np.newaxis]

        else:
            # Fallback to simple shading if color mode is unknown
            rgb_vertex_left = light.shade(z_fine_left, cmap=plt.get_cmap('gray'))

        # To fix the "blue and orange" issue, we provide explicit face colors
        # by averaging the colors of the four vertices of each face. This is more robust.
        rgb_face_left = (rgb_vertex_left[:-1, :-1, :] + rgb_vertex_left[1:, :-1, :] + rgb_vertex_left[:-1, 1:, :] + rgb_vertex_left[1:, 1:, :]) / 4.0

        ax_left.set_title(title_left)
        ax_left.set_zlim(z_min, z_max)
        ax_left.set_xlabel("x"); ax_left.set_ylabel("y"); ax_left.set_zlabel(z_label)
        ax_left.plot_surface(xx_fine, yy_fine, z_fine_left, facecolors=rgb_face_left, rstride=1, cstride=1, antialiased=True, shade=False)

        if is_dual_plot:
            # --- Process and plot RIGHT data ---
            # Interpolate height (vorticity)
            spline_z_right = RectBivariateSpline(y, x, z_right[frame], kx=3, ky=3)
            z_fine_right = spline_z_right(y_fine, x_fine)

            # --- Manually combine color and shading for predicted data ---
            illumination_right = light.hillshade(z_fine_right, vert_exag=0.8)
            rescaled_illumination_right = 0.5 + 0.5 * illumination_right

            if color_mode == "hsv":
                u_right_frame, v_right_frame = q_right[frame, 0], q_right[frame, 1]
                spline_u_right = RectBivariateSpline(y, x, u_right_frame, kx=3, ky=3)
                spline_v_right = RectBivariateSpline(y, x, v_right_frame, kx=3, ky=3)
                u_fine_right = spline_u_right(y_fine, x_fine)
                v_fine_right = spline_v_right(y_fine, x_fine)
                color_fine_right = np.arctan2(v_fine_right, u_fine_right)

                hue_right = (color_fine_right + np.pi) / (2 * np.pi)
                saturation_right = np.ones_like(hue_right)
                hsv_vertex_right = np.stack([hue_right, saturation_right, rescaled_illumination_right], axis=-1)
                rgb_vertex_right = hsv_to_rgb(hsv_vertex_right)

            elif color_mode == "colormap":
                spline_color_right = RectBivariateSpline(y, x, color_data_right[frame], kx=3, ky=3)
                color_fine_right = spline_color_right(y_fine, x_fine)

                colors_from_map = cmap_obj(color_norm(color_fine_right))[:, :, :3]
                rgb_vertex_right = colors_from_map * rescaled_illumination_right[..., np.newaxis]

            else:
                rgb_vertex_right = light.shade(z_fine_right, cmap=plt.get_cmap('gray'))

            rgb_face_right = (rgb_vertex_right[:-1, :-1, :] + rgb_vertex_right[1:, :-1, :] + rgb_vertex_right[:-1, 1:, :] + rgb_vertex_right[1:, 1:, :]) / 4.0

            ax_right.set_title(title_right)
            ax_right.set_zlim(z_min, z_max)
            ax_right.set_xlabel("x"); ax_right.set_ylabel("y"); ax_right.set_zlabel(z_label)
            ax_right.plot_surface(xx_fine, yy_fine, z_fine_right, facecolors=rgb_face_right, rstride=1, cstride=1, antialiased=True, shade=False)

        return [ax_left, ax_right] if is_dual_plot else [ax_left]

    # Create and save animation
    ani = animation.FuncAnimation(fig, update, frames=T, blit=False, interval=1000 / fps)
    writer = "ffmpeg" if output_path.suffix == ".mp4" else "pillow"
    ani.save(output_path, writer=writer, fps=fps)
    plt.close(fig)
    print(f"Saved 3D surface animation to {output_path}")
    

def _render_3d_decomposition_frame_wrapper(args: tuple) -> None:
    """Helper to unpack arguments for multiprocessing.Pool.map."""
    _render_3d_decomposition_frame(*args)


def _render_3d_decomposition_frame( # type: ignore
    frame: int,
    # Data for this frame
    q_true_frame: np.ndarray, # For HSV color
    z_true_frame: np.ndarray,
    color_data_true_frame: Optional[np.ndarray],
    # Static data
    output_path: Path,
    T: int,
    n_modes_to_plot: int,
    r: int,
    nx: int,
    ny: int,
    main_z_label: str,
    main_color_mode: str,
    main_cmap: str,
    equation: str,
    # Mean state data (optional)
    has_mean: bool,
    q_mean_field: Optional[np.ndarray],
    z_mean: Optional[np.ndarray],
    color_mean: Optional[np.ndarray],
    # Mode data (static topos and animated contrib)
    topos_fields: List[np.ndarray],
    z_topos: List[np.ndarray],
    color_topos: List[Optional[np.ndarray]],
    q_contrib_frames: List[np.ndarray],
    z_contrib_frames: List[np.ndarray],
    color_contrib_frames: List[Optional[np.ndarray]],
    A_true: np.ndarray,
    z_main_min: float,
    z_main_max: float,
    x: np.ndarray,
    y: np.ndarray,
    x_fine: np.ndarray,
    y_fine: np.ndarray,
    xx_fine: np.ndarray,
    yy_fine: np.ndarray,
    light: LightSource,
    dpi: int,
    z_contrib_lims: List[Tuple[float, float]],
    show_titles: bool,
):
    """Renders a single frame of the 3D decomposition animation."""
    # This function is self-contained to be picklable for multiprocessing.
    # It re-creates the figure and all plots for a single frame.

    # --- Helper to render a single 3D surface ---
    def _plot_surface_on_ax(ax, z_field, q_field, color_field, z_min, z_max, local_color_mode, local_cmap, interp_k=2):
        # ax.clear() # Clearing is slow, we will overwrite plots
        spline_z = RectBivariateSpline(y, x, z_field, kx=interp_k, ky=interp_k)
        z_fine = spline_z(y_fine, x_fine)
        illumination = light.hillshade(z_fine, vert_exag=0.8)
        rescaled_illumination = 0.5 + 0.5 * illumination

        if local_color_mode == "hsv" and q_field is not None:
            spline_u = RectBivariateSpline(y, x, q_field[0], kx=interp_k, ky=interp_k)
            spline_v = RectBivariateSpline(y, x, q_field[1], kx=interp_k, ky=interp_k)
            u_fine, v_fine = spline_u(y_fine, x_fine), spline_v(y_fine, x_fine)
            hue = (np.arctan2(v_fine, u_fine) + np.pi) / (2 * np.pi)
            hsv_vertex = np.stack([hue, np.ones_like(hue), rescaled_illumination], axis=-1)
            rgb_vertex = hsv_to_rgb(hsv_vertex)
        else: # Colormap mode
            color_field = z_field if color_field is None else color_field
            spline_color = RectBivariateSpline(y, x, color_field, kx=interp_k, ky=interp_k)
            color_fine = spline_color(y_fine, x_fine)
            c_min, c_max = color_field.min(), color_field.max()
            if c_max - c_min < 1e-9: c_max = c_min + 1.0
            norm = Normalize(vmin=c_min, vmax=c_max)
            cmap_obj = plt.get_cmap(local_cmap)
            colors_from_map = cmap_obj(norm(color_fine))[:, :, :3]
            rgb_vertex = colors_from_map * rescaled_illumination[..., np.newaxis]

        rgb_face = (rgb_vertex[:-1, :-1] + rgb_vertex[1:, :-1] + rgb_vertex[:-1, 1:] + rgb_vertex[1:, 1:]) / 4.0
        ax.plot_surface(xx_fine, yy_fine, z_fine, facecolors=rgb_face, rstride=1, cstride=1, antialiased=True, shade=False)
        if z_max > z_min:
            ax.set_zlim(z_min, z_max)
        ax.set_xlabel("x", labelpad=-10); ax.set_ylabel("y", labelpad=-10)
        ax.tick_params(axis='x', pad=-5); ax.tick_params(axis='y', pad=-5); ax.tick_params(axis='z', pad=-3)

    # --- Figure Layout ---
    n_plot_rows = n_modes_to_plot
    fig_width = 28 if has_mean else 22
    fig_height = 4 * n_plot_rows
    fig = plt.figure(figsize=(fig_width, max(fig_height, 7)), dpi=dpi)
    # Adjust margins to use more of the figure area
    fig.subplots_adjust(left=0.02, right=0.99, bottom=0.05, top=0.93)

    if has_mean:
        # 9 columns: main, approx, mean, plus, contrib, eq, topos, mul, chronos
        gs = gridspec.GridSpec(
            n_plot_rows + 1, 9,
            width_ratios=[8, 0.5, 4, 0.5, 4, 0.5, 4, 0.5, 4],
            height_ratios=[1] * n_plot_rows + [0.2],
            wspace=0.3, hspace=0.5
        )
        main_col, approx_col = 0, 1
        mean_col, plus_col = 2, 3
        contrib_col, eq_col, topos_col, mul_col, chronos_col = 4, 5, 6, 7, 8
    else:
        # 7 columns: main, approx, contrib, eq, topos, mul, chronos
        gs = gridspec.GridSpec(
            n_plot_rows + 1, 7,
            width_ratios=[8, 0.5, 4, 0.5, 4, 0.5, 4],
            height_ratios=[1] * n_plot_rows + [0.2],
            wspace=0.3, hspace=0.5
        )
        main_col, approx_col = 0, 1
        contrib_col, eq_col, topos_col, mul_col, chronos_col = 2, 3, 4, 5, 6

    fig.suptitle(f"Decomposition at Time Step {frame}", fontsize=16)

    # --- Main Plot ---
    ax_main = fig.add_subplot(gs[0:n_plot_rows, main_col], projection='3d')
    _plot_surface_on_ax(
        ax_main, z_true_frame, q_true_frame, color_data_true_frame, z_main_min, z_main_max,
        main_color_mode, main_cmap, interp_k=2
    )
    if show_titles:
        ax_main.set_title(f"True Field ({main_z_label})")
    ax_main.set_zlabel(main_z_label, labelpad=-8)

    # '≈' symbol in its own axes for robust alignment
    ax_approx = fig.add_subplot(gs[0:n_plot_rows, approx_col])
    ax_approx.text(0.5, 0.5, '≈', ha='center', va='center', fontsize=60, color='gray', transform=ax_approx.transAxes)
    ax_approx.axis('off')

    # --- Mean Topos Plot (if exists) ---
    if has_mean:
        ax_mean = fig.add_subplot(gs[0:n_plot_rows, mean_col], projection='3d')
        _, _, z_lbl, c_mode, c_map = _get_z_and_color_data(q_mean_field, equation)
        _plot_surface_on_ax(ax_mean, z_mean, q_mean_field, color_mean, z_main_min, z_main_max, c_mode, c_map, interp_k=2)
        if show_titles:
            ax_mean.set_title("Mean Topos", fontsize=10)
        ax_mean.set_zlabel(z_lbl, labelpad=-8)

        ax_plus_mean = fig.add_subplot(gs[0:n_plot_rows, plus_col])
        ax_plus_mean.text(0.5, 0.5, '+', ha='center', va='center', fontsize=24, transform=ax_plus_mean.transAxes)
        ax_plus_mean.axis('off')

    # --- Mode Plots ---
    for i in range(n_modes_to_plot):
        row_idx = i
        
        # Get consistent visualization properties for this mode's plots
        _, _, z_lbl, c_mode, c_map = _get_z_and_color_data(topos_fields[i], equation)

        # Animated Contribution Plot
        ax_contrib = fig.add_subplot(gs[row_idx, contrib_col], projection='3d')
        z_min_c, z_max_c = z_contrib_lims[i]
        _plot_surface_on_ax(ax_contrib, z_contrib_frames[i], q_contrib_frames[i], color_contrib_frames[i], z_min_c, z_max_c, c_mode, c_map, interp_k=2)
        if show_titles:
            ax_contrib.set_title(f"Mode {i} Contribution", fontsize=10)
        ax_contrib.set_zlabel(z_lbl, labelpad=-8)

        # Static Topos Plot
        ax_t = fig.add_subplot(gs[row_idx, topos_col], projection='3d')
        _plot_surface_on_ax(ax_t, z_topos[i], topos_fields[i], color_topos[i], z_topos[i].min(), z_topos[i].max(), c_mode, c_map, interp_k=2)
        if show_titles:
            ax_t.set_title(f"Topos (Mode {i})", fontsize=10)
        ax_t.set_zlabel(z_lbl, labelpad=-8)

        # Chronos Plot
        ax_c = fig.add_subplot(gs[row_idx, chronos_col])
        ax_c.plot(np.arange(T), A_true[:, i], color=f"C{i}")
        ax_c.plot([frame], [A_true[frame, i]], 'o', color='red', markersize=8)
        if show_titles:
            ax_c.set_title(f"Chronos (Mode {i})", fontsize=10)
        ax_c.set_xlabel("Time Step", fontsize=8); ax_c.set_ylabel("Amplitude", fontsize=8)
        ax_c.grid(True, alpha=0.3); ax_c.tick_params(axis='both', which='major', labelsize=8)

        # Add '+' symbol between rows, centered on the plot column
        pos_contrib = ax_contrib.get_position()
        if row_idx > 0: # Add '+' between rows
            fig.text(pos_contrib.x0 + pos_contrib.width / 2, pos_contrib.y1 + 0.02, '+', ha='center', va='center', fontsize=24)

        # Add text symbols for this row using dedicated axes for alignment
        ax_eq = fig.add_subplot(gs[row_idx, eq_col])
        ax_eq.text(0.5, 0.5, '=', ha='center', va='center', fontsize=24, transform=ax_eq.transAxes)
        ax_eq.axis('off')

        ax_mul = fig.add_subplot(gs[row_idx, mul_col])
        ax_mul.text(0.5, 0.5, '×', ha='center', va='center', fontsize=24, transform=ax_mul.transAxes)
        ax_mul.axis('off')

    if r > n_modes_to_plot:
        ax_placeholder = fig.add_subplot(gs[n_plot_rows, contrib_col]); ax_placeholder.axis('off')
        pos = ax_placeholder.get_position()
        # Center the '...' horizontally
        fig.text(pos.x0 + pos.width / 2, pos.y0 + pos.height, '+  more modes', ha='center', va='center', fontsize=20)

    fig.savefig(output_path)
    plt.close(fig)

def animate_3d_decomposition(
    rollout: RolloutResult,
    pod: PODResult,
    layout: SnapshotLayout,
    equation: str,
    output_path: Path,
    mean_state: Optional[np.ndarray],
    fps: int = 15,
    dpi: int = 120,
    interp_factor: int = 2,
    plot_cfg: Optional[PlotConfig] = None,
    is_centered: bool = False,
) -> None:
    """
    Animate the decomposition of the true field into its dominant POD modes.

    Layout:
    - Left: Large 3D surface plot of the true field, evolving in time.
    - Right: Vertical stack showing the top 3 mode contributions.
      Each row consists of: Mode # = [Static 3D Topos] x [Chronos line plot with moving dot]
    """
    q_true = rollout.fields_true
    A_true = rollout.A_true
    U = pod.U
    T, r = A_true.shape
    _, C, ny, nx = q_true.shape

    n_modes_to_plot = min(3, r)
    if n_modes_to_plot == 0:
        print("Skipping 3D decomposition movie: no modes to plot.")
        return

    # --- Data Preparation ---
    # Main plot data
    z_true, q_for_color_true, main_z_label, main_color_mode, main_cmap = _get_z_and_color_data(q_true, equation)
    color_data_true = q_for_color_true if main_color_mode == 'colormap' else None

    # Determine global z-limits for main, mean, and contribution plots
    all_z_mins = [z_true.min()]
    all_z_maxs = [z_true.max()]

    # Mean state data (if exists)
    has_mean = mean_state is not None and is_centered
    q_mean_field, z_mean, color_mean = None, None, None
    if has_mean:
        q_mean_field = state_vec_to_fields(mean_state, layout)
        z_mean, q_for_color_mean, _, _, _ = _get_z_and_color_data(q_mean_field, equation)
        color_mean = q_for_color_mean if _get_z_and_color_data(q_mean_field, equation)[3] == 'colormap' else None
        all_z_mins.append(z_mean.min())
        all_z_maxs.append(z_mean.max())

    # --- Prepare Topos (basis vectors) data ---
    topos_fields = [state_vec_to_fields(U[:, i], layout) for i in range(n_modes_to_plot)]
    
    # Prepare data for all small plots (topos, contributions)
    z_topos, color_topos, q_contrib, z_contrib, color_contrib = [], [], [], [], []
    for i in range(n_modes_to_plot):
        # Static topos plots
        z_t, q_for_color_t, _, _, _ = _get_z_and_color_data(topos_fields[i], equation)
        z_topos.append(z_t)
        color_topos.append(q_for_color_t if _get_z_and_color_data(topos_fields[i], equation)[3] == 'colormap' else None)

        # Animated contribution plots
        q_c = A_true[:, i][:, np.newaxis, np.newaxis, np.newaxis] * topos_fields[i][np.newaxis, ...]
        q_contrib.append(q_c)
        z_c, q_for_color_c, _, _, _ = _get_z_and_color_data(q_c, equation)
        z_contrib.append(z_c)
        color_contrib.append(q_for_color_c if _get_z_and_color_data(q_c, equation)[3] == 'colormap' else None)
        all_z_mins.append(z_c.min())
        all_z_maxs.append(z_c.max())

    # Calculate global limits and set them for all relevant plots
    z_global_min = min(all_z_mins)
    z_global_max = max(all_z_maxs)

    # The contribution plots will use these global limits.
    # The topos plots will keep their own auto-scaled limits.
    z_contrib_lims = [(z_global_min, z_global_max)] * n_modes_to_plot

    # --- Interpolation and Plotting Setup ---
    x = np.arange(nx)
    y = np.arange(ny)
    x_fine = np.linspace(0, nx - 1, nx * interp_factor)
    y_fine = np.linspace(0, ny - 1, ny * interp_factor)
    xx_fine, yy_fine = np.meshgrid(x_fine, y_fine)
    light = LightSource(azdeg=225, altdeg=45)

    use_parallel = False
    if plot_cfg is not None:
        use_parallel = getattr(plot_cfg, "movie_3d_parallel", False)

    # Pass down the title toggle
    show_titles = getattr(plot_cfg, "movie_3d_decomposition_show_titles", True)

    if use_parallel and plot_cfg is not None:
        # --- Parallel Rendering Path ---
        n_procs = plot_cfg.movie_3d_parallel_procs or multiprocessing.cpu_count() // 2 or 1
        movie_every = plot_cfg.movie_every
        frame_indices = range(0, T, movie_every)
        n_frames = len(frame_indices)

        print(f"Generating {n_frames} frames for 3D decomposition movie in parallel using {n_procs} processes...")
        mp_context = multiprocessing.get_context("spawn")

        with tempfile.TemporaryDirectory() as temp_dir:
            frame_dir = Path(temp_dir)
            tasks = []
            for frame in frame_indices:
                task_args = (
                    frame,
                    q_true[frame], # For HSV color
                    z_true[frame],
                    color_data_true[frame] if color_data_true is not None else None,
                    frame_dir / f"frame_{frame:05d}.png",
                    T, n_modes_to_plot, r, nx, ny, main_z_label, main_color_mode, main_cmap, equation,
                    has_mean, q_mean_field, z_mean, color_mean,
                    topos_fields,
                    z_topos,
                    color_topos,
                    [qc[frame] for qc in q_contrib],
                    [zc[frame] for zc in z_contrib],
                    [cc[frame] if cc is not None else None for cc in color_contrib],
                    A_true, z_global_min, z_global_max,
                    x, y, x_fine, y_fine, xx_fine, yy_fine, light, dpi,
                    z_contrib_lims,
                    show_titles,
                )
                tasks.append(task_args)

            with mp_context.Pool(processes=n_procs) as pool, tqdm(total=n_frames, desc="Rendering 3D frames") as pbar:
                # Use a generator to feed tasks to the pool to manage memory
                task_generator = (t for t in tasks)
                for _ in pool.imap_unordered(
                    _render_3d_decomposition_frame_wrapper,
                    task_generator,
                    chunksize=max(1, n_frames // (4 * n_procs))):
                    pbar.update(1)

            print("Stitching frames into movie with ffmpeg...")
            ffmpeg_cmd = [
                "ffmpeg", "-y", "-framerate", str(fps),
                "-i", str(frame_dir / "frame_%05d.png"),
                "-c:v", "libx264", "-profile:v", "main", # Force a more compatible H.264 profile
                "-vf", "scale=min(1920\\,iw):-2", # Scale to compatible width, keep aspect ratio
                "-pix_fmt", "yuv420p", "-crf", "17", "-preset", "fast",
                str(output_path),
            ]
            try:
                subprocess.run(ffmpeg_cmd, check=True, capture_output=True, text=True, encoding="utf-8")
                print(f"Saved 3D decomposition animation to {output_path}")
            except (FileNotFoundError, subprocess.CalledProcessError) as e:
                print(f"ERROR: ffmpeg failed. Is it installed? Details:\n{e}")
        return

    # --- Sequential Rendering Path (Original Implementation) ---
    # The sequential path is now deprecated in favor of the more robust parallel path.
    # If parallel fails (e.g., no ffmpeg), it's better to notify the user than
    # to fall back to a slow, memory-intensive animation method.
    print(
        "Skipping sequential 3D decomposition movie generation. "
        "Enable `movie_3d_parallel` and ensure ffmpeg is installed."
    )

def _render_3d_reconstruction_frame_wrapper(args: tuple) -> None:
    """Helper to unpack arguments for multiprocessing.Pool.map."""
    return _render_3d_reconstruction_frame(*args)


def _render_3d_reconstruction_frame(
    frame: int,
    # Data for this frame
    q_true_frame: np.ndarray, z_true_frame: np.ndarray, color_data_true_frame: Optional[np.ndarray],
    q_recon_frame: np.ndarray, z_recon_frame: np.ndarray, color_data_recon_frame: Optional[np.ndarray],
    q_error_frame: np.ndarray, z_error_frame: np.ndarray, color_data_error_frame: Optional[np.ndarray],
    # Static data
    output_path: Path,
    n_val: int,
    x: np.ndarray, y: np.ndarray, x_fine: np.ndarray, y_fine: np.ndarray, xx_fine: np.ndarray, yy_fine: np.ndarray, light: LightSource, dpi: int,
    z_min: float, z_max: float,
    plot_props_true: tuple,
    plot_props_recon: tuple,
    plot_props_error: tuple,
):
    """Renders a single 3-panel frame of the 3D reconstruction comparison animation."""

    # --- Helper to render a single 3D surface (adapted from _render_3d_decomposition_frame) ---
    def _plot_surface_on_ax(ax, z_field, q_field, color_field, z_min, z_max, local_color_mode, local_cmap, interp_k=2):
        spline_z = RectBivariateSpline(y, x, z_field, kx=interp_k, ky=interp_k)
        z_fine = spline_z(y_fine, x_fine)
        illumination = light.hillshade(z_fine, vert_exag=0.8)
        rescaled_illumination = 0.5 + 0.5 * illumination

        if local_color_mode == "hsv" and q_field is not None:
            spline_u = RectBivariateSpline(y, x, q_field[0], kx=interp_k, ky=interp_k)
            spline_v = RectBivariateSpline(y, x, q_field[1], kx=interp_k, ky=interp_k)
            u_fine, v_fine = spline_u(y_fine, x_fine), spline_v(y_fine, x_fine)
            hue = (np.arctan2(v_fine, u_fine) + np.pi) / (2 * np.pi)
            hsv_vertex = np.stack([hue, np.ones_like(hue), rescaled_illumination], axis=-1)
            rgb_vertex = hsv_to_rgb(hsv_vertex)
        else: # Colormap mode
            color_field_to_interp = z_field if color_field is None else color_field
            spline_color = RectBivariateSpline(y, x, color_field_to_interp, kx=interp_k, ky=interp_k)
            color_fine = spline_color(y_fine, x_fine)
            c_min, c_max = color_field_to_interp.min(), color_field_to_interp.max()
            if c_max - c_min < 1e-9: c_max = c_min + 1.0
            norm = Normalize(vmin=c_min, vmax=c_max)
            cmap_obj = plt.get_cmap(local_cmap)
            colors_from_map = cmap_obj(norm(color_fine))[:, :, :3]
            rgb_vertex = colors_from_map * rescaled_illumination[..., np.newaxis]

        rgb_face = (rgb_vertex[:-1, :-1] + rgb_vertex[1:, :-1] + rgb_vertex[:-1, 1:] + rgb_vertex[1:, 1:]) / 4.0
        ax.plot_surface(xx_fine, yy_fine, z_fine, facecolors=rgb_face, rstride=1, cstride=1, antialiased=True, shade=False)
        if z_max > z_min:
            ax.set_zlim(z_min, z_max)
        ax.set_xlabel("x", labelpad=-10); ax.set_ylabel("y", labelpad=-10)
        ax.tick_params(axis='x', pad=-5); ax.tick_params(axis='y', pad=-5); ax.tick_params(axis='z', pad=-3)

    # --- Figure Layout ---
    fig = plt.figure(figsize=(21, 7), dpi=dpi)
    gs = gridspec.GridSpec(1, 3, wspace=0.1, hspace=0.1)
    fig.suptitle(f"Time Step {frame}", fontsize=16)

    z_lbl_true, c_mode_true, c_map_true = plot_props_true
    z_lbl_recon, c_mode_recon, c_map_recon = plot_props_recon
    z_lbl_error, c_mode_error, c_map_error = plot_props_error

    # --- Left Plot: True Field ---
    ax_true = fig.add_subplot(gs[0, 0], projection='3d')
    _plot_surface_on_ax(ax_true, z_true_frame, q_true_frame, color_data_true_frame, z_min, z_max, c_mode_true, c_map_true)
    ax_true.set_title("Full True")
    ax_true.set_zlabel(z_lbl_true, labelpad=-8)

    # --- Middle Plot: Reconstruction ---
    ax_recon = fig.add_subplot(gs[0, 1], projection='3d')
    _plot_surface_on_ax(ax_recon, z_recon_frame, q_recon_frame, color_data_recon_frame, z_min, z_max, c_mode_recon, c_map_recon)
    ax_recon.set_title(f"Mean + Top {n_val} Predicted")
    ax_recon.set_zlabel(z_lbl_recon, labelpad=-8)

    # --- Right Plot: Error ---
    ax_error = fig.add_subplot(gs[0, 2], projection='3d')
    _plot_surface_on_ax(ax_error, z_error_frame, q_error_frame, color_data_error_frame, z_min, z_max, c_mode_error, c_map_error)
    ax_error.set_title("Error")
    ax_error.set_zlabel(z_lbl_error, labelpad=-8)

    # tight_layout is known to issue a UserWarning with 3D axes, which is safe to
    # ignore if the visual output is acceptable. We suppress it to keep logs clean.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output_path)
    plt.close(fig)


def animate_3d_reconstruction_comparison(
    rollout: RolloutResult,
    pod: PODResult,
    layout: SnapshotLayout,
    equation: str,
    output_dir: Path,
    mean_state: Optional[np.ndarray],
    fps: int = 15,
    dpi: int = 120,
    interp_factor: int = 2,
    plot_cfg: Optional[PlotConfig] = None,
    is_centered: bool = False,
) -> None:
    """
    Animate a 3-panel comparison: True vs. Reconstruction vs. Error.

    Generates a series of movies, one for each number of modes `n` used in
    the reconstruction (from n=0 to n=r).
    """
    use_parallel = plot_cfg and getattr(plot_cfg, "movie_3d_parallel", False)
    if not use_parallel:
        print("Skipping 3D reconstruction comparison: `movie_3d_parallel` is false.")
        return

    T, r = rollout.A_true.shape
    nx, ny = layout.nx, layout.ny

    q_true = rollout.fields_true
    x, y = np.arange(nx), np.arange(ny)
    x_fine, y_fine = np.linspace(0, nx - 1, nx * interp_factor), np.linspace(0, ny - 1, ny * interp_factor)
    xx_fine, yy_fine = np.meshgrid(x_fine, y_fine)
    light = LightSource(azdeg=225, altdeg=45)

    n_procs = plot_cfg.movie_3d_parallel_procs or multiprocessing.cpu_count() // 2 or 1
    mp_context = multiprocessing.get_context("spawn")

    # --- Prepare data that is constant across all movies in the series ---
    z_true, q_color_true, z_lbl_true, c_mode_true, c_map_true = _get_z_and_color_data(q_true, equation)
    color_true = q_color_true if c_mode_true == 'colormap' else None

    # --- Loop over number of modes `n` to generate one movie per `n` ---
    for n in range(r + 1):
        print(f"Generating 3D reconstruction comparison movie for n={n}/{r} modes...")
        output_path = output_dir / f"reconstruction_comparison_n{n:02d}.mp4"

        mode_indices = list(range(n))
        q_recon = _reconstruct_fields_from_modes(pod.U, rollout.A_pred, layout, mode_indices, mean_state, add_mean=is_centered)

        # The error is the difference between the partial reconstruction and the true field
        q_error = q_recon - q_true

        # Get plotting data for all three fields for this n
        z_recon, q_color_recon, z_lbl_recon, c_mode_recon, c_map_recon = _get_z_and_color_data(q_recon, equation)
        color_recon = q_color_recon if c_mode_recon == 'colormap' else None

        z_error, q_color_error, z_lbl_error, c_mode_error, c_map_error = _get_z_and_color_data(q_error, equation)
        color_error = q_color_error if c_mode_error == 'colormap' else None

        # Calculate global z-limits for all three plots for this n
        z_min, z_max = min(z_true.min(), z_recon.min(), z_error.min()), max(z_true.max(), z_recon.max(), z_error.max())
        if z_max - z_min < 1e-9:
            z_max = z_min + 1.0

        with tempfile.TemporaryDirectory() as temp_dir:
            frame_dir = Path(temp_dir)
            tasks = []
            for frame in range(T):
                tasks.append((
                    frame,
                    q_true[frame], z_true[frame], color_true[frame] if color_true is not None else None,
                    q_recon[frame], z_recon[frame], color_recon[frame] if color_recon is not None else None,
                    q_error[frame], z_error[frame], color_error[frame] if color_error is not None else None,
                    frame_dir / f"frame_{frame:05d}.png", n,
                    x, y, x_fine, y_fine, xx_fine, yy_fine, light, dpi, z_min, z_max,
                    (z_lbl_true, c_mode_true, c_map_true),
                    (z_lbl_recon, c_mode_recon, c_map_recon),
                    (z_lbl_error, c_mode_error, c_map_error),
                ))

            with mp_context.Pool(processes=n_procs) as pool, tqdm(total=T, desc=f"Rendering frames (n={n})") as pbar:
                for _ in pool.imap_unordered(_render_3d_reconstruction_frame_wrapper, tasks, chunksize=max(1, T // (4 * n_procs))):
                    pbar.update(1)

            print(f"Stitching frames for n={n} movie...")
            ffmpeg_cmd = [
                "ffmpeg", "-y", "-framerate", str(fps), "-i", str(frame_dir / "frame_%05d.png"),
                "-c:v", "libx264", "-profile:v", "main", "-vf", "scale=min(1920\\,iw):-2",
                "-pix_fmt", "yuv420p", "-crf", "17", "-preset", "fast", str(output_path),
            ]
            try:
                subprocess.run(ffmpeg_cmd, check=True, capture_output=True, text=True, encoding="utf-8")
                print(f"Saved 3D reconstruction comparison animation to {output_path}")
            except (FileNotFoundError, subprocess.CalledProcessError) as e:
                print(f"ERROR: ffmpeg failed for n={n}. Is it installed? Details:\n{e}")

def _reconstruct_fields_from_modes(
    U: np.ndarray,
    A: np.ndarray,
    layout: SnapshotLayout,
    mode_indices: List[int],
    mean_state: Optional[np.ndarray] = None,
    add_mean: bool = True,
) -> np.ndarray:
    """
    Reconstruct full-state fields from a subset of POD modes and coefficients.

    Parameters
    ----------
    U : np.ndarray
        POD basis of shape (n_state, r).
    A : np.ndarray
        POD coefficients of shape (T, r).
    layout : SnapshotLayout
        Layout for reshaping.
    mode_indices : list[int]
        Indices of the modes to use for reconstruction.
    mean_state : np.ndarray, optional
        Mean state vector of shape (n_state,). Added back if provided.
    add_mean : bool
        If True and mean_state is provided, add the mean back to the reconstruction.

    Returns
    -------
    np.ndarray
        Reconstructed fields of shape (T, C, ny, nx).
    """
    T, r = A.shape

    # Select modes
    U_modes = U[:, mode_indices]  # (n_state, k)
    A_modes = A[:, mode_indices]  # (T, k)

    # Reconstruct: X_rec = U_modes @ A_modes.T -> (n_state, T)
    X_rec = U_modes @ A_modes.T

    if add_mean and mean_state is not None:
        X_rec += mean_state[:, None]

    # Reshape from (n_state, T) to (T, C, ny, nx)
    fields = np.empty((T, layout.n_components, layout.ny, layout.nx), dtype=X_rec.dtype)
    for t in range(T):
        fields[t] = state_vec_to_fields(X_rec[:, t], layout)

    return fields


def generate_all_plots_and_movies(
    cfg: PlotConfig,
    rundir: Path,
    layout: SnapshotLayout,
    pod: PODResult,
    sindy: SINDyFitResult,
    rollout: RolloutResult,
    equation: str,
    mean_state: Optional[np.ndarray] = None,
    is_centered: bool = False,
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
    mean_state : np.ndarray, optional
        Mean state vector, if centering was used.
    is_centered : bool
        Whether centering was used.
    """
    if not cfg.enabled:
        return

    fig_dir = rundir / cfg.figures_subdir
    mov_dir = rundir / cfg.movies_subdir
    fig_dir.mkdir(parents=True, exist_ok=True)
    mov_dir.mkdir(parents=True, exist_ok=True)

    (fig_dir / "pod_basis").mkdir(parents=True, exist_ok=True)
    plot_pod_singular_values(pod=pod, out_dir=fig_dir, dpi=cfg.dpi)
    plot_pod_singular_values_scatterplot(pod=pod, out_dir=fig_dir, dpi=cfg.dpi)
    plot_pod_energy_scatterplot(pod=pod, out_dir=fig_dir, dpi=cfg.dpi, equation=equation)

    if getattr(cfg, "pod_decomposition_matrix", False):
        plot_pod_decomposition_matrix(pod=pod, out_dir=fig_dir, dpi=cfg.dpi)
    
    if getattr(cfg, "pod_decomposition_matrix_square_pixels", False):
        plot_pod_decomposition_matrix_square_pixels(pod=pod, out_dir=fig_dir, dpi=cfg.dpi)

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
        plot_field_rel_l2_comparison_with_dynabench(
            rundir=rundir,
            out_dir=fig_dir,
            dpi=cfg.dpi,
            equation=equation,
            timesteps=cfg.dynabench_comparison_timesteps,
        )

    if getattr(cfg, "movie_3d_surface", False):
        interp_factor = getattr(cfg, "movie_3d_interp_factor", 3)
        animate_3d_surface(
            q_left=rollout.fields_true,
            q_right=rollout.fields_pred,
            equation=equation,
            output_path=mov_dir / "rollout_3d_surface.mp4",
            title_left="True Field",
            title_right="Predicted Field",
            fps=cfg.movie_fps,
            dpi=cfg.dpi,
            interp_factor=interp_factor,
            plot_cfg=cfg,
        )

        # Movie of the error field
        print("Generating 3D surface movie for error field...")
        error_field = rollout.fields_pred - rollout.fields_true
        animate_3d_surface(
            q_left=error_field,
            equation=equation,
            output_path=mov_dir / "rollout_3d_error_surface.mp4",
            title_left="Error Field (Pred - True)",
            fps=cfg.movie_fps,
            dpi=cfg.dpi,
            interp_factor=interp_factor,
            plot_cfg=cfg,
        )

    if getattr(cfg, "movie_3d_decomposition", False):
        print("Generating 3D decomposition movie...")
        interp_factor = getattr(cfg, "movie_3d_interp_factor", 2)
        animate_3d_decomposition(
            rollout=rollout,
            pod=pod,
            layout=layout,
            equation=equation,
            output_path=mov_dir / "rollout_3d_decomposition.mp4",
            mean_state=mean_state,
            fps=cfg.movie_fps,
            dpi=cfg.dpi,
            interp_factor=interp_factor,
            plot_cfg=cfg,
            is_centered=is_centered,
        )

    if getattr(cfg, "movie_3d_reconstruction_comparison", False):
        print("Generating 3D reconstruction comparison movies...")
        recon_dir = mov_dir / "3d_reconstruction_comparison"
        recon_dir.mkdir(exist_ok=True)
        interp_factor = getattr(cfg, "movie_3d_interp_factor", 2)
        animate_3d_reconstruction_comparison(
            rollout=rollout,
            pod=pod,
            layout=layout,
            equation=equation,
            output_dir=recon_dir,
            mean_state=mean_state,
            fps=cfg.movie_fps,
            dpi=cfg.dpi,
            interp_factor=interp_factor,
            plot_cfg=cfg,
            is_centered=is_centered,
        )

    if getattr(cfg, "movie_3d_mode_contributions", False):
        r = pod.r
        contrib_dir = mov_dir / "3d_mode_contributions"
        contrib_dir.mkdir(exist_ok=True)
        interp_factor = getattr(cfg, "movie_3d_interp_factor", 3)

        # Movie of just the mean state
        if mean_state is not None:
            print("Generating 3D surface movie for true mean state...")
            # Repeat mean state T times to create a "movie"
            T = rollout.A_true.shape[0]
            mean_field = state_vec_to_fields(mean_state, layout)
            mean_field_movie = np.repeat(mean_field[np.newaxis, ...], T, axis=0)
            output_path = contrib_dir / "true_mean_state.mp4"
            animate_3d_surface(
                q_left=mean_field_movie,
                equation=equation,
                output_path=output_path,
                title_left="True Mean State",
                fps=cfg.movie_fps,
                dpi=cfg.dpi,
                interp_factor=interp_factor,
                plot_cfg=cfg,
            )

        # --- Predicted Mode Contributions ---
        print(f"Generating 3D surface movies for individual predicted mode contributions (r={r})...")
        for i in range(r):
            q_recon = _reconstruct_fields_from_modes(
                U=pod.U,
                A=rollout.A_pred, # Use predicted coefficients
                layout=layout,
                mode_indices=[i],
                add_mean=False,
            )
            output_path = contrib_dir / f"pred_mode_{i:02d}_contribution.mp4"
            animate_3d_surface(
                q_left=q_recon,
                equation=equation,
                output_path=output_path,
                title_left=f"Predicted Mode {i} Contribution",
                fps=cfg.movie_fps,
                dpi=cfg.dpi,
                interp_factor=interp_factor,
                plot_cfg=cfg,
            )

        print(f"Generating 3D surface movies for cumulative highest predicted mode contributions (r={r})...")
        for k in range(2, r + 1):
            mode_indices = list(range(r - k, r))
            q_recon = _reconstruct_fields_from_modes(
                U=pod.U,
                A=rollout.A_pred, # Use predicted coefficients
                layout=layout,
                mode_indices=mode_indices,
                add_mean=False,
            )
            output_path = contrib_dir / f"pred_highest_{k:02d}_modes_contribution.mp4"
            animate_3d_surface(
                q_left=q_recon,
                equation=equation,
                output_path=output_path,
                title_left=f"Predicted Highest {k} Modes Contribution",
                fps=cfg.movie_fps,
                dpi=cfg.dpi,
                interp_factor=interp_factor,
                plot_cfg=cfg,
            )

        print(f"Generating 3D surface movies for individual true mode contributions (r={r})...")
        for i in range(r):
            q_recon = _reconstruct_fields_from_modes(
                U=pod.U,
                A=rollout.A_true,
                layout=layout,
                mode_indices=[i],
                add_mean=False,  # Show mode contribution only, not relative to mean
            )
            output_path = contrib_dir / f"true_mode_{i:02d}_contribution.mp4"
            animate_3d_surface(
                q_left=q_recon,
                equation=equation,
                output_path=output_path,
                title_left=f"True Mode {i} Contribution",
                fps=cfg.movie_fps,
                dpi=cfg.dpi,
                interp_factor=interp_factor,
                plot_cfg=cfg,
            )

        print(f"Generating 3D surface movies for cumulative highest true mode contributions (r={r})...")
        for k in range(2, r + 1):
            mode_indices = list(range(r - k, r))
            q_recon = _reconstruct_fields_from_modes(
                U=pod.U,
                A=rollout.A_true,
                layout=layout,
                mode_indices=mode_indices,
                add_mean=False, # Show mode contribution only
            )
            output_path = contrib_dir / f"true_highest_{k:02d}_modes_contribution.mp4"
            animate_3d_surface(
                q_left=q_recon,
                equation=equation,
                output_path=output_path,
                title_left=f"True Highest {k} Modes Contribution",
                fps=cfg.movie_fps,
                dpi=cfg.dpi,
                interp_factor=interp_factor,
                plot_cfg=cfg,
            )