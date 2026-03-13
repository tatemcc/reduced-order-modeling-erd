"""Field/state helpers for the toy ERD plant.

This module defines initial conditions, profile helpers, and conversions between
FiPy variables and ``(N_r, N_phi)`` NumPy arrays used by the ROM pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import RunConfig
from .mesh import MeshBundle

try:
    from fipy import CellVariable
except Exception:  # pragma: no cover - runtime dependency
    CellVariable = None


@dataclass
class ERDState:
    """Live FiPy field variables for the toy ERD state.

    Attributes:
        n: Density-like field variable.
        omega: Vorticity-like field variable.
        psi: Streamfunction field variable.
    """

    n: object
    omega: object
    psi: object


def ring_equilibrium(r: np.ndarray, cfg: RunConfig) -> np.ndarray:
    """Evaluate the target ring profile ``n_eq(r)``.

    Args:
        r: Radial coordinates.
        cfg: Plant configuration containing ring parameters.

    Returns:
        Axisymmetric equilibrium profile sampled at ``r``.
    """

    rc = cfg.ring
    return rc.n_bg + rc.n_amp * np.exp(-((r - rc.r_star) ** 2) / (2.0 * rc.sigma_star**2))


def wall_loss_profile(r: np.ndarray, cfg: RunConfig) -> np.ndarray:
    """Evaluate wall-loss layer profile ``kappa_w(r)``.

    Args:
        r: Radial coordinates.
        cfg: Plant configuration containing wall-layer parameters.

    Returns:
        Nonnegative sink profile sampled at ``r``.
    """

    wc = cfg.wall
    return wc.kappa_0 * np.exp(-((cfg.domain.R_max - r) ** 2) / (2.0 * wc.delta_w**2))


def _reshape_to_grid(values: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """Reshape flattened cell data into ``(N_r, N_phi)`` grid form."""

    nr, nphi = shape
    return np.asarray(values).reshape(nr, nphi)


def variable_to_grid(var: object, shape: tuple[int, int]) -> np.ndarray:
    """Convert a FiPy ``CellVariable`` to a 2D grid array.

    Args:
        var: FiPy cell variable.
        shape: Grid shape ``(N_r, N_phi)``.

    Returns:
        NumPy array with shape ``shape``.
    """

    return _reshape_to_grid(np.asarray(var.value), shape)


def assign_grid(var: object, arr: np.ndarray) -> None:
    """Write a 2D grid array into a FiPy ``CellVariable``.

    Args:
        var: FiPy cell variable to update.
        arr: Grid array with shape ``(N_r, N_phi)``.

    Returns:
        None.
    """

    var.setValue(np.asarray(arr).reshape(-1))


def stack_state(n: np.ndarray, omega: np.ndarray) -> np.ndarray:
    """Stack ``n`` and ``omega`` grids into one full-order state vector.

    Args:
        n: Density field array.
        omega: Vorticity field array.

    Returns:
        Vector ``[vec(n), vec(omega)]``.
    """

    return np.concatenate([n.reshape(-1), omega.reshape(-1)], axis=0)


def unstack_state(x: np.ndarray, shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    """Recover ``(n, omega)`` grid arrays from a stacked state vector.

    Args:
        x: Stacked state vector ``[vec(n), vec(omega)]``.
        shape: Grid shape ``(N_r, N_phi)``.

    Returns:
        Tuple ``(n, omega)`` each with shape ``shape``.
    """

    nr, nphi = shape
    n_points = nr * nphi
    n = x[:n_points].reshape(shape)
    omega = x[n_points : 2 * n_points].reshape(shape)
    return n, omega


def make_initial_state(bundle: MeshBundle, cfg: RunConfig) -> ERDState:
    """Create initial FiPy variables using ring profile plus small random perturbations.

    Args:
        bundle: Mesh/geometric metadata for array shapes and coordinates.
        cfg: Plant run configuration including initialization amplitudes and seed.

    Returns:
        Initialized :class:`ERDState` containing ``n``, ``omega``, and ``psi``.
    """

    if CellVariable is None:
        raise RuntimeError("FiPy is required but not installed.")

    rng = np.random.default_rng(cfg.init.seed)

    n_eq_r = ring_equilibrium(bundle.r, cfg)
    n_eq = np.repeat(n_eq_r[:, None], bundle.shape[1], axis=1)

    xi_n = rng.standard_normal(bundle.shape)
    xi_n = xi_n - np.mean(xi_n)
    xi_w = rng.standard_normal(bundle.shape)
    xi_w = xi_w - np.mean(xi_w)

    phase1, phase2, phase5 = rng.uniform(0.0, 2.0 * np.pi, size=3)
    phi = bundle.phi[None, :]
    coherent = (
        cfg.init.mode1_amp * np.cos(phi + phase1)
        + cfg.init.mode2_amp * np.cos(2.0 * phi + phase2)
        + cfg.init.mode5_amp * np.cos(5.0 * phi + phase5)
    )
    ring_env = np.exp(-((bundle.r - cfg.ring.r_star) ** 2) / (2.0 * (1.15 * cfg.ring.sigma_star) ** 2))[:, None]

    n0 = n_eq * (1.0 + cfg.init.eps_n * xi_n + ring_env * coherent)
    n0 = np.clip(n0, 1e-6, None)
    omega0 = cfg.init.eps_omega * xi_w
    psi0 = np.zeros(bundle.shape, dtype=float)

    n_var = CellVariable(name="n", mesh=bundle.mesh, value=n0.reshape(-1))
    omega_var = CellVariable(name="omega", mesh=bundle.mesh, value=omega0.reshape(-1))
    psi_var = CellVariable(name="psi", mesh=bundle.mesh, value=psi0.reshape(-1))

    return ERDState(n=n_var, omega=omega_var, psi=psi_var)
