"""Mesh construction for the toy ERD annulus on a rectangularized ``(r, phi)`` grid.

The FiPy mesh uses ``x := phi`` and ``y := r`` so cell arrays can be reshaped to
``(N_r, N_phi)`` consistently throughout the pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import RunConfig

try:
    from fipy import Grid2D
except Exception:  # pragma: no cover - runtime dependency
    Grid2D = None


@dataclass(frozen=True)
class MeshBundle:
    """Discrete mesh and geometry arrays used by the plant and metrics.

    Attributes:
        mesh: FiPy 2D grid object with ``x := phi`` and ``y := r``.
        phi: 1D azimuthal cell-center coordinates.
        r: 1D radial cell-center coordinates.
        dphi: Uniform azimuthal spacing.
        dr: Uniform radial spacing.
        area_weights: Per-cell area weights ``r * dr * dphi`` with shape ``(N_r, N_phi)``.
        shape: Convenience tuple ``(N_r, N_phi)``.
    """

    mesh: object
    phi: np.ndarray
    r: np.ndarray
    dphi: float
    dr: float
    area_weights: np.ndarray
    shape: tuple[int, int]


def make_mesh(cfg: RunConfig) -> MeshBundle:
    """Create the rectangularized annular mesh used by the toy PDE system.

    Args:
        cfg: Plant run configuration containing domain resolution and bounds.

    Returns:
        A :class:`MeshBundle` with FiPy mesh and coordinate helpers.
    """

    if Grid2D is None:
        raise RuntimeError("FiPy is required but not installed.")

    nr = cfg.domain.N_r
    nphi = cfg.domain.N_phi

    dphi = 2.0 * np.pi / nphi
    dr = (cfg.domain.R_max - cfg.domain.R_min) / nr

    # Grid convention: x := phi, y := r
    mesh = Grid2D(dx=dphi, dy=dr, nx=nphi, ny=nr)

    phi = (np.arange(nphi) + 0.5) * dphi
    r = cfg.domain.R_min + (np.arange(nr) + 0.5) * dr

    area_weights = (r[:, None] * dr * dphi).astype(float)

    return MeshBundle(
        mesh=mesh,
        phi=phi,
        r=r,
        dphi=float(dphi),
        dr=float(dr),
        area_weights=area_weights,
        shape=(nr, nphi),
    )
