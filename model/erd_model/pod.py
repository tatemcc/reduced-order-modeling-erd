"""Proper Orthogonal Decomposition (POD) helpers for stacked ERD states."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass(frozen=True)
class PODResult:
    """POD decomposition outputs used by the ROM pipeline.

    Attributes:
        U: Retained POD basis vectors with shape ``(n_state, r)``.
        s: Full singular-value spectrum of the snapshot matrix.
        A: Reduced coefficients for all training samples with shape ``(r, M)``.
        r: Retained POD rank.
        energy_fraction: Fraction of total snapshot energy retained by rank ``r``.
    """

    U: np.ndarray
    s: np.ndarray
    A: np.ndarray
    r: int
    energy_fraction: float


def compute_pod(X: np.ndarray, rank: Optional[int]) -> PODResult:
    """Compute POD basis and reduced coefficients from a snapshot matrix.

    Args:
        X: Snapshot matrix with shape ``(n_state, M)``.
        rank: Requested retained rank. ``None`` defaults to 6 (clipped to range).

    Returns:
        :class:`PODResult` containing basis, coefficients, and energy stats.
    """

    X = np.asarray(X)
    if X.ndim != 2:
        raise ValueError("X must be 2D")

    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    if rank is None:
        rank = min(6, U.shape[1])
    r = int(max(1, min(rank, U.shape[1])))

    U_r = U[:, :r]
    s_r = s[:r]
    A = s_r[:, None] * Vt[:r, :]

    energy = s**2
    frac = float(np.sum(energy[:r]) / np.sum(energy)) if np.sum(energy) > 0 else 0.0
    return PODResult(U=U_r, s=s, A=A, r=r, energy_fraction=frac)


def reconstruct(U: np.ndarray, A: np.ndarray, mean_state: Optional[np.ndarray] = None) -> np.ndarray:
    """Reconstruct full-order states from POD basis and reduced coefficients.

    Args:
        U: POD basis matrix with shape ``(n_state, r)``.
        A: Reduced coefficients with shape ``(r, M)``.
        mean_state: Optional centering vector added back after reconstruction.

    Returns:
        Reconstructed snapshot matrix with shape ``(n_state, M)``.
    """

    X = U @ A
    if mean_state is not None:
        X = X + mean_state[:, None]
    return X
