"""
POD utilities based on SVD.

This module computes a POD basis from a snapshot matrix and returns:
- Truncated left singular vectors U_r
- Singular values s
- Time coefficients A such that X ≈ U_r @ A
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass(frozen=True)
class PODResult:
    """
    Result of a POD computation.

    Attributes
    ----------
    U : np.ndarray
        POD spatial modes of shape (n_state, r).
    s : np.ndarray
        Singular values of length min(n_state, M).
    A : np.ndarray
        POD coefficients of shape (r, M) satisfying X ≈ U @ A.
    r : int
        Selected rank.
    energy_fraction : float
        Fraction of snapshot energy captured by the selected rank.
    """

    U: np.ndarray
    s: np.ndarray
    A: np.ndarray
    r: int
    energy_fraction: float


def select_rank(
    s: np.ndarray,
    rank: Optional[int],
    energy_fraction: Optional[float],
) -> int:
    """
    Select POD rank from singular values.

    Parameters
    ----------
    s : np.ndarray
        Singular values.
    rank : int, optional
        Fixed rank selection.
    energy_fraction : float, optional
        Target cumulative energy fraction in (0, 1].

    Returns
    -------
    int
        Selected rank r.
    """
    if rank is None and energy_fraction is None:
        raise ValueError("Either rank or energy_fraction must be provided")

    if rank is not None and energy_fraction is not None:
        raise ValueError("Provide only one of rank or energy_fraction")

    if rank is not None:
        if rank <= 0:
            raise ValueError("rank must be positive")
        return min(rank, s.shape[0])

    ef = float(energy_fraction)
    if not (0.0 < ef <= 1.0):
        raise ValueError("energy_fraction must be in (0, 1]")

    energy = s**2
    cum = np.cumsum(energy)
    total = cum[-1]
    if total <= 0:
        return 1

    r = int(np.searchsorted(cum / total, ef, side="left") + 1)
    return max(1, min(r, s.shape[0]))


def compute_pod(
    X: np.ndarray,
    rank: Optional[int],
    energy_fraction: Optional[float],
) -> PODResult:
    """
    Compute POD via SVD on the snapshot matrix X.

    Convention:
        X has shape (n_state, M)
        X ≈ U_r @ A
    where:
        U_r has shape (n_state, r)
        A has shape (r, M)

    Parameters
    ----------
    X : np.ndarray
        Snapshot matrix of shape (n_state, M).
    rank : int, optional
        Fixed rank r.
    energy_fraction : float, optional
        Target cumulative energy fraction.

    Returns
    -------
    PODResult
        POD basis and coefficients.
    """
    X = np.asarray(X)
    if X.ndim != 2:
        raise ValueError("X must be 2D of shape (n_state, M)")

    U, s, Vt = np.linalg.svd(X, full_matrices=False)

    r = select_rank(s, rank=rank, energy_fraction=energy_fraction)
    U_r = U[:, :r]

    # X ≈ U_r @ (diag(s_r) @ Vt_r)
    s_r = s[:r]
    Vt_r = Vt[:r, :]
    A = s_r[:, None] * Vt_r

    energy = s**2
    captured = float(np.sum(energy[:r]))
    total = float(np.sum(energy))
    frac = 0.0 if total <= 0 else captured / total

    return PODResult(U=U_r, s=s, A=A, r=r, energy_fraction=frac)


def reconstruct_from_pod(
    U: np.ndarray,
    A: np.ndarray,
    mean_state: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Reconstruct snapshots from POD basis and coefficients.

    Parameters
    ----------
    U : np.ndarray
        POD basis of shape (n_state, r).
    A : np.ndarray
        Coefficients of shape (r, M).
    mean_state : np.ndarray, optional
        Mean state vector of shape (n_state,). Added back if provided.

    Returns
    -------
    np.ndarray
        Reconstructed snapshot matrix of shape (n_state, M).
    """
    X_hat = U @ A
    if mean_state is not None:
        X_hat = X_hat + mean_state[:, None]
    return X_hat
