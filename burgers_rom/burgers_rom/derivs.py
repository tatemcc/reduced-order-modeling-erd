"""
Time-derivative utilities for POD coefficient time series.

This module relies on PySINDy differentiation tools and provides:
- Finite-difference derivatives for coefficient trajectories
- Safe handling of multiple trajectories by not differencing across boundaries
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np
from pysindy.differentiation import FiniteDifference


@dataclass(frozen=True)
class DerivativeResult:
    """
    Container for coefficient derivatives.

    Attributes
    ----------
    A_used : np.ndarray
        Coefficients aligned with dA_dt. Shape (r, M_used).
    dA_dt : np.ndarray
        Time derivative of coefficients. Shape (r, M_used).
    used_indices : np.ndarray
        Indices into the original time axis used for A_used. Shape (M_used,).
    """

    A_used: np.ndarray
    dA_dt: np.ndarray
    used_indices: np.ndarray


def _segment_indices(
    segment_lengths: Sequence[int],
) -> List[np.ndarray]:
    """
    Build index arrays for each trajectory segment.

    Parameters
    ----------
    segment_lengths : sequence of int
        Length of each segment in number of time steps.

    Returns
    -------
    list of np.ndarray
        Each array contains the indices for one segment in concatenated time.
    """
    starts: List[int] = []
    s = 0
    for L in segment_lengths:
        starts.append(s)
        s += int(L)

    segs: List[np.ndarray] = []
    for start, L in zip(starts, segment_lengths, strict=True):
        segs.append(np.arange(start, start + int(L), dtype=int))
    return segs


def finite_difference_coeff_derivative(
    A: np.ndarray,
    dt: float,
    segment_lengths: Optional[Sequence[int]] = None,
    order: int = 2,
) -> DerivativeResult:
    """
    Compute dA_dt using PySINDy finite differences.

    Parameters
    ----------
    A : np.ndarray
        POD coefficients of shape (r, M).
    dt : float
        Time step between snapshots.
    segment_lengths : sequence of int, optional
        If provided, A is treated as concatenated segments and derivatives
        are computed within each segment only.
    order : int
        Finite difference order. Must be >= 1.

    Returns
    -------
    DerivativeResult
        dA_dt aligned with A_used, excluding invalid boundary points.
    """
    A = np.asarray(A)
    if A.ndim != 2:
        raise ValueError("A must have shape (r, M)")
    if dt <= 0:
        raise ValueError("dt must be positive")
    if order < 1:
        raise ValueError("order must be >= 1")

    r, M = A.shape
    fd = FiniteDifference(order=order)

    if segment_lengths is None:
        t = np.arange(M, dtype=float) * dt
        dA = fd._differentiate(A.T, t).T

        k = order
        used = np.arange(k, M - k, dtype=int)
        return DerivativeResult(A_used=A[:, used], dA_dt=dA[:, used], used_indices=used)

    segs = _segment_indices(segment_lengths)

    A_used_list: List[np.ndarray] = []
    dA_used_list: List[np.ndarray] = []
    used_idx_list: List[np.ndarray] = []

    k = order
    for seg in segs:
        if seg.size < 2 * k + 1:
            continue

        t_seg = np.arange(seg.size, dtype=float) * dt
        A_seg = A[:, seg]
        dA_seg = fd._differentiate(A_seg.T, t_seg).T

        used_local = np.arange(k, seg.size - k, dtype=int)
        used_global = seg[used_local]

        A_used_list.append(A[:, used_global])
        dA_used_list.append(dA_seg[:, used_local])
        used_idx_list.append(used_global)

    if not A_used_list:
        raise ValueError("No segment long enough for the requested finite difference order")

    A_used = np.concatenate(A_used_list, axis=1)
    dA_dt = np.concatenate(dA_used_list, axis=1)
    used_indices = np.concatenate(used_idx_list, axis=0)

    return DerivativeResult(A_used=A_used, dA_dt=dA_dt, used_indices=used_indices)
