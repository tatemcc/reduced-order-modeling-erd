"""Finite-difference derivative helpers for reduced trajectories."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence

import numpy as np


@dataclass(frozen=True)
class DerivativeResult:
    """Central-difference samples used to fit controlled SINDy.

    Attributes:
        A_used: Reduced states used for regression, shape ``(M_used, r)``.
        dA_dt: Central-difference derivatives, shape ``(M_used, r)``.
        U_used: Controls aligned with ``A_used``, shape ``(M_used, 5)``.
        used_indices: Original flattened sample indices kept after differencing.
    """

    A_used: np.ndarray
    dA_dt: np.ndarray
    U_used: np.ndarray
    used_indices: np.ndarray


def finite_difference_with_controls(
    A: np.ndarray,
    U_controls: np.ndarray,
    dt: float,
    segment_lengths: Sequence[int],
) -> DerivativeResult:
    """Compute central finite differences on each contiguous trajectory segment.

    Args:
        A: Reduced coefficients with shape ``(r, M)``.
        U_controls: Aligned controls with shape ``(M, 5)``.
        dt: Time step used between adjacent samples.
        segment_lengths: Lengths of contiguous segments in flattened ordering.

    Returns:
        :class:`DerivativeResult` containing derivative-ready samples.
    """

    A = np.asarray(A)
    U_controls = np.asarray(U_controls)

    if A.ndim != 2:
        raise ValueError("A must have shape (r, M)")
    if U_controls.ndim != 2:
        raise ValueError("U_controls must have shape (M, 5)")
    if A.shape[1] != U_controls.shape[0]:
        raise ValueError("A and controls length mismatch")

    start = 0
    a_used: List[np.ndarray] = []
    adot_used: List[np.ndarray] = []
    u_used: List[np.ndarray] = []
    idx_used: List[np.ndarray] = []

    for L in segment_lengths:
        end = start + int(L)
        if end - start < 3:
            start = end
            continue

        A_seg = A[:, start:end].T
        U_seg = U_controls[start:end]

        i = np.arange(1, A_seg.shape[0] - 1)
        dA = (A_seg[i + 1] - A_seg[i - 1]) / (2.0 * dt)

        a_used.append(A_seg[i])
        adot_used.append(dA)
        u_used.append(U_seg[i])
        idx_used.append(start + i)

        start = end

    if not a_used:
        raise ValueError("No valid derivative points (need segment length >= 3)")

    return DerivativeResult(
        A_used=np.vstack(a_used),
        dA_dt=np.vstack(adot_used),
        U_used=np.vstack(u_used),
        used_indices=np.concatenate(idx_used),
    )
