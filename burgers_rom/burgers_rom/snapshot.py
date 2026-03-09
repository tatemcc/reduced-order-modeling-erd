"""
Snapshot construction utilities for POD and reduced modeling.

This module defines:
- A consistent stacking convention for Burgers GRID fields
- Helpers to convert between field tensors and stacked vectors
- Construction of the snapshot matrix X used for POD
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class SnapshotLayout:
    """
    Shape metadata for mapping between field tensors and stacked vectors.

    Attributes
    ----------
    ny : int
        Grid size in y.
    nx : int
        Grid size in x.
    n_components : int
        Number of PDE components. Burgers uses 2.
    """

    ny: int
    nx: int
    n_components: int = 2

    @property
    def n_points(self) -> int:
        """Number of spatial grid points."""
        return self.ny * self.nx

    @property
    def n_state(self) -> int:
        """Length of stacked state vector."""
        return self.n_components * self.n_points


def fields_to_state_vec(
    fields: np.ndarray,
    layout: SnapshotLayout,
) -> np.ndarray:
    """
    Stack Burgers fields into a single state vector.

    Parameters
    ----------
    fields : np.ndarray
        Field array of shape (C, ny, nx).
    layout : SnapshotLayout
        Layout metadata.

    Returns
    -------
    np.ndarray
        Stacked state vector of shape (C*ny*nx,).
    """
    if fields.shape != (layout.n_components, layout.ny, layout.nx):
        raise ValueError(
            f"Expected fields shape {(layout.n_components, layout.ny, layout.nx)}, got {fields.shape}"
        )
    return fields.reshape(layout.n_state, order="C")


def state_vec_to_fields(
    state: np.ndarray,
    layout: SnapshotLayout,
) -> np.ndarray:
    """
    Unstack a state vector into Burgers fields.

    Parameters
    ----------
    state : np.ndarray
        State vector of shape (C*ny*nx,).
    layout : SnapshotLayout
        Layout metadata.

    Returns
    -------
    np.ndarray
        Field array of shape (C, ny, nx).
    """
    state = np.asarray(state)
    if state.ndim != 1 or state.shape[0] != layout.n_state:
        raise ValueError(f"Expected state shape ({layout.n_state},), got {state.shape}")
    return state.reshape((layout.n_components, layout.ny, layout.nx), order="C")


def build_snapshot_matrix(
    trajectories: np.ndarray,
    center: bool,
    time_stride: int = 1,
    time_limit: Optional[int] = None,
) -> Tuple[np.ndarray, SnapshotLayout, Optional[np.ndarray]]:
    """
    Construct the POD snapshot matrix X from one or more trajectories.

    The snapshot matrix is:
        X in R^{n_state x M}
    where columns are stacked state vectors at each sampled time.

    Parameters
    ----------
    trajectories : np.ndarray
        Array of shape (n_traj, T, C, ny, nx).
    center : bool
        If True, subtract the columnwise mean state (mean over time and trajectories).
    time_stride : int
        Subsampling factor along time before stacking.
    time_limit : int, optional
        Maximum number of time steps used per trajectory.

    Returns
    -------
    X : np.ndarray
        Snapshot matrix of shape (n_state, M).
    layout : SnapshotLayout
        Layout metadata.
    mean_state : np.ndarray or None
        Mean state vector of shape (n_state,) if center is True, else None.
    """
    traj = np.asarray(trajectories)
    if traj.ndim != 5:
        raise ValueError("Expected trajectories shape (n_traj, T, C, ny, nx)")

    n_traj, T, C, ny, nx = traj.shape
    if C != 2:
        raise ValueError(f"Expected Burgers components C=2, got C={C}")

    layout = SnapshotLayout(ny=ny, nx=nx, n_components=C)

    t_idx = np.arange(T)
    if time_stride > 1:
        t_idx = t_idx[::time_stride]
    if time_limit is not None:
        t_idx = t_idx[:time_limit]

    M = n_traj * len(t_idx)
    X = np.empty((layout.n_state, M), dtype=traj.dtype)

    col = 0
    for i in range(n_traj):
        for k in t_idx:
            X[:, col] = fields_to_state_vec(traj[i, k], layout)
            col += 1

    mean_state: Optional[np.ndarray] = None
    if center:
        mean_state = X.mean(axis=1)
        X = X - mean_state[:, None]

    return X, layout, mean_state
