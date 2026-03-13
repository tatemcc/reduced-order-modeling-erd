"""Snapshot stacking utilities for POD and SINDy training."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np

from .dataset_io import TrajectoryData


@dataclass(frozen=True)
class SnapshotLayout:
    """Layout metadata for stacked ERD states.

    Attributes:
        nr: Number of radial cells.
        nphi: Number of azimuthal cells.
        n_components: Number of stacked fields (``2`` for ``n`` and ``omega``).
    """

    nr: int
    nphi: int
    n_components: int = 2

    @property
    def n_points(self) -> int:
        """Return the number of cells per field, ``N_r * N_phi``."""

        return self.nr * self.nphi

    @property
    def n_state(self) -> int:
        """Return full stacked state dimension, ``n_components * n_points``."""

        return self.n_points * self.n_components


def fields_to_state_vec(fields: np.ndarray, layout: SnapshotLayout) -> np.ndarray:
    """Flatten one ``(n, omega)`` snapshot into a state vector.

    Args:
        fields: Snapshot with shape ``(2, N_r, N_phi)``.
        layout: Snapshot layout metadata.

    Returns:
        Flattened state vector ``x`` with length ``layout.n_state``.
    """

    if fields.shape != (layout.n_components, layout.nr, layout.nphi):
        raise ValueError(f"Unexpected fields shape: {fields.shape}")
    return fields.reshape(layout.n_state)


def state_vec_to_fields(state: np.ndarray, layout: SnapshotLayout) -> np.ndarray:
    """Reshape one stacked state vector back to field form.

    Args:
        state: State vector with length ``layout.n_state``.
        layout: Snapshot layout metadata.

    Returns:
        Field tensor with shape ``(2, N_r, N_phi)``.
    """

    return np.asarray(state).reshape(layout.n_components, layout.nr, layout.nphi)


def build_snapshot_matrix(
    trajectories: Sequence[TrajectoryData],
    center: bool,
    time_stride: int = 1,
    time_limit: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray, SnapshotLayout, Optional[np.ndarray], List[int], List[np.ndarray]]:
    """Assemble stacked snapshots and aligned controls from trajectories.

    Args:
        trajectories: Sequence of loaded trajectories from a manifest split.
        center: Whether to subtract global mean state from stacked snapshots.
        time_stride: Keep every ``time_stride`` snapshot in each trajectory.
        time_limit: Optional limit on snapshots retained per trajectory.

    Returns:
        Tuple ``(X, U, layout, mean_state, segment_lengths, times)`` where:
        - ``X`` has shape ``(n_state, M)``,
        - ``U`` has shape ``(M, 5)``,
        - ``segment_lengths`` identifies per-trajectory contiguous blocks,
        - ``times`` stores retained timestamp vectors per trajectory.
    """

    if not trajectories:
        raise ValueError("No trajectories provided")

    T = trajectories[0].fields.shape[0]
    c, nr, nphi = trajectories[0].fields.shape[1:]
    layout = SnapshotLayout(nr=nr, nphi=nphi, n_components=c)

    t_idx = np.arange(T)
    if time_stride > 1:
        t_idx = t_idx[::time_stride]
    if time_limit is not None:
        t_idx = t_idx[:time_limit]

    M = len(trajectories) * len(t_idx)
    X = np.empty((layout.n_state, M), dtype=float)
    U = np.empty((M, 5), dtype=float)

    seg_lengths: List[int] = []
    times: List[np.ndarray] = []

    col = 0
    for traj in trajectories:
        seg_lengths.append(len(t_idx))
        times.append(traj.t[t_idx])
        for k in t_idx:
            X[:, col] = fields_to_state_vec(traj.fields[k], layout)
            U[col] = traj.controls[k]
            col += 1

    mean_state = None
    if center:
        mean_state = X.mean(axis=1)
        X = X - mean_state[:, None]

    return X, U, layout, mean_state, seg_lengths, times
