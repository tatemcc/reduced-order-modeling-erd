"""Held-out ROM rollout helpers for model validation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .pod import reconstruct
from .sindy_model import rollout_controlled
from .snapshot import SnapshotLayout, state_vec_to_fields


@dataclass(frozen=True)
class RolloutResult:
    """Artifacts from one true-vs-predicted controlled rollout.

    Attributes:
        A_true: Ground-truth reduced coefficients.
        A_pred: Predicted reduced coefficients from SINDy rollout.
        q_true: Ground-truth full-order stacked states.
        q_pred: Predicted full-order stacked states.
        fields_true: Ground-truth fields shaped ``(T, 2, N_r, N_phi)``.
        fields_pred: Predicted fields shaped ``(T, 2, N_r, N_phi)``.
        U_used: Control sequence used during rollout.
    """

    A_true: np.ndarray
    A_pred: np.ndarray
    q_true: np.ndarray
    q_pred: np.ndarray
    fields_true: np.ndarray
    fields_pred: np.ndarray
    U_used: np.ndarray


def rollout_one(
    model: object,
    U_basis: np.ndarray,
    layout: SnapshotLayout,
    A_traj: np.ndarray,
    U_traj: np.ndarray,
    dt: float,
    horizon_steps: int,
    state_clip: float = 50.0,
    mean_state: Optional[np.ndarray] = None,
    start_idx: int = 0,
) -> RolloutResult:
    """Run one controlled rollout and reconstruct field predictions.

    Args:
        model: Fitted controlled model with ``predict(x, u)`` method.
        U_basis: POD basis matrix.
        layout: Layout metadata for stacked fields.
        A_traj: Ground-truth coefficient trajectory.
        U_traj: Control trajectory aligned with ``A_traj``.
        dt: Rollout integration step.
        horizon_steps: Maximum number of steps to forecast.
        state_clip: Absolute bound applied to predicted reduced coefficients.
        mean_state: Optional POD centering vector.
        start_idx: Trajectory index used as rollout initial condition.

    Returns:
        :class:`RolloutResult` containing reduced and lifted trajectories.
    """

    A_traj = np.asarray(A_traj)
    U_traj = np.asarray(U_traj)

    T = A_traj.shape[0]
    end_idx = min(T - 1, start_idx + horizon_steps)

    A_true = A_traj[start_idx : end_idx + 1]
    U_used = U_traj[start_idx : end_idx + 1]

    A_pred = rollout_controlled(model, a0=A_true[0], u_seq=U_used, dt=dt, state_clip=state_clip)

    q_true = reconstruct(U_basis, A_true.T, mean_state=mean_state).T
    q_pred = reconstruct(U_basis, A_pred.T, mean_state=mean_state).T

    fields_true = np.stack([state_vec_to_fields(x, layout) for x in q_true], axis=0)
    fields_pred = np.stack([state_vec_to_fields(x, layout) for x in q_pred], axis=0)

    return RolloutResult(
        A_true=A_true,
        A_pred=A_pred,
        q_true=q_true,
        q_pred=q_pred,
        fields_true=fields_true,
        fields_pred=fields_pred,
        U_used=U_used,
    )
