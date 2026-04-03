"""
Rollout utilities for forecasting in POD and physical space.

This module:
- Selects initial conditions and ground truth segments from POD coefficients
- Uses PySINDy model.simulate for coefficient rollouts
- Reconstructs stacked state vectors and field tensors
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pysindy as ps

from .pod import reconstruct_from_pod
from .snapshot import SnapshotLayout, state_vec_to_fields, fields_to_state_vec


@dataclass(frozen=True)
class RolloutResult:
    """
    Container for a single rollout result.

    Attributes
    ----------
    A_true : np.ndarray
        Ground truth coefficients of shape (T, r). Always present.
    q_true : np.ndarray
        Ground truth stacked states of shape (T, n_state). Always present.
    fields_true : np.ndarray
        Ground truth fields from the original high-fidelity dataset,
        of shape (T, C, ny, nx). Always present.
    A_pred : np.ndarray, optional
        Predicted coefficients of shape (T, r). None if SINDy is disabled.
    q_pred : np.ndarray, optional
        Predicted stacked states of shape (T, n_state). None if SINDy is disabled.
    fields_pred : np.ndarray
        Predicted fields of shape (T, C, ny, nx). None if SINDy is disabled.
    """

    A_true: np.ndarray
    q_true: np.ndarray
    fields_true: np.ndarray
    A_pred: Optional[np.ndarray] = None
    q_pred: Optional[np.ndarray] = None
    fields_pred: Optional[np.ndarray] = None


def _coeff_segment(
    A_traj: np.ndarray,
    start_idx: int,
    horizon_steps: int,
) -> np.ndarray:
    """
    Extract a coefficient segment of length horizon_steps + 1.

    Parameters
    ----------
    A_traj : np.ndarray
        Trajectory coefficients of shape (T, r).
    start_idx : int
        Start index for the segment.
    horizon_steps : int
        Number of forward steps.

    Returns
    -------
    np.ndarray
        Segment of shape (horizon_steps + 1, r).
    """
    T, _ = A_traj.shape
    end_idx = start_idx + horizon_steps
    if start_idx < 0 or end_idx >= T:
        raise ValueError("Requested rollout segment exceeds available trajectory length")
    return A_traj[start_idx : end_idx + 1]


def rollout_one(
    model: ps.SINDy,
    U: np.ndarray,
    layout: SnapshotLayout,
    full_field_traj: np.ndarray,
    A_traj: np.ndarray,
    dt: float,
    horizon_steps: int,
    mean_state: Optional[np.ndarray] = None,
    start_idx: int = 0,
) -> RolloutResult:
    """
    Roll out a learned SINDy model from a single trajectory initial condition.

    Parameters
    ----------
    model : ps.SINDy
        Fitted SINDy model in coefficient space.
    U : np.ndarray
        POD basis of shape (n_state, r).
    full_field_traj : np.ndarray
        The original, high-fidelity trajectory of shape (T, C, ny, nx).
    layout : SnapshotLayout
        Snapshot layout for unstacking fields.
    A_traj : np.ndarray
        Ground truth coefficients for one trajectory, shape (T, r).
    dt : float
        Time step between coefficient samples.
    horizon_steps : int
        Number of forward steps to predict.
    mean_state : np.ndarray, optional
        Mean state added back during reconstruction if POD was centered.
    start_idx : int
        Starting time index for the rollout.

    Returns
    -------
    RolloutResult
        Forecast results in coefficient, stacked-state, and field form.
    """
    A_traj = np.asarray(A_traj)
    if A_traj.ndim != 2:
        raise ValueError("A_traj must have shape (T, r)")
    if dt <= 0:
        raise ValueError("dt must be positive")
    if horizon_steps < 1:
        raise ValueError("horizon_steps must be >= 1")

    # Ground truth coefficients for the rollout segment
    A_true = _coeff_segment(A_traj, start_idx=start_idx, horizon_steps=horizon_steps)
    a0 = A_true[0]

    # Ground truth fields for the rollout segment (from original data)
    end_idx = start_idx + horizon_steps
    fields_true = full_field_traj[start_idx : end_idx + 1]
    q_true_mat = np.stack(
        [fields_to_state_vec(fields_true[i], layout) for i in range(fields_true.shape[0])], axis=0
    )

    # Predict coefficients using the SINDy model
    t = np.arange(horizon_steps + 1, dtype=float) * dt
    A_pred = model.simulate(a0, t=t, integrator="odeint")
    A_pred = np.asarray(A_pred)

    if A_pred.shape != A_true.shape:
        raise ValueError(f"simulate returned shape {A_pred.shape}, expected {A_true.shape}")

    # Reconstruct predicted fields from predicted coefficients
    q_pred_mat = reconstruct_from_pod(U, A_pred.T, mean_state=mean_state).T
    fields_pred = np.stack(
        [state_vec_to_fields(q_pred_mat[i], layout) for i in range(q_pred_mat.shape[0])], axis=0
    )

    return RolloutResult(
        A_true=A_true,
        q_true=q_true_mat,
        fields_true=fields_true,
        A_pred=A_pred,
        q_pred=q_pred_mat,
        fields_pred=fields_pred,
    )


def reshape_coeffs_by_trajectory(
    A: np.ndarray,
    n_trajectories: int,
    T_per_traj: int,
) -> np.ndarray:
    """
    Reshape concatenated coefficients into (n_traj, T, r).

    Parameters
    ----------
    A : np.ndarray
        Coefficients of shape (r, M) or (M, r).
    n_trajectories : int
        Number of trajectories.
    T_per_traj : int
        Number of time points per trajectory used in snapshots.

    Returns
    -------
    np.ndarray
        Coefficients of shape (n_traj, T, r).
    """
    A = np.asarray(A)
    if A.ndim != 2:
        raise ValueError("A must be 2D")

    if A.shape[0] == n_trajectories * T_per_traj:
        A_mt = A
    elif A.shape[1] == n_trajectories * T_per_traj:
        A_mt = A.T
    else:
        raise ValueError("A does not match provided (n_trajectories, T_per_traj) sizing")

    return A_mt.reshape(n_trajectories, T_per_traj, -1)
