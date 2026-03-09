"""
Metrics for ROM forecasting in coefficient and field space.

This module computes:
- Coefficient-space errors
- Field-space L2 and relative errors
- Discrete energy curves on a uniform grid
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np


@dataclass(frozen=True)
class ErrorCurves:
    """
    Time-resolved error curves.

    Attributes
    ----------
    coeff_mse : np.ndarray
        Mean squared error in coefficient space, shape (T,).
    field_l2 : np.ndarray
        L2 norm of field error, shape (T,).
    field_rel_l2 : np.ndarray
        Relative L2 error in field space, shape (T,).
    """

    coeff_mse: np.ndarray
    field_l2: np.ndarray
    field_rel_l2: np.ndarray


@dataclass(frozen=True)
class EnergyCurves:
    """
    Time-resolved energy curves.

    Energy definition for Burgers (2 components):
        E(t) = 0.5 * integral (u^2 + v^2) dA

    Attributes
    ----------
    energy_true : np.ndarray
        Energy of true fields, shape (T,).
    energy_pred : np.ndarray
        Energy of predicted fields, shape (T,).
    """

    energy_true: np.ndarray
    energy_pred: np.ndarray


def coefficient_mse(
    A_true: np.ndarray,
    A_pred: np.ndarray,
) -> np.ndarray:
    """
    Compute coefficient-space MSE per time step.

    Parameters
    ----------
    A_true : np.ndarray
        True coefficients, shape (T, r).
    A_pred : np.ndarray
        Predicted coefficients, shape (T, r).

    Returns
    -------
    np.ndarray
        MSE per time step, shape (T,).
    """
    A_true = np.asarray(A_true)
    A_pred = np.asarray(A_pred)
    if A_true.shape != A_pred.shape:
        raise ValueError("A_true and A_pred must have the same shape")
    return np.mean((A_true - A_pred) ** 2, axis=1)


def field_l2_errors(
    fields_true: np.ndarray,
    fields_pred: np.ndarray,
    dx: float,
    dy: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute L2 and relative L2 field errors per time step.

    Parameters
    ----------
    fields_true : np.ndarray
        True fields, shape (T, C, ny, nx).
    fields_pred : np.ndarray
        Predicted fields, shape (T, C, ny, nx).
    dx : float
        Grid spacing in x.
    dy : float
        Grid spacing in y.

    Returns
    -------
    (np.ndarray, np.ndarray)
        field_l2, field_rel_l2 each of shape (T,).
    """
    ft = np.asarray(fields_true)
    fp = np.asarray(fields_pred)
    if ft.shape != fp.shape:
        raise ValueError("fields_true and fields_pred must have the same shape")
    if ft.ndim != 4:
        raise ValueError("Expected fields shape (T, C, ny, nx)")
    if dx <= 0 or dy <= 0:
        raise ValueError("dx and dy must be positive")

    err = fp - ft

    wt = dx * dy
    err_sq = np.sum(err**2, axis=(1, 2, 3)) * wt
    true_sq = np.sum(ft**2, axis=(1, 2, 3)) * wt

    l2 = np.sqrt(err_sq)
    denom = np.sqrt(true_sq)
    rel = np.where(denom > 0, l2 / denom, 0.0)

    return l2, rel


def burgers_energy(
    fields: np.ndarray,
    dx: float,
    dy: float,
) -> np.ndarray:
    """
    Compute Burgers kinetic-energy-like quantity per time step.

    E(t) = 0.5 * integral (u^2 + v^2) dA

    Parameters
    ----------
    fields : np.ndarray
        Fields, shape (T, 2, ny, nx).
    dx : float
        Grid spacing in x.
    dy : float
        Grid spacing in y.

    Returns
    -------
    np.ndarray
        Energy curve, shape (T,).
    """
    f = np.asarray(fields)
    if f.ndim != 4 or f.shape[1] != 2:
        raise ValueError("Expected fields shape (T, 2, ny, nx)")
    if dx <= 0 or dy <= 0:
        raise ValueError("dx and dy must be positive")

    wt = dx * dy
    sq = np.sum(f**2, axis=1)
    e = 0.5 * np.sum(sq, axis=(1, 2)) * wt
    return e


def compute_curves(
    A_true: np.ndarray,
    A_pred: np.ndarray,
    fields_true: np.ndarray,
    fields_pred: np.ndarray,
    dx: float,
    dy: float,
) -> Tuple[ErrorCurves, EnergyCurves]:
    """
    Compute error and energy curves for one rollout.

    Parameters
    ----------
    A_true : np.ndarray
        True coefficients, shape (T, r).
    A_pred : np.ndarray
        Predicted coefficients, shape (T, r).
    fields_true : np.ndarray
        True fields, shape (T, 2, ny, nx).
    fields_pred : np.ndarray
        Predicted fields, shape (T, 2, ny, nx).
    dx : float
        Grid spacing in x.
    dy : float
        Grid spacing in y.

    Returns
    -------
    (ErrorCurves, EnergyCurves)
        Time-resolved metrics.
    """
    c_mse = coefficient_mse(A_true, A_pred)
    f_l2, f_rel = field_l2_errors(fields_true, fields_pred, dx=dx, dy=dy)

    e_true = burgers_energy(fields_true, dx=dx, dy=dy)
    e_pred = burgers_energy(fields_pred, dx=dx, dy=dy)

    return (
        ErrorCurves(coeff_mse=c_mse, field_l2=f_l2, field_rel_l2=f_rel),
        EnergyCurves(energy_true=e_true, energy_pred=e_pred),
    )


def summarize_aggregates(
    err: ErrorCurves,
    energy: EnergyCurves,
) -> Dict[str, float]:
    """
    Summarize rollout metrics into scalar aggregates.

    Parameters
    ----------
    err : ErrorCurves
        Error curves.
    energy : EnergyCurves
        Energy curves.

    Returns
    -------
    dict
        Aggregate metrics:
        - final_field_rel_l2: Relative L2 error of fields at the last step.
        - mean_field_rel_l2: Average relative L2 error over the trajectory.
        - mean_coeff_mse: Average MSE of POD coefficients over the trajectory.
        - final_energy_drift: Energy difference (pred - true) at the last step.
        - mean_abs_energy_drift: Average absolute energy difference over the trajectory.
    """
    final_rel = float(err.field_rel_l2[-1]) if err.field_rel_l2.size else 0.0
    mean_rel = float(np.mean(err.field_rel_l2)) if err.field_rel_l2.size else 0.0
    mean_coeff = float(np.mean(err.coeff_mse)) if err.coeff_mse.size else 0.0

    drift = energy.energy_pred - energy.energy_true
    final_drift = float(drift[-1]) if drift.size else 0.0
    mean_abs_drift = float(np.mean(np.abs(drift))) if drift.size else 0.0

    return {
        "final_field_rel_l2": final_rel,
        "mean_field_rel_l2": mean_rel,
        "mean_coeff_mse": mean_coeff,
        "final_energy_drift": final_drift,
        "mean_abs_energy_drift": mean_abs_drift,
    }
