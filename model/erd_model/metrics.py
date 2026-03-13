"""Validation metrics for ROM rollouts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np


@dataclass(frozen=True)
class ErrorCurves:
    """Time-series error curves from one rollout comparison.

    Attributes:
        coeff_mse: Mean-squared error in reduced coefficients per step.
        field_l2: Absolute L2 error in stacked fields per step.
        field_rel_l2: Relative L2 field error per step.
    """

    coeff_mse: np.ndarray
    field_l2: np.ndarray
    field_rel_l2: np.ndarray


@dataclass(frozen=True)
class EnergyCurves:
    """Pseudo-energy curves for true vs predicted fields.

    Attributes:
        energy_true: Pseudo-energy of true fields per step.
        energy_pred: Pseudo-energy of predicted fields per step.
    """

    energy_true: np.ndarray
    energy_pred: np.ndarray


def coefficient_mse(A_true: np.ndarray, A_pred: np.ndarray) -> np.ndarray:
    """Compute per-step mean-squared coefficient error.

    Args:
        A_true: Ground-truth reduced coefficients.
        A_pred: Predicted reduced coefficients.

    Returns:
        1D array of coefficient MSE values.
    """

    return np.mean((A_true - A_pred) ** 2, axis=1)


def field_l2_errors(fields_true: np.ndarray, fields_pred: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute absolute and relative L2 field errors.

    Args:
        fields_true: True stacked fields with shape ``(T, 2, N_r, N_phi)``.
        fields_pred: Predicted stacked fields with the same shape.

    Returns:
        Tuple ``(abs_l2, rel_l2)`` of per-step errors.
    """

    err = fields_pred - fields_true
    l2 = np.sqrt(np.sum(err**2, axis=(1, 2, 3)))
    denom = np.sqrt(np.sum(fields_true**2, axis=(1, 2, 3)))
    rel = np.where(denom > 0.0, l2 / denom, 0.0)
    return l2, rel


def pseudo_energy(fields: np.ndarray) -> np.ndarray:
    """Compute a simple quadratic field energy proxy.

    Args:
        fields: Stacked field tensor.

    Returns:
        1D pseudo-energy curve.
    """

    return 0.5 * np.sum(fields**2, axis=(1, 2, 3))


def compute_curves(
    A_true: np.ndarray,
    A_pred: np.ndarray,
    fields_true: np.ndarray,
    fields_pred: np.ndarray,
) -> Tuple[ErrorCurves, EnergyCurves]:
    """Compute error and energy curves for one rollout.

    Args:
        A_true: Ground-truth reduced coefficients.
        A_pred: Predicted reduced coefficients.
        fields_true: Ground-truth stacked fields.
        fields_pred: Predicted stacked fields.

    Returns:
        Tuple ``(error_curves, energy_curves)``.
    """

    coeff = coefficient_mse(A_true, A_pred)
    fl2, frel = field_l2_errors(fields_true, fields_pred)

    e_true = pseudo_energy(fields_true)
    e_pred = pseudo_energy(fields_pred)

    return (
        ErrorCurves(coeff_mse=coeff, field_l2=fl2, field_rel_l2=frel),
        EnergyCurves(energy_true=e_true, energy_pred=e_pred),
    )


def summarize_aggregates(err: ErrorCurves, energy: EnergyCurves) -> Dict[str, float]:
    """Summarize rollout curves into scalar report metrics.

    Args:
        err: Error curves for one rollout.
        energy: Energy curves for one rollout.

    Returns:
        Dictionary of aggregate statistics used in summaries.
    """

    drift = energy.energy_pred - energy.energy_true
    return {
        "final_field_rel_l2": float(err.field_rel_l2[-1]),
        "mean_field_rel_l2": float(np.mean(err.field_rel_l2)),
        "mean_coeff_mse": float(np.mean(err.coeff_mse)),
        "final_energy_drift": float(drift[-1]),
        "mean_abs_energy_drift": float(np.mean(np.abs(drift))),
    }
