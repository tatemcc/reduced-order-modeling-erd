"""
SINDy model fitting utilities for POD coefficient dynamics.

This module fits a sparse dynamical system in reduced coordinates:
    dA_dt = f(A)
using PySINDy.

Primary design choices:
- Use PySINDy PolynomialLibrary for linear/quadratic models
- Use STLSQ or SR3 optimizers
- Store enough metadata to reproduce runs and support rollout
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
import pysindy as ps


@dataclass(frozen=True)
class SINDyFitResult:
    """
    Container for a fitted SINDy model in reduced coordinates.

    Attributes
    ----------
    model : ps.SINDy
        Fitted PySINDy model.
    feature_names : list of str
        Names of library features in the same order as model coefficients.
    coefficient_matrix : np.ndarray
        Xi matrix of shape (n_targets, n_features).
    n_targets : int
        Number of state variables modeled (equals POD rank r).
    n_features : int
        Number of library features.
    """
    model: ps.SINDy
    feature_names: List[str]
    coefficient_matrix: np.ndarray
    n_targets: int
    n_features: int


def build_library(
    poly_order: int,
    include_bias: bool,
) -> ps.PolynomialLibrary:
    """
    Construct a polynomial feature library.

    Parameters
    ----------
    poly_order : int
        Maximum polynomial degree. Use 1 for linear, 2 for quadratic.
    include_bias : bool
        If True, include a constant feature.

    Returns
    -------
    ps.PolynomialLibrary
        Configured polynomial library.
    """
    if poly_order < 1:
        raise ValueError("poly_order must be >= 1")
    return ps.PolynomialLibrary(degree=poly_order, include_bias=include_bias)


def build_optimizer(
    optimizer_name: str,
    sparsity: float,
) -> ps.BaseOptimizer:
    """
    Construct a PySINDy optimizer.

    Parameters
    ----------
    optimizer_name : str
        'stlsq' or 'sr3'.
    sparsity : float
        Sparsity knob. Interpretation depends on optimizer:
        - STLSQ: threshold
        - SR3: threshold

    Returns
    -------
    ps.BaseOptimizer
        Configured optimizer.
    """
    if sparsity < 0:
        raise ValueError("sparsity must be nonnegative")

    if optimizer_name == "stlsq":
        return ps.STLSQ(threshold=sparsity)
    if optimizer_name == "sr3":
        return ps.SR3() # TODO fix

    raise ValueError(f"Unsupported optimizer: {optimizer_name}")


def fit_sindy_on_coeffs(
    A_used: np.ndarray,
    dA_dt: np.ndarray,
    dt: float,
    poly_order: int,
    include_bias: bool,
    optimizer_name: str,
    sparsity: float,
    feature_names: Optional[List[str]] = None,
) -> SINDyFitResult:
    """
    Fit a SINDy model for POD coefficient dynamics.

    Parameters
    ----------
    A_used : np.ndarray
        Coefficients aligned with derivatives. Shape (r, M_used).
    dA_dt : np.ndarray
        Coefficient derivatives. Shape (r, M_used).
    dt : float
        Time step (used by PySINDy and later by simulate).
    poly_order : int
        Polynomial library order (1 or 2 recommended).
    include_bias : bool
        Include constant feature in library.
    optimizer_name : str
        Optimizer name: 'stlsq' or 'sr3'.
    sparsity : float
        Sparsity knob (threshold).
    feature_names : list of str, optional
        Custom feature names for state variables. If None, defaults are used.

    Returns
    -------
    SINDyFitResult
        Fitted model and metadata.
    """
    A_used = np.asarray(A_used)
    dA_dt = np.asarray(dA_dt)

    if A_used.shape != dA_dt.shape:
        raise ValueError("A_used and dA_dt must have the same shape")
    if A_used.ndim != 2:
        raise ValueError("A_used must have shape (r, M_used)")
    if dt <= 0:
        raise ValueError("dt must be positive")

    r, _ = A_used.shape

    X = A_used.T
    Xdot = dA_dt.T

    lib = build_library(poly_order=poly_order, include_bias=include_bias)
    opt = build_optimizer(optimizer_name=optimizer_name, sparsity=sparsity)

    if feature_names is None:
        feature_names = [f"a{i}" for i in range(r)]

    model = ps.SINDy(
        feature_library=lib,
        optimizer=opt,
        differentiation_method=None # we did the derivatives ourselves
    )
    model.fit(X, t=dt, x_dot=Xdot)

    Xi = model.coefficients()
    names = model.get_feature_names()

    return SINDyFitResult(
        model=model,
        feature_names=names,
        coefficient_matrix=Xi,
        n_targets=Xi.shape[0],
        n_features=Xi.shape[1],
    )


def simulate_coeffs(
    model: ps.SINDy,
    a0: np.ndarray,
    dt: float,
    n_steps: int,
) -> np.ndarray:
    """
    Roll out POD coefficients using PySINDy model.simulate.

    Parameters
    ----------
    model : ps.SINDy
        Fitted SINDy model.
    a0 : np.ndarray
        Initial coefficient vector of shape (r,).
    dt : float
        Time step.
    n_steps : int
        Number of steps to simulate forward.

    Returns
    -------
    np.ndarray
        Simulated coefficients of shape (n_steps + 1, r).
    """
    a0 = np.asarray(a0).reshape(-1)
    if dt <= 0:
        raise ValueError("dt must be positive")
    if n_steps < 1:
        raise ValueError("n_steps must be >= 1")

    t = np.arange(n_steps + 1, dtype=float) * dt
    A_sim = model.simulate(a0, t=t)
    return np.asarray(A_sim)
