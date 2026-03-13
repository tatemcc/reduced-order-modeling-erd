"""Controlled SINDy fitting and rollout helpers for ERD reduced dynamics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence

import numpy as np
import pysindy as ps

from .config import TrappingSR3Config


@dataclass(frozen=True)
class SINDyControlFitResult:
    """Outputs from fitting ``a_dot = f(a, u)`` with PySINDy.

    Attributes:
        model: Fitted controlled model object with ``predict(x, u)`` API.
        feature_names: Human-readable feature names in library order.
        coefficient_matrix: Sparse coefficient matrix ``Xi``.
        fit_info: Dictionary with optimizer diagnostics.
    """

    model: object
    feature_names: List[str]
    coefficient_matrix: np.ndarray
    fit_info: dict[str, object]


class TrappingControlledHybridModel:
    """Hybrid controlled ROM with stabilized autonomous and sparse control parts.

    The model evaluates:
    ``a_dot = f_trap(a) + g_ctrl(a, u)``,
    where ``f_trap`` is a PySINDy model fitted with TrappingSR3 and
    ``g_ctrl`` is a sparse linear model in ``[u, a_i u_j]``.
    """

    def __init__(
        self,
        autonomous_model: ps.SINDy,
        control_coeffs: np.ndarray,
        control_feature_names: Sequence[str],
    ):
        """Store autonomous and control-residual submodels.

        Args:
            autonomous_model: TrappingSR3-fitted model for autonomous dynamics.
            control_coeffs: Coefficients for ``[u, a*u]`` features.
            control_feature_names: Feature labels matching ``control_coeffs`` columns.

        Returns:
            None.
        """

        self.autonomous_model = autonomous_model
        self.control_coeffs = np.asarray(control_coeffs, dtype=float)
        self.control_feature_names = list(control_feature_names)

    def _control_features(self, x: np.ndarray, u: np.ndarray) -> np.ndarray:
        """Build residual-control feature matrix ``[u, a_i u_j]``.

        Args:
            x: Reduced states with shape ``(M, r)``.
            u: Controls with shape ``(M, n_u)``.

        Returns:
            Feature matrix with shape ``(M, n_u + r*n_u)``.
        """

        bilinear = (x[:, :, None] * u[:, None, :]).reshape(x.shape[0], -1)
        return np.concatenate([u, bilinear], axis=1)

    def predict(self, x: np.ndarray, u: np.ndarray | None = None) -> np.ndarray:
        """Predict reduced derivatives for one or many samples.

        Args:
            x: Reduced states with shape ``(r,)`` or ``(M, r)``.
            u: Controls with shape ``(n_u,)`` or ``(M, n_u)``.

        Returns:
            Predicted derivatives with shape ``(M, r)``.
        """

        x_arr = np.asarray(x, dtype=float)
        if x_arr.ndim == 1:
            x_arr = x_arr.reshape(1, -1)

        if u is None:
            raise ValueError("Control input u is required for hybrid controlled predictions")
        u_arr = np.asarray(u, dtype=float)
        if u_arr.ndim == 1:
            u_arr = u_arr.reshape(1, -1)
        if u_arr.shape[0] == 1 and x_arr.shape[0] > 1:
            u_arr = np.repeat(u_arr, x_arr.shape[0], axis=0)
        if u_arr.shape[0] != x_arr.shape[0]:
            raise ValueError(f"x/u batch mismatch: {x_arr.shape[0]} vs {u_arr.shape[0]}")

        autonomous = np.asarray(self.autonomous_model.predict(x_arr), dtype=float)
        theta_u = self._control_features(x_arr, u_arr)
        controlled = theta_u @ self.control_coeffs.T
        return autonomous + controlled

    def get_feature_names(self) -> List[str]:
        """Return concatenated autonomous and control feature labels.

        Returns:
            Combined feature-name list.
        """

        return list(self.autonomous_model.get_feature_names()) + list(self.control_feature_names)

    def coefficients(self) -> np.ndarray:
        """Return the horizontally stacked coefficient matrix.

        Returns:
            Coefficient matrix with autonomous columns followed by control columns.
        """

        auto = np.asarray(self.autonomous_model.coefficients(), dtype=float)
        return np.hstack([auto, self.control_coeffs])


def build_control_library(n_state: int, n_control: int) -> ps.BaseFeatureLibrary:
    """Build the controlled SINDy library used in this project.

    The resulting library includes:
    - polynomial terms in reduced state ``a`` up to degree 2,
    - linear terms in controls ``u``,
    - bilinear terms ``a_i * u_j``.

    Args:
        n_state: Reduced-state dimension.
        n_control: Control dimension.

    Returns:
        PySINDy feature-library object.
    """

    x_idx = list(range(n_state))
    u_idx = list(range(n_state, n_state + n_control))

    lib_quad_x = ps.PolynomialLibrary(degree=2, include_bias=False)
    lib_lin_x = ps.PolynomialLibrary(degree=1, include_bias=False)
    lib_lin_u = ps.PolynomialLibrary(degree=1, include_bias=False)

    try:
        return ps.GeneralizedLibrary(
            libraries=[lib_quad_x, lib_lin_x, lib_lin_u],
            inputs_per_library=[x_idx, x_idx, u_idx],
            tensor_array=np.array([[0, 1, 1]], dtype=int),
        )
    except TypeError:
        # Compatibility fallback for alternate PySINDy signatures.
        return ps.GeneralizedLibrary(
            [lib_quad_x, lib_lin_x, lib_lin_u],
            inputs_per_library=[x_idx, x_idx, u_idx],
            tensor_array=np.array([[0, 1, 1]], dtype=int),
        )


def _build_control_residual_features(A_used: np.ndarray, U_used: np.ndarray) -> tuple[np.ndarray, List[str]]:
    """Build residual-control features ``[u, a_i u_j]`` and labels.

    Args:
        A_used: Reduced states with shape ``(M, r)``.
        U_used: Controls with shape ``(M, n_u)``.

    Returns:
        Tuple ``(Theta_u, names)``.
    """

    n_state = A_used.shape[1]
    n_control = U_used.shape[1]

    bilinear = (A_used[:, :, None] * U_used[:, None, :]).reshape(A_used.shape[0], -1)
    theta_u = np.concatenate([U_used, bilinear], axis=1)

    names: List[str] = [f"u{i}" for i in range(n_control)]
    for i in range(n_state):
        for j in range(n_control):
            names.append(f"x{i}*u{j}")
    return theta_u, names


def _fit_sparse_ridge(
    Theta: np.ndarray,
    Y: np.ndarray,
    alpha: float,
    threshold: float,
    max_iter: int,
) -> np.ndarray:
    """Fit sparse ridge regression with iterative hard thresholding.

    Args:
        Theta: Regression matrix.
        Y: Target matrix.
        alpha: Ridge regularization weight.
        threshold: Hard threshold applied to coefficients.
        max_iter: Number of sparse refit iterations.

    Returns:
        Coefficient matrix with shape ``(n_targets, n_features)``.
    """

    n_features = Theta.shape[1]
    reg = float(max(alpha, 0.0))
    gram = Theta.T @ Theta + reg * np.eye(n_features)
    rhs = Theta.T @ Y
    Xi = np.linalg.solve(gram, rhs).T

    if threshold <= 0.0:
        return Xi

    n_targets = Xi.shape[0]
    n_iter = int(max(1, max_iter))
    for _ in range(n_iter):
        active = np.abs(Xi) >= threshold
        Xi_new = np.zeros_like(Xi)
        for j in range(n_targets):
            idx = np.flatnonzero(active[j])
            if idx.size == 0:
                continue
            G = Theta[:, idx]
            Gj = G.T @ G + reg * np.eye(idx.size)
            bj = G.T @ Y[:, j]
            Xi_new[j, idx] = np.linalg.solve(Gj, bj)
        Xi = Xi_new
    return Xi


def _build_trapping_optimizer(
    n_state: int,
    threshold: float,
    max_iter: int,
    trapping: TrappingSR3Config,
) -> ps.TrappingSR3:
    """Construct a TrappingSR3 optimizer from project config values.

    Args:
        n_state: Reduced-state dimension.
        threshold: STLSQ threshold fallback for regularization scale.
        max_iter: STLSQ max-iter fallback when trapping value is invalid.
        trapping: Trapping optimizer configuration.

    Returns:
        Configured :class:`pysindy.TrappingSR3` instance.
    """

    reg_weight = trapping.reg_weight_lam if trapping.reg_weight_lam > 0.0 else max(float(threshold), 1.0e-6)
    trap_iter = trapping.max_iter if trapping.max_iter > 0 else int(max_iter)
    return ps.TrappingSR3(
        _n_tgts=n_state,
        _include_bias=False,
        method=trapping.method,
        reg_weight_lam=reg_weight,
        relax_coeff_nu=trapping.relax_coeff_nu,
        tol=trapping.tol,
        max_iter=trap_iter,
        eta=trapping.eta,
        alpha=trapping.stability_alpha,
        beta=trapping.stability_beta,
        normalize_columns=False,
        verbose=False,
        unbias=False,
    )


def fit_sindy_control(
    A_used: np.ndarray,
    dA_dt: np.ndarray,
    U_used: np.ndarray,
    dt: float,
    threshold: float,
    alpha: float,
    max_iter: int,
    optimizer: str = "stlsq",
    trapping: TrappingSR3Config | None = None,
) -> SINDyControlFitResult:
    """Fit controlled SINDy from precomputed reduced samples.

    Args:
        A_used: Regression states with shape ``(M, r)``.
        dA_dt: Target derivatives with shape ``(M, r)``.
        U_used: Controls aligned with states, shape ``(M, 5)``.
        dt: Sample spacing used for fit metadata.
        threshold: STLSQ sparsity threshold (or fallback scale for trapping path).
        alpha: STLSQ L2 regularization coefficient.
        max_iter: STLSQ iteration cap.
        optimizer: Optimizer mode (``stlsq`` or ``trappingsr3``).
        trapping: Optional trapping config for stabilized hybrid fitting.

    Returns:
        :class:`SINDyControlFitResult` with fitted model and coefficients.
    """

    A_used = np.asarray(A_used)
    dA_dt = np.asarray(dA_dt)
    U_used = np.asarray(U_used)

    if A_used.shape != dA_dt.shape:
        raise ValueError("A_used and dA_dt shape mismatch")

    n_state = A_used.shape[1]
    n_control = U_used.shape[1]
    mode = str(optimizer).strip().lower()

    if mode == "stlsq":
        library = build_control_library(n_state=n_state, n_control=n_control)
        sparse_opt = ps.STLSQ(threshold=threshold, alpha=alpha, max_iter=max_iter)
        model = ps.SINDy(feature_library=library, optimizer=sparse_opt, differentiation_method=None)
        model.fit(x=A_used, u=U_used, x_dot=dA_dt, t=dt)

        Xi = np.asarray(model.coefficients())
        names = model.get_feature_names()
        return SINDyControlFitResult(
            model=model,
            feature_names=names,
            coefficient_matrix=Xi,
            fit_info={
                "optimizer_used": "stlsq",
                "threshold": float(threshold),
                "alpha": float(alpha),
                "max_iter": int(max_iter),
            },
        )

    if mode in {"trappingsr3", "trapping_sr3", "trapping-sr3"}:
        trap_cfg = trapping or TrappingSR3Config()

        # TrappingSR3 in PySINDy 2.0 assumes a polynomial state library and is
        # not directly compatible with external control features. Fit a stable
        # autonomous model first, then regress control residual terms.
        autonomous_lib = ps.PolynomialLibrary(degree=2, include_bias=False)
        trap_opt = _build_trapping_optimizer(n_state=n_state, threshold=threshold, max_iter=max_iter, trapping=trap_cfg)
        autonomous_model = ps.SINDy(
            feature_library=autonomous_lib,
            optimizer=trap_opt,
            differentiation_method=None,
        )
        autonomous_model.fit(x=A_used, x_dot=dA_dt, t=dt)
        autonomous_pred = np.asarray(autonomous_model.predict(A_used), dtype=float)

        residual = dA_dt - autonomous_pred
        theta_u, control_names = _build_control_residual_features(A_used, U_used)
        Xi_u = _fit_sparse_ridge(
            Theta=theta_u,
            Y=residual,
            alpha=trap_cfg.control_alpha,
            threshold=trap_cfg.control_threshold,
            max_iter=trap_cfg.control_max_iter,
        )

        model = TrappingControlledHybridModel(
            autonomous_model=autonomous_model,
            control_coeffs=Xi_u,
            control_feature_names=control_names,
        )
        Xi = np.hstack([np.asarray(autonomous_model.coefficients(), dtype=float), Xi_u])
        names = model.get_feature_names()
        return SINDyControlFitResult(
            model=model,
            feature_names=names,
            coefficient_matrix=Xi,
            fit_info={
                "optimizer_used": "trappingsr3_hybrid",
                "threshold": float(threshold),
                "alpha": float(alpha),
                "max_iter": int(max_iter),
                "control_threshold": float(trap_cfg.control_threshold),
                "control_alpha": float(trap_cfg.control_alpha),
                "control_max_iter": int(trap_cfg.control_max_iter),
            },
        )

    raise ValueError(f"Unsupported SINDy optimizer: {optimizer!r}")


def rollout_controlled(
    model: object,
    a0: np.ndarray,
    u_seq: np.ndarray,
    dt: float,
    state_clip: float = 50.0,
) -> np.ndarray:
    """Roll out a controlled ROM trajectory with finite-value guards.

    Args:
        model: Fitted controlled model with ``predict(x, u)`` method.
        a0: Initial reduced state vector with shape ``(r,)``.
        u_seq: Control sequence with shape ``(T, n_u)``.
        dt: Integration step used in forward Euler rollout.
        state_clip: Absolute clip bound applied to predicted reduced states.

    Returns:
        Predicted reduced trajectory with shape ``(T, r)``.
    """

    a = np.asarray(a0, dtype=float).reshape(-1)
    u_seq = np.asarray(u_seq, dtype=float)
    n_steps = int(u_seq.shape[0])
    if n_steps <= 0:
        raise ValueError("u_seq must have at least one step")

    out = np.empty((n_steps, a.size), dtype=float)
    out[0] = a

    for k in range(n_steps - 1):
        u_k = u_seq[k : k + 1]
        try:
            adot = np.asarray(model.predict(a.reshape(1, -1), u=u_k), dtype=float).reshape(-1)
        except Exception:
            out[k + 1 :] = a
            return out

        if not np.all(np.isfinite(adot)):
            out[k + 1 :] = a
            return out

        a_next = a + dt * adot
        if state_clip > 0.0:
            a_next = np.clip(a_next, -state_clip, state_clip)

        if not np.all(np.isfinite(a_next)):
            out[k + 1 :] = a
            return out

        a = a_next
        out[k + 1] = a

    return out
