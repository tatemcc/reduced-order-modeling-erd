"""Metric definitions for the toy ERD plant.

All metrics are evaluated from discrete fields using the exact formulas adopted
for the toy ROM+MPC demonstration, including area weights ``r * dr * dphi``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List

import numpy as np

from .config import RunConfig


@dataclass(frozen=True)
class MetricState:
    """Per-step scalar metrics used for analysis and control.

    Attributes:
        E_wob: Combined wobble energy ``E1 + E2``.
        E1: Azimuthal mode-1 weighted energy of ``n``.
        E2: Azimuthal mode-2 weighted energy of ``n``.
        L_w: Wall-loss rate.
        sigma_r: Ring thickness (standard deviation in radius).
        r_mean: Mass-weighted mean radius.
        M: Axisymmetric mass proxy.
        E_u: Flow kinetic-energy proxy from streamfunction velocity.
        P_ctrl: Control power proxy from current actuator values.
    """

    E_wob: float
    E1: float
    E2: float
    L_w: float
    sigma_r: float
    r_mean: float
    M: float
    E_u: float
    P_ctrl: float


def _mode_energy(n: np.ndarray, r: np.ndarray, dr: float, sigma_w: float, r_star: float, m: int) -> float:
    """Compute weighted radial energy of azimuthal Fourier mode ``m``."""

    n_hat = np.fft.fft(n, axis=1) / n.shape[1]
    radial_weight = r * np.exp(-((r - r_star) ** 2) / (2.0 * sigma_w**2))
    e = 0.5 * np.sum(radial_weight * (np.abs(n_hat[:, m]) ** 2)) * dr
    return float(e)


def compute_metric_state(
    cfg: RunConfig,
    n: np.ndarray,
    kappa_w_r: np.ndarray,
    area_weights: np.ndarray,
    u_r: np.ndarray,
    u_phi: np.ndarray,
    u: np.ndarray,
    r: np.ndarray,
    dr: float,
) -> MetricState:
    """Evaluate all scalar metrics for one plant step.

    Args:
        cfg: Run configuration containing metric constants.
        n: Density field ``n(r, phi)``.
        kappa_w_r: Radial wall-loss profile.
        area_weights: Cell area weights ``r * dr * dphi``.
        u_r: Radial velocity component.
        u_phi: Azimuthal velocity component.
        u: Applied control vector ``[u0, ..., u4]``.
        r: Radial coordinate vector.
        dr: Radial grid spacing.

    Returns:
        :class:`MetricState` with all per-step scalar quantities.
    """

    mc = cfg.metrics

    e1 = _mode_energy(n, r, dr, mc.sigma_w, cfg.ring.r_star, m=1)
    e2 = _mode_energy(n, r, dr, mc.sigma_w, cfg.ring.r_star, m=2)
    e_wob = e1 + e2

    kappa_w = np.repeat(kappa_w_r[:, None], n.shape[1], axis=1)
    l_w = float(np.sum(kappa_w * n * area_weights))

    nbar = np.mean(n, axis=1)
    mass = float(np.sum(nbar * r) * dr)
    mass_safe = max(mass, 1e-12)

    r_mean = float(np.sum(r * nbar * r) * dr / mass_safe)
    sigma_r2 = float(np.sum(((r - r_mean) ** 2) * nbar * r) * dr / mass_safe)

    e_u = float(0.5 * np.sum((u_r**2 + u_phi**2) * area_weights))

    u = np.asarray(u).reshape(5)
    p_ctrl = float(
        u[0] ** 2
        + mc.lambda_1 * (u[1] ** 2 + u[2] ** 2)
        + mc.lambda_2 * (u[3] ** 2 + u[4] ** 2)
    )

    return MetricState(
        E_wob=e_wob,
        E1=e1,
        E2=e2,
        L_w=l_w,
        sigma_r=np.sqrt(max(sigma_r2, 0.0)),
        r_mean=r_mean,
        M=mass,
        E_u=e_u,
        P_ctrl=p_ctrl,
    )


def add_window_efficiency(curves: Dict[str, List[float]], cfg: RunConfig) -> None:
    """Append windowed efficiency ``eta_T`` to metric curves in place.

    Args:
        curves: Mutable metric curves dictionary produced during stepping.
        cfg: Run configuration containing efficiency window constants.

    Returns:
        None.
    """

    eps = cfg.metrics.eta_eps
    k = max(1, int(cfg.metrics.eta_window_steps))
    eu = np.asarray(curves["E_u"], dtype=float)
    pctrl = np.asarray(curves["P_ctrl"], dtype=float)

    eta = np.zeros_like(eu)
    n = eu.size
    for i in range(n):
        j1 = i
        j2 = min(n, i + k)
        eta[i] = float(np.sum(eu[j1:j2]) / (eps + np.sum(pctrl[j1:j2])))

    curves["eta_T"] = eta.tolist()


def add_degradation_metrics(curves: Dict[str, List[float]], cfg: RunConfig) -> None:
    """Append cumulative and excess degradation metrics in place.

    Added keys:
    - ``L_w_cum``: Cumulative wall-loss integral.
    - ``E_wob_excess``: Positive excess above initial wobble energy.
    - ``sigma_r_excess``: Positive excess above initial ring thickness.
    - ``badness_score``: Weighted combination of the three degradation terms.

    Args:
        curves: Mutable metric-curve dictionary produced during stepping.
        cfg: Run configuration containing degradation metric weights.

    Returns:
        None.
    """

    t = np.asarray(curves.get("t", []), dtype=float)
    e_wob = np.asarray(curves.get("E_wob", []), dtype=float)
    l_w = np.asarray(curves.get("L_w", []), dtype=float)
    sigma_r = np.asarray(curves.get("sigma_r", []), dtype=float)

    if t.size == 0:
        curves["L_w_cum"] = []
        curves["E_wob_excess"] = []
        curves["sigma_r_excess"] = []
        curves["badness_score"] = []
        return

    dt = np.diff(np.concatenate([[0.0], t]))
    dt = np.clip(dt, 0.0, None)

    l_w_cum = np.cumsum(l_w * dt)
    e_wob_excess = np.maximum(e_wob - e_wob[0], 0.0)
    sigma_r_excess = np.maximum(sigma_r - sigma_r[0], 0.0)

    mc = cfg.metrics
    badness = mc.wE * e_wob_excess + mc.wL * l_w_cum + mc.wS * sigma_r_excess

    curves["L_w_cum"] = l_w_cum.tolist()
    curves["E_wob_excess"] = e_wob_excess.tolist()
    curves["sigma_r_excess"] = sigma_r_excess.tolist()
    curves["badness_score"] = badness.tolist()


def dominant_spectrum(n: np.ndarray) -> np.ndarray:
    """Compute radial-averaged azimuthal power spectrum of ``n``.

    Args:
        n: Density field at one time step.

    Returns:
        1D power spectrum indexed by azimuthal mode number ``m``.
    """

    n_hat = np.fft.fft(n, axis=1) / n.shape[1]
    power = np.mean(np.abs(n_hat) ** 2, axis=0)
    return power


def compute_variation_diagnostics(
    cfg: RunConfig,
    n_snapshots: np.ndarray,
    r: np.ndarray,
    t_snap: np.ndarray,
    m_nonax_max: int = 24,
) -> Dict[str, Any]:
    """Compute variation diagnostics used for visual-activity acceptance checks.

    Diagnostics include:
    - non-axisymmetric energy ratio ``E_nonax / E0``,
    - low/mid/high azimuthal band power ratios in the ring band,
    - frame-to-frame relative field change ``Delta``.

    Args:
        cfg: Run configuration used to define ring-centered weights.
        n_snapshots: Density snapshots with shape ``(T, N_r, N_phi)``.
        r: Radial coordinate vector of length ``N_r``.
        t_snap: Snapshot time vector of length ``T``.
        m_nonax_max: Maximum azimuthal mode index included in ``E_nonax``.

    Returns:
        Dictionary containing per-time curves and summary scalars.
    """

    n_data = np.asarray(n_snapshots, dtype=float)
    r = np.asarray(r, dtype=float)
    t = np.asarray(t_snap, dtype=float)
    if n_data.ndim != 3 or n_data.shape[0] < 2:
        return {
            "ratio_nonax_over_axisym": [],
            "band_ratio_mid_to_low": [],
            "band_ratio_high_to_low": [],
            "delta_frame_l2_rel": [],
            "summary": {},
        }

    dr = float(np.mean(np.diff(r))) if r.size > 1 else 1.0
    n_hat = np.fft.fft(n_data, axis=2) / n_data.shape[2]

    w = r * np.exp(-((r - cfg.ring.r_star) ** 2) / (2.0 * cfg.metrics.sigma_w**2))
    mmax = min(int(m_nonax_max), n_data.shape[2] // 2 - 1)
    e0 = 0.5 * np.sum(w[None, :] * np.abs(n_hat[:, :, 0]) ** 2, axis=1) * dr
    enon = 0.5 * np.sum(w[None, :, None] * np.abs(n_hat[:, :, 1 : mmax + 1]) ** 2, axis=(1, 2)) * dr
    ratio_nonax = enon / np.maximum(e0, 1e-14)

    ring_band = np.where(np.abs(r - cfg.ring.r_star) <= cfg.ring.sigma_star)[0]
    if ring_band.size == 0:
        ring_band = np.arange(r.size)
    power = np.mean(np.abs(n_hat[:, ring_band, :]) ** 2, axis=1)

    low = np.sum(power[:, 1:4], axis=1)
    mid = np.sum(power[:, 6:13], axis=1)
    high = np.sum(power[:, 15 : min(25, power.shape[1])], axis=1)
    mid_to_low = mid / np.maximum(low, 1e-14)
    high_to_low = high / np.maximum(low, 1e-14)

    diff_num = np.linalg.norm(n_data[1:] - n_data[:-1], axis=(1, 2))
    diff_den = np.maximum(np.linalg.norm(n_data[:-1], axis=(1, 2)), 1e-14)
    delta = diff_num / diff_den

    first_twenty = max(1, int(0.2 * ratio_nonax.size))
    summary = {
        "ratio_nonax_over_axisym_first20_max": float(np.max(ratio_nonax[:first_twenty])),
        "ratio_nonax_over_axisym_first20_std": float(np.std(ratio_nonax[:first_twenty])),
        "ratio_nonax_over_axisym_global_max": float(np.max(ratio_nonax)),
        "band_ratio_mid_to_low_max": float(np.max(mid_to_low)),
        "band_ratio_high_to_low_max": float(np.max(high_to_low)),
        "delta_frame_l2_rel_median": float(np.median(delta)),
        "delta_frame_l2_rel_mean": float(np.mean(delta)),
    }

    return {
        "t": t.tolist(),
        "ratio_nonax_over_axisym": ratio_nonax.tolist(),
        "band_ratio_mid_to_low": mid_to_low.tolist(),
        "band_ratio_high_to_low": high_to_low.tolist(),
        "delta_frame_l2_rel_t": t[1:].tolist(),
        "delta_frame_l2_rel": delta.tolist(),
        "summary": summary,
    }
