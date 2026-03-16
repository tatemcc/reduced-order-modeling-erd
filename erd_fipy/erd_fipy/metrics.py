"""Metric definitions for the annular mHW ERD plant.

All metrics are evaluated from discrete fields using the run-folder contracts
consumed by the ROM and MPC pipelines.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List

import numpy as np

from .config import RunConfig
from .fields import ring_equilibrium


@dataclass(frozen=True)
class MetricState:
    """Per-step scalar metrics used for analysis and control.

    Attributes:
        J_prof: Profile-deviation energy of the zonal mean relative to ``n_eq``.
        E_low: Weighted low-mode energy summed over azimuthal modes ``m=1..4``.
        L_w: Wall-loss rate.
        sigma_r: Ring thickness (standard deviation in radius).
        r_mean: Mass-weighted mean radius.
        M: Axisymmetric mass proxy.
        E_u: Flow kinetic-energy proxy from streamfunction velocity.
        P_ctrl: Control power proxy from current actuator values.
    """

    J_prof: float
    E_low: float
    L_w: float
    sigma_r: float
    r_mean: float
    M: float
    E_u: float
    P_ctrl: float


def _low_mode_energy(n: np.ndarray, r: np.ndarray, dr: float, sigma_w: float, r_star: float, m_c: int = 4) -> float:
    """Compute weighted low-mode energy summed over azimuthal modes ``1..m_c``."""

    n_hat = np.fft.fft(n, axis=1) / n.shape[1]
    radial_weight = r * np.exp(-((r - r_star) ** 2) / (2.0 * sigma_w**2))
    e = 0.5 * np.sum(radial_weight[:, None] * np.abs(n_hat[:, 1 : m_c + 1]) ** 2) * dr
    return float(e)


def _diff_r_rows(data: np.ndarray, dr: float) -> np.ndarray:
    """Differentiate a ``(T, N_r)`` array along the radial axis.

    Args:
        data: Time-by-radius array.
        dr: Radial grid spacing.

    Returns:
        Array with the same shape as ``data`` containing centered radial
        derivatives with one-sided boundary stencils.
    """

    out = np.zeros_like(data)
    out[:, 1:-1] = (data[:, 2:] - data[:, :-2]) / (2.0 * dr)
    out[:, 0] = (data[:, 1] - data[:, 0]) / dr
    out[:, -1] = (data[:, -1] - data[:, -2]) / dr
    return out


def _diff_r_snapshots(data: np.ndarray, dr: float) -> np.ndarray:
    """Differentiate a ``(T, N_r, N_phi)`` array along the radial axis.

    Args:
        data: Time-by-radius-by-angle array.
        dr: Radial grid spacing.

    Returns:
        Array with the same shape as ``data`` containing centered radial
        derivatives with one-sided boundary stencils.
    """

    out = np.zeros_like(data)
    out[:, 1:-1, :] = (data[:, 2:, :] - data[:, :-2, :]) / (2.0 * dr)
    out[:, 0, :] = (data[:, 1, :] - data[:, 0, :]) / dr
    out[:, -1, :] = (data[:, -1, :] - data[:, -2, :]) / dr
    return out


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
    n_eq_r = ring_equilibrium(r, cfg)
    nbar = np.mean(n, axis=1)

    j_prof = 0.5 * np.sum(r * (nbar - n_eq_r) ** 2) * dr
    e_low = _low_mode_energy(n, r, dr, mc.sigma_w, cfg.ring.r_star, m_c=4)

    kappa_w = np.repeat(kappa_w_r[:, None], n.shape[1], axis=1)
    l_w = float(np.sum(kappa_w * n * area_weights))

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
        J_prof=float(j_prof),
        E_low=e_low,
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
    - ``J_prof_excess``: Positive excess above initial profile-deviation energy.
    - ``E_low_excess``: Positive excess above initial low-mode energy.
    - ``sigma_r_excess``: Positive excess above initial ring thickness.
    - ``badness_score``: Weighted combination of degradation terms.

    Args:
        curves: Mutable metric-curve dictionary produced during stepping.
        cfg: Run configuration containing degradation metric weights.

    Returns:
        None.
    """

    t = np.asarray(curves.get("t", []), dtype=float)
    j_prof = np.asarray(curves.get("J_prof", []), dtype=float)
    e_low = np.asarray(curves.get("E_low", []), dtype=float)
    l_w = np.asarray(curves.get("L_w", []), dtype=float)
    sigma_r = np.asarray(curves.get("sigma_r", []), dtype=float)

    if t.size == 0:
        curves["L_w_cum"] = []
        curves["J_prof_excess"] = []
        curves["E_low_excess"] = []
        curves["sigma_r_excess"] = []
        curves["badness_score"] = []
        return

    dt = np.diff(np.concatenate([[0.0], t]))
    dt = np.clip(dt, 0.0, None)

    l_w_cum = np.cumsum(l_w * dt)
    j_prof_excess = np.maximum(j_prof - j_prof[0], 0.0)
    e_low_excess = np.maximum(e_low - e_low[0], 0.0)
    sigma_r_excess = np.maximum(sigma_r - sigma_r[0], 0.0)

    mc = cfg.metrics
    badness = mc.wJ * j_prof_excess + mc.wE * e_low_excess + mc.wL * l_w_cum + mc.wS * sigma_r_excess

    curves["L_w_cum"] = l_w_cum.tolist()
    curves["J_prof_excess"] = j_prof_excess.tolist()
    curves["E_low_excess"] = e_low_excess.tolist()
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

    positive_power = np.clip(power[:, 1 : power.shape[1] // 2], 0.0, None)
    positive_sum = np.maximum(np.sum(positive_power, axis=1, keepdims=True), 1e-14)
    p_mode = positive_power / positive_sum
    entropy = -np.sum(p_mode * np.log(np.maximum(p_mode, 1e-14)), axis=1)
    if positive_power.shape[1] > 1:
        entropy /= np.log(float(positive_power.shape[1]))
    else:
        entropy[:] = 0.0

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
        "spectral_entropy_mean": float(np.mean(entropy)),
        "spectral_entropy_peak": float(np.max(entropy)),
        "delta_frame_l2_rel_median": float(np.median(delta)),
        "delta_frame_l2_rel_mean": float(np.mean(delta)),
    }

    return {
        "t": t.tolist(),
        "ratio_nonax_over_axisym": ratio_nonax.tolist(),
        "band_ratio_mid_to_low": mid_to_low.tolist(),
        "band_ratio_high_to_low": high_to_low.tolist(),
        "spectral_entropy": entropy.tolist(),
        "delta_frame_l2_rel_t": t[1:].tolist(),
        "delta_frame_l2_rel": delta.tolist(),
        "summary": summary,
    }


def compute_transport_diagnostics(
    cfg: RunConfig,
    n_snapshots: np.ndarray,
    omega_snapshots: np.ndarray,
    u_mag_snapshots: np.ndarray | None,
    u_r_snapshots: np.ndarray | None,
    u_phi_snapshots: np.ndarray | None,
    r: np.ndarray,
    t_snap: np.ndarray,
) -> Dict[str, Any]:
    """Compute transport-oriented diagnostics for qualitative tuning.

    These diagnostics complement the scalar degradation metrics by measuring
    motion, shear balance, threshold-to-instability margin, and burst onset.

    Args:
        cfg: Run configuration used to define the ring band.
        n_snapshots: Density snapshots with shape ``(T, N_r, N_phi)``.
        omega_snapshots: Vorticity snapshots with shape ``(T, N_r, N_phi)``.
        u_mag_snapshots: Optional velocity-magnitude snapshots with shape
            ``(T, N_r, N_phi)``.
        u_r_snapshots: Optional radial-velocity snapshots with shape
            ``(T, N_r, N_phi)``.
        u_phi_snapshots: Optional azimuthal-velocity snapshots with shape
            ``(T, N_r, N_phi)``.
        r: Radial coordinate vector of length ``N_r``.
        t_snap: Snapshot time vector of length ``T``.

    Returns:
        Dictionary containing time traces and summary scalars.
    """

    n_data = np.asarray(n_snapshots, dtype=float)
    omega_data = np.asarray(omega_snapshots, dtype=float)
    u_data = None if u_mag_snapshots is None else np.asarray(u_mag_snapshots, dtype=float)
    u_r_data = None if u_r_snapshots is None else np.asarray(u_r_snapshots, dtype=float)
    u_phi_data = None if u_phi_snapshots is None else np.asarray(u_phi_snapshots, dtype=float)
    r = np.asarray(r, dtype=float)
    t = np.asarray(t_snap, dtype=float)

    if n_data.ndim != 3 or n_data.shape[0] < 2:
        return {"t": t.tolist(), "summary": {}}

    dt = np.diff(t)
    dt_safe = np.maximum(dt, 1.0e-12)
    dphi = float(2.0 * np.pi / n_data.shape[2]) if n_data.shape[2] > 0 else 1.0

    def _rel_delta(data: np.ndarray) -> np.ndarray:
        prev = data[:-1].reshape(data.shape[0] - 1, -1)
        nxt = data[1:].reshape(data.shape[0] - 1, -1)
        numer = np.linalg.norm(nxt - prev, axis=1)
        denom = np.maximum(np.linalg.norm(prev, axis=1), 1.0e-12)
        return numer / denom

    delta_n = _rel_delta(n_data)
    delta_omega = _rel_delta(omega_data)

    ring_mask = np.abs(r - cfg.ring.r_star) <= max(cfg.ring.sigma_star, 1.0e-12)
    if not np.any(ring_mask):
        ring_mask = np.ones_like(r, dtype=bool)

    if u_data is not None and u_data.shape == n_data.shape:
        u_ring = u_data[:, ring_mask, :]
        u_ring_mean = np.mean(u_ring, axis=(1, 2))
        u_ring_p95 = np.percentile(u_ring, 95.0, axis=(1, 2))
    else:
        u_ring_mean = np.zeros(n_data.shape[0], dtype=float)
        u_ring_p95 = np.zeros(n_data.shape[0], dtype=float)

    if u_r_data is not None and u_r_data.shape == n_data.shape:
        u_r_ring = np.abs(u_r_data[:, ring_mask, :])
        u_r_ring_abs_mean = np.mean(u_r_ring, axis=(1, 2))
        u_r_ring_p95 = np.percentile(u_r_ring, 95.0, axis=(1, 2))
        gamma_r_ring = n_data[:, ring_mask, :] * u_r_data[:, ring_mask, :]
        gamma_r_ring_signed_mean = np.mean(gamma_r_ring, axis=(1, 2))
        gamma_r_ring_abs_q90 = np.percentile(np.abs(gamma_r_ring), 90.0, axis=(1, 2))
    else:
        u_r_ring_abs_mean = np.zeros(n_data.shape[0], dtype=float)
        u_r_ring_p95 = np.zeros(n_data.shape[0], dtype=float)
        gamma_r_ring_signed_mean = np.zeros(n_data.shape[0], dtype=float)
        gamma_r_ring_abs_q90 = np.zeros(n_data.shape[0], dtype=float)

    if u_phi_data is not None and u_phi_data.shape == n_data.shape:
        u_phi_pos_frac = np.mean(u_phi_data[:, ring_mask, :] > 0.0, axis=(1, 2))
        u_phi_neg_frac = np.mean(u_phi_data[:, ring_mask, :] < 0.0, axis=(1, 2))
        u_phi_ring_signed = np.mean(u_phi_data[:, ring_mask, :], axis=(1, 2))
        u_phi_ring = np.abs(u_phi_data[:, ring_mask, :])
        u_phi_ring_abs_mean = np.mean(u_phi_ring, axis=(1, 2))
        u_phi_ring_p95 = np.percentile(u_phi_ring, 95.0, axis=(1, 2))
        u_phi_zonal = np.mean(u_phi_data, axis=2)
    else:
        u_phi_pos_frac = np.zeros(n_data.shape[0], dtype=float)
        u_phi_neg_frac = np.zeros(n_data.shape[0], dtype=float)
        u_phi_ring_signed = np.zeros(n_data.shape[0], dtype=float)
        u_phi_ring_abs_mean = np.zeros(n_data.shape[0], dtype=float)
        u_phi_ring_p95 = np.zeros(n_data.shape[0], dtype=float)
        u_phi_zonal = np.zeros((n_data.shape[0], n_data.shape[1]), dtype=float)

    u_phi_over_u_r = u_phi_ring_abs_mean / np.maximum(u_r_ring_abs_mean, 1.0e-12)
    u_phi_sign = np.sign(u_phi_ring_signed)
    u_phi_sign[np.abs(u_phi_ring_signed) < 1.0e-10] = 0.0
    u_phi_sign_changes = int(np.count_nonzero(np.diff(u_phi_sign) != 0.0))
    u_phi_signed_bias = float(np.mean(u_phi_ring_signed) / np.maximum(np.mean(u_phi_ring_abs_mean), 1.0e-12))

    n_ring = n_data[:, ring_mask, :]
    n_hat = np.fft.fft(n_ring, axis=2) / n_ring.shape[2]
    m1 = np.mean(n_hat[:, :, 1], axis=1)
    m1_amp = np.abs(m1)
    m1_phase = np.unwrap(np.angle(m1))
    m1_phase_rate = np.diff(m1_phase) / dt_safe if m1_phase.size > 1 else np.zeros(0, dtype=float)

    nbar = np.mean(n_data, axis=2)
    mass = np.sum(nbar * r[None, :], axis=1)
    mass_safe = np.maximum(mass, 1.0e-12)
    r_mean = np.sum((r[None, :] ** 2) * nbar, axis=1) / mass_safe

    dr = float(np.mean(np.diff(r))) if r.size > 1 else 1.0
    grad_ref = max(cfg.ring.n_amp / max(cfg.ring.sigma_star, 1.0e-12), 1.0e-12)
    grad_ratio = -_diff_r_rows(nbar, dr) / grad_ref
    grad_abs = np.abs(grad_ratio)
    shear_ref = max(float(getattr(cfg.pde, "shear_ref", 1.0)), 1.0e-8)
    shear_norm = np.abs(_diff_r_rows(u_phi_zonal, dr)) / shear_ref
    crit = float(getattr(cfg.pde, "critical_gradient_ratio", 0.0))
    shear_gain = float(getattr(cfg.pde, "shear_suppression_gain", 0.0))
    width = max(float(getattr(cfg.pde, "threshold_width", 0.08)), 1.0e-8)
    instability_margin = grad_abs - crit - shear_gain * shear_norm
    instability_activation = 0.5 * (
        instability_margin + np.sqrt(instability_margin * instability_margin + width * width)
    )
    phase_drive = float(getattr(cfg.pde, "phase_advection_gain", 0.0)) * instability_activation * grad_ratio
    coupling_boost = 1.0 + float(getattr(cfg.pde, "supercritical_coupling_gain", 0.0)) * instability_activation
    feedback_gain = float(getattr(cfg.pde, "supercritical_feedback_gain", 0.0)) * instability_activation
    transport_boost = 1.0 + float(getattr(cfg.pde, "supercritical_transport_gain", 0.0)) * instability_activation

    ring_grad_ratio = np.mean(grad_ratio[:, ring_mask], axis=1)
    ring_grad_abs = np.mean(grad_abs[:, ring_mask], axis=1)
    ring_shear_norm = np.mean(shear_norm[:, ring_mask], axis=1)
    ring_instability_margin = np.mean(instability_margin[:, ring_mask], axis=1)
    ring_instability_activation = np.mean(instability_activation[:, ring_mask], axis=1)
    ring_phase_drive_signed = np.mean(phase_drive[:, ring_mask], axis=1)
    ring_phase_drive_abs = np.mean(np.abs(phase_drive[:, ring_mask]), axis=1)
    ring_coupling_boost = np.mean(coupling_boost[:, ring_mask], axis=1)
    ring_feedback_gain = np.mean(feedback_gain[:, ring_mask], axis=1)
    ring_transport_boost = np.mean(transport_boost[:, ring_mask], axis=1)
    supercritical_fraction = np.mean(instability_margin[:, ring_mask] > 0.0, axis=1)

    n_tilde = n_data - nbar[:, :, None]
    omega_tilde = omega_data - np.mean(omega_data, axis=2, keepdims=True)
    n_ref = max(cfg.ring.n_amp, 1.0e-8)
    omega_ref = max(12.0, 1.0)
    omega_grad_ref = max(omega_ref / max(cfg.ring.sigma_star, 1.0e-12), 1.0e-8)
    dn_tilde_dphi = (np.roll(n_tilde, -1, axis=2) - np.roll(n_tilde, 1, axis=2)) / (2.0 * dphi)
    domega_tilde_dphi = (np.roll(omega_tilde, -1, axis=2) - np.roll(omega_tilde, 1, axis=2)) / (2.0 * dphi)
    dn_tilde_dr = _diff_r_snapshots(n_tilde, dr)
    domega_tilde_dr = _diff_r_snapshots(omega_tilde, dr)

    packet_structure_phi = (
        0.55 * (omega_tilde / omega_ref)
        + 0.30 * (dn_tilde_dphi / n_ref)
        + 0.15 * (domega_tilde_dphi / omega_ref)
    )
    packet_structure_r = (
        0.55 * (dn_tilde_dr / grad_ref)
        - 0.25 * (domega_tilde_dphi / omega_ref)
        + 0.20 * (domega_tilde_dr / omega_grad_ref)
    )
    wave_burst_gain = (
        1.0
        + float(getattr(cfg.pde, "wave_burst_speed_gain", 0.0)) * instability_activation[:, :, None]
        + 0.35 * np.tanh(np.abs(omega_tilde) / omega_ref)
    )
    wave_phi_speed = phase_drive[:, :, None] + float(getattr(cfg.pde, "wave_packet_phi_gain", 0.0)) * wave_burst_gain * np.tanh(packet_structure_phi)
    wave_r_speed = float(getattr(cfg.pde, "wave_packet_r_gain", 0.0)) * wave_burst_gain * np.tanh(packet_structure_r)
    wave_packet_speed = np.sqrt(wave_phi_speed**2 + wave_r_speed**2)
    ring_wave_phi_speed_signed = np.mean(wave_phi_speed[:, ring_mask, :], axis=(1, 2))
    ring_wave_phi_speed_abs = np.mean(np.abs(wave_phi_speed[:, ring_mask, :]), axis=(1, 2))
    ring_wave_r_speed_signed = np.mean(wave_r_speed[:, ring_mask, :], axis=(1, 2))
    ring_wave_r_speed_abs = np.mean(np.abs(wave_r_speed[:, ring_mask, :]), axis=(1, 2))
    ring_wave_packet_speed = np.mean(wave_packet_speed[:, ring_mask, :], axis=(1, 2))
    ring_wave_burst_gain = np.mean(wave_burst_gain[:, ring_mask, :], axis=(1, 2))

    wave_phi_sign = np.sign(ring_wave_phi_speed_signed)
    wave_phi_sign[np.abs(ring_wave_phi_speed_signed) < 1.0e-10] = 0.0
    wave_phi_sign_changes = int(np.count_nonzero(np.diff(wave_phi_sign) != 0.0))

    r_off = (r - cfg.ring.r_star) / max(cfg.ring.sigma_star, 1.0e-12)
    confinement_shape = (-r_off * np.exp(-0.5 * (r_off / 1.6) ** 2))[None, :, None]
    confinement_release = 1.0 / (
        1.0 + float(getattr(cfg.pde, "radial_release_gain", 0.0)) * instability_activation[:, :, None]
    )
    confinement_v = float(getattr(cfg.pde, "radial_confinement_gain", 0.0)) * confinement_release * confinement_shape
    ring_confinement_abs = np.mean(np.abs(confinement_v[:, ring_mask, :]), axis=(1, 2))
    ring_confinement_signed = np.mean(confinement_v[:, ring_mask, :], axis=(1, 2))
    ring_confinement_release = np.mean(confinement_release[:, ring_mask, :], axis=(1, 2))

    zonal_restore_gain = float(getattr(cfg.pde, "zonal_profile_restore_gain", 0.0))
    zonal_release_gain = float(getattr(cfg.pde, "zonal_profile_release_gain", 0.0))
    zonal_damage_release_gain = float(getattr(cfg.pde, "zonal_profile_damage_release_gain", 0.0))
    zonal_damage = np.abs(nbar - ring_equilibrium(r, cfg)[None, :]) / max(cfg.ring.n_amp, 1.0e-8)
    zonal_damage = np.clip(zonal_damage, 0.0, 3.0)
    zonal_restore_release = 1.0 / (
        1.0 + zonal_release_gain * instability_activation + zonal_damage_release_gain * zonal_damage
    )
    zonal_restore = zonal_restore_gain * zonal_restore_release * (ring_equilibrium(r, cfg)[None, :] - nbar)
    ring_zonal_restore_abs = np.mean(np.abs(zonal_restore[:, ring_mask]), axis=1)
    ring_zonal_restore_signed = np.mean(zonal_restore[:, ring_mask], axis=1)
    ring_zonal_restore_release = np.mean(zonal_restore_release[:, ring_mask], axis=1)

    mass_r = np.maximum(nbar, 0.0) * r[None, :]
    total_mass_r = np.maximum(np.sum(mass_r, axis=1), 1.0e-12)
    edge_width = max(2.0 * dr, 0.6 * cfg.ring.sigma_star)
    inner_mask = r <= (r[0] + edge_width)
    outer_mask = r >= (r[-1] - edge_width)
    core_mask = np.abs(r - cfg.ring.r_star) <= max(1.2 * cfg.ring.sigma_star, edge_width)
    inner_edge_mass_frac = np.sum(mass_r[:, inner_mask], axis=1) / total_mass_r
    outer_edge_mass_frac = np.sum(mass_r[:, outer_mask], axis=1) / total_mass_r
    ring_core_mass_frac = np.sum(mass_r[:, core_mask], axis=1) / total_mass_r

    row_abs_corr = np.zeros(n_data.shape[0], dtype=float)
    row_adj_corr = np.zeros(n_data.shape[0], dtype=float)
    n_tilde_rms = np.zeros(n_data.shape[0], dtype=float)
    for idx, snap in enumerate(n_data):
        snap_phi = snap - np.mean(snap, axis=1, keepdims=True)
        n_tilde_rms[idx] = float(np.sqrt(np.mean(snap_phi**2)))
        row = snap_phi / np.maximum(np.linalg.norm(snap_phi, axis=1, keepdims=True), 1.0e-12)
        corr = row @ row.T
        iu = np.triu_indices_from(corr, k=1)
        row_abs_corr[idx] = float(np.mean(np.abs(corr[iu]))) if iu[0].size else 0.0
        if row.shape[0] > 1:
            row_adj_corr[idx] = float(np.mean(np.sum(row[:-1] * row[1:], axis=1)))

    n_third = max(1, delta_n.size // 3)
    delta_growth_ratio = float(np.median(delta_n[-n_third:]) / np.maximum(np.median(delta_n[:n_third]), 1.0e-12))
    m1_growth_factor = float(m1_amp[-1] / np.maximum(m1_amp[0], 1.0e-12))
    u_ring_accel = np.diff(u_ring_mean) / dt_safe if u_ring_mean.size > 1 else np.zeros(0, dtype=float)
    peak_time_frac_j = float(t[np.argmax(np.sum((nbar - ring_equilibrium(r, cfg)[None, :]) ** 2, axis=1))] / max(t[-1], 1.0e-12))
    peak_time_frac_e = float(t[np.argmax(m1_amp)] / max(t[-1], 1.0e-12))
    last_third = slice(max(0, (2 * n_data.shape[0]) // 3), n_data.shape[0])

    flux_ref = float(np.median(gamma_r_ring_abs_q90[: max(1, gamma_r_ring_abs_q90.size // 3)]))
    flux_ref = max(flux_ref, 1.0e-12)
    flux_burst_threshold = 2.5 * flux_ref
    flux_burst_mask = gamma_r_ring_abs_q90 > flux_burst_threshold
    flux_burst_count = int(np.count_nonzero(np.logical_and(flux_burst_mask[1:], ~flux_burst_mask[:-1])))

    instability_active_mask = ring_instability_activation > 0.05
    instability_burst_count = int(
        np.count_nonzero(np.logical_and(instability_active_mask[1:], ~instability_active_mask[:-1]))
    )

    summary = {
        "delta_n_rel_median": float(np.median(delta_n)),
        "delta_n_rel_max": float(np.max(delta_n)),
        "delta_omega_rel_median": float(np.median(delta_omega)),
        "delta_omega_rel_max": float(np.max(delta_omega)),
        "u_ring_mean_peak": float(np.max(u_ring_mean)),
        "u_ring_mean_final": float(u_ring_mean[-1]),
        "u_ring_p95_peak": float(np.max(u_ring_p95)),
        "u_r_ring_abs_mean_peak": float(np.max(u_r_ring_abs_mean)),
        "u_r_ring_abs_mean_final": float(u_r_ring_abs_mean[-1]),
        "u_r_ring_p95_peak": float(np.max(u_r_ring_p95)),
        "gamma_r_ring_signed_peak": float(np.max(np.abs(gamma_r_ring_signed_mean))),
        "gamma_r_ring_abs_q90_peak": float(np.max(gamma_r_ring_abs_q90)),
        "gamma_r_ring_abs_q90_final": float(gamma_r_ring_abs_q90[-1]),
        "gamma_r_burst_factor": float(np.max(gamma_r_ring_abs_q90) / flux_ref),
        "gamma_r_burst_count": flux_burst_count,
        "u_phi_ring_abs_mean_peak": float(np.max(u_phi_ring_abs_mean)),
        "u_phi_ring_abs_mean_final": float(u_phi_ring_abs_mean[-1]),
        "u_phi_ring_p95_peak": float(np.max(u_phi_ring_p95)),
        "u_phi_ring_signed_peak": float(np.max(np.abs(u_phi_ring_signed))),
        "u_phi_over_u_r_peak_ratio": float(np.max(u_phi_over_u_r)),
        "u_phi_over_u_r_mean_ratio": float(np.mean(u_phi_over_u_r)),
        "u_phi_over_u_r_final_ratio": float(u_phi_over_u_r[-1]),
        "u_phi_sign_changes": u_phi_sign_changes,
        "u_phi_signed_bias": u_phi_signed_bias,
        "u_phi_pos_frac_mean": float(np.mean(u_phi_pos_frac)),
        "u_phi_sign_mix_peak": float(np.max(np.minimum(u_phi_pos_frac, u_phi_neg_frac))),
        "u_phi_sign_balance_mean": float(np.mean(np.abs(u_phi_pos_frac - u_phi_neg_frac))),
        "m1_amp_peak": float(np.max(m1_amp)),
        "m1_amp_final": float(m1_amp[-1]),
        "m1_phase_drift_total": float(abs(m1_phase[-1] - m1_phase[0])),
        "m1_phase_rate_rms": float(np.sqrt(np.mean(m1_phase_rate**2))) if m1_phase_rate.size else 0.0,
        "u_ring_accel_rms": float(np.sqrt(np.mean(u_ring_accel**2))) if u_ring_accel.size else 0.0,
        "u_ring_accel_peak": float(np.max(np.abs(u_ring_accel))) if u_ring_accel.size else 0.0,
        "r_mean_shift": float(r_mean[-1] - r_mean[0]),
        "mass_rel_drift": float((mass[-1] - mass[0]) / max(abs(mass[0]), 1.0e-12)),
        "mass_rel_span": float((np.max(mass) - np.min(mass)) / max(abs(mass[0]), 1.0e-12)),
        "n_tilde_rms_peak": float(np.max(n_tilde_rms)),
        "n_tilde_rms_final": float(n_tilde_rms[-1]),
        "n_tilde_retention_final_over_peak": float(n_tilde_rms[-1] / max(np.max(n_tilde_rms), 1.0e-12)),
        "row_abs_corr_mean": float(np.mean(row_abs_corr)),
        "row_abs_corr_final": float(row_abs_corr[-1]),
        "row_adjacent_corr_mean": float(np.mean(row_adj_corr)),
        "row_adjacent_corr_final": float(row_adj_corr[-1]),
        "delta_n_growth_ratio": delta_growth_ratio,
        "m1_growth_factor": m1_growth_factor,
        "ring_grad_ratio_peak": float(np.max(ring_grad_ratio)),
        "ring_grad_ratio_final": float(ring_grad_ratio[-1]),
        "ring_grad_abs_peak": float(np.max(ring_grad_abs)),
        "ring_grad_abs_final": float(ring_grad_abs[-1]),
        "ring_shear_norm_peak": float(np.max(ring_shear_norm)),
        "ring_shear_norm_final": float(ring_shear_norm[-1]),
        "ring_instability_margin_peak": float(np.max(ring_instability_margin)),
        "ring_instability_margin_final": float(ring_instability_margin[-1]),
        "ring_instability_activation_peak": float(np.max(ring_instability_activation)),
        "ring_instability_activation_final": float(ring_instability_activation[-1]),
        "ring_phase_drive_abs_peak": float(np.max(ring_phase_drive_abs)),
        "ring_phase_drive_abs_final": float(ring_phase_drive_abs[-1]),
        "ring_phase_drive_signed_peak": float(np.max(np.abs(ring_phase_drive_signed))),
        "ring_phase_drive_signed_final": float(ring_phase_drive_signed[-1]),
        "ring_wave_phi_speed_abs_peak": float(np.max(ring_wave_phi_speed_abs)),
        "ring_wave_phi_speed_abs_final": float(ring_wave_phi_speed_abs[-1]),
        "ring_wave_phi_speed_signed_peak": float(np.max(np.abs(ring_wave_phi_speed_signed))),
        "ring_wave_phi_speed_signed_final": float(ring_wave_phi_speed_signed[-1]),
        "ring_wave_r_speed_abs_peak": float(np.max(ring_wave_r_speed_abs)),
        "ring_wave_r_speed_abs_final": float(ring_wave_r_speed_abs[-1]),
        "ring_wave_r_speed_signed_peak": float(np.max(np.abs(ring_wave_r_speed_signed))),
        "ring_wave_r_speed_signed_final": float(ring_wave_r_speed_signed[-1]),
        "ring_wave_packet_speed_peak": float(np.max(ring_wave_packet_speed)),
        "ring_wave_packet_speed_final": float(ring_wave_packet_speed[-1]),
        "ring_wave_burst_gain_peak": float(np.max(ring_wave_burst_gain)),
        "ring_wave_burst_gain_final": float(ring_wave_burst_gain[-1]),
        "wave_phi_sign_changes": wave_phi_sign_changes,
        "ring_confinement_abs_peak": float(np.max(ring_confinement_abs)),
        "ring_confinement_abs_final": float(ring_confinement_abs[-1]),
        "ring_confinement_signed_peak": float(np.max(np.abs(ring_confinement_signed))),
        "ring_confinement_signed_final": float(ring_confinement_signed[-1]),
        "ring_confinement_release_peak": float(np.max(ring_confinement_release)),
        "ring_confinement_release_final": float(ring_confinement_release[-1]),
        "ring_zonal_restore_abs_peak": float(np.max(ring_zonal_restore_abs)),
        "ring_zonal_restore_abs_final": float(ring_zonal_restore_abs[-1]),
        "ring_zonal_restore_signed_peak": float(np.max(np.abs(ring_zonal_restore_signed))),
        "ring_zonal_restore_signed_final": float(ring_zonal_restore_signed[-1]),
        "ring_zonal_restore_release_peak": float(np.max(ring_zonal_restore_release)),
        "ring_zonal_restore_release_final": float(ring_zonal_restore_release[-1]),
        "ring_coupling_boost_peak": float(np.max(ring_coupling_boost)),
        "ring_coupling_boost_final": float(ring_coupling_boost[-1]),
        "ring_feedback_gain_peak": float(np.max(ring_feedback_gain)),
        "ring_feedback_gain_final": float(ring_feedback_gain[-1]),
        "ring_transport_boost_peak": float(np.max(ring_transport_boost)),
        "ring_transport_boost_final": float(ring_transport_boost[-1]),
        "inner_edge_mass_frac_peak": float(np.max(inner_edge_mass_frac)),
        "inner_edge_mass_frac_final": float(inner_edge_mass_frac[-1]),
        "outer_edge_mass_frac_peak": float(np.max(outer_edge_mass_frac)),
        "outer_edge_mass_frac_final": float(outer_edge_mass_frac[-1]),
        "ring_core_mass_frac_peak": float(np.max(ring_core_mass_frac)),
        "ring_core_mass_frac_final": float(ring_core_mass_frac[-1]),
        "J_peak_time_fraction": peak_time_frac_j,
        "E_peak_time_fraction": peak_time_frac_e,
        "inner_edge_mass_frac_late_mean": float(np.mean(inner_edge_mass_frac[last_third])),
        "outer_edge_mass_frac_late_mean": float(np.mean(outer_edge_mass_frac[last_third])),
        "ring_core_mass_frac_late_mean": float(np.mean(ring_core_mass_frac[last_third])),
        "supercritical_fraction_peak": float(np.max(supercritical_fraction)),
        "supercritical_fraction_final": float(supercritical_fraction[-1]),
        "supercritical_fraction_mean": float(np.mean(supercritical_fraction)),
        "instability_burst_count": instability_burst_count,
    }

    return {
        "t": t.tolist(),
        "delta_t": t[1:].tolist(),
        "delta_n_rel": delta_n.tolist(),
        "delta_omega_rel": delta_omega.tolist(),
        "u_ring_mean": u_ring_mean.tolist(),
        "u_ring_p95": u_ring_p95.tolist(),
        "u_r_ring_abs_mean": u_r_ring_abs_mean.tolist(),
        "u_r_ring_p95": u_r_ring_p95.tolist(),
        "gamma_r_ring_signed_mean": gamma_r_ring_signed_mean.tolist(),
        "gamma_r_ring_abs_q90": gamma_r_ring_abs_q90.tolist(),
        "u_phi_pos_frac": u_phi_pos_frac.tolist(),
        "u_phi_neg_frac": u_phi_neg_frac.tolist(),
        "u_phi_ring_signed": u_phi_ring_signed.tolist(),
        "u_phi_ring_abs_mean": u_phi_ring_abs_mean.tolist(),
        "u_phi_ring_p95": u_phi_ring_p95.tolist(),
        "u_phi_over_u_r": u_phi_over_u_r.tolist(),
        "m1_amp": m1_amp.tolist(),
        "m1_phase": m1_phase.tolist(),
        "m1_phase_rate_t": t[1:].tolist(),
        "m1_phase_rate": m1_phase_rate.tolist(),
        "u_ring_accel_t": t[1:].tolist(),
        "u_ring_accel": u_ring_accel.tolist(),
        "r_mean": r_mean.tolist(),
        "mass": mass.tolist(),
        "n_tilde_rms": n_tilde_rms.tolist(),
        "row_abs_corr": row_abs_corr.tolist(),
        "row_adjacent_corr": row_adj_corr.tolist(),
        "ring_grad_ratio": ring_grad_ratio.tolist(),
        "ring_grad_abs": ring_grad_abs.tolist(),
        "ring_shear_norm": ring_shear_norm.tolist(),
        "ring_instability_margin": ring_instability_margin.tolist(),
        "ring_instability_activation": ring_instability_activation.tolist(),
        "ring_phase_drive_signed": ring_phase_drive_signed.tolist(),
        "ring_phase_drive_abs": ring_phase_drive_abs.tolist(),
        "ring_wave_phi_speed_signed": ring_wave_phi_speed_signed.tolist(),
        "ring_wave_phi_speed_abs": ring_wave_phi_speed_abs.tolist(),
        "ring_wave_r_speed_signed": ring_wave_r_speed_signed.tolist(),
        "ring_wave_r_speed_abs": ring_wave_r_speed_abs.tolist(),
        "ring_wave_packet_speed": ring_wave_packet_speed.tolist(),
        "ring_wave_burst_gain": ring_wave_burst_gain.tolist(),
        "ring_confinement_signed": ring_confinement_signed.tolist(),
        "ring_confinement_abs": ring_confinement_abs.tolist(),
        "ring_confinement_release": ring_confinement_release.tolist(),
        "ring_zonal_restore_signed": ring_zonal_restore_signed.tolist(),
        "ring_zonal_restore_abs": ring_zonal_restore_abs.tolist(),
        "ring_zonal_restore_release": ring_zonal_restore_release.tolist(),
        "ring_coupling_boost": ring_coupling_boost.tolist(),
        "ring_feedback_gain": ring_feedback_gain.tolist(),
        "ring_transport_boost": ring_transport_boost.tolist(),
        "inner_edge_mass_frac": inner_edge_mass_frac.tolist(),
        "outer_edge_mass_frac": outer_edge_mass_frac.tolist(),
        "ring_core_mass_frac": ring_core_mass_frac.tolist(),
        "supercritical_fraction": supercritical_fraction.tolist(),
        "summary": summary,
    }
