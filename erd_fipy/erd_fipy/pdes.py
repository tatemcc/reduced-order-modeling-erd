"""Numerical operators for the annular modified Hasegawa-Wakatani ERD plant.

The implementation uses rectangularized finite differences on a uniform
``(r, phi)`` grid with:
- periodic differences in ``phi``,
- Neumann radial boundaries for ``n`` and ``omega``,
- Dirichlet radial boundaries for ``psi`` in the Poisson solve,
- an Arakawa-type bracket for nonlinear transport,
- RK4 time stepping for the coupled ``(n, omega)`` dynamics.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import DisturbanceModeConfig, RunConfig
from .fields import ring_equilibrium, wall_loss_profile
from .mesh import MeshBundle


@dataclass(frozen=True)
class BasisFields:
    """Precomputed radial/azimuthal forcing basis arrays.

    Attributes:
        b0: Axisymmetric drive basis.
        b1: ``m=1`` cosine basis.
        b2: ``m=1`` sine basis.
        b3: ``m=2`` cosine basis.
        b4: ``m=2`` sine basis.
        g1: Radial envelope for ``m=1`` terms.
        g2: Radial envelope for ``m=2`` terms.
    """

    b0: np.ndarray
    b1: np.ndarray
    b2: np.ndarray
    b3: np.ndarray
    b4: np.ndarray
    g1: np.ndarray
    g2: np.ndarray


@dataclass(frozen=True)
class StepDiagnostics:
    """Diagnostics produced during one PDE update.

    Attributes:
        forcing: Total vorticity forcing field used in the step.
        u_r: Radial velocity component used for density update.
        u_phi: Azimuthal velocity component used for density update.
    """

    forcing: np.ndarray
    u_r: np.ndarray
    u_phi: np.ndarray


def _diff_phi_periodic(f: np.ndarray, dphi: float) -> np.ndarray:
    """Centered first derivative in ``phi`` with periodic wrapping."""

    return (np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1)) / (2.0 * dphi)


def _d2_phi_periodic(f: np.ndarray, dphi: float) -> np.ndarray:
    """Centered second derivative in ``phi`` with periodic wrapping."""

    return (np.roll(f, -1, axis=1) - 2.0 * f + np.roll(f, 1, axis=1)) / (dphi**2)


def _diff_r_centered(f: np.ndarray, dr: float) -> np.ndarray:
    """Centered first derivative in ``r`` with one-sided boundary stencils."""

    out = np.zeros_like(f)
    out[1:-1] = (f[2:] - f[:-2]) / (2.0 * dr)
    out[0] = (f[1] - f[0]) / dr
    out[-1] = (f[-1] - f[-2]) / dr
    return out


def _d2_r_neumann(f: np.ndarray, dr: float) -> np.ndarray:
    """Second derivative in ``r`` consistent with radial Neumann boundaries."""

    out = np.zeros_like(f)
    out[1:-1] = (f[2:] - 2.0 * f[1:-1] + f[:-2]) / (dr**2)
    out[0] = 2.0 * (f[1] - f[0]) / (dr**2)
    out[-1] = 2.0 * (f[-2] - f[-1]) / (dr**2)
    return out


def _enforce_neumann_radial(f: np.ndarray) -> None:
    """Project field values onto zero radial-gradient boundaries in place."""

    f[0] = f[1]
    f[-1] = f[-2]


def _gauss(r: np.ndarray, center: float, sigma: float) -> np.ndarray:
    """Evaluate a radial Gaussian envelope."""

    return np.exp(-((r - center) ** 2) / (2.0 * sigma**2))


def _smooth_positive(x: np.ndarray, width: float) -> np.ndarray:
    """Smooth approximation of ``max(x, 0)`` with width ``width``.

    Args:
        x: Input array.
        width: Transition width controlling how sharp the activation is.

    Returns:
        Array that is approximately zero for negative ``x`` and linear for
        sufficiently positive ``x``.
    """

    width = max(float(width), 1.0e-8)
    return 0.5 * (x + np.sqrt(x * x + width * width))


def _periodic_angle_delta(phi: np.ndarray, center: float) -> np.ndarray:
    """Return wrapped angular distance ``phi - center`` in ``[-pi, pi)``."""

    return ((phi - center + np.pi) % (2.0 * np.pi)) - np.pi


def _mode_amplitude(mode: DisturbanceModeConfig, t: float) -> float:
    """Evaluate deterministic multi-sine disturbance amplitude at time ``t``."""

    total = 0.0
    for amp, freq, phase in zip(mode.amplitudes, mode.frequencies, mode.phases, strict=True):
        total += amp * np.sin(2.0 * np.pi * freq * t + phase)
    return float(total)


def _mode_phase(mode: DisturbanceModeConfig, t: float) -> float:
    """Evaluate deterministic carrier phase for one disturbance mode."""

    return float(
        mode.phase_offset
        + mode.phase_rate * t
        + mode.phase_mod_amp * np.sin(2.0 * np.pi * mode.phase_mod_freq * t)
    )


def _lowpass_phi_modes(f: np.ndarray, m_cut: int) -> np.ndarray:
    """Retain only low azimuthal Fourier modes up to ``m_cut``.

    Args:
        f: Real field with shape ``(N_r, N_phi)``.
        m_cut: Highest retained nonnegative mode index.

    Returns:
        Field reconstructed from the retained low modes.
    """

    if m_cut <= 0:
        return np.array(f, copy=True)
    nphi = f.shape[1]
    keep = int(min(max(m_cut, 0), nphi // 2))
    hat = np.fft.rfft(f, axis=1)
    if keep + 1 < hat.shape[1]:
        hat[:, keep + 1 :] = 0.0
    return np.fft.irfft(hat, n=nphi, axis=1)


class ToyERDOperators:
    """Numerical operator bundle for one-step updates of ``(n, omega, psi)``."""

    def __init__(self, cfg: RunConfig, bundle: MeshBundle):
        """Precompute static profiles and Poisson operators.

        Args:
            cfg: Plant run configuration.
            bundle: Mesh and geometry metadata.

        Returns:
            None.
        """

        self.cfg = cfg
        self.bundle = bundle
        self.nr, self.nphi = bundle.shape

        self.n_eq_r = ring_equilibrium(bundle.r, cfg)
        self.n_eq = np.repeat(self.n_eq_r[:, None], self.nphi, axis=1)
        self.kappa_w_r = wall_loss_profile(bundle.r, cfg)
        self.kappa_w = np.repeat(self.kappa_w_r[:, None], self.nphi, axis=1)
        self.kappa_r = _gauss(bundle.r, cfg.ring.r_star, cfg.pde.sigma_kappa) * cfg.pde.kappa_0
        self.kappa = np.repeat(self.kappa_r[:, None], self.nphi, axis=1)
        self.S_r = _gauss(bundle.r, cfg.ring.r_star, cfg.pde.sigma_S) * cfg.pde.S_0
        self.S = np.repeat(self.S_r[:, None], self.nphi, axis=1)
        self._g_noise = _gauss(bundle.r, cfg.ring.r_star, 0.9 * cfg.ring.sigma_star)[:, None]
        self._area = bundle.area_weights
        self._mass_target = float(np.sum(self.n_eq * self._area))
        self._source_int = float(np.sum(self.S * self._area))
        self._n_cap = max(4.0, 4.0 * float(cfg.ring.n_bg + cfg.ring.n_amp))
        self._omega_cap = 80.0
        proj_profile = 0.7 * self.S + 0.3 * np.maximum(self.n_eq - cfg.ring.n_bg, 0.0)
        proj_norm = float(np.sum(proj_profile * self._area))
        self._mass_proj_profile = proj_profile / max(proj_norm, 1.0e-12)

        self.basis = self._build_basis()
        self._poisson_inv = self._build_poisson_inverse()
        self._disturbance_jitter = self._build_disturbance_jitter()

    def _build_basis(self) -> BasisFields:
        """Build and cache control/drive basis fields ``b0..b4``."""

        r = self.bundle.r
        phi = self.bundle.phi
        fcfg = self.cfg.forcing
        rcfg = self.cfg.ring

        g0 = _gauss(r, rcfg.r_star, fcfg.sigma_0)
        g1 = _gauss(r, rcfg.r_star, fcfg.sigma_1)
        g2 = _gauss(r, rcfg.r_star, fcfg.sigma_2)

        cos1 = np.cos(phi)
        sin1 = np.sin(phi)
        cos2 = np.cos(2.0 * phi)
        sin2 = np.sin(2.0 * phi)

        b0 = np.repeat(g0[:, None], self.nphi, axis=1)
        b1 = g1[:, None] * cos1[None, :]
        b2 = g1[:, None] * sin1[None, :]
        b3 = g2[:, None] * cos2[None, :]
        b4 = g2[:, None] * sin2[None, :]

        return BasisFields(b0=b0, b1=b1, b2=b2, b3=b3, b4=b4, g1=g1[:, None], g2=g2[:, None])

    def _build_disturbance_jitter(self) -> dict[str, np.ndarray | float]:
        """Build deterministic per-trajectory jitter for disturbance schedules.

        Jitter is derived from ``cfg.init.seed`` so trajectories remain
        reproducible while still differing across run seeds.

        Args:
            None.

        Returns:
            Dictionary of amplitude/frequency/phase perturbations for both
            disturbance modes.
        """

        rng = np.random.default_rng(self.cfg.init.seed + 7919)
        frac = float(max(0.0, self.cfg.disturbance.mode_jitter_amp_frac))
        ff = float(max(0.0, self.cfg.disturbance.mode_jitter_freq_frac))
        ph = float(max(0.0, self.cfg.disturbance.mode_jitter_phase))

        m1_len = len(self.cfg.disturbance.mode1.amplitudes)
        m2_len = len(self.cfg.disturbance.mode2.amplitudes)
        ms_len = len(self.cfg.disturbance.multiscale_modes)

        m_min = int(self.cfg.disturbance.noise_m_min)
        m_max = int(self.cfg.disturbance.noise_m_max)
        if m_max < m_min:
            m_min, m_max = m_max, m_min
        noise_modes = np.arange(max(1, m_min), max(1, m_max) + 1, dtype=int)
        noise_len = int(noise_modes.size)

        noise_refresh_dt = float(max(self.cfg.time.dt, self.cfg.disturbance.noise_refresh_dt))
        noise_t = np.arange(0.0, self.cfg.time.T_final + noise_refresh_dt, noise_refresh_dt, dtype=float)
        if noise_t.size < 2:
            noise_t = np.array([0.0, self.cfg.time.T_final], dtype=float)
        corr_time = float(max(noise_refresh_dt, self.cfg.disturbance.noise_corr_time))
        rho = float(np.exp(-noise_refresh_dt / corr_time))

        noise_a = np.zeros((noise_len, noise_t.size), dtype=float)
        noise_b = np.zeros((noise_len, noise_t.size), dtype=float)
        if noise_len > 0:
            noise_a[:, 0] = rng.normal(0.0, 1.0, size=noise_len)
            noise_b[:, 0] = rng.normal(0.0, 1.0, size=noise_len)
            sigma = np.sqrt(max(1.0 - rho**2, 1.0e-12))
            for k in range(1, noise_t.size):
                noise_a[:, k] = rho * noise_a[:, k - 1] + sigma * rng.normal(0.0, 1.0, size=noise_len)
                noise_b[:, k] = rho * noise_b[:, k - 1] + sigma * rng.normal(0.0, 1.0, size=noise_len)
            noise_a -= noise_a.mean(axis=1, keepdims=True)
            noise_b -= noise_b.mean(axis=1, keepdims=True)
            noise_a /= np.maximum(noise_a.std(axis=1, keepdims=True), 1.0e-12)
            noise_b /= np.maximum(noise_b.std(axis=1, keepdims=True), 1.0e-12)

        noise_mode_decay = float(max(0.0, self.cfg.disturbance.noise_mode_decay))
        if noise_len > 0:
            ref_mode = float(max(1, int(noise_modes.min())))
            noise_weight = (ref_mode / noise_modes.astype(float)) ** noise_mode_decay
            noise_weight = noise_weight / np.maximum(np.linalg.norm(noise_weight), 1.0e-12)
        else:
            noise_weight = np.zeros(0, dtype=float)

        phase_strength = float(max(0.0, self.cfg.disturbance.phase_drift_strength))
        phase_corr = float(max(noise_refresh_dt, self.cfg.disturbance.phase_drift_corr_time))
        rho_phase = float(np.exp(-noise_refresh_dt / phase_corr))
        sigma_phase = np.sqrt(max(1.0 - rho_phase**2, 1.0e-12))

        carrier_v1 = np.zeros(noise_t.size, dtype=float)
        carrier_v2 = np.zeros(noise_t.size, dtype=float)
        if phase_strength > 0.0:
            for k in range(1, noise_t.size):
                carrier_v1[k] = rho_phase * carrier_v1[k - 1] + sigma_phase * rng.normal(0.0, phase_strength)
                carrier_v2[k] = rho_phase * carrier_v2[k - 1] + sigma_phase * rng.normal(0.0, phase_strength)
        carrier_phase1 = np.cumsum(carrier_v1) * noise_refresh_dt
        carrier_phase2 = np.cumsum(carrier_v2) * noise_refresh_dt

        ms_phase_drift = np.zeros((ms_len, noise_t.size), dtype=float)
        if phase_strength > 0.0 and ms_len > 0:
            for m_idx in range(ms_len):
                vel = np.zeros(noise_t.size, dtype=float)
                local_strength = phase_strength / np.sqrt(1.0 + m_idx)
                for k in range(1, noise_t.size):
                    vel[k] = rho_phase * vel[k - 1] + sigma_phase * rng.normal(0.0, local_strength)
                ms_phase_drift[m_idx] = np.cumsum(vel) * noise_refresh_dt

        eddy_count = int(max(0, self.cfg.disturbance.eddy_count))
        eddy_phi = np.zeros((eddy_count, noise_t.size), dtype=float)
        eddy_r = np.zeros((eddy_count, noise_t.size), dtype=float)
        eddy_amp = np.zeros((eddy_count, noise_t.size), dtype=float)
        eddy_theta = np.zeros((eddy_count, noise_t.size), dtype=float)
        eddy_tilt = np.zeros(eddy_count, dtype=float)
        if eddy_count > 0 and float(self.cfg.disturbance.eddy_strength) > 0.0:
            eddy_corr = float(max(noise_refresh_dt, self.cfg.disturbance.eddy_corr_time))
            rho_eddy = float(np.exp(-noise_refresh_dt / eddy_corr))
            sigma_eddy = np.sqrt(max(1.0 - rho_eddy**2, 1.0e-12))
            speed_std = float(max(1.0e-8, self.cfg.disturbance.eddy_speed_std))
            radial_jitter = float(max(1.0e-8, self.cfg.disturbance.eddy_radial_jitter))

            phi_pos = rng.uniform(0.0, 2.0 * np.pi, size=eddy_count)
            phi_vel = rng.normal(0.0, speed_std, size=eddy_count)
            r_off = rng.normal(0.0, 0.4 * radial_jitter, size=eddy_count)
            amp_state = rng.normal(0.0, 1.0, size=eddy_count)
            theta = rng.uniform(0.0, 2.0 * np.pi, size=eddy_count)
            theta_vel = rng.normal(0.0, 0.75 * speed_std, size=eddy_count)
            eddy_tilt = rng.choice([-1.0, 1.0], size=eddy_count)
            phi_vel -= np.mean(phi_vel)
            r_off -= np.mean(r_off)
            amp_state -= np.mean(amp_state)
            theta_vel -= np.mean(theta_vel)

            eddy_phi[:, 0] = phi_pos
            eddy_r[:, 0] = np.clip(self.cfg.ring.r_star + r_off, self.bundle.r[0], self.bundle.r[-1])
            eddy_amp[:, 0] = amp_state
            eddy_theta[:, 0] = theta

            for k in range(1, noise_t.size):
                phi_vel = rho_eddy * phi_vel + sigma_eddy * rng.normal(0.0, speed_std, size=eddy_count)
                amp_state = rho_eddy * amp_state + sigma_eddy * rng.normal(0.0, 1.0, size=eddy_count)
                r_off = rho_eddy * r_off + sigma_eddy * rng.normal(0.0, radial_jitter, size=eddy_count)
                theta_vel = rho_eddy * theta_vel + sigma_eddy * rng.normal(0.0, 0.75 * speed_std, size=eddy_count)
                phi_vel -= np.mean(phi_vel)
                r_off -= np.mean(r_off)
                amp_state -= np.mean(amp_state)
                theta_vel -= np.mean(theta_vel)
                phi_pos = (phi_pos + phi_vel * noise_refresh_dt) % (2.0 * np.pi)
                theta = (theta + theta_vel * noise_refresh_dt) % (2.0 * np.pi)

                eddy_phi[:, k] = phi_pos
                eddy_r[:, k] = np.clip(self.cfg.ring.r_star + r_off, self.bundle.r[0], self.bundle.r[-1])
                eddy_amp[:, k] = amp_state
                eddy_theta[:, k] = theta
        else:
            eddy_theta = np.zeros((0, noise_t.size), dtype=float)
            eddy_tilt = np.zeros(0, dtype=float)

        return {
            "m1_amp": 1.0 + frac * rng.uniform(-1.0, 1.0, size=m1_len),
            "m2_amp": 1.0 + frac * rng.uniform(-1.0, 1.0, size=m2_len),
            "m1_freq": 1.0 + ff * rng.uniform(-1.0, 1.0, size=m1_len),
            "m2_freq": 1.0 + ff * rng.uniform(-1.0, 1.0, size=m2_len),
            "m1_phase": ph * rng.uniform(-1.0, 1.0, size=m1_len),
            "m2_phase": ph * rng.uniform(-1.0, 1.0, size=m2_len),
            "carrier_m1": ph * float(rng.uniform(-1.0, 1.0)),
            "carrier_m2": ph * float(rng.uniform(-1.0, 1.0)),
            "ms_phase": 2.0 * np.pi * rng.uniform(0.0, 1.0, size=ms_len),
            "ms_amp": 1.0 + frac * rng.uniform(-1.0, 1.0, size=ms_len),
            "noise_modes": noise_modes,
            "noise_t": noise_t,
            "noise_a": noise_a,
            "noise_b": noise_b,
            "noise_weight": noise_weight,
            "carrier_phase1": carrier_phase1,
            "carrier_phase2": carrier_phase2,
            "ms_phase_drift": ms_phase_drift,
            "eddy_phi": eddy_phi,
            "eddy_r": eddy_r,
            "eddy_amp": eddy_amp,
            "eddy_theta": eddy_theta,
            "eddy_tilt": eddy_tilt,
        }

    def _interp_noise_state(self, t: float) -> tuple[np.ndarray, np.ndarray]:
        """Interpolate seeded stochastic disturbance coefficients at time ``t``.

        Args:
            t: Simulation time.

        Returns:
            Tuple of cosine and sine modal coefficients for the configured
            stochastic disturbance band.
        """

        noise_t = np.asarray(self._disturbance_jitter["noise_t"], dtype=float)
        noise_a = np.asarray(self._disturbance_jitter["noise_a"], dtype=float)
        noise_b = np.asarray(self._disturbance_jitter["noise_b"], dtype=float)
        if noise_a.size == 0 or noise_b.size == 0:
            return np.zeros(0, dtype=float), np.zeros(0, dtype=float)
        coeff_a = np.array([np.interp(t, noise_t, row) for row in noise_a], dtype=float)
        coeff_b = np.array([np.interp(t, noise_t, row) for row in noise_b], dtype=float)
        return coeff_a, coeff_b

    def _interp_rows(self, t: float, rows: np.ndarray) -> np.ndarray:
        """Interpolate a 2D ``(n_series, n_time)`` array at time ``t``."""

        noise_t = np.asarray(self._disturbance_jitter["noise_t"], dtype=float)
        if rows.size == 0:
            return np.zeros(0, dtype=float)
        return np.array([np.interp(t, noise_t, row) for row in rows], dtype=float)

    def _mass_feedback_source(self, n_eval: np.ndarray) -> np.ndarray:
        """Compute a dynamic fueling field that keeps ring mass near target.

        The fueling profile retains the configured radial shape ``S(r)``, but
        its amplitude is adjusted each stage so wall loss is replenished and the
        total density inventory does not drift unrealistically over the run.
        A fixed baseline fueling floor is always present so the plant can
        actually build and hold a near-threshold profile instead of only
        servoing to the instantaneous wall sink.
        """

        p = self.cfg.pde
        wall_loss = float(np.sum(self.kappa_w * n_eval * self._area))
        mass_now = float(np.sum(n_eval * self._area))
        mass_err = (self._mass_target - mass_now) / max(self._mass_target, 1.0e-12)
        source_scale = (
            p.source_floor_scale
            + p.source_balance_gain * wall_loss / max(self._source_int, 1.0e-12)
            + p.source_mass_gain * np.clip(mass_err, -0.12, 0.12)
        )
        source_scale = float(np.clip(source_scale, 0.0, p.source_scale_max))
        return source_scale * self.S

    def _phase_advection_speed(
        self,
        grad_signed: np.ndarray,
        activation: np.ndarray,
    ) -> np.ndarray:
        """Compute a supercritical diamagnetic phase-advection speed.

        Args:
            grad_signed: Signed zonal-gradient ratio field.
            activation: Smooth supercritical activation field.

        Returns:
            Signed azimuthal phase-speed field that accelerates unstable
            nonzonal structures once the ring becomes supercritical.
        """

        gain = float(getattr(self.cfg.pde, "phase_advection_gain", 0.0))
        if gain == 0.0:
            return np.zeros_like(activation)
        return gain * activation * grad_signed

    def _wave_packet_speeds(
        self,
        n_tilde: np.ndarray,
        omega_tilde: np.ndarray,
        activation: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Build local packet speeds for the ``wave_landau`` branch.

        The goal is to make propagation emerge from the evolving state rather
        than from an externally scripted carrier:
        - azimuthal packet direction follows local density/vorticity structure,
        - radial packet speed follows cross-gradient structure,
        - both speeds increase when the local instability turns on.

        Args:
            n_tilde: Nonzonal density component at the current RK stage.
            omega_tilde: Nonzonal vorticity component at the current RK stage.
            activation: Smooth threshold activation field.

        Returns:
            Tuple ``(c_phi, c_r, burst_gain)`` of local azimuthal speed,
            radial speed, and multiplicative burst amplification.
        """

        p = self.cfg.pde
        phi_gain = float(getattr(p, "wave_packet_phi_gain", 0.0))
        r_gain = float(getattr(p, "wave_packet_r_gain", 0.0))
        burst_gain = float(getattr(p, "wave_burst_speed_gain", 0.0))
        if phi_gain == 0.0 and r_gain == 0.0 and burst_gain == 0.0:
            zeros = np.zeros_like(activation)
            ones = np.ones_like(activation)
            return zeros, zeros, ones

        n_ref = max(self.cfg.ring.n_amp, 1.0e-8)
        grad_ref = max(self.cfg.ring.n_amp / max(self.cfg.ring.sigma_star, 1.0e-12), 1.0e-8)
        omega_ref = max(0.15 * self._omega_cap, 1.0)
        omega_grad_ref = max(omega_ref / max(self.cfg.ring.sigma_star, 1.0e-12), 1.0e-8)

        dn_dphi = _diff_phi_periodic(n_tilde, self.bundle.dphi) / n_ref
        dn_dr = _diff_r_centered(n_tilde, self.bundle.dr) / grad_ref
        domega_dphi = _diff_phi_periodic(omega_tilde, self.bundle.dphi) / omega_ref
        domega_dr = _diff_r_centered(omega_tilde, self.bundle.dr) / omega_grad_ref

        structure_phi = 0.55 * (omega_tilde / omega_ref) + 0.30 * dn_dphi + 0.15 * domega_dphi
        structure_r = 0.55 * dn_dr - 0.25 * domega_dphi + 0.20 * domega_dr
        local_burst = 1.0 + burst_gain * activation + 0.35 * np.tanh(np.abs(omega_tilde) / omega_ref)
        activation_gate = activation / (0.20 + activation)
        phi_gate = 0.35 + 0.65 * activation_gate
        radial_gate = activation_gate

        c_phi = phi_gain * phi_gate * local_burst * np.tanh(structure_phi)
        c_r = r_gain * radial_gate * local_burst * np.tanh(structure_r)
        return c_phi, c_r, local_burst

    def _radial_confinement_flux(
        self,
        n_eval: np.ndarray,
        activation: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Build a state-aware inward density confinement flux.

        The ring should be coherent and fast in the marginal regime rather than
        instantly exploding to the radial boundaries. This term acts like a
        weak inward pinch around ``r_star`` that relaxes as the local
        instability activation grows.

        Args:
            n_eval: Density field at the current RK stage.
            activation: Smooth instability activation field.

        Returns:
            Tuple ``(flux_divergence, v_conf)`` where ``flux_divergence`` is
            added to the density equation and ``v_conf`` is the local
            confinement velocity used for diagnostics.
        """

        gain = float(getattr(self.cfg.pde, "radial_confinement_gain", 0.0))
        if gain == 0.0:
            zeros = np.zeros_like(n_eval)
            return zeros, zeros

        release_gain = float(getattr(self.cfg.pde, "radial_release_gain", 0.0))
        damage_release_gain = float(getattr(self.cfg.pde, "radial_damage_release_gain", 0.0))
        r_off = (self.bundle.r - self.cfg.ring.r_star) / max(self.cfg.ring.sigma_star, 1.0e-12)
        envelope = np.exp(-0.5 * (r_off / 1.6) ** 2)
        shape = (-r_off * envelope)[:, None]
        degradation = np.abs(n_eval - self.n_eq) / max(self.cfg.ring.n_amp, 1.0e-8)
        degradation = np.clip(degradation, 0.0, 3.0)
        release = 1.0 / (1.0 + release_gain * activation + damage_release_gain * degradation)
        v_conf = gain * release * shape
        v_conf[0, :] = 0.0
        v_conf[-1, :] = 0.0

        flux = -_diff_r_centered(v_conf * n_eval, self.bundle.dr)
        flux -= float(np.sum(flux * self._area) / np.sum(self._area))
        return flux, v_conf

    def _zonal_profile_restore(
        self,
        n_eval: np.ndarray,
        activation: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Build an axisymmetric operating-point support term for ``n``.

        The active wave branch should be marginally stable in the sense that
        the mean ring persists and propagating nonzonal packets live on top of
        it. A local inward pinch on the full 2D density field was too blunt:
        it could suppress the very wave structures we want and bias the run
        toward inward collapse. This term instead acts only on the zonal mean,
        nudging it toward ``n_eq(r)`` while releasing smoothly once the local
        instability or accumulated zonal damage becomes significant.

        Args:
            n_eval: Density field at the current RK stage.
            activation: Smooth instability activation field.

        Returns:
            Tuple ``(support, release)`` where ``support`` is an axisymmetric
            density source term and ``release`` is the zonal support-release
            factor used for diagnostics.
        """

        gain = float(getattr(self.cfg.pde, "zonal_profile_restore_gain", 0.0))
        if gain == 0.0:
            zeros = np.zeros_like(n_eval)
            return zeros, zeros

        zonal_n = self._zonal_mean(n_eval)
        zonal_activation = self._zonal_mean(activation)
        release_gain = float(getattr(self.cfg.pde, "zonal_profile_release_gain", 0.0))
        damage_release_gain = float(getattr(self.cfg.pde, "zonal_profile_damage_release_gain", 0.0))

        zonal_damage = np.abs(zonal_n - self.n_eq_r[:, None]) / max(self.cfg.ring.n_amp, 1.0e-8)
        zonal_damage = np.clip(zonal_damage, 0.0, 3.0)
        release = 1.0 / (1.0 + release_gain * zonal_activation + damage_release_gain * zonal_damage)

        support = gain * release * (self.n_eq_r[:, None] - zonal_n)
        support = np.repeat(support, self.nphi, axis=1)
        release_2d = np.repeat(release, self.nphi, axis=1)
        return support, release_2d

    def _collapse_guard_terms(
        self,
        n_eval: np.ndarray,
        omega_eval: np.ndarray,
        omega_packet: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, float]:
        """Build targeted safeguards against late annular collapse.

        The good historical ``wave_landau`` runs fail in a very specific way:
        density eventually spills across the full radial strip, touches the
        walls, and then the mass-control machinery plus strong zonal support
        rebuild an unrealistically smooth axisymmetric ring. The remedy should
        therefore be equally specific. This helper:

        1. Detects when edge occupancy is growing or ring-core occupancy is
           collapsing.
        2. Reconstructs a mass-conserving annular target that preserves the
           current azimuthal column masses, so blob structure in ``phi`` is not
           erased.
        3. Applies edge-localized high-mode vorticity damping only during that
           collapse regime.

        Args:
            n_eval: Density field at the current RK stage.
            omega_eval: Vorticity field at the current RK stage.
            omega_packet: Low-/packet-mode vorticity field used to separate
                edge-localized high-mode blow-up from the coherent packet
                dynamics that should remain active.

        Returns:
            Tuple ``(n_guard, omega_guard, guard_activation)`` where:
            - ``n_guard`` is a mass-conserving annular support term for
              density,
            - ``omega_guard`` is targeted vorticity damping,
            - ``guard_activation`` is a scalar collapse-guard activation level.
        """

        p = self.cfg.pde
        support_gain = float(getattr(p, "collapse_guard_gain", 0.0))
        omega_damp_gain = float(getattr(p, "collapse_omega_damp", 0.0))
        if support_gain == 0.0 and omega_damp_gain == 0.0:
            zeros = np.zeros_like(n_eval)
            return zeros, zeros, 0.0

        total_mass = float(np.sum(np.maximum(n_eval, 0.0) * self._area))
        if not np.isfinite(total_mass) or total_mass <= 1.0e-12:
            zeros = np.zeros_like(n_eval)
            return zeros, zeros, 0.0

        r = self.bundle.r
        dr = self.bundle.dr
        span = max(float(self.cfg.domain.R_max - self.cfg.domain.R_min), 1.0e-12)
        edge_width = max(2.0 * dr, 0.10 * span, 0.75 * self.cfg.ring.sigma_star)
        core_width = max(1.35 * self.cfg.ring.sigma_star, edge_width)

        inner_mask = r <= (r[0] + edge_width)
        outer_mask = r >= (r[-1] - edge_width)
        core_mask = np.abs(r - self.cfg.ring.r_star) <= core_width

        edge_mask = np.logical_or(inner_mask, outer_mask)
        edge_frac = float(np.sum(np.maximum(n_eval, 0.0) * self._area * edge_mask[:, None]) / total_mass)
        core_frac = float(np.sum(np.maximum(n_eval, 0.0) * self._area * core_mask[:, None]) / total_mass)

        edge_thresh = float(getattr(p, "collapse_edge_threshold", 0.18))
        core_thresh = float(getattr(p, "collapse_core_threshold", 0.42))
        edge_excess = float(_smooth_positive(edge_frac - edge_thresh, 0.01))
        core_deficit = float(_smooth_positive(core_thresh - core_frac, 0.01))
        raw_activation = edge_excess + 0.85 * core_deficit
        guard_activation = float(raw_activation / (0.03 + raw_activation)) if raw_activation > 0.0 else 0.0
        if guard_activation <= 0.0:
            zeros = np.zeros_like(n_eval)
            return zeros, zeros, 0.0

        line_weight = r * dr
        column_line_mass = np.sum(np.maximum(n_eval, 0.0) * line_weight[:, None], axis=0)
        guard_sigma = max(1.15 * self.cfg.ring.sigma_star, 1.6 * dr)
        annulus_profile = np.exp(-((r - self.cfg.ring.r_star) ** 2) / (2.0 * guard_sigma**2))
        annulus_profile /= max(float(np.sum(annulus_profile * line_weight)), 1.0e-12)
        annular_target = annulus_profile[:, None] * column_line_mass[None, :]
        density_guard = support_gain * guard_activation * (annular_target - n_eval)

        edge_env_inner = np.exp(-((r - (r[0] + 0.55 * edge_width)) ** 2) / (2.0 * edge_width**2))
        edge_env_outer = np.exp(-((r - (r[-1] - 0.55 * edge_width)) ** 2) / (2.0 * edge_width**2))
        edge_env = np.clip(edge_env_inner + edge_env_outer, 0.0, 1.0)[:, None]
        omega_high = omega_eval - omega_packet
        omega_guard = omega_damp_gain * guard_activation * edge_env * (0.35 * omega_eval + 0.65 * omega_high)
        return density_guard, omega_guard, guard_activation

    def _project_density_mass(self, n_field: np.ndarray) -> np.ndarray:
        """Relax ``n`` toward the target inventory while preserving ``n_tilde``.

        The earlier pointwise projection logic could erase the nonzonal density
        structure exactly when the active wave branch became interesting. The
        current implementation instead:

        1. Repairs positivity/caps by rescaling the nonzonal component
           ``n_tilde`` against the zonal mean rather than clipping cells
           independently.
        2. Applies any mass correction to the zonal part only, leaving
           ``n_tilde`` unchanged.

        This is the right numerical behavior for the intended digital twin:
        mass control should keep the annular operating point in range without
        deleting the coherent blobs that the ROM and controller are meant to
        predict and regulate.

        Args:
            n_field: Density field produced by one explicit RK stage update.

        Returns:
            Density field projected back onto the configured target mass.
        """

        projected = self._repair_density_shape(n_field)
        mass_now = float(np.sum(projected * self._area))
        if not np.isfinite(mass_now) or mass_now <= 1.0e-12:
            return projected

        delta_mass = self._mass_target - mass_now
        rel_err = abs(delta_mass) / max(self._mass_target, 1.0e-12)
        if rel_err <= 0.01:
            return projected

        zonal = self._zonal_mean(projected)
        tilde = projected - zonal
        zonal_profile = np.maximum(zonal[:, 0], 1.0e-8)
        zonal_profile = 0.55 * zonal_profile + 0.45 * np.maximum(self.n_eq_r, 1.0e-8)
        zonal_profile_2d = np.repeat(zonal_profile[:, None], self.nphi, axis=1)
        zonal_profile_2d /= max(float(np.sum(zonal_profile_2d * self._area)), 1.0e-12)
        relax = float(np.clip(1.8 * rel_err, 0.08, 0.35))
        zonal += relax * delta_mass * zonal_profile_2d[:, :1]
        projected = zonal + tilde
        return self._repair_density_shape(projected)

    def _repair_density_shape(self, n_field: np.ndarray) -> np.ndarray:
        """Enforce density bounds by rescaling ``n_tilde`` instead of clipping.

        Args:
            n_field: Candidate density field.

        Returns:
            Density field that satisfies the configured lower/upper bounds while
            preserving as much of the nonzonal structure as possible.
        """

        floor = 1.0e-8
        cap = float(self._n_cap)
        repaired = np.array(n_field, copy=True)
        _enforce_neumann_radial(repaired)

        zonal = np.maximum(self._zonal_mean(repaired), 5.0 * floor)
        tilde = repaired - zonal
        zonal_2d = np.repeat(zonal, self.nphi, axis=1)

        alpha = 1.0
        neg_mask = tilde < 0.0
        if np.any(neg_mask):
            alpha_low = np.min((zonal_2d[neg_mask] - floor) / np.maximum(-tilde[neg_mask], 1.0e-12))
            alpha = min(alpha, float(alpha_low))
        pos_mask = tilde > 0.0
        if np.any(pos_mask):
            alpha_high = np.min((cap - zonal_2d[pos_mask]) / np.maximum(tilde[pos_mask], 1.0e-12))
            alpha = min(alpha, float(alpha_high))
        alpha = float(np.clip(alpha, 0.0, 1.0))

        repaired = zonal_2d + alpha * tilde
        repaired = np.clip(repaired, floor, cap)
        _enforce_neumann_radial(repaired)
        return repaired

    def _limit_density_rhs(self, n_eval: np.ndarray, dn_dt: np.ndarray) -> np.ndarray:
        """Bound density time derivatives so one step cannot destroy structure.

        The active branch should be able to generate strong post-critical
        density transport, but an explicit RK stage that tries to remove nearly
        all local density in one shot causes the field to hit the positivity
        floor and then be rebuilt as an axisymmetric ring by the mass-control
        machinery. Limiting the *rate* instead of the full update preserves
        nonzonal structure while preventing that shock/reset failure mode.

        Args:
            n_eval: Density field at the current RK stage.
            dn_dt: Proposed density right-hand side.

        Returns:
            Rate-limited density right-hand side.
        """

        dt_ref = max(float(self.cfg.time.dt), 1.0e-12)
        floor = 1.0e-8
        cap = float(self._n_cap)
        protect_ref = np.maximum(0.18 * self.n_eq, 5.0e-3)
        scarcity = np.clip((protect_ref - n_eval) / protect_ref, 0.0, 1.0)
        cap_fill = np.clip((n_eval - 0.82 * cap) / max(0.18 * cap, 1.0e-8), 0.0, 1.0)

        max_drop_frac = 0.32 - 0.20 * scarcity
        max_rise_frac = 0.90 - 0.30 * cap_fill

        lower = -max_drop_frac * np.maximum(n_eval - floor, 0.0) / dt_ref
        upper = max_rise_frac * np.maximum(cap - n_eval, 0.0) / dt_ref
        limited = np.clip(dn_dt, lower, upper)
        return limited

    def _limit_omega_rhs(self, omega_eval: np.ndarray, domega_dt: np.ndarray) -> np.ndarray:
        """Bound vorticity time derivatives to avoid explosive one-step bursts.

        The active branch should still accelerate when it becomes unstable, but
        once the explicit vorticity update starts growing by order-one factors
        in a single RK stage, the density field is forced into the same
        pathological reset path observed in the failing runs. This limiter keeps
        the post-critical regime violent without allowing unbounded stage-wise
        amplification.

        Args:
            omega_eval: Vorticity field at the current RK stage.
            domega_dt: Proposed vorticity right-hand side.

        Returns:
            Rate-limited vorticity right-hand side.
        """

        dt_ref = max(float(self.cfg.time.dt), 1.0e-12)
        local_scale = np.maximum(np.abs(omega_eval), 1.0)
        omega_ref = max(0.15 * self._omega_cap, 1.0)
        blowup = np.clip((np.abs(omega_eval) / omega_ref - 0.9) / 1.2, 0.0, 1.0)
        max_frac = 1.40 - 0.85 * blowup
        limit = max_frac * local_scale / dt_ref
        limited = np.clip(domega_dt, -limit, limit)
        return limited

    def _build_poisson_inverse(self) -> np.ndarray:
        """Precompute per-Fourier-mode radial inverses for ``∇²ψ = -ω``."""

        n_in = self.nr - 2
        invs = np.zeros((self.nphi, n_in, n_in), dtype=np.complex128)

        for m in range(self.nphi):
            lam_phi = -4.0 * np.sin(np.pi * m / self.nphi) ** 2 / (self.bundle.dphi**2)
            a = np.full(n_in - 1, 1.0 / (self.bundle.dr**2), dtype=np.complex128)
            b = np.full(n_in, -2.0 / (self.bundle.dr**2) + lam_phi, dtype=np.complex128)
            c = np.full(n_in - 1, 1.0 / (self.bundle.dr**2), dtype=np.complex128)

            mat = np.zeros((n_in, n_in), dtype=np.complex128)
            np.fill_diagonal(mat, b)
            np.fill_diagonal(mat[1:], a)
            np.fill_diagonal(mat[:, 1:], c)
            invs[m] = np.linalg.inv(mat)
        return invs

    def _disturbance_stage(self, t: float) -> tuple[float, bool]:
        """Return stage amplitude scale and final-stage flag.

        Args:
            t: Simulation time.

        Returns:
            Tuple ``(scale, is_final_stage)`` for the configured three-stage
            disturbance schedule.
        """

        dc = self.cfg.disturbance
        t_final = max(self.cfg.time.T_final, 1e-12)
        frac = float(np.clip(t / t_final, 0.0, 1.0))
        warm = float(np.clip(dc.warmup_fraction, 0.0, 1.0))
        excite = float(np.clip(dc.excite_fraction, 0.0, 1.0 - warm))

        if frac < warm:
            return float(dc.warmup_scale), False
        if frac < warm + excite:
            return float(dc.excite_scale), False
        return float(dc.hold_scale), True

    def solve_poisson(self, omega: np.ndarray) -> np.ndarray:
        """Solve ``∇²ψ = -ω`` with radial Dirichlet boundaries on ``ψ``.

        Args:
            omega: Vorticity field ``omega(r, phi)``.

        Returns:
            Streamfunction field ``psi(r, phi)`` satisfying ``psi=0`` on radial walls.
        """

        omega_hat = np.fft.fft(omega, axis=1)
        rhs_hat = -omega_hat[1:-1, :]

        psi_hat = np.zeros_like(omega_hat, dtype=np.complex128)
        for m in range(self.nphi):
            psi_hat[1:-1, m] = self._poisson_inv[m] @ rhs_hat[:, m]

        psi = np.fft.ifft(psi_hat, axis=1).real
        psi[0, :] = 0.0
        psi[-1, :] = 0.0
        return psi

    def velocity(self, psi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Compute rectangularized velocity from streamfunction.

        Args:
            psi: Streamfunction field.

        Returns:
            Tuple ``(u_r, u_phi) = (-dpsi/dphi, dpsi/dr)``.
        """

        u_r = -_diff_phi_periodic(psi, self.bundle.dphi)
        u_phi = _diff_r_centered(psi, self.bundle.dr)
        return u_r, u_phi

    def _zonal_mean(self, f: np.ndarray) -> np.ndarray:
        """Compute the zonal mean ``<f>(r)`` as a broadcastable column array."""

        return np.mean(f, axis=1, keepdims=True)

    def _tilde(self, f: np.ndarray) -> np.ndarray:
        """Compute the nonzonal component ``f - <f>``."""

        return f - self._zonal_mean(f)

    def _laplacian(self, f: np.ndarray) -> np.ndarray:
        """Apply the rectangularized Laplacian with Neumann/periodic BCs."""

        return _d2_r_neumann(f, self.bundle.dr) + _d2_phi_periodic(f, self.bundle.dphi)

    def _instability_state(
        self,
        n_eval: np.ndarray,
        phi: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Compute zonal-gradient / zonal-shear instability diagnostics.

        The near-threshold branch is driven by the evolving zonal density
        gradient and suppressed by zonal ExB shear. The returned fields are
        dimensionless and broadcast to the full ``(r, phi)`` grid.

        Args:
            n_eval: Density field at the current RK stage.
            phi: Diagnostic potential field at the current RK stage.

        Returns:
            Tuple ``(grad_signed, grad_mag, shear_norm, margin, activation)`` where:
            - ``grad_signed`` is the signed zonal density-gradient normalized by
              ``n_amp / sigma_star``,
            - ``grad_mag`` is the corresponding absolute-gradient magnitude,
            - ``shear_norm`` is the zonal ExB shear normalized by
              ``pde.shear_ref``,
            - ``margin`` is the signed instability margin based on
              ``grad_mag``,
            - ``activation`` is a smooth positive part of ``margin``.
        """

        p = self.cfg.pde
        zonal_n = self._zonal_mean(n_eval)
        grad_ref = max(self.cfg.ring.n_amp / max(self.cfg.ring.sigma_star, 1.0e-12), 1.0e-12)
        grad_signed_1d = -_diff_r_centered(zonal_n, self.bundle.dr) / grad_ref
        grad_signed = np.repeat(grad_signed_1d, self.nphi, axis=1)
        grad_mag = np.abs(grad_signed)

        _, u_phi = self.velocity(phi)
        zonal_u_phi = self._zonal_mean(u_phi)
        shear_ref = max(float(getattr(p, "shear_ref", 1.0)), 1.0e-8)
        shear_norm_1d = np.abs(_diff_r_centered(zonal_u_phi, self.bundle.dr)) / shear_ref
        shear_norm = np.repeat(shear_norm_1d, self.nphi, axis=1)

        crit = float(getattr(p, "critical_gradient_ratio", 0.0))
        shear_gain = float(getattr(p, "shear_suppression_gain", 0.0))
        margin = grad_mag - crit - shear_gain * shear_norm
        activation = _smooth_positive(margin, float(getattr(p, "threshold_width", 0.08)))
        return grad_signed, grad_mag, shear_norm, margin, activation

    def _drift_drive_profile(
        self,
        n_eval: np.ndarray,
        phi: np.ndarray | None = None,
        activation: np.ndarray | None = None,
    ) -> np.ndarray:
        """Build the effective drift-wave drive profile for the density equation.

        The baseline annular mHW term uses the prescribed profile ``kappa(r)``.
        To reduce radially columnar structures, an additional self-consistent
        contribution based on the evolving zonal density gradient can be
        activated. That allows the inner and outer flanks of the ring to drift
        differently as the profile deforms.

        Args:
            n_eval: Density field at the current RK stage.
            phi: Optional potential field used to compute thresholded
                supercritical drive activation.
            activation: Optional precomputed instability activation field.

        Returns:
            Broadcasted ``(N_r, N_phi)`` drift-drive profile.
        """

        p = self.cfg.pde
        drive = np.array(self.kappa, copy=True)
        grad_gain = float(p.gradient_drive_gain)
        threshold_enabled = (
            float(getattr(p, "critical_gradient_ratio", 0.0)) > 0.0
            or float(getattr(p, "shear_suppression_gain", 0.0)) != 0.0
            or float(getattr(p, "shear_damping_gain", 0.0)) != 0.0
        )
        if grad_gain == 0.0 and not threshold_enabled:
            return drive

        zonal_n = self._zonal_mean(n_eval)
        grad_n = -_diff_r_centered(zonal_n, self.bundle.dr)
        grad_ref = max(self.cfg.ring.n_amp / max(self.cfg.ring.sigma_star, 1.0e-12), 1.0e-12)
        if threshold_enabled:
            if activation is None:
                if phi is None:
                    raise ValueError("Thresholded drift drive requires phi or activation")
                grad_signed, _, _, _, activation = self._instability_state(n_eval, phi)
            else:
                grad_signed, _, _, _, _ = self._instability_state(n_eval, phi)
            drive += self.cfg.pde.kappa_0 * grad_gain * activation * grad_signed
        elif grad_gain != 0.0:
            grad_term = grad_gain * (grad_n / grad_ref)
            drive += np.repeat(grad_term, self.nphi, axis=1)
        return np.clip(drive, -2.5 * self.cfg.pde.kappa_0, 2.5 * self.cfg.pde.kappa_0)

    def _hyperdiff(self, f: np.ndarray, nu: float, p: int) -> np.ndarray:
        """Apply dissipative order-``2p`` regularization to ``f``.

        ``p=1`` gives standard Laplacian diffusion and ``p=2`` gives
        dissipative biharmonic regularization.
        """

        lap = self._laplacian(f)
        if p == 1:
            return nu * lap
        if p == 2:
            return -nu * self._laplacian(lap)
        raise ValueError(f"Unsupported hyper_p={p}; expected 1 or 2")

    def _flux_balance_omega_source(self, n_eval: np.ndarray, phi: np.ndarray) -> np.ndarray:
        """Compute a flux-balanced Hasegawa-Wakatani-style vorticity source.

        The source follows the structure of the auxiliary vorticity equation in
        flux-balanced / bHW-style models: radial density gradients couple to
        fluctuating radial transport and can produce bursts when local flux
        channels amplify. This is the most direct way to encode
        ``marginally-stable until transport runs away`` while preserving the
        two-field ``(n, omega)`` contract.

        Args:
            n_eval: Density field at the current RK stage.
            phi: Diagnostic potential at the current RK stage.

        Returns:
            Zero-mean vorticity source field.
        """

        gain = float(getattr(self.cfg.pde, "flux_balance_omega_gain", 0.0))
        if gain == 0.0:
            return np.zeros_like(n_eval)

        n_tilde = self._tilde(n_eval)
        zonal_n = self._zonal_mean(n_eval)
        grad_zonal_n = _diff_r_centered(zonal_n, self.bundle.dr)
        grad_zonal_n_2d = np.repeat(grad_zonal_n, self.nphi, axis=1)

        u_r = -_diff_phi_periodic(phi, self.bundle.dphi)
        u_r_tilde = self._tilde(u_r)
        flux_nonlinear = _diff_r_centered(u_r_tilde * n_tilde, self.bundle.dr)

        source = gain * (grad_zonal_n_2d * u_r_tilde - flux_nonlinear)
        source -= float(np.sum(source * self._area) / np.sum(self._area))
        return source

    def _arakawa_bracket(self, phi: np.ndarray, f: np.ndarray, dr: float, dphi: float) -> np.ndarray:
        """Compute an Arakawa-type Jacobian approximation of ``{phi, f}``."""

        phi_pad = np.pad(phi, ((1, 1), (0, 0)), mode="edge")
        f_pad = np.pad(f, ((1, 1), (0, 0)), mode="edge")

        pc = phi_pad[1:-1]
        pp = phi_pad[2:]
        pm = phi_pad[:-2]
        fc = f_pad[1:-1]
        fp = f_pad[2:]
        fm = f_pad[:-2]

        pc_p = np.roll(pc, -1, axis=1)
        pc_m = np.roll(pc, 1, axis=1)
        pp_p = np.roll(pp, -1, axis=1)
        pp_m = np.roll(pp, 1, axis=1)
        pm_p = np.roll(pm, -1, axis=1)
        pm_m = np.roll(pm, 1, axis=1)

        fc_p = np.roll(fc, -1, axis=1)
        fc_m = np.roll(fc, 1, axis=1)
        fp_p = np.roll(fp, -1, axis=1)
        fp_m = np.roll(fp, 1, axis=1)
        fm_p = np.roll(fm, -1, axis=1)
        fm_m = np.roll(fm, 1, axis=1)

        coeff = 1.0 / (12.0 * dr * dphi)
        j1 = ((pp - pm) * (fc_p - fc_m) - (pc_p - pc_m) * (fp - fm)) * (3.0 * coeff)
        j2 = (
            pp * (fp_p - fp_m)
            - pm * (fm_p - fm_m)
            - pc_p * (pp_p - pm_p)
            + pc_m * (pp_m - pm_m)
        ) * coeff
        j3 = (
            pp_p * (fc_p - fp)
            - pm_m * (fm - fc_m)
            - pm_p * (fc_p - fm)
            + pp_m * (fp - fc_m)
        ) * coeff
        return (j1 + j2 + j3) / 3.0

    def disturbance(self, t: float) -> np.ndarray:
        """Evaluate deterministic disturbance forcing field at time ``t``.

        Args:
            t: Simulation time.

        Returns:
            Disturbance contribution to vorticity forcing.
        """

        phi = self.bundle.phi[None, :]

        m1 = self.cfg.disturbance.mode1
        m2 = self.cfg.disturbance.mode2

        scale, is_final_stage = self._disturbance_stage(t)

        d1 = 0.0
        for amp, freq, phase, a_j, f_j, p_j in zip(
            m1.amplitudes,
            m1.frequencies,
            m1.phases,
            self._disturbance_jitter["m1_amp"],
            self._disturbance_jitter["m1_freq"],
            self._disturbance_jitter["m1_phase"],
            strict=True,
        ):
            d1 += (amp * a_j) * np.sin(2.0 * np.pi * (freq * f_j) * t + phase + p_j)

        d2 = 0.0
        for amp, freq, phase, a_j, f_j, p_j in zip(
            m2.amplitudes,
            m2.frequencies,
            m2.phases,
            self._disturbance_jitter["m2_amp"],
            self._disturbance_jitter["m2_freq"],
            self._disturbance_jitter["m2_phase"],
            strict=True,
        ):
            d2 += (amp * a_j) * np.sin(2.0 * np.pi * (freq * f_j) * t + phase + p_j)

        p1 = (
            _mode_phase(m1, t)
            + float(self._disturbance_jitter["carrier_m1"])
            + float(np.interp(t, self._disturbance_jitter["noise_t"], self._disturbance_jitter["carrier_phase1"]))
        )
        p2 = (
            _mode_phase(m2, t)
            + float(self._disturbance_jitter["carrier_m2"])
            + float(np.interp(t, self._disturbance_jitter["noise_t"], self._disturbance_jitter["carrier_phase2"]))
        )

        if is_final_stage:
            boost = float(self.cfg.disturbance.phase_c_boost_freq)
            p1 += 0.18 * np.sin(2.0 * np.pi * boost * t)
            p2 += 0.16 * np.sin(2.0 * np.pi * (1.2 * boost) * t)

        term1 = scale * d1 * self.basis.g1 * np.cos(phi - p1)
        term2 = scale * d2 * self.basis.g2 * np.cos(2.0 * phi - p2)
        dc = self.cfg.disturbance

        term_ms = np.zeros_like(term1)
        if dc.multiscale_modes:
            n_ms = min(
                len(dc.multiscale_modes),
                len(dc.multiscale_amplitudes),
                len(dc.multiscale_frequencies),
                len(self._disturbance_jitter["ms_phase"]),
                len(self._disturbance_jitter["ms_amp"]),
            )
            phase_drift = self._interp_rows(t, np.asarray(self._disturbance_jitter["ms_phase_drift"], dtype=float))
            for idx, (m, amp, freq, phase0, amp_j) in enumerate(zip(
                dc.multiscale_modes[:n_ms],
                dc.multiscale_amplitudes[:n_ms],
                dc.multiscale_frequencies[:n_ms],
                self._disturbance_jitter["ms_phase"][:n_ms],
                self._disturbance_jitter["ms_amp"][:n_ms],
            )):
                phase_offset = float(phase_drift[idx]) if idx < phase_drift.size else 0.0
                phase_t = (
                    2.0 * np.pi * freq * t
                    + phase0
                    + phase_offset
                )
                if is_final_stage:
                    phase_t += 0.10 * np.sin(2.0 * np.pi * dc.phase_c_boost_freq * t)
                term_ms += (amp * amp_j) * self.basis.g2 * np.cos(float(m) * phi - phase_t)
            term_ms *= float(scale * dc.multiscale_strength)

        term_noise = np.zeros_like(term1)
        noise_strength = float(max(0.0, dc.noise_strength))
        if noise_strength > 0.0 and self._disturbance_jitter["noise_modes"].size > 0:
            noise_a_t, noise_b_t = self._interp_noise_state(t)
            for m, a_t, b_t, w_m in zip(
                self._disturbance_jitter["noise_modes"],
                noise_a_t,
                noise_b_t,
                self._disturbance_jitter["noise_weight"],
                strict=True,
            ):
                modal = a_t * np.cos(float(m) * phi) + b_t * np.sin(float(m) * phi)
                term_noise += float(w_m) * modal

            term_noise = (noise_strength * scale) * self._g_noise * term_noise

        term_eddy = np.zeros_like(term1)
        if dc.eddy_strength > 0.0 and np.asarray(self._disturbance_jitter["eddy_amp"]).size > 0:
            eddy_phi = self._interp_rows(t, np.asarray(self._disturbance_jitter["eddy_phi"], dtype=float))
            eddy_r = self._interp_rows(t, np.asarray(self._disturbance_jitter["eddy_r"], dtype=float))
            eddy_amp = self._interp_rows(t, np.asarray(self._disturbance_jitter["eddy_amp"], dtype=float))
            eddy_theta = self._interp_rows(t, np.asarray(self._disturbance_jitter["eddy_theta"], dtype=float))
            eddy_tilt = np.asarray(self._disturbance_jitter["eddy_tilt"], dtype=float)
            sigma_r = max(dc.eddy_sigma_r, 1.0e-8)
            sigma_phi = max(dc.eddy_sigma_phi, 1.0e-8)
            r_col = self.bundle.r[:, None]
            phi_grid = self.bundle.phi[None, :]
            geometry = str(getattr(dc, "eddy_geometry", "fixed_pair")).strip().lower()
            phi_sep_base = 0.75 * sigma_phi
            r_sep_base = 0.45 * sigma_r
            for idx, (phi_c, r_c, amp_c, theta_c) in enumerate(zip(eddy_phi, eddy_r, eddy_amp, eddy_theta, strict=True)):
                if geometry == "rotating_pair":
                    theta_val = float(theta_c)
                    r_sep = 0.42 * sigma_r * np.cos(theta_val)
                    phi_sep = 0.85 * sigma_phi * np.sin(theta_val)
                else:
                    tilt = float(eddy_tilt[idx]) if idx < eddy_tilt.size else 1.0
                    r_sep = tilt * r_sep_base
                    phi_sep = phi_sep_base

                r_plus = float(np.clip(r_c + r_sep, self.bundle.r[0], self.bundle.r[-1]))
                r_minus = float(np.clip(r_c - r_sep, self.bundle.r[0], self.bundle.r[-1]))
                phi_plus = float(phi_c + phi_sep)
                phi_minus = float(phi_c - phi_sep)

                env_plus = np.exp(-((r_col - r_plus) ** 2) / (2.0 * sigma_r**2)) * np.exp(
                    -(_periodic_angle_delta(phi_grid, phi_plus) ** 2) / (2.0 * sigma_phi**2)
                )
                env_minus = np.exp(-((r_col - r_minus) ** 2) / (2.0 * sigma_r**2)) * np.exp(
                    -(_periodic_angle_delta(phi_grid, phi_minus) ** 2) / (2.0 * sigma_phi**2)
                )
                term_eddy += float(amp_c) * (env_plus - env_minus)
            term_eddy *= float(scale * dc.eddy_strength)

        total = term1 + term2 + term_ms + term_noise + term_eddy
        total -= float(np.sum(total * self._area) / np.sum(self._area))
        return total

    def forcing(self, u: np.ndarray, t: float) -> np.ndarray:
        """Assemble control-actuated vorticity forcing ``f_omega(r,phi,t;u)``.

        Args:
            u: Control vector ``[u0, ..., u4]``.
            t: Simulation time.

        Returns:
            Control-actuated forcing field for the vorticity equation.
        """

        _ = t
        u = np.asarray(u, dtype=float).reshape(5)
        drive = u[0] * self.basis.b0
        ctrl = u[1] * self.basis.b1 + u[2] * self.basis.b2 + u[3] * self.basis.b3 + u[4] * self.basis.b4
        return drive + ctrl

    def _rhs(
        self,
        n: np.ndarray,
        omega: np.ndarray,
        t: float,
        u: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Evaluate the mHW right-hand side for one RK stage."""

        if str(getattr(self.cfg.pde, "dynamics_model", "mhw")).strip().lower() == "wave_landau":
            return self._rhs_wave_landau(n, omega, t, u)

        p = self.cfg.pde
        n_eval = np.array(n, copy=True)
        omega_eval = np.array(omega, copy=True)
        n_eval = np.clip(n_eval, 1e-8, self._n_cap)
        omega_eval = np.clip(omega_eval, -self._omega_cap, self._omega_cap)
        _enforce_neumann_radial(n_eval)
        _enforce_neumann_radial(omega_eval)

        phi = self.solve_poisson(omega_eval)
        total_forcing = self.forcing(u, t) + self.disturbance(t)
        threshold_enabled = (
            float(getattr(p, "critical_gradient_ratio", 0.0)) > 0.0
            or float(getattr(p, "shear_suppression_gain", 0.0)) != 0.0
            or float(getattr(p, "shear_damping_gain", 0.0)) != 0.0
        )
        grad_signed, grad_mag, shear_norm, margin, activation = self._instability_state(n_eval, phi)
        if not threshold_enabled:
            activation = np.ones_like(n_eval)
        coupling_gain = 1.0 + float(getattr(p, "supercritical_coupling_gain", 0.0)) * activation
        coupling = p.C * coupling_gain * (self._tilde(phi) - self._tilde(n_eval))
        drift_drive = self._drift_drive_profile(n_eval, phi=phi, activation=activation)
        dn_dphi = _diff_phi_periodic(n_eval, self.bundle.dphi)
        density_curvature = float(p.curvature_omega_gain) * activation * drift_drive * dn_dphi
        grad_ref_r = max(self.cfg.ring.n_amp / max(self.cfg.ring.sigma_star, 1.0e-12), 1.0e-12)
        grad_ref_phi = max(self.cfg.ring.n_amp, 1.0e-12)
        dn_dr_tilde = _diff_r_centered(self._tilde(n_eval), self.bundle.dr)
        density_baroclinic = (
            float(p.baroclinic_omega_gain)
            * activation
            * drift_drive
            * (dn_dr_tilde / grad_ref_r)
            * (dn_dphi / grad_ref_phi)
        )
        density_baroclinic -= float(np.sum(density_baroclinic * self._area) / np.sum(self._area))
        flux_balance_source = activation * self._flux_balance_omega_source(n_eval, phi)
        landau_phi = float(getattr(p, "landau_phi_gain", 0.0)) * self._tilde(phi)
        shear_damping = float(getattr(p, "shear_damping_gain", 0.0)) * shear_norm
        phase_speed = self._phase_advection_speed(grad_signed, activation)
        phase_adv_n = phase_speed * _diff_phi_periodic(self._tilde(n_eval), self.bundle.dphi)
        phase_adv_omega = phase_speed * _diff_phi_periodic(self._tilde(omega_eval), self.bundle.dphi)
        feedback_gain = float(getattr(p, "supercritical_feedback_gain", 0.0)) * activation
        n_tilde = self._tilde(n_eval)
        omega_tilde = self._tilde(omega_eval)
        n_ref = max(self.cfg.ring.n_amp, 1.0e-8)
        omega_ref = max(0.15 * self._omega_cap, 1.0)
        density_feedback = 0.35 * feedback_gain * (n_tilde - 0.35 * (n_tilde**3) / (n_ref**2))
        omega_feedback = 1.00 * feedback_gain * (omega_tilde - 0.30 * (omega_tilde**3) / (omega_ref**2))
        transport_gain = 1.0 + float(getattr(p, "supercritical_transport_gain", 0.0)) * activation
        bracket_n = transport_gain * self._arakawa_bracket(phi, n_eval, self.bundle.dr, self.bundle.dphi)
        bracket_omega = transport_gain * self._arakawa_bracket(phi, omega_eval, self.bundle.dr, self.bundle.dphi)

        dn_dt = (
            coupling
            - bracket_n
            - drift_drive * _diff_phi_periodic(phi, self.bundle.dphi)
            - phase_adv_n
            + density_feedback
            + self._hyperdiff(n_eval, p.nu_n, int(p.hyper_p))
            + self._mass_feedback_source(n_eval)
            - self.kappa_w * n_eval
            - shear_damping * self._tilde(n_eval)
        )
        domega_dt = (
            coupling
            - bracket_omega
            + self._hyperdiff(omega_eval, p.nu_omega, int(p.hyper_p))
            - p.gamma_omega * omega_eval
            - phase_adv_omega
            - 0.75 * shear_damping * self._tilde(omega_eval)
            + density_curvature
            + density_baroclinic
            + flux_balance_source
            + omega_feedback
            + landau_phi
            + total_forcing
        )
        return dn_dt, domega_dt, phi, total_forcing

    def _rhs_wave_landau(
        self,
        n: np.ndarray,
        omega: np.ndarray,
        t: float,
        u: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Evaluate the alternative near-threshold wave/instability branch.

        This branch is intentionally more direct than the mHW closure. It uses
        propagating ``phi``-coupled wave terms together with a thresholded
        Landau growth/saturation law in nonzonal vorticity so the ring can stay
        coherent and fast below threshold, then break once seeded
        perturbations become supercritical.
        """

        p = self.cfg.pde
        n_eval = np.array(n, copy=True)
        omega_eval = np.array(omega, copy=True)
        n_eval = np.clip(n_eval, 1e-8, self._n_cap)
        omega_eval = np.clip(omega_eval, -self._omega_cap, self._omega_cap)
        _enforce_neumann_radial(n_eval)
        _enforce_neumann_radial(omega_eval)

        phi = self.solve_poisson(omega_eval)
        total_forcing = self.forcing(u, t) + self.disturbance(t)
        threshold_enabled = (
            float(getattr(p, "critical_gradient_ratio", 0.0)) > 0.0
            or float(getattr(p, "shear_suppression_gain", 0.0)) != 0.0
            or float(getattr(p, "shear_damping_gain", 0.0)) != 0.0
        )
        grad_signed, _, shear_norm, _, activation = self._instability_state(n_eval, phi)
        if not threshold_enabled:
            activation = np.ones_like(n_eval)

        n_tilde = self._tilde(n_eval)
        omega_tilde = self._tilde(omega_eval)
        phi_tilde = self._tilde(phi)
        packet_cut = int(max(0, getattr(p, "packet_mode_cut", 0)))
        n_packet = self._tilde(_lowpass_phi_modes(n_tilde, packet_cut))
        omega_packet = self._tilde(_lowpass_phi_modes(omega_tilde, packet_cut))
        phi_packet = self._tilde(_lowpass_phi_modes(phi_tilde, packet_cut))
        omega_high = omega_tilde - omega_packet

        transport_gain = 1.0 + float(getattr(p, "supercritical_transport_gain", 0.0)) * activation
        bracket_n = transport_gain * self._arakawa_bracket(phi, n_eval, self.bundle.dr, self.bundle.dphi)
        bracket_omega = transport_gain * self._arakawa_bracket(phi, omega_eval, self.bundle.dr, self.bundle.dphi)

        phase_speed = self._phase_advection_speed(grad_signed, activation)
        dn_tilde_dphi = _diff_phi_periodic(n_packet, self.bundle.dphi)
        domega_tilde_dphi = _diff_phi_periodic(omega_packet, self.bundle.dphi)
        dn_tilde_dr = _diff_r_centered(n_packet, self.bundle.dr)
        domega_tilde_dr = _diff_r_centered(omega_packet, self.bundle.dr)
        phase_adv_n = phase_speed * dn_tilde_dphi
        phase_adv_omega = phase_speed * domega_tilde_dphi

        packet_phi_speed, packet_r_speed, packet_burst = self._wave_packet_speeds(n_packet, omega_packet, activation)
        confinement_flux, _ = self._radial_confinement_flux(n_eval, activation)
        zonal_restore, _ = self._zonal_profile_restore(n_eval, activation)
        collapse_guard_n, collapse_guard_omega, _ = self._collapse_guard_terms(n_eval, omega_eval, omega_packet)

        drift_drive = self._drift_drive_profile(n_eval, phi=phi, activation=activation)
        wave_n = (
            float(getattr(p, "wave_density_coupling_gain", 0.0))
            * (packet_phi_speed * domega_tilde_dphi + packet_r_speed * domega_tilde_dr)
            - drift_drive * _diff_phi_periodic(phi, self.bundle.dphi)
        )
        wave_omega = float(getattr(p, "wave_vorticity_coupling_gain", 0.0)) * (
            packet_phi_speed * dn_tilde_dphi + packet_r_speed * dn_tilde_dr
        )
        wave_n -= self._zonal_mean(wave_n)
        wave_omega -= float(np.sum(wave_omega * self._area) / np.sum(self._area))

        coupling_gain = 1.0 + float(getattr(p, "supercritical_coupling_gain", 0.0)) * activation
        coupling = p.C * coupling_gain * (phi_tilde - n_tilde)

        shear_damping = float(getattr(p, "shear_damping_gain", 0.0)) * shear_norm

        n_ref = max(self.cfg.ring.n_amp, 1.0e-8)
        omega_ref = max(0.15 * self._omega_cap, 1.0)
        density_growth_gain = float(getattr(p, "density_landau_gain", 0.0)) * activation
        omega_growth_gain = float(getattr(p, "omega_landau_gain", 0.0)) * activation
        density_growth = density_growth_gain * (
            n_packet - float(getattr(p, "density_landau_sat", 0.30)) * (n_packet**3) / (n_ref**2)
        )
        omega_growth = omega_growth_gain * (
            omega_packet - float(getattr(p, "omega_landau_sat", 0.30)) * (omega_packet**3) / (omega_ref**2)
        )
        omega_density_driver = 0.85 * omega_packet + 0.15 * omega_tilde
        density_vortex = (
            float(getattr(p, "density_vortex_gain", 0.0))
            * (0.20 + 0.85 * activation)
            * packet_burst
            * n_ref
            * np.tanh(omega_density_driver / omega_ref)
        )
        density_vortex -= self._zonal_mean(density_vortex)
        density_packet_drive = (
            float(getattr(p, "density_packet_drive_gain", 0.0))
            * self.S
            * packet_burst
            * (0.12 + activation)
            * np.tanh((omega_packet + 0.35 * phi_packet) / omega_ref)
        )
        density_packet_drive -= self._zonal_mean(density_packet_drive)
        packet_gap = _smooth_positive(
            np.abs(omega_density_driver) / omega_ref - 0.45 * np.abs(n_tilde) / n_ref,
            0.06,
        )
        density_retention = (
            float(getattr(p, "density_retention_gain", 0.0))
            * n_ref
            * packet_burst
            * (0.15 + activation)
            * packet_gap
            * np.tanh((omega_density_driver + 0.35 * phi_packet) / omega_ref)
        )
        density_retention -= self._zonal_mean(density_retention)

        feedback_gain = float(getattr(p, "supercritical_feedback_gain", 0.0)) * activation
        density_feedback = 0.30 * feedback_gain * phi_tilde
        omega_feedback = 0.45 * feedback_gain * n_tilde
        landau_phi = float(getattr(p, "landau_phi_gain", 0.0)) * phi_tilde
        omega_high_damp = float(getattr(p, "omega_high_damp", 0.0)) * omega_high

        dn_dt = (
            coupling
            - bracket_n
            - phase_adv_n
            + wave_n
            + density_growth
            + density_vortex
            + density_packet_drive
            + density_retention
            + density_feedback
            + confinement_flux
            + zonal_restore
            + collapse_guard_n
            + self._hyperdiff(n_eval, p.nu_n, int(p.hyper_p))
            + self._mass_feedback_source(n_eval)
            - self.kappa_w * n_eval
            - shear_damping * n_tilde
        )
        dn_dt = self._limit_density_rhs(n_eval, dn_dt)
        domega_dt = (
            coupling
            - bracket_omega
            - phase_adv_omega
            + wave_omega
            + omega_growth
            + omega_feedback
            + landau_phi
            + self._hyperdiff(omega_eval, p.nu_omega, int(p.hyper_p))
            - p.gamma_omega * omega_eval
            - omega_high_damp
            - collapse_guard_omega
            - 0.50 * shear_damping * omega_tilde
            + total_forcing
        )
        domega_dt = self._limit_omega_rhs(omega_eval, domega_dt)
        return dn_dt, domega_dt, phi, total_forcing

    def step(
        self,
        n: np.ndarray,
        omega: np.ndarray,
        u: np.ndarray,
        t: float,
        dt: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, StepDiagnostics]:
        """Advance ``(n, omega, psi)`` by one RK4 time step.

        Args:
            n: Current density field.
            omega: Current vorticity field.
            u: Control vector at the current step.
            t: Current simulation time.
            dt: Time-step size.

        Returns:
            Tuple ``(n_next, omega_next, psi_next, diagnostics)``.
        """

        u = np.asarray(u, dtype=float).reshape(5)

        k1_n, k1_o, _, _ = self._rhs(n, omega, t, u)
        k2_n, k2_o, _, _ = self._rhs(
            n + 0.5 * dt * k1_n,
            omega + 0.5 * dt * k1_o,
            t + 0.5 * dt,
            u,
        )
        k3_n, k3_o, _, _ = self._rhs(
            n + 0.5 * dt * k2_n,
            omega + 0.5 * dt * k2_o,
            t + 0.5 * dt,
            u,
        )
        k4_n, k4_o, _, _ = self._rhs(
            n + dt * k3_n,
            omega + dt * k3_o,
            t + dt,
            u,
        )

        n_next = n + (dt / 6.0) * (k1_n + 2.0 * k2_n + 2.0 * k3_n + k4_n)
        omega_next = omega + (dt / 6.0) * (k1_o + 2.0 * k2_o + 2.0 * k3_o + k4_o)
        n_next = self._project_density_mass(n_next)
        omega_next = np.clip(omega_next, -self._omega_cap, self._omega_cap)
        _enforce_neumann_radial(omega_next)

        psi_next = self.solve_poisson(omega_next)
        u_r_next, u_phi_next = self.velocity(psi_next)
        f_total = self.forcing(u, t + dt) + self.disturbance(t + dt)

        diag = StepDiagnostics(forcing=f_total, u_r=u_r_next, u_phi=u_phi_next)
        return n_next, omega_next, psi_next, diag
