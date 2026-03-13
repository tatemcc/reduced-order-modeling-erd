"""Numerical operators for the toy ERD ``(n, omega, psi)`` system.

The implementation intentionally uses simple rectangularized finite differences
on a uniform ``(r, phi)`` grid:
- periodic differences in ``phi``,
- Neumann radial boundaries for ``n`` and ``omega``,
- Dirichlet radial boundaries for ``psi`` in the Poisson solve.
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


class ToyERDOperators:
    """Numerical operator bundle for one-step updates of ``(n, omega, psi)``.

    This class precomputes forcing bases and FFT-by-mode Poisson inverses so the
    simulation loop can call a lightweight :meth:`step` at each time index.
    """

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
        self._g_noise = _gauss(bundle.r, cfg.ring.r_star, 0.65 * cfg.ring.sigma_star)[:, None]

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
            "noise_cos": rng.normal(0.0, 1.0, size=noise_len),
            "noise_sin": rng.normal(0.0, 1.0, size=noise_len),
            "noise_phase": 2.0 * np.pi * rng.uniform(0.0, 1.0, size=noise_len),
        }

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

        p1 = _mode_phase(m1, t) + float(self._disturbance_jitter["carrier_m1"])
        p2 = _mode_phase(m2, t) + float(self._disturbance_jitter["carrier_m2"])

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
            for m, amp, freq, phase0, amp_j in zip(
                dc.multiscale_modes[:n_ms],
                dc.multiscale_amplitudes[:n_ms],
                dc.multiscale_frequencies[:n_ms],
                self._disturbance_jitter["ms_phase"][:n_ms],
                self._disturbance_jitter["ms_amp"][:n_ms],
            ):
                phase_t = 2.0 * np.pi * freq * t + phase0
                if is_final_stage:
                    phase_t += 0.10 * np.sin(2.0 * np.pi * dc.phase_c_boost_freq * t)
                term_ms += (amp * amp_j) * self.basis.g2 * np.cos(float(m) * phi - phase_t)
            term_ms *= float(scale * dc.multiscale_strength)

        term_noise = np.zeros_like(term1)
        noise_strength = float(max(0.0, dc.noise_strength))
        if noise_strength > 0.0 and self._disturbance_jitter["noise_modes"].size > 0:
            for m, c_m, s_m, p_m in zip(
                self._disturbance_jitter["noise_modes"],
                self._disturbance_jitter["noise_cos"],
                self._disturbance_jitter["noise_sin"],
                self._disturbance_jitter["noise_phase"],
                strict=True,
            ):
                freq_t = dc.noise_freq_base + 0.025 * float(m)
                osc_t = np.sin(2.0 * np.pi * freq_t * t + p_m)
                modal = c_m * np.cos(float(m) * phi) + s_m * np.sin(float(m) * phi)
                term_noise += osc_t * modal

            norm = np.sqrt(float(self._disturbance_jitter["noise_modes"].size))
            term_noise = (noise_strength * scale / max(norm, 1.0)) * self._g_noise * term_noise

        return term1 + term2 + term_ms + term_noise

    def forcing(self, u: np.ndarray, t: float) -> np.ndarray:
        """Assemble total vorticity forcing ``f_drive + f_ctrl + f_dist``.

        Args:
            u: Control vector ``[u0, ..., u4]``.
            t: Simulation time.

        Returns:
            Full forcing field for the vorticity equation.
        """

        u = np.asarray(u, dtype=float).reshape(5)
        drive = u[0] * self.basis.b0
        ctrl = u[1] * self.basis.b1 + u[2] * self.basis.b2 + u[3] * self.basis.b3 + u[4] * self.basis.b4
        return drive + ctrl + self.disturbance(t)

    def step(
        self,
        n: np.ndarray,
        omega: np.ndarray,
        u: np.ndarray,
        t: float,
        dt: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, StepDiagnostics]:
        """Advance ``(n, omega, psi)`` by one explicit time step.

        Step order:
        1) Solve Poisson for ``psi`` from current ``omega``.
        2) Build velocity and vorticity forcing.
        3) Update ``omega`` with advection, diffusion, damping, forcing.
        4) Recompute ``psi`` and velocity from updated ``omega``.
        5) Update ``n`` with advection, anisotropic diffusion, relaxation, wall sink.

        Args:
            n: Current density field.
            omega: Current vorticity field.
            u: Control vector at the current step.
            t: Current simulation time.
            dt: Time-step size.

        Returns:
            Tuple ``(n_next, omega_next, psi_next, diagnostics)``.
        """

        p = self.cfg.pde

        psi = self.solve_poisson(omega)
        u_r, u_phi = self.velocity(psi)
        f_omega = self.forcing(u, t)

        domega_dt = (
            -u_r * _diff_r_centered(omega, self.bundle.dr)
            -u_phi * _diff_phi_periodic(omega, self.bundle.dphi)
            + p.nu * (_d2_r_neumann(omega, self.bundle.dr) + _d2_phi_periodic(omega, self.bundle.dphi))
            - p.gamma * omega
            + p.beta_instab * self.basis.b0 * (n - self.n_eq)
            + f_omega
        )

        omega_next = omega + dt * domega_dt
        _enforce_neumann_radial(omega_next)

        psi_next = self.solve_poisson(omega_next)
        u_r_next, u_phi_next = self.velocity(psi_next)

        dn_dt = (
            -u_r_next * _diff_r_centered(n, self.bundle.dr)
            -u_phi_next * _diff_phi_periodic(n, self.bundle.dphi)
            + p.D_r * _d2_r_neumann(n, self.bundle.dr)
            + p.D_phi * _d2_phi_periodic(n, self.bundle.dphi)
            + p.alpha * (self.n_eq - n)
            - self.kappa_w * n
        )

        n_next = n + dt * dn_dt
        n_next = np.clip(n_next, 1e-8, None)
        _enforce_neumann_radial(n_next)

        diag = StepDiagnostics(forcing=f_omega, u_r=u_r_next, u_phi=u_phi_next)
        return n_next, omega_next, psi_next, diag
