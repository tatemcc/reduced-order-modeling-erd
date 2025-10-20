"""
closures.py
-----------
Constitutive relations and source terms for the ERD FiPy plant.

These are intentionally compact. All knobs live in `ClosureParams` so you can
tune without code changes. Units are noted in each docstring.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from .config import gas as GAS, icbc as ICBC

# Physical constants
e = 1.602176634e-19         # C
me = 9.10938356e-31         # kg
kB = 1.380649e-23           # J/K
mi_Xe = 2.1801714e-25       # ~131.3 amu in kg (Xenon, default gas)


@dataclass
class ClosureParams:
    # Electron-neutral collision frequency: νc ≈ ν0 * p_Pa / sqrt(Tg)
    # Very rough scaling; replace when you have a preferred cross-section model.
    nu0_HzPa: float = 5.0e6

    # Ambipolar diffusion: Da ≈ Da0 * Te[eV] / p[Torr]
    Da0_m2ps_per_eV_over_Torr: float = 0.05

    # Ionization (Townsend-like): S = alpha0 * p[Torr] * exp( - Ethr * p / |E| ) * ne
    alpha0_1ps_per_Torr: float = 5.0e2
    Ethr_over_Vpm_per_Torr: float = 2.0e3

    # Electron thermal conductivity (lumped): kappa ≈ k0 * Te[eV]  [W/m/K]
    k0_WpmK_per_eV: float = 0.5

    # Cooling: Qloss = c1 * ne * Te + c2 * ne  [W/m^3]
    c1_Wpm3peV: float = 1.0e-17
    c2_Wpm3: float = 5.0e-18


params = ClosureParams()


def nu_c_Hz(p_Torr: float, Tgas_K: float) -> float:
    """
    Electron-neutral collision frequency [Hz].
    νc ≈ ν0 * p_Pa / Tgas_K
    """
    p_Pa = p_Torr * 133.322368
    return params.nu0_HzPa * p_Pa / max(Tgas_K, 1.0)


def sigma_Spm(ne_m3: float, Te_eV: float, p_Torr: float, Tgas_K: float) -> float:
    """
    Ohmic conductivity [S/m]: σ ≈ e^2 ne / (me νc).
    Ignores temperature dependence of νc beyond Tgas; adequate for ERD prototype.
    """
    nu = nu_c_Hz(p_Torr, Tgas_K)
    return (e * e * max(ne_m3, 0.0)) / (me * max(nu, 1.0))


def sigma_uniform_equiv(ne_avg: float, Te_avg: float) -> float:
    """
    Convenience wrapper for a uniform effective σ̄ based on average (ne, Te)
    under current gas conditions from config.
    """
    return sigma_Spm(ne_avg, Te_avg, GAS.p_Torr, GAS.Tgas_K)


def Da_m2ps(Te_eV: float, p_Torr: float) -> float:
    """
    Ambipolar diffusion coefficient [m^2/s].
    Simple scaling Da ~ Te / p.
    """
    return params.Da0_m2ps_per_eV_over_Torr * max(Te_eV, 0.01) / max(p_Torr, 0.1)


def S_ion_Hz(E_eff_Vpm: float, p_Torr: float) -> float:
    """
    Townsend-like ionization frequency coefficient [1/s] multiplying ne.
    S = α0 * p * exp( - Ethr * p / |E| )
    """
    x = params.Ethr_over_Vpm_per_Torr * p_Torr / max(abs(E_eff_Vpm), 1.0)
    x = np.clip(x, -100.0, 100.0)
    return params.alpha0_1ps_per_Torr * p_Torr * np.exp(-x)


def Q_ohmic_Wpm3(sigma_Spm_val: float, E_Vpm: float) -> float:
    """
    Ohmic heating density [W/m^3]: QΩ = Re(σ) * |E|^2 (σ assumed real here).
    """
    return max(sigma_Spm_val, 0.0) * (abs(E_Vpm) ** 2)


def Q_loss_Wpm3(ne_m3: float, Te_eV: float) -> float:
    """
    Lumped cooling [W/m^3]: radiation + inelastic + misc.
    """
    return params.c1_Wpm3peV * max(ne_m3, 0.0) * max(Te_eV, 0.0) + params.c2_Wpm3 * max(ne_m3, 0.0)


def v_loss_mps(Te_eV: float) -> float:
    """
    Wall loss velocity scale [m/s], used in Robin flux BCs:
        -Da * ∂n/∂n_wall = ne * v_loss
    We use ion thermal speed scaling ~ sqrt(kB * Te[K] / mi), with a tunable multiplier.
    """
    Te_K = max(Te_eV, 0.0) * e / kB
    return float(np.sqrt(max(Te_K, 1.0) * kB / mi_Xe) * ICBC.vloss_scale)


def tau_wall_s(R_m: float, Da_m2ps_val: float) -> float:
    """
    Lumped wall-loss time constant [s] for the slowest radial diffusion mode:
        τ ≈ R^2 / (π^2 * Da)
    Suitable for capturing overall depletion rate when using a bulk loss term.
    """
    return (R_m ** 2) / (np.pi ** 2 * max(Da_m2ps_val, 1e-9))
