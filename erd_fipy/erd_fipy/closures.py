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


def nu_c_Hz(p_Torr: float, Tgas_K: float):
    """Vectorized collision frequency [Hz]."""
    import numpy as _np
    p_Pa = p_Torr * 133.322368  # scalar
    T = max(Tgas_K, 1.0)        # scalar
    return params.nu0_HzPa * p_Pa / T

def sigma_Spm(ne_m3, Te_eV, p_Torr: float, Tgas_K: float):
    """Vectorized Ohmic conductivity σ = e^2 ne / (me νc)."""
    import numpy as _np
    ne = _np.asarray(ne_m3, dtype=float)
    nu = nu_c_Hz(p_Torr, Tgas_K)  # scalar
    return (e * e * _np.maximum(ne, 0.0)) / (me * max(nu, 1.0))

def sigma_uniform_equiv(ne_avg, Te_avg):
    """Uniform effective σ̄ from average (ne,Te). Accepts scalars/arrays; uses means."""
    import numpy as _np
    ne = float(_np.mean(ne_avg))
    Te = float(_np.mean(Te_avg))
    return sigma_Spm(ne, Te, GAS.p_Torr, GAS.Tgas_K)

def Da_m2ps(Te_eV, p_Torr: float):
    """Vectorized ambipolar diffusion Da ~ Te / p."""
    import numpy as _np
    Te = _np.asarray(Te_eV, dtype=float)
    return params.Da0_m2ps_per_eV_over_Torr * _np.maximum(Te, 0.01) / max(p_Torr, 0.1)

def S_ion_Hz(E_eff_Vpm, p_Torr: float):
    """Vectorized Townsend-like ionization coefficient (per ne)."""
    import numpy as _np
    E = _np.asarray(E_eff_Vpm, dtype=float)
    denom = _np.maximum(_np.abs(E), 1.0)
    x = params.Ethr_over_Vpm_per_Torr * p_Torr / denom
    x = _np.clip(x, -100.0, 100.0)
    return params.alpha0_1ps_per_Torr * p_Torr * _np.exp(-x)

def Q_ohmic_Wpm3(sigma_Spm_val, E_Vpm):
    """Vectorized Ohmic heating QΩ = σ |E|^2."""
    import numpy as _np
    sig = _np.asarray(sigma_Spm_val, dtype=float)
    E = _np.asarray(E_Vpm, dtype=float)
    return _np.maximum(sig, 0.0) * (_np.abs(E) ** 2)

def Q_loss_Wpm3(ne_m3, Te_eV):
    """Vectorized lumped cooling Qloss = c1 ne Te + c2 ne."""
    import numpy as _np
    ne = _np.asarray(ne_m3, dtype=float)
    Te = _np.asarray(Te_eV, dtype=float)
    return params.c1_Wpm3peV * _np.maximum(ne, 0.0) * _np.maximum(Te, 0.0) + params.c2_Wpm3 * _np.maximum(ne, 0.0)


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
