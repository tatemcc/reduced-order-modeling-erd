"""
fields.py
---------
Surrogate backend for azimuthal E-field and axial B-field in a cylindrical ERD.

API:
    compute_fields(r_m, ne_m3, Te_eV, sigma_profile=None) -> FieldProfiles

- Default: uses a uniform-σ "Bessel-like" surrogate (no FiPy coupling), fast and
  phase-stepping friendly.
- If `sigma_profile` is provided (callable or array), we derive an effective
  σ̄ from it for the skin-depth scaling while still returning Eφ(r) explicitly.
  This lets you explore arbitrary σ without rewriting the PDEs.

Later, we can swap this module for a variable-σ 1D(r) ODE solver (same API).
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from .config import geometry as G, rf as RF
from .closures import sigma_uniform_equiv


@dataclass
class FieldProfiles:
    r_m: np.ndarray        # radial cell centers [m]
    Ephi_Vpm: np.ndarray   # azimuthal electric field [V/m] (magnitude profile)
    Bz_T: np.ndarray       # axial magnetic field [T] (magnitude profile)


def _parabolic_Ephi(r: np.ndarray, R: float, E0: float) -> np.ndarray:
    """
    Smooth non-singular Eφ profile: finite at axis, vanishes at wall.
    """
    x = np.clip(r / max(R, 1e-9), 0.0, 1.0)
    return E0 * (1.0 - x**2)


def _effective_sigma(ne_m3, Te_eV, sigma_profile, r_m) -> float:
    """
    Return an effective uniform conductivity σ̄ for skin-depth scaling.
    Priority:
      1) If sigma_profile is an array-like matching r_m -> average (area-weighted).
      2) If sigma_profile is callable(r) -> sample & average.
      3) Else: use closures.sigma_uniform_equiv(ne_avg, Te_avg).
    """
    # Case (1) array-like
    if sigma_profile is not None and hasattr(sigma_profile, "__array__"):
        sig = np.asarray(sigma_profile, dtype=float)
        if sig.shape != r_m.shape:
            raise ValueError("sigma_profile array must have same shape as r_m")
        # area weighting in cylinder: weight ~ r
        w = np.maximum(r_m, 1e-9)
        return float(np.sum(sig * w) / np.sum(w))

    # Case (2) callable
    if callable(sigma_profile):
        sig = np.asarray([float(sigma_profile(r)) for r in r_m])
        w = np.maximum(r_m, 1e-9)
        return float(np.sum(sig * w) / np.sum(w))

    # Case (3) fallback to uniform σ from average ne, Te
    ne_avg = float(np.mean(ne_m3)) if np.ndim(ne_m3) else float(ne_m3)
    Te_avg = float(np.mean(Te_eV)) if np.ndim(Te_eV) else float(Te_eV)
    return sigma_uniform_equiv(ne_avg, Te_avg)


def compute_fields(
    r_m: np.ndarray,
    ne_m3: np.ndarray | float,
    Te_eV: np.ndarray | float,
    sigma_profile: np.ndarray | callable | None = None,
) -> FieldProfiles:
    """
    Compute Eφ(r) and Bz(r) quickly from a surrogate.

    Parameters
    ----------
    r_m : (Nr,) ndarray
        Radial cell-center coordinates [m].
    ne_m3, Te_eV : scalar or arrays
        Electron density and temperature (used for σ; arrays will be averaged).
    sigma_profile : array-like or callable, optional
        If array-like of shape (Nr,), interpreted as σ(r) [S/m].
        If callable, evaluated as σ(r) over r_m.
        If None, a uniform σ̄ is computed from (ne,Te) via closures.

    Returns
    -------
    FieldProfiles
        r_m, Ephi(r), Bz(r) magnitude profiles.
    """
    R = G.R_cm * 1e-2
    mu0 = 4e-7 * np.pi
    omega = 2.0 * np.pi * RF.freq_Hz
    E0 = float(RF.E0_Vpm)

    # Eφ profile (magnitude); phase handled by controller/scheduler elsewhere
    Ephi = _parabolic_Ephi(r_m, R, E0)

    # Effective σ̄ for skin-depth-like attenuation of Bz
    sigma_bar = max(_effective_sigma(ne_m3, Te_eV, sigma_profile, r_m), 1e-9)

    # Characteristic attenuation wavenumber (quasi-static, conduction-dominated)
    # k ≈ sqrt(ω μ0 σ / 2)  → skin depth δ ~ sqrt(2/(ω μ0 σ))
    k = np.sqrt(omega * mu0 * sigma_bar / 2.0)

    # Scale Bz so that characteristic |∂t B| ~ |∇×E|
    # Use |B| ~ |E| / ω near axis, then apply radial attenuation exp(-k r)
    B_scale = max(E0, 1e-12) / max(omega, 1e-12)
    Bz = B_scale * np.exp(-k * r_m)

    return FieldProfiles(r_m=r_m, Ephi_Vpm=Ephi, Bz_T=Bz)
