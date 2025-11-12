"""
pdes.py
-------
FiPy PDE definitions for electron density n_e and optional electron temperature T_e.

- Axisymmetric r–z grid (we map mesh.x -> r, mesh.y -> z).
- Uses a fast fields backend to provide Eφ(r) and Bz(r) per control step.
- Density equation includes ambipolar diffusion, Townsend-like ionization, and
  a lumped wall-loss term via τ_wall.
- Energy equation (toggle in config) includes thermal diffusion, Ohmic heating,
  and lumped cooling.

Note: We begin with simple (homogeneous Neumann) boundary behavior from FiPy.
Robin-type wall flux can be added later if needed; for v1 we emulate wall losses
with the τ_wall bulk sink, which captures the correct time scale.
"""

from __future__ import annotations
import numpy as np

try:
    from fipy import CellVariable, DiffusionTerm, TransientTerm, ImplicitSourceTerm
    from fipy.tools import numerix
except ImportError as e:
    CellVariable = None

from .config import geometry as G, rf as RF, gas as GAS, icbc as ICBC
from .closures import (
    Da_m2ps,
    S_ion_Hz,
    Q_ohmic_Wpm3,
    Q_loss_Wpm3,
    sigma_Spm,
    tau_wall_s,
)


def build_state_vars(mesh):
    """
    Create FiPy CellVariables for n_e and T_e with configurable initial conditions.
    Returns (ne, Te).
    """
    if CellVariable is None:
        raise RuntimeError("FiPy is not installed. Please install FiPy to run the PDE solver.")

    Rm = G.R_cm * 1e-2
    # Parabolic initial ne profile in r, uniform in z
    ne_init = ICBC.ne0_m3 * (1.0 - (mesh.x / max(Rm, 1e-12)) ** 2)
    ne_init = np.clip(ne_init, 0.0, None)

    ne = CellVariable(name="ne", mesh=mesh, value=ne_init)
    Te = CellVariable(name="Te", mesh=mesh, value=ICBC.Te0_eV)
    return ne, Te


def build_equations(mesh, ne, Te, fields_module, sigma_profile=None):
    """
    Build FiPy equations for n_e and optionally T_e.

    Parameters
    ----------
    mesh : FiPy mesh
    ne, Te : CellVariable
    fields_module : module with compute_fields(r, ne, Te, sigma_profile=None) -> profiles
    sigma_profile : None | array-like | callable
        If provided, passed to the fields backend to define σ(r) for effective skin scaling.

    Returns
    -------
    eqs : list of FiPy terms (equations)
    profiles : FieldProfiles (r, Eφ(r), Bz(r)) for the current step
    """
    if CellVariable is None:
        raise RuntimeError("FiPy is not installed. Please install FiPy to run the PDE solver.")

    Rm = G.R_cm * 1e-2

    # Radial coordinates for cells (FiPy Grid2D: x->r, y->z)
    r_cells = numerix.array(mesh.x)

    # Update field surrogate (use current averages for uniform-σ fallback)
    profiles = fields_module.compute_fields(
        r_m=np.unique(r_cells),
        ne_m3=ne.value.mean(),
        Te_eV=Te.value.mean(),
        sigma_profile=sigma_profile,
    )

    # Interpolate Eφ(r) to cell centers
    Ephi_cells = np.interp(r_cells, profiles.r_m, profiles.Ephi_Vpm)

    # -------- Electron density equation --------
    # Ambipolar diffusion coefficient Da(Te, p) as a CellVariable for spatial variation
    Da_array = Da_m2ps(np.maximum(Te.value, 0.05), GAS.p_Torr)   # numpy array [m^2/s]
    Da_cells = CellVariable(name="Da", mesh=mesh, value=Da_array)

    # --- Ionization - loss coefficient (per-ne) as a CellVariable ---
    # S_coeff: numpy array [1/s] from |Eφ|
    # loss_coeff: scalar [1/s] from τ_wall
    S_coeff = S_ion_Hz(np.maximum(np.abs(Ephi_cells), 1.0), GAS.p_Torr)  # numpy array
    Da_mean = float(np.maximum(Te.value.mean(), 0.05))
    Da_mean = Da_m2ps(Da_mean, GAS.p_Torr)
    tauw = tau_wall_s(Rm, Da_mean)
    loss_coeff = 1.0 / max(tauw, 1e-9)

    reaction_cells = CellVariable(name="S_minus_loss", mesh=mesh, value=S_coeff - loss_coeff)

    # --- Fully implicit ne-equation ---
    ne_eq = (
        TransientTerm(var=ne)
        == DiffusionTerm(coeff=Da_cells, var=ne)
        + ImplicitSourceTerm(coeff=reaction_cells, var=ne)
    )

    eqs = [ne_eq]

    # -------- Electron energy equation (optional) --------
    if RF.use_energy_eq:
        # Heat capacity (lumped): (3/2) n_e  in "FiPy units"
        heat_capacity_arr = 1.5 * np.maximum(ne.value, 1e10)
        heat_capacity = CellVariable(name="hc", mesh=mesh, value=heat_capacity_arr)

        # Thermal conductivity κ ~ k0 * Te; keep scalar & modest to avoid stiffness
        kappa_arr = 0.5 * np.maximum(Te.value, 0.1)  # W/m/K (placeholder scaling)
        kappa = CellVariable(name="kappa", mesh=mesh, value=kappa_arr)

        # Ohmic heating density: QΩ = σ |E|^2  (σ from local (ne,Te))
        sigma_cells_arr = sigma_Spm(
            ne_m3=np.maximum(ne.value, 1e10),
            Te_eV=np.maximum(Te.value, 0.1),
            p_Torr=GAS.p_Torr,
            Tgas_K=GAS.Tgas_K,
        )
        Qohm = Q_ohmic_Wpm3(sigma_cells_arr, np.maximum(np.abs(Ephi_cells), 1.0))

        # Lumped cooling
        Qloss = Q_loss_Wpm3(np.maximum(ne.value, 1e10), np.maximum(Te.value, 0.1))
        heat_source = CellVariable(
            name="Q_net",
            mesh=mesh,
            value=Qohm - Qloss,
        )

        Te_eq = (
            TransientTerm(coeff=heat_capacity, var=Te)
            == DiffusionTerm(coeff=kappa, var=Te)
            + heat_source  # explicit source term is OK
        )
        eqs.append(Te_eq)

    return eqs, profiles
