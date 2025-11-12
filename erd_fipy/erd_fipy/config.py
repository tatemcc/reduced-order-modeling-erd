"""
Configuration for ERD FiPy plant.
All parameters live here so you can tune without touching solver code.
"""

from dataclasses import dataclass

from datetime import datetime
now = datetime.now()
timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")

@dataclass
class GeometryConfig:
    R_cm: float = 6.0        # cylinder radius [cm]
    H_cm: float = 12.0       # cylinder height [cm]
    Nr: int = 256            # radial cells
    Nz: int = 512            # axial cells

@dataclass
class RFConfig:
    freq_Hz: float = 13.56e6         # RF frequency
    n_phase: int = 32                # phase samples per RF cycle
    E0_Vpm: float = 200.0            # default amplitude [V/m]
    phase_deg: float = 0.0           # default phase [deg]
    use_energy_eq: bool = True       # toggle electron energy equation

@dataclass
class GasConfig:
    gas: str = "Xe"
    p_Torr: float = 15.0             # neutral pressure [Torr]
    Tgas_K: float = 300.0            # gas temperature [K]

@dataclass
class InitBCConfig:
    ne0_m3: float = 1e15             # initial ne center value [1/m^3]
    Te0_eV: float = 2.0              # initial Te [eV]
    wall_bc_kind: str = "robin"      # "robin" or "neumann"
    # Robin coefficient scaling for wall loss flux: -Da * dn/dr = ne * v_loss
    vloss_scale: float = 1.0         # multiplier on thermal speed used for v_loss

@dataclass
class TimeConfig:
    dt_s: float = 5e-6               # slow-time PDE step [s]
    n_steps: int = 2              # total steps (e.g., 10 ms) (should be 2000)
    save_every: int = 50             # write HDF5 every N steps

@dataclass
class OutputConfig:
    outdir: str = "outputs"
    run_name: str = f"run_{timestamp}"

@dataclass
class ModelToggles:
    use_uniform_sigma_for_fields: bool = True  # Bessel-like surrogate vs variable-sigma backend
    allow_custom_sigma_profile: bool = True    # allow user-supplied σ profile

# Instantiate defaults
geometry = GeometryConfig()
rf       = RFConfig()
gas      = GasConfig()
icbc     = InitBCConfig()
time     = TimeConfig()
output   = OutputConfig()
toggles  = ModelToggles()
