"""
Controller/plant stepping: set inputs -> update fields -> advance PDEs -> write outputs.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np, os
from time import time as pytime

from .config import rf as RF, time as TCFG, output as OUT
from .mesh import make_mesh
from . import fields as fields_backend
from .pdes import build_state_vars, build_equations
from .io import H5Writer

@dataclass
class ControlInputs:
    E0_Vpm: float
    phase_deg: float
    freq_Hz: float

def run_sim(seed: int = 0, sigma_profile=None):
    """
    Runs the ERD plant for time.n_steps slow-time steps.
    Returns path to the HDF5 results file.
    """
    np.random.seed(seed)

    mesh, r, z = make_mesh()
    ne, Te = build_state_vars(mesh)
    # only make directory if it doesn't exist
    if not os.path.exists(OUT.outdir):
        os.makedirs(OUT.outdir)
    # output names are unique, so no overwrite check needed

    writer = H5Writer(OUT.outdir, OUT.run_name, r, z)

    # print save info
    print(f"Running simulation for {TCFG.n_steps} steps, saving every {TCFG.save_every} steps at absolute directory: {OUT.outdir}")
    ti = pytime()
    for k in range(TCFG.n_steps):
        # (Optional) vary control here if you want closed-loop later
        u = ControlInputs(E0_Vpm=RF.E0_Vpm, phase_deg=RF.phase_deg, freq_Hz=RF.freq_Hz)

        # Build equations with current fields (uses ne, Te averages internally)
        eqs, prof = build_equations(mesh, ne, Te, fields_backend, sigma_profile=sigma_profile)

        # Advance one slow-time step
        dt = TCFG.dt_s
        for eq in eqs:
            eq.sweep(dt=dt)

        # Save periodically (or on final step)
        if (k % TCFG.save_every) == 0 or (k == TCFG.n_steps - 1):
            print(f"Step {k+1}/{TCFG.n_steps}: writing output...\tElapsed: {pytime() - ti:.2f} s\tStep dt: {dt:.2e} s")
            t = (k + 1) * dt
            writer.write_snapshot(t, prof.Bz_T, ne.value.copy(), Te.value.copy(), u)

    writer.close()
    return os.path.join(OUT.outdir, f"{OUT.run_name}.h5")
