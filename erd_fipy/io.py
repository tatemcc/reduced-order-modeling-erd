"""
HDF5 output writer for snapshots compatible with downstream ROM/SINDy ingestion.
"""
from __future__ import annotations
import os, h5py, numpy as np

class H5Writer:
    def __init__(self, outdir: str, run_name: str, r: np.ndarray, z: np.ndarray):
        self.path = os.path.join(outdir, f"{run_name}.h5")
        self._h = h5py.File(self.path, "w")
        # Coordinates
        self._h.create_dataset("/coords/r", data=r)
        self._h.create_dataset("/coords/z", data=z)
        # Inputs log
        self._inputs_grp = self._h.create_group("/inputs")
        # Time axis
        self._t = self._h.create_dataset("/time/t", shape=(0,), maxshape=(None,), dtype="f8")
        # Fields (extendable along time dim)
        self._Bz = self._h.create_dataset("/fields/Bz", shape=(0, r.size), maxshape=(None, r.size), dtype="f8")
        self._ne = self._h.create_dataset("/fields/ne", shape=(0, z.size, r.size), maxshape=(None, z.size, r.size), dtype="f8")
        self._Te = self._h.create_dataset("/fields/Te", shape=(0, z.size, r.size), maxshape=(None, z.size, r.size), dtype="f8")
        self._i = 0

    def _append(self, ds, arr):
        ds.resize((self._i + 1,) + ds.shape[1:])
        ds[self._i] = arr

    def write_snapshot(self, t: float, Bz_r: np.ndarray, ne_zr: np.ndarray, Te_zr: np.ndarray, inputs):
        self._append(self._t, np.array(t))
        self._append(self._Bz, Bz_r)
        self._append(self._ne, ne_zr.reshape(self._ne.shape[1:]))
        self._append(self._Te, Te_zr.reshape(self._Te.shape[1:]))
        # Log inputs for this snapshot
        g = self._inputs_grp.create_group(f"{self._i}")
        g.create_dataset("E0_Vpm", data=float(inputs.E0_Vpm))
        g.create_dataset("phase_deg", data=float(inputs.phase_deg))
        g.create_dataset("freq_Hz", data=float(inputs.freq_Hz))
        self._i += 1

    def close(self):
        if self._h is not None:
            self._h.close()
            self._h = None
