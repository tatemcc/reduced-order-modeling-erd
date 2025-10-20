"""
Mesh creation. We use a uniform r–z grid (axisymmetric interpretation).
FiPy's Grid2D works fine; we treat its x-coordinate as radius (r) and y as axial (z).
"""

import numpy as np
from .config import geometry as G

try:
    from fipy import Grid2D
except ImportError:
    Grid2D = None


def make_mesh():
    """
    Returns (mesh, r_coords, z_coords)
    Units: r, z in meters.
    """
    if Grid2D is None:
        raise RuntimeError("FiPy is not installed. Please install FiPy to run the PDE solver.")

    Rm = G.R_cm * 1e-2
    Hm = G.H_cm * 1e-2
    dr = Rm / G.Nr
    dz = Hm / G.Nz

    mesh = Grid2D(dx=dr, dy=dz, nx=G.Nr, ny=G.Nz)

    # coordinate arrays (cell centers)
    r = (np.arange(G.Nr) + 0.5) * dr
    z = (np.arange(G.Nz) + 0.5) * dz
    return mesh, r, z

make_mesh()