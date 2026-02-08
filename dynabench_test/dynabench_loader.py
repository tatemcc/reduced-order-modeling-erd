"""
Dynabench Dataset Loader

Summary of Data Structure and Simulation Parameters (based on Dulny et al., 2023):

1. Data Structure (.h5 files):
   - The dataset is split into sharded HDF5 files (e.g., `*_0_499.h5`).
   - 'data': Simulation field values. Shape: (N_sims, Time, Channels, X, Y).
     Example (Low Res): (500, 201, 1, 15, 15).
   - 'points': Spatial coordinates. Shape: (N_sims, X, Y, 2) for grid data.

2. Simulation Setup:
   - Spatial Domain: Ω = [0, 1] × [0, 1].
   - Boundary Conditions: Periodic in x and y.
   - Temporal Domain: T = [0, 200].
   - Discretization: Original simulations computed on a 64x64 grid (Δx = Δy = 0.0156).

3. Data Processing:
   - Temporal Resolution: Saved with Δt = 1, yielding 201 observations per simulation.
   - Initialization: Zeros, uniform/normal noise, or sum of Gaussian curves.
   - Normalization: Data is normalized to ensure consistent value ranges.
   - Sampling: Measurements are sampled from the domain (uniformly for cloud, grid for grid)
     and bilinearly interpolated from the original high-resolution simulation.
"""

import os
import glob
import h5py
import numpy as np


def _parse_file_indices(filename):
    """Parses start and end indices from a filename like '..._0_499.h5'."""
    try:
        parts = os.path.splitext(filename)[0].split("_")
        start_idx = int(parts[-2])
        end_idx = int(parts[-1])
        return start_idx, end_idx
    except (ValueError, IndexError):
        return None


def load_dynabench_data(
    equation="advection",
    structure="grid",
    resolution="low",
    split="train",
    base_path="data",
    sim_index=0,
):
    """
    Loads a specific simulation from the Dynabench dataset.

    This function is designed to replace DynabenchIterator and DynabenchSimulationIterator.
    It manually locates and loads data from the split HDF5 files (e.g. *_0_499.h5).

    Args:
        equation (str): The equation name (e.g., "advection").
        structure (str): The data structure ("grid" or "cloud").
        resolution (str): The resolution ("low", "medium", "high").
        split (str): The dataset split ("train", "val", "test").
        base_path (str): The root directory of the dataset.
        sim_index (int): The global index of the simulation to load.

    Returns:
        tuple: (u_data, points)
            u_data (np.ndarray): Simulation data of shape (Time, Channels, X, Y).
            points (np.ndarray): Spatial points.
    """
    print(f"      [Loader] Searching for {split} data ({structure}) in {base_path}...")

    # Construct search pattern for split files
    # Format: {equation}_{split}_{structure}_{resolution}_{start}_{end}.h5
    file_pattern = f"{equation}_{split}_{structure}_{resolution}_*.h5"
    search_path = os.path.join(base_path, equation, structure, resolution, file_pattern)
    files = sorted(glob.glob(search_path))

    if not files:
        raise FileNotFoundError(f"No .h5 files found at: {search_path}")

    target_file = None
    local_index = 0

    # Attempt to find the file containing the requested sim_index
    for fpath in files:
        fname = os.path.basename(fpath)
        indices = _parse_file_indices(fname)
        if indices:
            start_idx, end_idx = indices
            if start_idx <= sim_index <= end_idx:
                target_file = fpath
                local_index = sim_index - start_idx
                break

    # Fallback if parsing fails or index not found
    if target_file is None:
        print(f"      [Loader] Warning: Could not locate sim_index={sim_index} in filenames.")
        print(f"      [Loader] Defaulting to first file: {files[0]}")
        target_file = files[0]
        local_index = 0

    print(
        f"      [Loader] Loading file: {target_file} (Global Index: {sim_index}, Local Index: {local_index})"
    )

    with h5py.File(target_file, "r") as f:
        # Ensure index is valid within the file
        if local_index >= f["data"].shape[0]:
            raise IndexError(
                f"Local index {local_index} out of bounds for file {target_file} with {f['data'].shape[0]} simulations."
            )

        u_data = f["data"][local_index]

        if "points" in f:
            points = f["points"][local_index]
        else:
            # Generate implicit grid
            nx, ny = u_data.shape[-2], u_data.shape[-1]
            x = np.linspace(0, 1, nx)
            y = np.linspace(0, 1, ny)
            X, Y = np.meshgrid(x, y, indexing="ij")
            points = np.stack([X, Y], axis=-1)

    return u_data, points
