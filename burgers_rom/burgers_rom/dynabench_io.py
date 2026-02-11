"""
I/O utilities for loading DynaBench GRID data for the Burgers equation.

This module:
- Constructs DynabenchIterator objects (optionally downloading via Dynabench)
- Extracts trajectories into a standardized tensor layout
  suitable for POD and SINDy
"""

from __future__ import annotations

from typing import Iterable, List, Optional, Tuple

import numpy as np
# import torch
from dynabench.dataset import DynabenchIterator


def build_iterator(
    split: str,
    structure: str,
    resolution: str,
    base_path: str,
    lookback: int,
    rollout: int,
    squeeze_lookback_dim: bool,
    download: bool,
) -> DynabenchIterator:
    """
    Construct a DynabenchIterator for Burgers data.

    Parameters
    ----------
    split : str
        Dataset split: 'train', 'val', or 'test'.
    structure : str
        Data structure type. Must be 'grid'.
    resolution : str
        Resolution name: 'low', 'medium', or 'high'.
    base_path : str
        Filesystem location where data are stored (passed to Dynabench).
    lookback : int
        Number of past steps provided by the iterator.
    rollout : int
        Number of future steps provided by the iterator.
    squeeze_lookback_dim : bool
        Whether to remove the lookback dimension when lookback == 1.
    download : bool
        If True, let Dynabench download missing data to base_path.

    Returns
    -------
    DynabenchIterator
        Iterator yielding Burgers samples.
    """
    return DynabenchIterator(
        split=split,
        equation="burgers",
        structure=structure,
        resolution=resolution,
        base_path=base_path,
        lookback=lookback,
        rollout=rollout,
        squeeze_lookback_dim=squeeze_lookback_dim,
        download=download,
    )


def load_trajectories(
    iterator: DynabenchIterator,
    n_trajectories: int,
    trajectory_ids: Optional[Iterable[int]] = None,
    time_stride: int = 1,
    time_limit: Optional[int] = None,
) -> np.ndarray:
    """
    Load full trajectories from a DynabenchIterator into a dense tensor.

    Output layout is:
        U[t, c, y, x]

    where:
        t = time index
        c = component index (2 for Burgers)
        y, x = spatial grid coordinates

    Parameters
    ----------
    iterator : DynabenchIterator
        Iterator yielding Burgers samples.
    n_trajectories : int
        Number of trajectories to load.
    trajectory_ids : iterable of int, optional
        Explicit trajectory indices to load.
    time_stride : int
        Subsampling factor along the time axis.
    time_limit : int, optional
        Maximum number of time steps to keep per trajectory.

    Returns
    -------
    np.ndarray
        Array of shape (n_traj, T, C, Ny, Nx).
    """
    trajectories: List[np.ndarray] = []

    if trajectory_ids is None:
        trajectory_ids = range(n_trajectories)

    for idx in trajectory_ids:
        sample = iterator[idx]
        # DynabenchIterator returns a DataItem with .x holding the window
        data = sample.x

        # if isinstance(data, torch.Tensor):
        #     data = data.cpu().numpy()
        # else:
        #     data = np.asarray(data)
        data = np.asarray(data)

        if data.ndim == 3:
            data = np.expand_dims(data, axis=0)  # ensure time axis exists

        if time_stride > 1:
            data = data[::time_stride]

        if time_limit is not None:
            data = data[:time_limit]

        trajectories.append(data)

        if len(trajectories) >= n_trajectories:
            break

    return np.stack(trajectories, axis=0)


def infer_grid_spacing(
    resolution: str,
) -> Tuple[float, float]:
    """
    Infer grid spacing (dx, dy) from the DynaBench GRID resolution.

    Parameters
    ----------
    resolution : str
        Resolution name: 'low', 'medium', or 'high'.

    Returns
    -------
    (float, float)
        Spatial grid spacings (dx, dy).
    """
    if resolution == "low":
        nx = ny = 15
    elif resolution == "medium":
        nx = ny = 22
    elif resolution == "high":
        nx = ny = 30
    else:
        raise ValueError(f"Unsupported resolution: {resolution}")

    dx = 1.0 / (nx - 1)
    dy = 1.0 / (ny - 1)

    return dx, dy
