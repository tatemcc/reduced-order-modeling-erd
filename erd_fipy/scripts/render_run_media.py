"""Render GIF diagnostics from an existing ERD run folder.

This script is intentionally narrower than the full run finalizer in
``erd_fipy.io``. It postprocesses an already-completed run directory into
rectangular and polar GIFs without reopening the whole simulation pipeline.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import h5py
import imageio.v2 as imageio
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

if __package__ in (None, ""):
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)

from erd_fipy.io import _polar_mesh


STYLE = {
    "n": ("n", "viridis"),
    "omega": ("omega", "RdBu_r"),
    "u_mag": ("|u|", "magma"),
    "u_r": ("u_r", "RdBu_r"),
    "u_phi": ("u_phi", "RdBu_r"),
}


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for run-folder GIF rendering.

    Args:
        None.

    Returns:
        Configured argument parser.
    """

    parser = argparse.ArgumentParser(description="Render GIFs from an existing ERD run folder")
    parser.add_argument("--run-dir", required=True, type=str, help="Existing run directory containing fields/snapshots.h5")
    parser.add_argument(
        "--field",
        action="append",
        choices=sorted(STYLE),
        help="Render only the selected field(s); defaults to all available fields.",
    )
    parser.add_argument("--fps", type=int, default=15, help="GIF frame rate.")
    return parser


def _canvas_to_rgb(fig: plt.Figure) -> np.ndarray:
    """Convert a drawn Matplotlib canvas to an RGB frame.

    Args:
        fig: Figure with a fully drawn canvas.

    Returns:
        RGB image array.
    """

    fig.canvas.draw()
    width, height = fig.canvas.get_width_height()
    return np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(height, width, 4)[..., :3].copy()


def _frame_imshow(field: np.ndarray, phi: np.ndarray, r: np.ndarray, label: str, cmap: str, vmin: float, vmax: float) -> np.ndarray:
    """Render one rectangular field frame.

    Args:
        field: Two-dimensional field on the ``(r, phi)`` grid.
        phi: Azimuthal coordinate array.
        r: Radial coordinate array.
        label: Colorbar label.
        cmap: Matplotlib colormap name.
        vmin: Fixed lower color bound.
        vmax: Fixed upper color bound.

    Returns:
        RGB frame array.
    """

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    im = ax.imshow(
        field,
        origin="lower",
        aspect="auto",
        extent=(float(phi[0]), float(phi[-1]), float(r[0]), float(r[-1])),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_xlabel("phi")
    ax.set_ylabel("r")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(label)
    fig.tight_layout()
    frame = _canvas_to_rgb(fig)
    plt.close(fig)
    return frame


def _frame_polar(field: np.ndarray, phi: np.ndarray, r: np.ndarray, label: str, cmap: str, vmin: float, vmax: float) -> np.ndarray:
    """Render one annular polar frame.

    Args:
        field: Two-dimensional field on the ``(r, phi)`` grid.
        phi: Azimuthal coordinate array.
        r: Radial coordinate array.
        label: Colorbar label.
        cmap: Matplotlib colormap name.
        vmin: Fixed lower color bound.
        vmax: Fixed upper color bound.

    Returns:
        RGB frame array.
    """

    x, y = _polar_mesh(phi, r)
    fig, ax = plt.subplots(figsize=(5.8, 5.8))
    pcm = ax.pcolormesh(x, y, field, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_aspect("equal")
    ax.set_axis_off()
    cbar = fig.colorbar(pcm, ax=ax, shrink=0.82)
    cbar.set_label(label)
    fig.tight_layout()
    frame = _canvas_to_rgb(fig)
    plt.close(fig)
    return frame


def _render_field_stack(data: np.ndarray, phi: np.ndarray, r: np.ndarray, run_dir: Path, name: str, fps: int) -> None:
    """Render rectangular and polar GIFs for one stored field stack.

    Args:
        data: Snapshot tensor with shape ``(n_frames, N_r, N_phi)``.
        phi: Azimuthal coordinate array.
        r: Radial coordinate array.
        run_dir: Existing run directory.
        name: Dataset / output stem.
        fps: GIF frame rate.

    Returns:
        None.
    """

    label, cmap = STYLE[name]
    vmin = float(np.nanpercentile(data, 2.0))
    vmax = float(np.nanpercentile(data, 98.0))
    if np.isclose(vmin, vmax):
        delta = max(1.0e-6, abs(vmax) * 0.1 + 1.0e-6)
        vmin -= delta
        vmax += delta

    rect_frames = [_frame_imshow(frame, phi, r, label, cmap, vmin, vmax) for frame in data]
    polar_frames = [_frame_polar(frame, phi, r, label, cmap, vmin, vmax) for frame in data]

    movies_dir = run_dir / "movies"
    movies_dir.mkdir(parents=True, exist_ok=True)
    imageio.mimsave(movies_dir / f"{name}.gif", rect_frames, fps=fps)
    imageio.mimsave(movies_dir / f"{name}_polar.gif", polar_frames, fps=fps)


def main() -> None:
    """Load one run folder and render GIFs for available field datasets.

    Args:
        None.

    Returns:
        None.
    """

    args = build_parser().parse_args()
    run_dir = Path(args.run_dir).resolve()
    h5_path = run_dir / "fields" / "snapshots.h5"
    if not h5_path.is_file():
        raise FileNotFoundError(f"Missing snapshot file: {h5_path}")

    with h5py.File(h5_path, "r") as f:
        r = np.asarray(f["r"][:], dtype=float)
        phi = np.asarray(f["phi"][:], dtype=float)
        available = [name for name in STYLE if name in f]
        selected = args.field or available
        for name in selected:
            if name not in f:
                continue
            _render_field_stack(np.asarray(f[name][:], dtype=float), phi, r, run_dir, name, int(args.fps))

    print(run_dir)


if __name__ == "__main__":
    main()
