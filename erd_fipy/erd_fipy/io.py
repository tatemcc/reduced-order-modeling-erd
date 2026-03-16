"""Run-folder artifact writer for ERD plant simulations.

This module creates self-contained run directories and writes all required
outputs: config, snapshots, controls, metric curves, plots, and GIF movies.
"""

from __future__ import annotations

from datetime import datetime
import json
import logging
import os
from pathlib import Path
from typing import Dict, List, Optional

import h5py
import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np

from .config import RunConfig, save_run_config
from .metrics import compute_transport_diagnostics, compute_variation_diagnostics, dominant_spectrum


SUBDIRS = ("fields", "controls", "metrics", "plots", "movies", "model", "logs", "notes")


def make_run_dir(outputs_root: str, tag: str, run_dir: Optional[str] = None) -> Path:
    """Create a timestamped run directory with standard subfolders.

    Args:
        outputs_root: Root directory that stores all run folders.
        tag: Human-readable run tag appended to the timestamp.
        run_dir: Optional explicit run directory path. When provided, the
            timestamp/tag scheme is bypassed.

    Returns:
        Absolute path to the created run directory.
    """

    if run_dir is not None:
        rd = Path(run_dir).resolve()
    else:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base = Path(outputs_root).resolve() / f"{stamp}_{tag}"
        rd = base
        suffix = 1
        while rd.exists():
            rd = Path(f"{base}_{suffix:02d}")
            suffix += 1

    rd.mkdir(parents=True, exist_ok=True)
    for name in SUBDIRS:
        (rd / name).mkdir(parents=True, exist_ok=True)
    return rd


def _center_to_edges(centers: np.ndarray, periodic: bool = False) -> np.ndarray:
    """Convert 1D cell centers into cell-edge coordinates.

    Args:
        centers: Monotone center coordinates.
        periodic: Whether the coordinate wraps periodically.

    Returns:
        Edge coordinates with length ``len(centers) + 1``.
    """

    centers = np.asarray(centers, dtype=float)
    if centers.size == 1:
        width = 1.0
        return np.array([centers[0] - 0.5 * width, centers[0] + 0.5 * width], dtype=float)

    delta = float(np.mean(np.diff(centers)))
    edges = np.empty(centers.size + 1, dtype=float)
    edges[:-1] = centers - 0.5 * delta
    edges[-1] = centers[-1] + 0.5 * delta
    if periodic:
        return edges
    return edges


def _polar_mesh(phi: np.ndarray, r: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Build Cartesian ``(x, y)`` edge grids for annular polar plots.

    Args:
        phi: Azimuthal cell centers.
        r: Radial cell centers.

    Returns:
        Tuple ``(x_edges, y_edges)`` with shape ``(N_r+1, N_phi+1)``.
    """

    phi_edges = _center_to_edges(phi, periodic=True)
    r_edges = _center_to_edges(r, periodic=False)
    rr, pp = np.meshgrid(r_edges, phi_edges, indexing="ij")
    x = rr * np.cos(pp)
    y = rr * np.sin(pp)
    return x, y


class RunWriter:
    """Write a complete artifact set for one ERD plant simulation.

    The writer keeps HDF5 files open during stepping and flushes/finalizes all
    outputs in :meth:`finalize`.
    """

    def __init__(self, cfg: RunConfig, run_dir: Optional[str] = None):
        """Initialize run directory, loggers, and HDF5 datasets.

        Args:
            cfg: Plant run configuration.
            run_dir: Optional explicit run directory path override.

        Returns:
            None.
        """

        self.cfg = cfg
        self.run_dir = make_run_dir(cfg.output.outputs_root, cfg.output.tag, run_dir=run_dir or cfg.output.run_dir)

        self._logger = logging.getLogger(f"erd_run_{self.run_dir.name}")
        self._logger.setLevel(logging.INFO)
        self._logger.handlers.clear()
        self._logger.propagate = False

        file_handler = logging.FileHandler(self.run_dir / "logs" / "run.log")
        file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        self._logger.addHandler(file_handler)

        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(logging.Formatter("%(message)s"))
        self._logger.addHandler(stream_handler)

        self._fields_h5 = h5py.File(self.run_dir / "fields" / "snapshots.h5", "w")
        self._ctrl_h5 = h5py.File(self.run_dir / "controls" / "control_timeseries.h5", "w")

        self._t_snap = self._fields_h5.create_dataset("t_snap", shape=(0,), maxshape=(None,), dtype="f8")
        self._n_ds = None
        self._omega_ds = None
        self._u_mag_ds = None
        self._u_r_ds = None
        self._u_phi_ds = None
        self._fields_h5.create_dataset("r", data=np.array([], dtype=float))
        self._fields_h5.create_dataset("phi", data=np.array([], dtype=float))

        self._t_ctrl = self._ctrl_h5.create_dataset("t_ctrl", shape=(0,), maxshape=(None,), dtype="f8")
        self._u_ds = self._ctrl_h5.create_dataset("u", shape=(0, 5), maxshape=(None, 5), dtype="f8")
        self._ctrl_h5.create_dataset("channel_names", data=np.array([b"u0", b"u1", b"u2", b"u3", b"u4"]))

        self.curves: Dict[str, List[float]] = {
            "t": [],
            "J_prof": [],
            "E_low": [],
            "L_w": [],
            "sigma_r": [],
            "r_mean": [],
            "M": [],
            "E_u": [],
            "P_ctrl": [],
        }

    @property
    def logger(self) -> logging.Logger:
        """Return the run-specific logger instance."""

        return self._logger

    def log(self, message: str) -> None:
        """Write one info message to ``logs/run.log`` and stdout.

        Args:
            message: Log line content.

        Returns:
            None.
        """

        self._logger.info(message)

    def write_grid(self, r: np.ndarray, phi: np.ndarray) -> None:
        """Persist radial and azimuthal coordinates in ``fields/snapshots.h5``.

        Args:
            r: Radial cell-center array.
            phi: Azimuthal cell-center array.

        Returns:
            None.
        """

        del self._fields_h5["r"]
        del self._fields_h5["phi"]
        self._fields_h5.create_dataset("r", data=np.asarray(r, dtype=float))
        self._fields_h5.create_dataset("phi", data=np.asarray(phi, dtype=float))

    def _ensure_field_datasets(
        self,
        shape: tuple[int, int],
        include_u_mag: bool = False,
        include_velocity_components: bool = False,
    ) -> None:
        """Create variable-sized field datasets on first snapshot append.

        Args:
            shape: Grid shape ``(N_r, N_phi)``.
            include_u_mag: If ``True``, ensure ``u_mag`` snapshot dataset exists.
            include_velocity_components: If ``True``, ensure ``u_r`` and
                ``u_phi`` snapshot datasets exist.

        Returns:
            None.
        """

        nr, nphi = shape
        if self._n_ds is None:
            self._n_ds = self._fields_h5.create_dataset(
                "n", shape=(0, nr, nphi), maxshape=(None, nr, nphi), dtype="f8"
            )
        if self._omega_ds is None:
            self._omega_ds = self._fields_h5.create_dataset(
                "omega", shape=(0, nr, nphi), maxshape=(None, nr, nphi), dtype="f8"
            )
        if include_u_mag and self._u_mag_ds is None:
            self._u_mag_ds = self._fields_h5.create_dataset(
                "u_mag", shape=(0, nr, nphi), maxshape=(None, nr, nphi), dtype="f8"
            )
        if include_velocity_components and self._u_r_ds is None:
            self._u_r_ds = self._fields_h5.create_dataset(
                "u_r", shape=(0, nr, nphi), maxshape=(None, nr, nphi), dtype="f8"
            )
        if include_velocity_components and self._u_phi_ds is None:
            self._u_phi_ds = self._fields_h5.create_dataset(
                "u_phi", shape=(0, nr, nphi), maxshape=(None, nr, nphi), dtype="f8"
            )

    @staticmethod
    def _append(ds: h5py.Dataset, arr: np.ndarray) -> None:
        """Append one row/frame to a resizable HDF5 dataset."""

        idx = ds.shape[0]
        ds.resize((idx + 1,) + ds.shape[1:])
        ds[idx] = arr

    def append_snapshot(
        self,
        t: float,
        n: np.ndarray,
        omega: np.ndarray,
        u_mag: np.ndarray | None = None,
        u_r: np.ndarray | None = None,
        u_phi: np.ndarray | None = None,
    ) -> None:
        """Append one field snapshot to ``fields/snapshots.h5``.

        Args:
            t: Snapshot timestamp.
            n: Density field array.
            omega: Vorticity field array.
            u_mag: Optional velocity-magnitude field for movie generation.
            u_r: Optional radial velocity field for diagnostics and movies.
            u_phi: Optional azimuthal velocity field for diagnostics and movies.

        Returns:
            None.
        """

        self._ensure_field_datasets(
            n.shape,
            include_u_mag=u_mag is not None,
            include_velocity_components=(u_r is not None) or (u_phi is not None),
        )
        self._append(self._t_snap, np.array(float(t)))
        self._append(self._n_ds, np.asarray(n, dtype=float))
        self._append(self._omega_ds, np.asarray(omega, dtype=float))
        if u_mag is not None and self._u_mag_ds is not None:
            self._append(self._u_mag_ds, np.asarray(u_mag, dtype=float))
        if u_r is not None and self._u_r_ds is not None:
            self._append(self._u_r_ds, np.asarray(u_r, dtype=float))
        if u_phi is not None and self._u_phi_ds is not None:
            self._append(self._u_phi_ds, np.asarray(u_phi, dtype=float))

    def append_control(self, t: float, u: np.ndarray) -> None:
        """Append one control sample to ``controls/control_timeseries.h5``.

        Args:
            t: Control timestamp.
            u: Control vector ``[u0, ..., u4]``.

        Returns:
            None.
        """

        self._append(self._t_ctrl, np.array(float(t)))
        self._append(self._u_ds, np.asarray(u, dtype=float).reshape(5))

    def append_metrics(self, t: float, metric: Dict[str, float]) -> None:
        """Append one metric sample to in-memory curve arrays.

        Args:
            t: Metric timestamp.
            metric: Per-step metric dictionary from :class:`ERDPlant`.

        Returns:
            None.
        """

        self.curves["t"].append(float(t))
        for k in self.curves:
            if k == "t":
                continue
            self.curves[k].append(float(metric[k]))

    def write_config(self) -> None:
        """Serialize ``config.yaml`` for this run."""

        save_run_config(self.run_dir / "config.yaml", self.cfg)

    def write_curves(self, curves: Dict[str, List[float]]) -> None:
        """Write metric time series to ``metrics/curves.json``.

        Args:
            curves: Metric curves dictionary.

        Returns:
            None.
        """

        with (self.run_dir / "metrics" / "curves.json").open("w", encoding="utf-8") as f:
            json.dump(curves, f, indent=2)

    def write_summary(self, summary: Dict[str, float]) -> None:
        """Write scalar summary metrics to ``summary.json``.

        Args:
            summary: Aggregate metric dictionary.

        Returns:
            None.
        """

        with (self.run_dir / "summary.json").open("w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)

    def _plot_metrics(self, curves: Dict[str, List[float]]) -> None:
        """Render multi-panel metric history figure."""

        t = np.asarray(curves["t"])
        fig, axes = plt.subplots(2, 3, figsize=(14, 7))
        axes = axes.reshape(-1)

        axes[0].plot(t, curves["J_prof"], color="C0")
        axes[0].set_title("J_prof")

        axes[1].plot(t, curves["E_low"], color="C1")
        axes[1].set_title("E_low")

        axes[2].plot(t, curves["L_w"], color="C2")
        axes[2].set_title("L_w")

        axes[3].plot(t, curves["sigma_r"], label="sigma_r")
        axes[3].plot(t, curves["r_mean"], label="r_mean", alpha=0.8)
        axes[3].set_title("Ring Thickness")
        axes[3].legend()

        axes[4].plot(t, curves["E_u"], label="E_u")
        axes[4].plot(t, curves["P_ctrl"], label="P_ctrl")
        if "eta_T" in curves:
            axes[4].plot(t, curves["eta_T"], label="eta_T")
        axes[4].set_title("Flow/Actuation")
        axes[4].legend()

        axes[5].plot(t, curves.get("L_w_cum", []), label="L_w_cum")
        axes[5].plot(t, curves.get("J_prof_excess", []), label="J_prof_excess")
        axes[5].plot(t, curves.get("E_low_excess", []), label="E_low_excess")
        axes[5].plot(t, curves.get("sigma_r_excess", []), label="sigma_r_excess")
        axes[5].plot(t, curves.get("badness_score", []), label="badness_score", alpha=0.7)
        axes[5].set_title("Degradation Terms")
        axes[5].legend()

        for ax in axes:
            ax.set_xlabel("t")
            ax.grid(alpha=0.3)

        fig.tight_layout()
        fig.savefig(self.run_dir / "plots" / "metrics_vs_time.png", dpi=150)
        plt.close(fig)

    def _plot_controls(self) -> None:
        """Render control-channel history figure."""

        t = np.asarray(self._t_ctrl)
        u = np.asarray(self._u_ds)
        fig = plt.figure(figsize=(10, 4.5))
        for j in range(5):
            plt.plot(t, u[:, j], label=f"u{j}")
        plt.grid(alpha=0.3)
        plt.xlabel("t")
        plt.ylabel("control")
        plt.title("Control Channels")
        plt.legend(ncol=5)
        plt.tight_layout()
        plt.savefig(self.run_dir / "plots" / "controls_vs_time.png", dpi=150)
        plt.close(fig)

    def _plot_spectrum(self) -> None:
        """Render final-snapshot azimuthal spectrum figure for ``n``."""

        n_last = np.asarray(self._n_ds[-1])
        spec = dominant_spectrum(n_last)

        fig = plt.figure(figsize=(7, 4))
        plt.semilogy(spec[: max(16, min(40, spec.shape[0]))], marker="o")
        plt.grid(alpha=0.3)
        plt.xlabel("Fourier mode m")
        plt.ylabel("radial-avg |n_hat_m|^2")
        plt.title("Azimuthal Spectrum (final snapshot)")
        plt.tight_layout()
        plt.savefig(self.run_dir / "plots" / "spectrum_n.png", dpi=150)
        plt.close(fig)

    def _write_field_movie(
        self,
        data: np.ndarray,
        phi: np.ndarray,
        r: np.ndarray,
        out_name: str,
        label: str,
        cmap: str,
    ) -> None:
        """Write one GIF movie from a stored field time-series.

        Args:
            data: Field snapshots with shape ``(T, N_r, N_phi)``.
            phi: Azimuthal coordinate vector.
            r: Radial coordinate vector.
            out_name: Output GIF filename.
            label: Colorbar label / title prefix.
            cmap: Matplotlib colormap name.

        Returns:
            None.
        """

        if data.size == 0:
            return

        vmin = float(np.percentile(data, 1))
        vmax = float(np.percentile(data, 99))
        if np.isclose(vmax, vmin):
            vmax = vmin + 1.0

        frames: List[np.ndarray] = []
        step = max(1, int(data.shape[0] / 180))

        for k in range(0, data.shape[0], step):
            fig = plt.figure(figsize=(6, 4.2))
            plt.imshow(
                data[k],
                origin="lower",
                aspect="auto",
                extent=[phi[0], phi[-1], r[0], r[-1]],
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
            )
            plt.colorbar(label=label)
            plt.xlabel("phi")
            plt.ylabel("r")
            plt.title(f"{label}(r,phi), frame {k}")
            plt.tight_layout()

            fig.canvas.draw()
            frame = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
            frames.append(frame)
            plt.close(fig)

        imageio.mimsave(self.run_dir / "movies" / out_name, frames, fps=15)

    def _write_polar_field_movie(
        self,
        data: np.ndarray,
        phi: np.ndarray,
        r: np.ndarray,
        out_name: str,
        label: str,
        cmap: str,
    ) -> None:
        """Write one GIF movie using an annular polar remap for visualization.

        Args:
            data: Field snapshots with shape ``(T, N_r, N_phi)``.
            phi: Azimuthal coordinate vector.
            r: Radial coordinate vector.
            out_name: Output GIF filename.
            label: Colorbar label / title prefix.
            cmap: Matplotlib colormap name.

        Returns:
            None.
        """

        if data.size == 0:
            return

        x, y = _polar_mesh(phi, r)
        vmax_abs = float(np.percentile(np.abs(data), 99))
        if "RdBu" in cmap:
            vmin = -max(vmax_abs, 1.0e-12)
            vmax = max(vmax_abs, 1.0e-12)
        else:
            vmin = float(np.percentile(data, 1))
            vmax = float(np.percentile(data, 99))
            if np.isclose(vmax, vmin):
                vmax = vmin + 1.0

        frames: List[np.ndarray] = []
        step = max(1, int(data.shape[0] / 180))
        lim = 1.05 * float(np.max(r))

        for k in range(0, data.shape[0], step):
            fig, ax = plt.subplots(figsize=(5.4, 5.4))
            mesh = ax.pcolormesh(x, y, data[k], shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
            fig.colorbar(mesh, ax=ax, label=label, shrink=0.85)
            ax.set_aspect("equal")
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_title(f"{label}(x,y), frame {k}")
            fig.tight_layout()

            fig.canvas.draw()
            frame = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
            frames.append(frame)
            plt.close(fig)

        imageio.mimsave(self.run_dir / "movies" / out_name, frames, fps=15)

    def _write_movie(self) -> None:
        """Write GIF movies from stored ``n``, ``omega``, and velocity fields."""

        n_data = np.asarray(self._n_ds)
        omega_data = np.asarray(self._omega_ds)
        phi = np.asarray(self._fields_h5["phi"])
        r = np.asarray(self._fields_h5["r"])
        self._write_field_movie(n_data, phi, r, out_name="n.gif", label="n", cmap="viridis")
        self._write_polar_field_movie(n_data, phi, r, out_name="n_polar.gif", label="n", cmap="viridis")
        self._write_field_movie(omega_data, phi, r, out_name="omega.gif", label="omega", cmap="RdBu_r")
        self._write_polar_field_movie(omega_data, phi, r, out_name="omega_polar.gif", label="omega", cmap="RdBu_r")
        if self._u_mag_ds is not None:
            u_mag = np.asarray(self._u_mag_ds)
            self._write_field_movie(u_mag, phi, r, out_name="u_mag.gif", label="|u|", cmap="magma")
            self._write_polar_field_movie(u_mag, phi, r, out_name="u_mag_polar.gif", label="|u|", cmap="magma")
        if self._u_r_ds is not None:
            u_r = np.asarray(self._u_r_ds)
            self._write_field_movie(u_r, phi, r, out_name="u_r.gif", label="u_r", cmap="RdBu_r")
            self._write_polar_field_movie(u_r, phi, r, out_name="u_r_polar.gif", label="u_r", cmap="RdBu_r")
        if self._u_phi_ds is not None:
            u_phi = np.asarray(self._u_phi_ds)
            self._write_field_movie(u_phi, phi, r, out_name="u_phi.gif", label="u_phi", cmap="RdBu_r")
            self._write_polar_field_movie(u_phi, phi, r, out_name="u_phi_polar.gif", label="u_phi", cmap="RdBu_r")

    def _write_contact_sheet(
        self,
        data: np.ndarray,
        phi: np.ndarray,
        r: np.ndarray,
        out_name: str,
        label: str,
        cmap: str,
    ) -> None:
        """Write a fixed-scale contact sheet for one field movie.

        Args:
            data: Field snapshots with shape ``(T, N_r, N_phi)``.
            phi: Azimuthal coordinate vector.
            r: Radial coordinate vector.
            out_name: Output PNG filename.
            label: Plot title prefix.
            cmap: Matplotlib colormap name.

        Returns:
            None.
        """

        if data.size == 0:
            return

        n_frames = data.shape[0]
        picks = np.unique(np.clip(np.round(np.linspace(0, n_frames - 1, 8)).astype(int), 0, n_frames - 1))
        vmin = float(np.percentile(data, 1))
        vmax = float(np.percentile(data, 99))
        if np.isclose(vmax, vmin):
            vmax = vmin + 1.0

        fig, axes = plt.subplots(4, 2, figsize=(10, 12))
        axes = axes.reshape(-1)
        for ax, idx in zip(axes, picks, strict=False):
            im = ax.imshow(
                data[idx],
                origin="lower",
                aspect="auto",
                extent=[phi[0], phi[-1], r[0], r[-1]],
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
            )
            ax.set_title(f"{label}, t={float(self._t_snap[idx]):.3f}")
            ax.set_xlabel("phi")
            ax.set_ylabel("r")
        for ax in axes[len(picks) :]:
            ax.axis("off")

        fig.colorbar(im, ax=axes.tolist(), shrink=0.95, label=label)
        fig.tight_layout()
        fig.savefig(self.run_dir / "plots" / out_name, dpi=150)
        plt.close(fig)

    def _write_polar_contact_sheet(
        self,
        data: np.ndarray,
        phi: np.ndarray,
        r: np.ndarray,
        out_name: str,
        label: str,
        cmap: str,
    ) -> None:
        """Write a polar-mapped contact sheet for one field.

        Args:
            data: Field snapshots with shape ``(T, N_r, N_phi)``.
            phi: Azimuthal coordinate vector.
            r: Radial coordinate vector.
            out_name: Output PNG filename.
            label: Plot title prefix.
            cmap: Matplotlib colormap name.

        Returns:
            None.
        """

        if data.size == 0:
            return

        n_frames = data.shape[0]
        picks = np.unique(np.clip(np.round(np.linspace(0, n_frames - 1, 8)).astype(int), 0, n_frames - 1))
        x, y = _polar_mesh(phi, r)
        vmax_abs = float(np.percentile(np.abs(data), 99))
        if "RdBu" in cmap:
            vmin = -max(vmax_abs, 1.0e-12)
            vmax = max(vmax_abs, 1.0e-12)
        else:
            vmin = float(np.percentile(data, 1))
            vmax = float(np.percentile(data, 99))
            if np.isclose(vmax, vmin):
                vmax = vmin + 1.0

        lim = 1.05 * float(np.max(r))
        fig, axes = plt.subplots(4, 2, figsize=(10, 12))
        axes = axes.reshape(-1)
        for ax, idx in zip(axes, picks, strict=False):
            im = ax.pcolormesh(x, y, data[idx], shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_aspect("equal")
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_title(f"{label}, t={float(self._t_snap[idx]):.3f}")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
        for ax in axes[len(picks) :]:
            ax.axis("off")

        fig.colorbar(im, ax=axes.tolist(), shrink=0.95, label=label)
        fig.tight_layout()
        fig.savefig(self.run_dir / "plots" / out_name, dpi=150)
        plt.close(fig)

    def _write_contact_sheets(self) -> None:
        """Write contact sheets for density, vorticity, and velocity fields."""

        phi = np.asarray(self._fields_h5["phi"])
        r = np.asarray(self._fields_h5["r"])
        self._write_contact_sheet(np.asarray(self._n_ds), phi, r, "n_contact_sheet.png", "n", "viridis")
        self._write_polar_contact_sheet(
            np.asarray(self._n_ds),
            phi,
            r,
            "n_polar_contact_sheet.png",
            "n",
            "viridis",
        )
        self._write_contact_sheet(
            np.asarray(self._omega_ds),
            phi,
            r,
            "omega_contact_sheet.png",
            "omega",
            "RdBu_r",
        )
        self._write_polar_contact_sheet(
            np.asarray(self._omega_ds),
            phi,
            r,
            "omega_polar_contact_sheet.png",
            "omega",
            "RdBu_r",
        )
        if self._u_mag_ds is not None:
            self._write_contact_sheet(
                np.asarray(self._u_mag_ds),
                phi,
                r,
                "u_mag_contact_sheet.png",
                "|u|",
                "magma",
            )
            self._write_polar_contact_sheet(
                np.asarray(self._u_mag_ds),
                phi,
                r,
                "u_mag_polar_contact_sheet.png",
                "|u|",
                "magma",
            )
        if self._u_r_ds is not None:
            self._write_contact_sheet(
                np.asarray(self._u_r_ds),
                phi,
                r,
                "u_r_contact_sheet.png",
                "u_r",
                "RdBu_r",
            )
            self._write_polar_contact_sheet(
                np.asarray(self._u_r_ds),
                phi,
                r,
                "u_r_polar_contact_sheet.png",
                "u_r",
                "RdBu_r",
            )
        if self._u_phi_ds is not None:
            self._write_contact_sheet(
                np.asarray(self._u_phi_ds),
                phi,
                r,
                "u_phi_contact_sheet.png",
                "u_phi",
                "RdBu_r",
            )
            self._write_polar_contact_sheet(
                np.asarray(self._u_phi_ds),
                phi,
                r,
                "u_phi_polar_contact_sheet.png",
                "u_phi",
                "RdBu_r",
            )

    def _write_transport_diagnostics(self) -> None:
        """Compute and render transport diagnostics aimed at visual tuning."""

        n_data = np.asarray(self._n_ds)
        omega_data = np.asarray(self._omega_ds)
        u_mag = None if self._u_mag_ds is None else np.asarray(self._u_mag_ds)
        u_r = None if self._u_r_ds is None else np.asarray(self._u_r_ds)
        u_phi = None if self._u_phi_ds is None else np.asarray(self._u_phi_ds)
        r = np.asarray(self._fields_h5["r"])
        t_snap = np.asarray(self._t_snap)

        diag = compute_transport_diagnostics(
            cfg=self.cfg,
            n_snapshots=n_data,
            omega_snapshots=omega_data,
            u_mag_snapshots=u_mag,
            u_r_snapshots=u_r,
            u_phi_snapshots=u_phi,
            r=r,
            t_snap=t_snap,
        )
        with (self.run_dir / "metrics" / "transport_diagnostics.json").open("w", encoding="utf-8") as f:
            json.dump(diag, f, indent=2)

        t = np.asarray(diag.get("t", []), dtype=float)
        td = np.asarray(diag.get("delta_t", []), dtype=float)
        delta_n = np.asarray(diag.get("delta_n_rel", []), dtype=float)
        delta_omega = np.asarray(diag.get("delta_omega_rel", []), dtype=float)
        u_ring_mean = np.asarray(diag.get("u_ring_mean", []), dtype=float)
        u_ring_p95 = np.asarray(diag.get("u_ring_p95", []), dtype=float)
        u_r_ring_abs_mean = np.asarray(diag.get("u_r_ring_abs_mean", []), dtype=float)
        gamma_r_ring_signed_mean = np.asarray(diag.get("gamma_r_ring_signed_mean", []), dtype=float)
        gamma_r_ring_abs_q90 = np.asarray(diag.get("gamma_r_ring_abs_q90", []), dtype=float)
        u_phi_pos_frac = np.asarray(diag.get("u_phi_pos_frac", []), dtype=float)
        u_phi_ring_signed = np.asarray(diag.get("u_phi_ring_signed", []), dtype=float)
        u_phi_ring_abs_mean = np.asarray(diag.get("u_phi_ring_abs_mean", []), dtype=float)
        u_phi_over_u_r = np.asarray(diag.get("u_phi_over_u_r", []), dtype=float)
        m1_amp = np.asarray(diag.get("m1_amp", []), dtype=float)
        m1_phase = np.asarray(diag.get("m1_phase", []), dtype=float)
        m1_phase_rate_t = np.asarray(diag.get("m1_phase_rate_t", []), dtype=float)
        m1_phase_rate = np.asarray(diag.get("m1_phase_rate", []), dtype=float)
        u_ring_accel_t = np.asarray(diag.get("u_ring_accel_t", []), dtype=float)
        u_ring_accel = np.asarray(diag.get("u_ring_accel", []), dtype=float)
        r_mean = np.asarray(diag.get("r_mean", []), dtype=float)
        ring_grad_ratio = np.asarray(diag.get("ring_grad_ratio", []), dtype=float)
        ring_grad_abs = np.asarray(diag.get("ring_grad_abs", []), dtype=float)
        ring_shear_norm = np.asarray(diag.get("ring_shear_norm", []), dtype=float)
        ring_instability_margin = np.asarray(diag.get("ring_instability_margin", []), dtype=float)
        ring_instability_activation = np.asarray(diag.get("ring_instability_activation", []), dtype=float)
        ring_phase_drive_signed = np.asarray(diag.get("ring_phase_drive_signed", []), dtype=float)
        ring_phase_drive_abs = np.asarray(diag.get("ring_phase_drive_abs", []), dtype=float)
        ring_wave_phi_speed_signed = np.asarray(diag.get("ring_wave_phi_speed_signed", []), dtype=float)
        ring_wave_phi_speed_abs = np.asarray(diag.get("ring_wave_phi_speed_abs", []), dtype=float)
        ring_wave_r_speed_signed = np.asarray(diag.get("ring_wave_r_speed_signed", []), dtype=float)
        ring_wave_r_speed_abs = np.asarray(diag.get("ring_wave_r_speed_abs", []), dtype=float)
        ring_wave_packet_speed = np.asarray(diag.get("ring_wave_packet_speed", []), dtype=float)
        ring_wave_burst_gain = np.asarray(diag.get("ring_wave_burst_gain", []), dtype=float)
        ring_confinement_signed = np.asarray(diag.get("ring_confinement_signed", []), dtype=float)
        ring_confinement_abs = np.asarray(diag.get("ring_confinement_abs", []), dtype=float)
        ring_confinement_release = np.asarray(diag.get("ring_confinement_release", []), dtype=float)
        ring_zonal_restore_signed = np.asarray(diag.get("ring_zonal_restore_signed", []), dtype=float)
        ring_zonal_restore_abs = np.asarray(diag.get("ring_zonal_restore_abs", []), dtype=float)
        ring_zonal_restore_release = np.asarray(diag.get("ring_zonal_restore_release", []), dtype=float)
        ring_coupling_boost = np.asarray(diag.get("ring_coupling_boost", []), dtype=float)
        ring_feedback_gain = np.asarray(diag.get("ring_feedback_gain", []), dtype=float)
        ring_transport_boost = np.asarray(diag.get("ring_transport_boost", []), dtype=float)
        inner_edge_mass_frac = np.asarray(diag.get("inner_edge_mass_frac", []), dtype=float)
        outer_edge_mass_frac = np.asarray(diag.get("outer_edge_mass_frac", []), dtype=float)
        ring_core_mass_frac = np.asarray(diag.get("ring_core_mass_frac", []), dtype=float)
        supercritical_fraction = np.asarray(diag.get("supercritical_fraction", []), dtype=float)
        summary = diag.get("summary", {})

        fig, axes = plt.subplots(3, 2, figsize=(12, 10))
        axes = axes.reshape(-1)

        axes[0].plot(td, delta_n, label="n Δ")
        axes[0].plot(td, delta_omega, label="omega Δ")
        axes[0].set_title("Frame-to-Frame Relative Change")
        axes[0].legend()

        axes[1].plot(t, u_ring_mean, label="ring mean |u|")
        axes[1].plot(t, u_ring_p95, label="ring p95 |u|")
        axes[1].plot(t, u_r_ring_abs_mean, label="ring mean |u_r|")
        axes[1].plot(t, u_phi_ring_abs_mean, label="ring mean |u_phi|")
        axes[1].plot(t, u_phi_ring_signed, label="ring mean u_phi", alpha=0.7)
        axes[1].set_title("Ring-Band Velocity Split")
        axes[1].legend()

        axes[2].plot(t, m1_amp, label="|m1|")
        axes[2].set_title("Ring-Band m=1 Amplitude")

        axes[3].plot(t, m1_phase, label="phase")
        if gamma_r_ring_abs_q90.size == t.size:
            axes[3].plot(t, gamma_r_ring_abs_q90, label="q90 |n u_r|")
        if gamma_r_ring_signed_mean.size == t.size:
            axes[3].plot(t, gamma_r_ring_signed_mean, label="mean n u_r", alpha=0.7)
        axes[3].set_title("Phase and Radial Flux Bursts")
        axes[3].legend()

        axes[4].plot(m1_phase_rate_t, m1_phase_rate, label="dphase/dt")
        axes[4].plot(t, u_phi_over_u_r, label="|u_phi| / |u_r|")
        if u_phi_pos_frac.size == t.size:
            axes[4].plot(t, u_phi_pos_frac, label="u_phi > 0 frac", alpha=0.7)
        axes[4].set_title("Phase Rate and Transport Ratio")
        axes[4].legend()

        if ring_instability_margin.size == t.size:
            axes[5].plot(t, ring_grad_ratio, label="grad ratio")
            axes[5].plot(t, ring_grad_abs, label="|grad|", alpha=0.7)
            axes[5].plot(t, ring_shear_norm, label="shear norm")
            axes[5].plot(t, ring_instability_margin, label="margin")
            axes[5].plot(t, ring_instability_activation, label="activation")
            axes[5].plot(t, ring_phase_drive_abs, label="|phase drive|", alpha=0.8)
            axes[5].plot(t, ring_phase_drive_signed, label="phase drive", alpha=0.6)
            axes[5].plot(t, ring_wave_phi_speed_abs, label="|wave c_phi|", alpha=0.8)
            axes[5].plot(t, ring_wave_phi_speed_signed, label="wave c_phi", alpha=0.6)
            axes[5].plot(t, ring_wave_r_speed_abs, label="|wave c_r|", alpha=0.8)
            axes[5].plot(t, ring_wave_packet_speed, label="|wave c|", alpha=0.7)
            axes[5].plot(t, ring_wave_burst_gain, label="wave burst", alpha=0.7)
            axes[5].plot(t, ring_confinement_abs, label="|confine|", alpha=0.7)
            axes[5].plot(t, ring_confinement_release, label="confine release", alpha=0.7)
            axes[5].plot(t, ring_zonal_restore_abs, label="|zonal restore|", alpha=0.8)
            axes[5].plot(t, ring_zonal_restore_release, label="restore release", alpha=0.8)
            axes[5].plot(t, ring_coupling_boost, label="coupling boost", alpha=0.8)
            axes[5].plot(t, ring_feedback_gain, label="feedback gain", alpha=0.8)
            axes[5].plot(t, ring_transport_boost, label="transport boost", alpha=0.8)
            axes[5].plot(t, supercritical_fraction, label="supercritical frac", alpha=0.8)
            axes[5].plot(t, inner_edge_mass_frac, label="inner edge mass", alpha=0.7)
            axes[5].plot(t, outer_edge_mass_frac, label="outer edge mass", alpha=0.7)
            axes[5].plot(t, ring_core_mass_frac, label="ring core mass", alpha=0.7)
            if u_ring_accel.size == u_ring_accel_t.size:
                axes[5].plot(u_ring_accel_t, u_ring_accel, label="d|u|/dt", alpha=0.7)
            axes[5].set_title("Threshold and Burst State")
            axes[5].legend()
        else:
            axes[5].axis("off")

        fig2, ax_text = plt.subplots(figsize=(8.5, 8.5))
        ax_text.axis("off")
        lines = [
            f"delta_n_median: {summary.get('delta_n_rel_median', 0.0):.4e}",
            f"delta_omega_median: {summary.get('delta_omega_rel_median', 0.0):.4e}",
            f"u_ring_mean_peak: {summary.get('u_ring_mean_peak', 0.0):.4e}",
            f"u_ring_p95_peak: {summary.get('u_ring_p95_peak', 0.0):.4e}",
            f"u_r_ring_abs_mean_peak: {summary.get('u_r_ring_abs_mean_peak', 0.0):.4e}",
            f"gamma_r_ring_abs_q90_peak: {summary.get('gamma_r_ring_abs_q90_peak', 0.0):.4e}",
            f"gamma_r_ring_signed_peak: {summary.get('gamma_r_ring_signed_peak', 0.0):.4e}",
            f"gamma_r_burst_factor: {summary.get('gamma_r_burst_factor', 0.0):.4e}",
            f"gamma_r_burst_count: {int(summary.get('gamma_r_burst_count', 0))}",
            f"u_phi_ring_abs_mean_peak: {summary.get('u_phi_ring_abs_mean_peak', 0.0):.4e}",
            f"u_phi_ring_signed_peak: {summary.get('u_phi_ring_signed_peak', 0.0):.4e}",
            f"u_phi_over_u_r_mean_ratio: {summary.get('u_phi_over_u_r_mean_ratio', 0.0):.4e}",
            f"u_phi_over_u_r_peak_ratio: {summary.get('u_phi_over_u_r_peak_ratio', 0.0):.4e}",
            f"u_phi_sign_changes: {int(summary.get('u_phi_sign_changes', 0))}",
            f"u_phi_signed_bias: {summary.get('u_phi_signed_bias', 0.0):.4e}",
            f"u_phi_pos_frac_mean: {summary.get('u_phi_pos_frac_mean', 0.0):.4e}",
            f"u_phi_sign_mix_peak: {summary.get('u_phi_sign_mix_peak', 0.0):.4e}",
            f"u_phi_sign_balance_mean: {summary.get('u_phi_sign_balance_mean', 0.0):.4e}",
            f"m1_amp_peak: {summary.get('m1_amp_peak', 0.0):.4e}",
            f"m1_amp_final: {summary.get('m1_amp_final', 0.0):.4e}",
            f"m1_phase_drift_total: {summary.get('m1_phase_drift_total', 0.0):.4e}",
            f"m1_phase_rate_rms: {summary.get('m1_phase_rate_rms', 0.0):.4e}",
            f"u_ring_accel_peak: {summary.get('u_ring_accel_peak', 0.0):.4e}",
            f"u_ring_accel_rms: {summary.get('u_ring_accel_rms', 0.0):.4e}",
            f"r_mean_shift: {summary.get('r_mean_shift', 0.0):.4e}",
            f"mass_rel_drift: {summary.get('mass_rel_drift', 0.0):.4e}",
            f"mass_rel_span: {summary.get('mass_rel_span', 0.0):.4e}",
            f"n_tilde_rms_peak: {summary.get('n_tilde_rms_peak', 0.0):.4e}",
            f"n_tilde_rms_final: {summary.get('n_tilde_rms_final', 0.0):.4e}",
            f"n_tilde_retention_final_over_peak: {summary.get('n_tilde_retention_final_over_peak', 0.0):.4e}",
            f"row_abs_corr_mean: {summary.get('row_abs_corr_mean', 0.0):.4e}",
            f"row_adjacent_corr_mean: {summary.get('row_adjacent_corr_mean', 0.0):.4e}",
            f"delta_n_growth_ratio: {summary.get('delta_n_growth_ratio', 0.0):.4e}",
            f"m1_growth_factor: {summary.get('m1_growth_factor', 0.0):.4e}",
            f"ring_grad_ratio_peak: {summary.get('ring_grad_ratio_peak', 0.0):.4e}",
            f"ring_grad_abs_peak: {summary.get('ring_grad_abs_peak', 0.0):.4e}",
            f"ring_shear_norm_peak: {summary.get('ring_shear_norm_peak', 0.0):.4e}",
            f"ring_margin_peak: {summary.get('ring_instability_margin_peak', 0.0):.4e}",
            f"ring_activation_peak: {summary.get('ring_instability_activation_peak', 0.0):.4e}",
            f"ring_phase_drive_abs_peak: {summary.get('ring_phase_drive_abs_peak', 0.0):.4e}",
            f"ring_wave_phi_speed_abs_peak: {summary.get('ring_wave_phi_speed_abs_peak', 0.0):.4e}",
            f"ring_wave_r_speed_abs_peak: {summary.get('ring_wave_r_speed_abs_peak', 0.0):.4e}",
            f"ring_wave_packet_speed_peak: {summary.get('ring_wave_packet_speed_peak', 0.0):.4e}",
            f"ring_wave_burst_gain_peak: {summary.get('ring_wave_burst_gain_peak', 0.0):.4e}",
            f"wave_phi_sign_changes: {int(summary.get('wave_phi_sign_changes', 0))}",
            f"ring_confinement_abs_peak: {summary.get('ring_confinement_abs_peak', 0.0):.4e}",
            f"ring_confinement_release_final: {summary.get('ring_confinement_release_final', 0.0):.4e}",
            f"ring_zonal_restore_abs_peak: {summary.get('ring_zonal_restore_abs_peak', 0.0):.4e}",
            f"ring_zonal_restore_release_final: {summary.get('ring_zonal_restore_release_final', 0.0):.4e}",
            f"ring_coupling_boost_peak: {summary.get('ring_coupling_boost_peak', 0.0):.4e}",
            f"ring_feedback_gain_peak: {summary.get('ring_feedback_gain_peak', 0.0):.4e}",
            f"ring_transport_boost_peak: {summary.get('ring_transport_boost_peak', 0.0):.4e}",
            f"inner_edge_mass_frac_peak: {summary.get('inner_edge_mass_frac_peak', 0.0):.4e}",
            f"inner_edge_mass_frac_final: {summary.get('inner_edge_mass_frac_final', 0.0):.4e}",
            f"outer_edge_mass_frac_peak: {summary.get('outer_edge_mass_frac_peak', 0.0):.4e}",
            f"outer_edge_mass_frac_final: {summary.get('outer_edge_mass_frac_final', 0.0):.4e}",
            f"ring_core_mass_frac_peak: {summary.get('ring_core_mass_frac_peak', 0.0):.4e}",
            f"ring_core_mass_frac_final: {summary.get('ring_core_mass_frac_final', 0.0):.4e}",
            f"J_peak_time_fraction: {summary.get('J_peak_time_fraction', 0.0):.4e}",
            f"E_peak_time_fraction: {summary.get('E_peak_time_fraction', 0.0):.4e}",
            f"supercritical_fraction_peak: {summary.get('supercritical_fraction_peak', 0.0):.4e}",
            f"instability_burst_count: {int(summary.get('instability_burst_count', 0))}",
        ]
        ax_text.text(0.0, 1.0, "\n".join(lines), va="top", family="monospace")

        for ax in axes:
            ax.grid(alpha=0.3)
            ax.set_xlabel("t")

        fig.tight_layout()
        fig.savefig(self.run_dir / "plots" / "transport_diagnostics.png", dpi=150)
        plt.close(fig)
        fig2.tight_layout()
        fig2.savefig(self.run_dir / "plots" / "transport_diagnostics_summary.png", dpi=150)
        plt.close(fig2)

    def _write_variation_note(self) -> None:
        """Write a brief markdown note describing variation-driving settings.

        Args:
            None.

        Returns:
            None.
        """

        cfg = self.cfg
        lines = [
            "# Variation Note",
            "",
            "This run uses an annular modified Hasegawa-Wakatani core with multiscale disturbance forcing.",
            "",
            "## Key Levers",
            f"- time: dt={cfg.time.dt}, T_final={cfg.time.T_final}, snapshot_every={cfg.time.snapshot_every}",
            f"- pde: dynamics_model={cfg.pde.dynamics_model}, C={cfg.pde.C}, kappa_0={cfg.pde.kappa_0}, sigma_kappa={cfg.pde.sigma_kappa}, gradient_drive_gain={cfg.pde.gradient_drive_gain}, critical_gradient_ratio={cfg.pde.critical_gradient_ratio}, shear_suppression_gain={cfg.pde.shear_suppression_gain}, shear_ref={cfg.pde.shear_ref}, threshold_width={cfg.pde.threshold_width}, shear_damping_gain={cfg.pde.shear_damping_gain}, curvature_omega_gain={cfg.pde.curvature_omega_gain}, baroclinic_omega_gain={cfg.pde.baroclinic_omega_gain}, flux_balance_omega_gain={cfg.pde.flux_balance_omega_gain}, phase_advection_gain={cfg.pde.phase_advection_gain}, supercritical_coupling_gain={cfg.pde.supercritical_coupling_gain}, supercritical_feedback_gain={cfg.pde.supercritical_feedback_gain}, supercritical_transport_gain={cfg.pde.supercritical_transport_gain}, wave_density_coupling_gain={cfg.pde.wave_density_coupling_gain}, wave_vorticity_coupling_gain={cfg.pde.wave_vorticity_coupling_gain}, wave_packet_phi_gain={cfg.pde.wave_packet_phi_gain}, wave_packet_r_gain={cfg.pde.wave_packet_r_gain}, wave_burst_speed_gain={cfg.pde.wave_burst_speed_gain}, packet_mode_cut={cfg.pde.packet_mode_cut}, omega_high_damp={cfg.pde.omega_high_damp}, omega_landau_gain={cfg.pde.omega_landau_gain}, density_landau_gain={cfg.pde.density_landau_gain}, density_vortex_gain={cfg.pde.density_vortex_gain}, density_packet_drive_gain={cfg.pde.density_packet_drive_gain}, density_retention_gain={cfg.pde.density_retention_gain}, collapse_guard_gain={cfg.pde.collapse_guard_gain}, collapse_edge_threshold={cfg.pde.collapse_edge_threshold}, collapse_core_threshold={cfg.pde.collapse_core_threshold}, collapse_omega_damp={cfg.pde.collapse_omega_damp}, radial_confinement_gain={cfg.pde.radial_confinement_gain}, radial_release_gain={cfg.pde.radial_release_gain}, radial_damage_release_gain={cfg.pde.radial_damage_release_gain}, zonal_profile_restore_gain={cfg.pde.zonal_profile_restore_gain}, zonal_profile_release_gain={cfg.pde.zonal_profile_release_gain}, zonal_profile_damage_release_gain={cfg.pde.zonal_profile_damage_release_gain}, landau_phi_gain={cfg.pde.landau_phi_gain}, S_0={cfg.pde.S_0}, sigma_S={cfg.pde.sigma_S}",
            f"- dissipation: nu_n={cfg.pde.nu_n}, nu_omega={cfg.pde.nu_omega}, hyper_p={cfg.pde.hyper_p}, gamma_omega={cfg.pde.gamma_omega}",
            f"- source feedback: floor={cfg.pde.source_floor_scale}, balance_gain={cfg.pde.source_balance_gain}, mass_gain={cfg.pde.source_mass_gain}, max_scale={cfg.pde.source_scale_max}",
            f"- forcing: drive_u0_base={cfg.forcing.drive_u0_base}, u_bounds={list(cfg.forcing.u_bounds)}",
            f"- disturbance scales: warmup={cfg.disturbance.warmup_scale}, excite={cfg.disturbance.excite_scale}, hold={cfg.disturbance.hold_scale}",
            f"- multiscale modes: {list(cfg.disturbance.multiscale_modes)}",
            f"- pseudo-random band forcing: m=[{cfg.disturbance.noise_m_min}, {cfg.disturbance.noise_m_max}], strength={cfg.disturbance.noise_strength}, phase_drift={cfg.disturbance.phase_drift_strength}",
            f"- eddy packets: count={cfg.disturbance.eddy_count}, strength={cfg.disturbance.eddy_strength}, sigma_r={cfg.disturbance.eddy_sigma_r}, sigma_phi={cfg.disturbance.eddy_sigma_phi}",
            "",
            "Expected effect: a near-threshold ring where zonal gradients feed drift activity, zonal shear suppresses breakup, and stochastic perturbations trigger bursts only after the instability margin turns positive.",
        ]
        (self.run_dir / "notes" / "variation_note.md").write_text("\n".join(lines), encoding="utf-8")

    def _write_variation_diagnostics(self) -> None:
        """Compute and persist variation diagnostics from saved snapshots.

        The method writes:
        - ``metrics/variation_diagnostics.json``
        - ``plots/variation_diagnostics.png``

        It optionally compares temporal activity against a baseline run when
        ``ERD_VARIATION_BASELINE_RUN`` points to another run directory. If no
        external baseline is supplied, the configured
        ``metrics.temporal_activity_baseline`` is used so single-run demo
        executions can still satisfy the acceptance contract.

        Args:
            None.

        Returns:
            None.
        """

        n_data = np.asarray(self._n_ds)
        r = np.asarray(self._fields_h5["r"])
        t_snap = np.asarray(self._t_snap)
        diag = compute_variation_diagnostics(cfg=self.cfg, n_snapshots=n_data, r=r, t_snap=t_snap)

        baseline_run = os.environ.get("ERD_VARIATION_BASELINE_RUN", "").strip()
        baseline_info: Dict[str, float | str] = {}
        if baseline_run:
            try:
                with h5py.File(Path(baseline_run) / "fields" / "snapshots.h5", "r") as hf_base:
                    n_base = np.asarray(hf_base["n"])
                    r_base = np.asarray(hf_base["r"])
                    t_base = np.asarray(hf_base["t_snap"])
                base_diag = compute_variation_diagnostics(
                    cfg=self.cfg,
                    n_snapshots=n_base,
                    r=r_base,
                    t_snap=t_base,
                )
                base_med = float(base_diag["summary"].get("delta_frame_l2_rel_median", 0.0))
                cur_med = float(diag["summary"].get("delta_frame_l2_rel_median", 0.0))
                factor = cur_med / max(base_med, 1e-14)
                baseline_info = {
                    "baseline_run_dir": baseline_run,
                    "baseline_delta_median": base_med,
                    "delta_median_vs_baseline_factor": factor,
                }
                diag["summary"]["delta_median_vs_baseline_factor"] = factor
            except Exception as exc:
                baseline_info = {
                    "baseline_run_dir": baseline_run,
                    "baseline_error": str(exc),
                }
        else:
            base_med = float(self.cfg.metrics.temporal_activity_baseline)
            cur_med = float(diag["summary"].get("delta_frame_l2_rel_median", 0.0))
            factor = cur_med / max(base_med, 1e-14)
            baseline_info = {
                "baseline_source": "config.metrics.temporal_activity_baseline",
                "baseline_delta_median": base_med,
                "delta_median_vs_baseline_factor": factor,
            }
            diag["summary"]["delta_median_vs_baseline_factor"] = factor

        s = diag.get("summary", {})
        acceptance = {
            "criterion1_nonaxis_growth": bool(
                s.get("ratio_nonax_over_axisym_first20_max", 0.0) >= 0.2
                and s.get("ratio_nonax_over_axisym_first20_std", 0.0) > 1e-4
            ),
            "criterion2_mid_band": bool(s.get("band_ratio_mid_to_low_max", 0.0) >= 0.05),
            "criterion2_high_band": bool(s.get("band_ratio_high_to_low_max", 0.0) >= 0.05),
            "criterion3_temporal_activity": bool(s.get("delta_median_vs_baseline_factor", 0.0) >= 5.0),
        }
        acceptance["all_passed"] = bool(all(acceptance.values()))

        payload: Dict[str, object] = {
            "diagnostics": diag,
            "baseline": baseline_info,
            "acceptance": acceptance,
        }
        with (self.run_dir / "metrics" / "variation_diagnostics.json").open("w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)

        t = np.asarray(diag.get("t", []), dtype=float)
        ratio = np.asarray(diag.get("ratio_nonax_over_axisym", []), dtype=float)
        mid = np.asarray(diag.get("band_ratio_mid_to_low", []), dtype=float)
        high = np.asarray(diag.get("band_ratio_high_to_low", []), dtype=float)
        td = np.asarray(diag.get("delta_frame_l2_rel_t", []), dtype=float)
        delta = np.asarray(diag.get("delta_frame_l2_rel", []), dtype=float)

        fig, axes = plt.subplots(2, 2, figsize=(11, 7))
        axes = axes.reshape(-1)

        axes[0].plot(t, ratio, label="E_nonax / E0")
        axes[0].axhline(0.2, color="r", linestyle="--", label="target 0.2")
        axes[0].set_title("Non-axisymmetry Ratio")
        axes[0].legend()

        axes[1].plot(t, mid, label="mid / low")
        axes[1].plot(t, high, label="high / low")
        axes[1].axhline(0.05, color="r", linestyle="--", label="target 0.05")
        axes[1].set_title("Spectral Band Ratios")
        axes[1].legend()

        axes[2].plot(td, delta, label="frame Δ")
        if "baseline_delta_median" in baseline_info:
            baseline_target = 5.0 * float(baseline_info["baseline_delta_median"])
            axes[2].axhline(baseline_target, color="r", linestyle="--", label="5x baseline")
        axes[2].set_title("Temporal Activity")
        axes[2].legend()

        axes[3].axis("off")
        lines = [
            f"criterion1_nonaxis_growth: {acceptance['criterion1_nonaxis_growth']}",
            f"criterion2_mid_band: {acceptance['criterion2_mid_band']}",
            f"criterion2_high_band: {acceptance['criterion2_high_band']}",
            f"criterion3_temporal_activity: {acceptance['criterion3_temporal_activity']}",
            f"all_passed: {acceptance['all_passed']}",
        ]
        if "delta_median_vs_baseline_factor" in s:
            lines.append(f"delta_vs_baseline_factor: {s['delta_median_vs_baseline_factor']:.3f}")
        axes[3].text(0.0, 1.0, "\n".join(lines), va="top", family="monospace")

        for ax in axes[:3]:
            ax.set_xlabel("t")
            ax.grid(alpha=0.3)

        fig.tight_layout()
        fig.savefig(self.run_dir / "plots" / "variation_diagnostics.png", dpi=150)
        plt.close(fig)

    def finalize(self, curves: Dict[str, List[float]], summary: Dict[str, float]) -> None:
        """Finalize run artifacts and close all open resources.

        Args:
            curves: Full metric curves dictionary.
            summary: Scalar summary metrics.

        Returns:
            None.
        """

        self.write_config()
        self.write_curves(curves)
        self.write_summary(summary)

        skip_media = os.environ.get("ERD_SKIP_MEDIA", "").strip() == "1"
        if not skip_media:
            self._plot_metrics(curves)
            self._plot_controls()
            self._plot_spectrum()
            self._write_movie()
        else:
            self.log("Skipping plot/movie rendering because ERD_SKIP_MEDIA=1")

        self._write_variation_note()
        self._write_variation_diagnostics()
        self._write_transport_diagnostics()
        self._write_contact_sheets()

        self._fields_h5.flush()
        self._ctrl_h5.flush()
        self._fields_h5.close()
        self._ctrl_h5.close()

        self.log(f"Run artifacts written: {self.run_dir}")
