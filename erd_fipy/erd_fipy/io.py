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
from .metrics import compute_variation_diagnostics, dominant_spectrum


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
        self._fields_h5.create_dataset("r", data=np.array([], dtype=float))
        self._fields_h5.create_dataset("phi", data=np.array([], dtype=float))

        self._t_ctrl = self._ctrl_h5.create_dataset("t_ctrl", shape=(0,), maxshape=(None,), dtype="f8")
        self._u_ds = self._ctrl_h5.create_dataset("u", shape=(0, 5), maxshape=(None, 5), dtype="f8")
        self._ctrl_h5.create_dataset("channel_names", data=np.array([b"u0", b"u1", b"u2", b"u3", b"u4"]))

        self.curves: Dict[str, List[float]] = {
            "t": [],
            "E_wob": [],
            "E1": [],
            "E2": [],
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

    def _ensure_field_datasets(self, shape: tuple[int, int], include_u_mag: bool = False) -> None:
        """Create variable-sized field datasets on first snapshot append.

        Args:
            shape: Grid shape ``(N_r, N_phi)``.
            include_u_mag: If ``True``, ensure ``u_mag`` snapshot dataset exists.

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

    @staticmethod
    def _append(ds: h5py.Dataset, arr: np.ndarray) -> None:
        """Append one row/frame to a resizable HDF5 dataset."""

        idx = ds.shape[0]
        ds.resize((idx + 1,) + ds.shape[1:])
        ds[idx] = arr

    def append_snapshot(self, t: float, n: np.ndarray, omega: np.ndarray, u_mag: np.ndarray | None = None) -> None:
        """Append one field snapshot to ``fields/snapshots.h5``.

        Args:
            t: Snapshot timestamp.
            n: Density field array.
            omega: Vorticity field array.
            u_mag: Optional velocity-magnitude field for movie generation.

        Returns:
            None.
        """

        self._ensure_field_datasets(n.shape, include_u_mag=u_mag is not None)
        self._append(self._t_snap, np.array(float(t)))
        self._append(self._n_ds, np.asarray(n, dtype=float))
        self._append(self._omega_ds, np.asarray(omega, dtype=float))
        if u_mag is not None and self._u_mag_ds is not None:
            self._append(self._u_mag_ds, np.asarray(u_mag, dtype=float))

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

        axes[0].plot(t, curves["E_wob"], label="E_wob")
        axes[0].plot(t, curves["E1"], label="E1", alpha=0.7)
        axes[0].plot(t, curves["E2"], label="E2", alpha=0.7)
        axes[0].set_title("Wobble Energies")
        axes[0].legend()

        axes[1].plot(t, curves["L_w"], color="C1")
        axes[1].set_title("Boundary Loss L_w")

        axes[2].plot(t, curves["sigma_r"], label="sigma_r")
        axes[2].plot(t, curves["r_mean"], label="r_mean", alpha=0.8)
        axes[2].set_title("Ring Coherence")
        axes[2].legend()

        axes[3].plot(t, curves["E_u"], label="E_u")
        axes[3].plot(t, curves["P_ctrl"], label="P_ctrl")
        if "eta_T" in curves:
            axes[3].plot(t, curves["eta_T"], label="eta_T")
        axes[3].set_title("Flow/Actuation")
        axes[3].legend()

        axes[4].plot(t, curves.get("L_w_cum", []), label="L_w_cum")
        axes[4].plot(t, curves.get("E_wob_excess", []), label="E_wob_excess")
        axes[4].plot(t, curves.get("sigma_r_excess", []), label="sigma_r_excess")
        axes[4].set_title("Degradation Terms")
        axes[4].legend()

        axes[5].plot(t, curves.get("badness_score", []), color="C3")
        axes[5].set_title("badness_score")

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

    def _write_movie(self) -> None:
        """Write GIF movies from stored ``n``, ``omega``, and ``u_mag`` fields."""

        n_data = np.asarray(self._n_ds)
        omega_data = np.asarray(self._omega_ds)
        phi = np.asarray(self._fields_h5["phi"])
        r = np.asarray(self._fields_h5["r"])
        self._write_field_movie(n_data, phi, r, out_name="n.gif", label="n", cmap="viridis")
        self._write_field_movie(omega_data, phi, r, out_name="omega.gif", label="omega", cmap="RdBu_r")
        if self._u_mag_ds is not None:
            u_mag = np.asarray(self._u_mag_ds)
            self._write_field_movie(u_mag, phi, r, out_name="u_mag.gif", label="|u|", cmap="magma")

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
            "This run uses stronger non-axisymmetric transport and multi-scale forcing to increase visible dynamics.",
            "",
            "## Key Levers",
            f"- time: dt={cfg.time.dt}, T_final={cfg.time.T_final}, snapshot_every={cfg.time.snapshot_every}",
            f"- pde: nu={cfg.pde.nu}, gamma={cfg.pde.gamma}, D_r={cfg.pde.D_r}, D_phi={cfg.pde.D_phi}, alpha={cfg.pde.alpha}, beta_instab={cfg.pde.beta_instab}",
            f"- forcing: drive_u0_base={cfg.forcing.drive_u0_base}, u_bounds={list(cfg.forcing.u_bounds)}",
            f"- disturbance scales: warmup={cfg.disturbance.warmup_scale}, excite={cfg.disturbance.excite_scale}, hold={cfg.disturbance.hold_scale}",
            f"- multiscale modes: {list(cfg.disturbance.multiscale_modes)}",
            f"- pseudo-random band forcing: m=[{cfg.disturbance.noise_m_min}, {cfg.disturbance.noise_m_max}], strength={cfg.disturbance.noise_strength}",
            "",
            "Expected effect: faster axisymmetry breaking, richer azimuthal spectrum, and larger frame-to-frame changes.",
        ]
        (self.run_dir / "notes" / "variation_note.md").write_text("\n".join(lines), encoding="utf-8")

    def _write_variation_diagnostics(self) -> None:
        """Compute and persist variation diagnostics from saved snapshots.

        The method writes:
        - ``metrics/variation_diagnostics.json``
        - ``plots/variation_diagnostics.png``

        It optionally compares temporal activity against a baseline run when
        ``ERD_VARIATION_BASELINE_RUN`` points to another run directory.

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

        self._fields_h5.flush()
        self._ctrl_h5.flush()
        self._fields_h5.close()
        self._ctrl_h5.close()

        self.log(f"Run artifacts written: {self.run_dir}")
