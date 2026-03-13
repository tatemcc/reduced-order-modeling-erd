"""Dataset loading helpers for ERD run folders and manifest files.

The manifest is the single source of truth for mapping train/test splits to
concrete run directories produced by ``erd_fipy/scripts/generate_training_runs.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import h5py
import numpy as np
import yaml


@dataclass(frozen=True)
class TrajectoryData:
    """Model-ready trajectory payload.

    Attributes:
        run_dir: Source run directory path.
        t: Snapshot time vector of length ``T``.
        fields: Field tensor with shape ``(T, 2, N_r, N_phi)``.
        controls: Controls aligned to snapshots, shape ``(T, 5)``.
        dt: Mean snapshot spacing.
    """

    run_dir: str
    t: np.ndarray
    fields: np.ndarray  # (T, 2, Nr, Nphi)
    controls: np.ndarray  # (T, 5)
    dt: float


def _candidate_manifests(search_dir: Path) -> List[Path]:
    """List manifest files in newest-first order from a directory.

    Args:
        search_dir: Directory to scan.

    Returns:
        Candidate manifest paths sorted by descending modification time.
    """

    if not search_dir.exists() or not search_dir.is_dir():
        return []
    cands = list(search_dir.glob("*manifest*.yaml"))
    cands.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return cands


def _resolve_manifest_path(path: str | Path) -> Path:
    """Resolve a manifest path from an explicit file or directory.

    Accepted forms:
    1) Direct file path to a manifest.
    2) Directory path containing one or more ``*manifest*.yaml`` files.
       The newest manifest in that directory is selected.

    Args:
        path: Path-like value provided by the model config.

    Returns:
        Resolved absolute manifest file path.

    Raises:
        ValueError: If ``path`` is empty or no manifest can be found.
    """

    raw = str(path).strip()
    if not raw:
        raise ValueError(
            "data.manifest_path is required. "
            "Pass an explicit manifest file path or a directory containing one."
        )

    p = Path(raw).expanduser().resolve()
    if p.is_file():
        return p
    if p.is_dir():
        cands = _candidate_manifests(p)
        if cands:
            return cands[0].resolve()
        raise ValueError(f"No manifest files found in directory: {p}")

    raise ValueError(f"Manifest path does not exist: {p}")


def load_manifest(path: str | Path) -> Dict:
    """Load a dataset manifest YAML.

    Args:
        path: Manifest file path or directory path.

    Returns:
        Parsed manifest dictionary containing split metadata and lineage fields.
    """

    manifest_path = _resolve_manifest_path(path)
    with manifest_path.open("r", encoding="utf-8") as f:
        payload = yaml.safe_load(f) or {}
    if "splits" not in payload:
        raise ValueError(f"manifest missing 'splits': {manifest_path}")
    payload["_resolved_manifest_path"] = str(manifest_path)
    payload.setdefault("dataset_id", manifest_path.stem.replace("_manifest", ""))
    return payload


def _align_controls(t_target: np.ndarray, t_ctrl: np.ndarray, u_ctrl: np.ndarray) -> np.ndarray:
    """Align controls to snapshot timestamps using zero-order hold.

    The plant writes controls as piecewise-constant values indexed by control
    update time. For ROM identification, linear interpolation is intentionally
    avoided because it smooths abrupt actuator changes.

    Args:
        t_target: Snapshot timestamp vector of shape ``(T,)``.
        t_ctrl: Control timestamp vector of shape ``(K,)``.
        u_ctrl: Raw controls with shape ``(K, n_u)``.

    Returns:
        Snapshot-aligned control matrix of shape ``(T, n_u)``.
    """

    t_target = np.asarray(t_target, dtype=float).reshape(-1)
    t_ctrl = np.asarray(t_ctrl, dtype=float).reshape(-1)
    u_ctrl = np.asarray(u_ctrl, dtype=float)

    if t_ctrl.size == 0:
        raise ValueError("Control time series is empty")
    if u_ctrl.ndim != 2:
        raise ValueError(f"Expected control array with shape (K, n_u), got {u_ctrl.shape}")
    if u_ctrl.shape[0] != t_ctrl.size:
        raise ValueError(
            f"Control timestamp/value length mismatch: t_ctrl={t_ctrl.size}, u_ctrl={u_ctrl.shape[0]}"
        )
    if np.any(np.diff(t_ctrl) < 0.0):
        raise ValueError("Control timestamps must be nondecreasing")

    idx = np.searchsorted(t_ctrl, t_target, side="right") - 1
    idx = np.clip(idx, 0, t_ctrl.size - 1)
    return u_ctrl[idx]


def load_trajectory(run_dir: str | Path) -> TrajectoryData:
    """Load one trajectory from a plant run folder.

    Args:
        run_dir: Path to a run directory with ``fields`` and ``controls`` HDF5 files.

    Returns:
        Structured trajectory data with controls aligned to snapshot timestamps.
    """

    run = Path(run_dir)
    fields_path = run / "fields" / "snapshots.h5"
    ctrl_path = run / "controls" / "control_timeseries.h5"

    with h5py.File(fields_path, "r") as hf:
        t = np.asarray(hf["t_snap"])
        n = np.asarray(hf["n"])
        omega = np.asarray(hf["omega"])

    with h5py.File(ctrl_path, "r") as hc:
        t_ctrl = np.asarray(hc["t_ctrl"])
        u_ctrl = np.asarray(hc["u"])

    if t.size < 3:
        raise ValueError(f"Trajectory too short in {run}")

    controls = _align_controls(t, t_ctrl, u_ctrl)
    fields = np.stack([n, omega], axis=1)
    dt = float(np.mean(np.diff(t)))
    return TrajectoryData(run_dir=str(run), t=t, fields=fields, controls=controls, dt=dt)


def load_split(
    manifest: Dict,
    split_name: str,
    max_trajectories: Optional[int] = None,
) -> List[TrajectoryData]:
    """Load all trajectories from one manifest split.

    Args:
        manifest: Parsed manifest dictionary.
        split_name: Split key to load (for example ``train``).
        max_trajectories: Optional cap on number of trajectories to load.

    Returns:
        List of trajectories trimmed to a common temporal length.
    """

    entries = list(manifest["splits"].get(split_name, []))
    if max_trajectories is not None:
        entries = entries[:max_trajectories]

    traj = [load_trajectory(entry["run_dir"]) for entry in entries]
    if not traj:
        return []

    nmin = min(t.fields.shape[0] for t in traj)
    nr = traj[0].fields.shape[2]
    nphi = traj[0].fields.shape[3]

    out: List[TrajectoryData] = []
    for t in traj:
        if t.fields.shape[2] != nr or t.fields.shape[3] != nphi:
            raise ValueError("Inconsistent grid shape across trajectories")
        out.append(
            TrajectoryData(
                run_dir=t.run_dir,
                t=t.t[:nmin],
                fields=t.fields[:nmin],
                controls=t.controls[:nmin],
                dt=t.dt,
            )
        )
    return out
