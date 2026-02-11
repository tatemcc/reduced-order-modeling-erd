"""
Filesystem and path utilities for the Burgers ROM + SINDy pipeline.

This module centralizes:
- Repository root discovery
- Deterministic run ID generation from configurations
- Construction of run output directories
- Creation of standard subdirectory structure for artifacts
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Dict, Optional


def find_repo_root(start: Optional[Path] = None) -> Path:
    """
    Locate the repository root directory.

    The repository root is defined as the nearest parent directory
    (including the starting directory) containing a pyproject.toml file.

    Parameters
    ----------
    start : pathlib.Path, optional
        Directory to start searching from. Defaults to this file's
        location.

    Returns
    -------
    pathlib.Path
        Path to the repository root directory.
    """
    here = (start or Path(__file__)).resolve()
    for parent in (here, *here.parents):
        if (parent / "pyproject.toml").is_file():
            return parent
    raise RuntimeError("Could not locate repo root (pyproject.toml not found).")


def _to_jsonable(obj: Any) -> Any:
    """
    Convert an object to a JSON-serializable representation.

    This function recursively converts dataclasses, dictionaries,
    lists, and tuples into plain Python types suitable for JSON
    serialization.

    Parameters
    ----------
    obj : Any
        Object to convert.

    Returns
    -------
    Any
        JSON-serializable representation of the object.
    """
    if is_dataclass(obj):
        return asdict(obj)
    if isinstance(obj, dict):
        return {k: _to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(v) for v in obj]
    return obj


def config_hash(config_dict: Dict[str, Any]) -> str:
    """
    Compute a deterministic short hash for a configuration dictionary.

    The hash is computed from a canonical JSON encoding of the
    configuration with sorted keys.

    Parameters
    ----------
    config_dict : dict
        Dictionary representation of a configuration.

    Returns
    -------
    str
        Short hexadecimal hash string identifying the configuration.
    """
    payload = json.dumps(
        config_dict, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:12]


def run_dir(
    outputs_dir: str,
    equation: str,
    run_id: str,
) -> Path:
    """
    Construct the path to a run's outputs, nested
    underneath the configured output directory.

    Parameters
    ----------
    outputs_dirn : str
        Path of the outputs directory.
    equation : str
        Equation name (e.g., 'burgers').
    run_id : str
        Unique identifier for the run.

    Returns
    -------
    pathlib.Path
        Path to the run's output directory.
    """
    return Path(outputs_dir) / equation / run_id


def ensure_run_subdirs(rundir: Path) -> Dict[str, Path]:
    """
    Create the standard subdirectory structure for a run.

    Subdirectories created:
    - pod
    - sindy
    - rollout
    - metrics
    - figures
    - movies

    Parameters
    ----------
    rundir : pathlib.Path
        Path to the run's root directory.

    Returns
    -------
    dict
        Mapping from subdirectory names to their paths.
    """
    sub = {
        "pod": rundir / "pod",
        "sindy": rundir / "sindy",
        "rollout": rundir / "rollout",
        "metrics": rundir / "metrics",
        "figures": rundir / "figures",
        "movies": rundir / "movies",
    }
    rundir.mkdir(parents=True, exist_ok=True)
    for p in sub.values():
        p.mkdir(parents=True, exist_ok=True)
    return sub
