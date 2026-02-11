"""
Artifact serialization utilities for the Burgers ROM + SINDy pipeline.

This module provides small, explicit helpers for saving run outputs
to disk in a structured and reproducible way. It does not impose
any directory structure by itself; it only writes to paths provided
by the caller.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

import numpy as np
import yaml

from .paths import _to_jsonable


def save_config_yaml(path: Path | str, config: Any) -> None:
    """
    Save a configuration object to a YAML file.

    The configuration is first converted to a JSON-serializable
    dictionary (including nested dataclasses).

    Parameters
    ----------
    path : pathlib.Path or str
        Destination path for the YAML file.
    config : Any
        Configuration object, typically a RunConfig instance.
    """
    data = _to_jsonable(config)
    path = Path(path).resolve()
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False)


def save_json(path: Path, data: Dict[str, Any]) -> None:
    """
    Save a dictionary to a JSON file.

    Parameters
    ----------
    path : pathlib.Path
        Destination path for the JSON file.
    data : dict
        Dictionary containing JSON-serializable values.
    """
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)


def save_npy(path: Path, arr: np.ndarray) -> None:
    """
    Save a NumPy array to disk in .npy format.

    Parameters
    ----------
    path : pathlib.Path
        Destination path for the array.
    arr : np.ndarray
        Array to be saved.
    """
    np.save(path, arr)
