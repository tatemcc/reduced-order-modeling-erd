"""Lightweight artifact serialization helpers for model runs."""

from __future__ import annotations

import json
import pickle
from pathlib import Path
from typing import Any, Dict

import numpy as np
import yaml


def save_json(path: Path, data: Dict[str, Any]) -> None:
    """Write a dictionary to JSON.

    Args:
        path: Destination file path.
        data: JSON-serializable mapping.

    Returns:
        None.
    """

    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)


def save_yaml(path: Path, data: Dict[str, Any]) -> None:
    """Write a dictionary to YAML.

    Args:
        path: Destination file path.
        data: YAML-serializable mapping.

    Returns:
        None.
    """

    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False)


def save_npy(path: Path, arr: np.ndarray) -> None:
    """Persist an array using NumPy ``.npy`` format.

    Args:
        path: Destination file path.
        arr: NumPy array to store.

    Returns:
        None.
    """

    np.save(path, arr)


def save_pickle(path: Path, obj: Any) -> None:
    """Serialize an object using Python pickle.

    Args:
        path: Destination file path.
        obj: Object to serialize.

    Returns:
        None.
    """

    with path.open("wb") as f:
        pickle.dump(obj, f)


def load_pickle(path: Path) -> Any:
    """Load a pickled object from disk.

    Args:
        path: Pickle file path.

    Returns:
        Deserialized Python object.
    """

    with path.open("rb") as f:
        return pickle.load(f)
