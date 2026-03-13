"""Path helpers for model-run artifact directories."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Dict


SUBDIRS = {
    "fields": "fields",
    "controls": "controls",
    "metrics": "metrics",
    "plots": "plots",
    "movies": "movies",
    "model": "model",
    "logs": "logs",
}


def make_run_dir(outputs_root: str, tag: str) -> Path:
    """Create a timestamped model run directory and standard subfolders.

    Args:
        outputs_root: Parent path for all run directories.
        tag: Run tag appended to the timestamp.

    Returns:
        Absolute path to the created run directory.
    """

    root = Path(outputs_root).resolve()
    root.mkdir(parents=True, exist_ok=True)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base = root / f"{stamp}_{tag}"
    run_dir = base
    idx = 1
    while run_dir.exists():
        run_dir = Path(f"{base}_{idx:02d}")
        idx += 1

    run_dir.mkdir(parents=True, exist_ok=False)
    for name in SUBDIRS.values():
        (run_dir / name).mkdir(parents=True, exist_ok=True)
    return run_dir


def ensure_subdirs(run_dir: Path) -> Dict[str, Path]:
    """Ensure all standard subdirectories exist under ``run_dir``.

    Args:
        run_dir: Existing or new run directory.

    Returns:
        Mapping from logical subdir keys to absolute paths.
    """

    out: Dict[str, Path] = {}
    run_dir.mkdir(parents=True, exist_ok=True)
    for key, name in SUBDIRS.items():
        p = run_dir / name
        p.mkdir(parents=True, exist_ok=True)
        out[key] = p
    return out
