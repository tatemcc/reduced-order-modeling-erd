"""Configuration schema and YAML helpers for closed-loop ERD control runs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

import yaml



def _repo_root(start: Optional[Path] = None) -> Path:
    """Locate the ERD-capstone workspace root.

    Args:
        start: Optional path used as the search anchor.

    Returns:
        Workspace root path containing ``erd_fipy``, ``model``, and ``control``.
    """

    here = (start or Path(__file__)).resolve()
    for parent in (here, *here.parents):
        if all((parent / name).is_dir() for name in ("erd_fipy", "model", "control")):
            return parent
    for parent in (here, *here.parents):
        if (parent / "pyproject.toml").is_file():
            return parent
    raise RuntimeError("Could not locate repository root")


REPO_ROOT = _repo_root()
DEFAULT_OUTPUTS_ROOT = (REPO_ROOT / "outputs" / "runs").as_posix()


@dataclass(frozen=True)
class ControlRunConfig:
    """Top-level closed-loop orchestration config.

    Attributes:
        model_run_dir: Path to a completed ERD model run folder.
        plant_config_path: Optional explicit plant YAML config path.
        preset: Plant preset name used when no explicit plant config is supplied.
        seed: Random seed for MPC shooting reproducibility.
        outputs_root: Parent directory for canonical closed-loop run folders.
        tag: Suffix tag for the timestamped closed-loop run folder.
    """

    model_run_dir: str
    plant_config_path: Optional[str] = None
    preset: str = "smoke"
    seed: int = 0
    outputs_root: str = DEFAULT_OUTPUTS_ROOT
    tag: str = "auto"



def load_control_config(path: str | Path) -> ControlRunConfig:
    """Load closed-loop configuration from YAML.

    Args:
        path: Path to YAML config file.

    Returns:
        Parsed :class:`ControlRunConfig`.
    """

    p = Path(path)
    with p.open("r", encoding="utf-8") as f:
        payload = yaml.safe_load(f) or {}

    cfg = ControlRunConfig(**payload)
    if not str(cfg.model_run_dir).strip() or str(cfg.model_run_dir).startswith("REQUIRED_"):
        raise ValueError("model_run_dir must be set to a completed model run directory")
    return cfg



def to_dict(cfg: ControlRunConfig) -> Dict[str, Any]:
    """Convert a control config object into a serializable dictionary.

    Args:
        cfg: Closed-loop run configuration.

    Returns:
        Dictionary representation of ``cfg``.
    """

    return {
        "model_run_dir": cfg.model_run_dir,
        "plant_config_path": cfg.plant_config_path,
        "preset": cfg.preset,
        "seed": cfg.seed,
        "outputs_root": cfg.outputs_root,
        "tag": cfg.tag,
    }
