"""
Run a single Burgers ROM + SINDy experiment and write artifacts.

This script:
- Loads a RunConfig from YAML or uses defaults
- Runs the full pipeline
- Writes outputs under outputs/<equation>/<run_id>/
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, Optional

import yaml

from burgers_rom.config import (
    DataConfig,
    DerivConfig,
    PlotConfig,
    PODConfig,
    RunConfig,
    SINDyConfig,
)
from burgers_rom.run_pipeline import run


def _load_yaml(path: Path) -> Dict[str, Any]:
    """
    Load a YAML file into a Python dictionary.

    Parameters
    ----------
    path : Path
        Path to the YAML file.

    Returns
    -------
    dict
        Parsed YAML content.
    """
    with path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError("Config YAML must parse to a dictionary")
    return data


# Using dictionary unpacking to populate config objects with fallbacks instead of raising errors
def load_run_config(path: Optional[Path]) -> RunConfig:
    """
    Load RunConfig from YAML or return default RunConfig if path is None.
    Parameters
    ----------
    path : Path, optional
        YAML config path.

    Returns
    -------
    RunConfig
        Run configuration.
    """
    if path is None:
        return RunConfig()

    cfg = _load_yaml(path)

    # Load sub-configs using unpacking to respect dataclass defaults
    data = DataConfig(**cfg.get("data", {}))
    pod = PODConfig(**cfg.get("pod", {}))
    deriv = DerivConfig(**cfg.get("deriv", {}))
    sindy = SINDyConfig(**cfg.get("sindy", {}))
    plots = PlotConfig(**cfg.get("plots", {}))

    # Filter top-level arguments to avoid passing sub-config dicts to RunConfig
    run_kwargs = {k: v for k, v in cfg.items() if k not in ["data", "pod", "deriv", "sindy", "plots"]}

    return RunConfig(data=data, pod=pod, deriv=deriv, sindy=sindy, plots=plots, **run_kwargs)


def build_argparser() -> argparse.ArgumentParser:
    """
    Build CLI argument parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for run_one.
    """
    p = argparse.ArgumentParser(description="Run Burgers ROM + SINDy pipeline once.")
    p.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to YAML config. If omitted, uses defaults.",
    )
    p.add_argument(
        "--no-write",
        action="store_true",
        help="Do not write artifacts to disk.",
    )
    return p


def main() -> None:
    """
    CLI entry point.
    """
    parser = build_argparser()
    args = parser.parse_args()

    cfg_path = Path(args.config).resolve() if args.config is not None else None
    cfg = load_run_config(cfg_path)

    result = run(cfg, write_artifacts=not args.no_write)

    print(f"run_id: {result.run_id}")
    print(f"rundir: {result.rundir}")
    for k, v in result.aggregates.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
