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
    PODConfig,
    RolloutConfig,
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


def _maybe_get(d: Dict[str, Any], key: str, default: Any = None, raise_if_missing: bool = True) -> Any:
    """
    Return d[key] if present else default.

    Parameters
    ----------
    d : dict
        Dictionary to read.
    key : str
        Key to lookup.
    default : Any
        Fallback value.

    Returns
    -------
    Any
        Value.
    """
    if key in d:
        return d[key]
    elif not raise_if_missing:
        print(f"Couldn't find {key} in config, using default={default}")
        return default
    else:
        raise ValueError(f"Couldn't find {key} in config")


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

    data_d = cfg.get("data", {})
    pod_d = cfg.get("pod", {})
    deriv_d = cfg.get("deriv", {})
    sindy_d = cfg.get("sindy", {})
    rollout_d = cfg.get("rollout", {})

    data = DataConfig(
        equation=_maybe_get(data_d, "equation"),
        structure=_maybe_get(data_d, "structure"),
        resolution=_maybe_get(data_d, "resolution"),
        split=_maybe_get(data_d, "split", "train"),
        n_trajectories=_maybe_get(data_d, "n_trajectories"),
        trajectory_ids=_maybe_get(data_d, "trajectory_ids"),
        data_path=_maybe_get(data_d, "data_path"),
        download_if_missing=_maybe_get(data_d, "download_if_missing"),
        lookback=_maybe_get(data_d, "lookback"),
        rollout=_maybe_get(data_d, "rollout"),
        squeeze_lookback_dim=_maybe_get(data_d, "squeeze_lookback_dim"),
        time_stride=_maybe_get(data_d, "time_stride"),
        time_limit=_maybe_get(data_d, "time_limit"),
    )

    pod = PODConfig(
        rank=_maybe_get(pod_d, "rank"),
        energy_fraction=_maybe_get(pod_d, "energy_fraction"),
        center=_maybe_get(pod_d, "center"),
    )

    deriv = DerivConfig(
        method=_maybe_get(deriv_d, "method"),
        scheme=_maybe_get(deriv_d, "scheme"),
    )

    sindy = SINDyConfig(
        poly_order=_maybe_get(sindy_d, "poly_order"),
        include_bias=_maybe_get(sindy_d, "include_bias"),
        optimizer=_maybe_get(sindy_d, "optimizer"),
        sparsity=_maybe_get(sindy_d, "sparsity"),
        constrain_energy=_maybe_get(sindy_d, "constrain_energy"),
    )

    rollout = RolloutConfig(
        horizon_steps=_maybe_get(rollout_d, "horizon_steps"),
    )

    return RunConfig(
        name=_maybe_get(cfg, "name"),
        seed=_maybe_get(cfg, "seed"),
        data=data,
        pod=pod,
        deriv=deriv,
        sindy=sindy,
        rollout=rollout,
        outputs_dir=_maybe_get(cfg, "outputs_dir"),
    )


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
