"""CLI entrypoint for a single open-loop ERD plant demo run."""

from __future__ import annotations

import argparse
from dataclasses import replace
import sys
from pathlib import Path

if __package__ in (None, ""):
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)

from erd_fipy import ERDPlant, load_run_config


def build_parser() -> argparse.ArgumentParser:
    """Build CLI parser for open-loop demo runs.

    Returns:
        Configured argument parser.
    """

    p = argparse.ArgumentParser(description="Run toy ERD open-loop demo")
    p.add_argument("--config", type=str, default=None, help="Optional plant YAML config")
    p.add_argument("--preset", type=str, default="smoke", choices=["smoke", "report"])
    p.add_argument("--tag", type=str, default="open_loop_demo")
    return p


def main() -> None:
    """Run one open-loop simulation and print run outputs."""

    args = build_parser().parse_args()

    cfg = load_run_config(args.config, preset=args.preset)
    cfg = replace(cfg, output=replace(cfg.output, tag=args.tag, preset=args.preset))

    plant = ERDPlant(cfg)
    result = plant.run(write_artifacts=True)
    print(f"run_dir: {result.run_dir}")
    print(f"summary: {result.summary}")


if __name__ == "__main__":
    main()
