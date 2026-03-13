"""CLI entrypoint for baseline-vs-MPC closed-loop ERD runs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

if __package__ in (None, ""):
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)

from erd_control import load_control_config, run_closed_loop


def build_parser() -> argparse.ArgumentParser:
    """Build CLI arguments for closed-loop execution.

    Returns:
        Parser configured with required control config argument.
    """

    p = argparse.ArgumentParser(description="Run closed-loop ERD control")
    p.add_argument("--config", type=str, required=True, help="Control YAML config path")
    return p


def main() -> None:
    """Load control config, run baseline and MPC stages, print run folder."""

    args = build_parser().parse_args()
    cfg = load_control_config(args.config)
    run_dir = run_closed_loop(cfg)
    print(f"run_dir: {run_dir}")


if __name__ == "__main__":
    main()
