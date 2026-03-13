"""CLI entrypoint for ERD POD + controlled SINDy model fitting."""

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

from erd_model import load_run_config, run


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI argument parser.

    Returns:
        Parser configured for model pipeline execution.
    """

    p = argparse.ArgumentParser(description="Run ERD POD+SINDy model pipeline")
    p.add_argument("--config", type=str, required=True, help="Model YAML config path")
    return p


def main() -> None:
    """Load configuration, run the pipeline, and print artifact locations."""

    args = build_parser().parse_args()
    cfg = load_run_config(args.config)
    res = run(cfg)
    print(f"run_dir: {res.run_dir}")
    print(f"aggregates: {res.aggregates}")


if __name__ == "__main__":
    main()
