"""
CLI driver for the ERD FiPy simulator.

Supports being executed either as ``python -m erd_fipy.scripts.run_demo`` or
directly via ``python erd_fipy/scripts/run_demo.py`` (or ``uv run``).
"""
from __future__ import annotations

import sys
from pathlib import Path

# Allow running the script directly (`python erd_fipy/scripts/run_demo.py`)
# by ensuring the project root (which contains the `erd_fipy` package) is on sys.path.
if __package__ in (None, ""):
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)


from erd_fipy.stepping import run_sim

if __name__ == "__main__":
    outpath = run_sim(seed=0)
    print("Wrote:", outpath)
