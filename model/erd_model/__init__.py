"""Public package API for ERD model fitting.

This module keeps imports lightweight by avoiding eager loading of the
full pipeline stack (PySINDy, plotting, etc.) at package import time.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .config import RunConfig, load_run_config

if TYPE_CHECKING:
    from .run_pipeline import ModelRunResult as ModelRunResult


def run(cfg: RunConfig):
    """Execute the model fitting pipeline.

    Args:
        cfg: Parsed run configuration.

    Returns:
        Model run result bundle from :mod:`erd_model.run_pipeline`.
    """

    from .run_pipeline import run as _run

    return _run(cfg)


def __getattr__(name: str):
    """Lazily resolve heavyweight attributes on first access.

    Args:
        name: Attribute name requested from the package namespace.

    Returns:
        Resolved attribute object.
    """

    if name == "ModelRunResult":
        from .run_pipeline import ModelRunResult as _ModelRunResult

        return _ModelRunResult
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ["RunConfig", "load_run_config", "ModelRunResult", "run"]
