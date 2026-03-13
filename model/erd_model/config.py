"""Configuration schema and YAML loading for ERD ROM/SINDy fitting."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Optional

import yaml



def _repo_root(start: Optional[Path] = None) -> Path:
    """Locate the ERD-capstone workspace root directory.

    Args:
        start: Optional path to start searching from.

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
class DataConfig:
    """Data-loading and split selection options.

    Attributes:
        manifest_path: Path to the generated training manifest YAML.
        train_split: Split key used for training trajectories.
        val_split: Split key used for validation trajectories.
        test_split: Split key used for test trajectories.
        max_train: Optional cap on number of training trajectories.
        max_val: Optional cap on number of validation trajectories.
        max_test: Optional cap on number of test trajectories.
        time_stride: Temporal stride for snapshot subsampling.
        time_limit: Optional max number of time points per trajectory.
    """

    manifest_path: str = ""
    train_split: str = "train"
    val_split: str = "val"
    test_split: str = "test"
    max_train: Optional[int] = None
    max_val: Optional[int] = None
    max_test: Optional[int] = None
    time_stride: int = 1
    time_limit: Optional[int] = None


@dataclass(frozen=True)
class PODConfig:
    """POD basis construction options.

    Attributes:
        rank: Requested POD rank (later clipped to allowed range).
        center: Whether to subtract the global mean state before SVD.
    """

    rank: int = 6
    center: bool = True


@dataclass(frozen=True)
class DerivConfig:
    """Finite-difference derivative options.

    Attributes:
        order: Finite-difference order (currently central difference footprint).
    """

    order: int = 2


@dataclass(frozen=True)
class TrappingSR3Config:
    """TrappingSR3 settings for stabilized autonomous reduced dynamics.

    Attributes:
        method: Trapping objective mode (``local`` or ``global``).
        reg_weight_lam: Sparse regularization weight for TrappingSR3.
        relax_coeff_nu: SR3 relaxation coefficient.
        tol: Optimizer stopping tolerance.
        max_iter: Maximum TrappingSR3 iterations.
        eta: Weight for the stability term ``||Pw - A||^2``.
        stability_alpha: Weight for the ``||Qijk||`` local-stability term.
        stability_beta: Weight for the ``||Qijk + Qjik + Qkij||`` term.
        control_threshold: Sparsity threshold for residual control terms.
        control_alpha: Ridge regularization for residual control fit.
        control_max_iter: Iterations for residual control sparse refit.
    """

    method: str = "local"
    reg_weight_lam: float = 0.01
    relax_coeff_nu: float = 1.0
    tol: float = 1e-6
    max_iter: int = 160
    eta: float = 10.0
    stability_alpha: float = 1e-3
    stability_beta: float = 1e-3
    control_threshold: float = 0.02
    control_alpha: float = 1e-4
    control_max_iter: int = 4


@dataclass(frozen=True)
class SINDyConfig:
    """Controlled SINDy optimizer settings.

    Attributes:
        optimizer: Sparse optimizer name (``stlsq`` or ``trappingsr3``).
        threshold: Sparsity threshold used by STLSQ.
        alpha: L2 regularization coefficient.
        max_iter: Maximum STLSQ iterations.
        backoff_factor: Multiplicative sparsity-threshold reduction on retries.
        backoff_retries: Maximum retry count for sparsity backoff.
        fallback_optimizer: Optimizer to try if primary fitting fails repeatedly.
        validation_spike_threshold: Held-out mean relative error threshold that
            triggers a backoff retry.
        trapping: Trapping-SR3 settings used when ``optimizer=trappingsr3``.
    """

    optimizer: str = "stlsq"
    threshold: float = 0.05
    alpha: float = 0.0
    max_iter: int = 80
    backoff_factor: float = 0.1
    backoff_retries: int = 3
    fallback_optimizer: str = "stlsq"
    validation_spike_threshold: float = 0.60
    trapping: TrappingSR3Config = TrappingSR3Config()


@dataclass(frozen=True)
class AcceptanceConfig:
    """Model acceptance gates for strict ROM quality checks.

    Attributes:
        min_pod_energy: Minimum retained POD energy required at chosen rank.
        max_rank: Maximum allowed POD rank after clipping.
        max_mean_field_rel_l2: Maximum acceptable held-out mean relative field
            error.
    """

    min_pod_energy: float = 0.92
    max_rank: int = 8
    max_mean_field_rel_l2: float = 0.20


@dataclass(frozen=True)
class RolloutConfig:
    """Validation rollout settings.

    Attributes:
        horizon_steps: Maximum forecast horizon for held-out trajectory checks.
    """

    horizon_steps: int = 120


@dataclass(frozen=True)
class PlotConfig:
    """Plot and GIF generation settings.

    Attributes:
        dpi: Figure resolution in dots-per-inch.
        movie_fps: GIF frame rate for rollout animations.
    """

    dpi: int = 150
    movie_fps: int = 15


@dataclass(frozen=True)
class OutputConfig:
    """Run-folder output naming.

    Attributes:
        outputs_root: Parent directory for model run folders.
        tag: Tag appended to timestamped run folder names.
    """

    outputs_root: str = DEFAULT_OUTPUTS_ROOT
    tag: str = "auto"


@dataclass(frozen=True)
class RunConfig:
    """Top-level model-fitting configuration.

    Attributes:
        name: Human-readable run name.
        seed: Random seed used for deterministic components.
        data: Data loading and split options.
        pod: POD settings.
        deriv: Derivative estimation settings.
        sindy: Sparse-regression settings.
        acceptance: Model acceptance gates.
        rollout: Held-out rollout settings.
        plots: Plot/movie generation settings.
        output: Output directory and tag settings.
    """

    name: str = "erd_model_fit"
    seed: int = 0
    data: DataConfig = DataConfig()
    pod: PODConfig = PODConfig()
    deriv: DerivConfig = DerivConfig()
    sindy: SINDyConfig = SINDyConfig()
    acceptance: AcceptanceConfig = AcceptanceConfig()
    rollout: RolloutConfig = RolloutConfig()
    plots: PlotConfig = PlotConfig()
    output: OutputConfig = OutputConfig()

    def to_dict(self) -> Dict[str, Any]:
        """Convert the config tree into plain serializable dictionaries.

        Args:
            None.

        Returns:
            Nested dictionary representation of this run configuration.
        """

        return asdict(self)



def load_run_config(path: str | Path) -> RunConfig:
    """Load model pipeline configuration from YAML.

    Args:
        path: Path to YAML config file.

    Returns:
        Parsed and typed :class:`RunConfig`.
    """

    p = Path(path)
    with p.open("r", encoding="utf-8") as f:
        payload = yaml.safe_load(f) or {}

    sindy_payload = payload.get("sindy", {})
    trapping_payload = sindy_payload.get("trapping", {})

    cfg = RunConfig(
        name=payload.get("name", RunConfig().name),
        seed=payload.get("seed", RunConfig().seed),
        data=DataConfig(**payload.get("data", {})),
        pod=PODConfig(**payload.get("pod", {})),
        deriv=DerivConfig(**payload.get("deriv", {})),
        sindy=SINDyConfig(
            optimizer=sindy_payload.get("optimizer", SINDyConfig().optimizer),
            threshold=float(sindy_payload.get("threshold", SINDyConfig().threshold)),
            alpha=float(sindy_payload.get("alpha", SINDyConfig().alpha)),
            max_iter=int(sindy_payload.get("max_iter", SINDyConfig().max_iter)),
            backoff_factor=float(sindy_payload.get("backoff_factor", SINDyConfig().backoff_factor)),
            backoff_retries=int(sindy_payload.get("backoff_retries", SINDyConfig().backoff_retries)),
            fallback_optimizer=str(sindy_payload.get("fallback_optimizer", SINDyConfig().fallback_optimizer)),
            validation_spike_threshold=float(
                sindy_payload.get("validation_spike_threshold", SINDyConfig().validation_spike_threshold)
            ),
            trapping=TrappingSR3Config(**trapping_payload),
        ),
        acceptance=AcceptanceConfig(**payload.get("acceptance", {})),
        rollout=RolloutConfig(**payload.get("rollout", {})),
        plots=PlotConfig(**payload.get("plots", {})),
        output=OutputConfig(**payload.get("output", {})),
    )
    manifest_raw = str(cfg.data.manifest_path).strip()
    if not manifest_raw or manifest_raw.startswith("REQUIRED_"):
        raise ValueError("data.manifest_path must be set to a dataset manifest path")
    return cfg
