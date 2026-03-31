"""
Configuration definitions for the Burgers ROM + SINDy pipeline.

This module defines all configuration objects as dataclasses.
Each configuration group controls one stage of the pipeline:
- data loading
- POD
- time differentiation
- SINDy regression
- rollout forecasting

Configurations are immutable and serializable.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, field
from typing import Literal, Optional, Sequence, Dict, Any
from .paths import find_repo_root

# ---------------------------
# Constants
# ---------------------------
REPO_ROOT = find_repo_root()
DATA_DIR = REPO_ROOT / "data"
OUTPUTS_DIR = REPO_ROOT / "outputs"

# ---------------------------
# Type aliases
# ---------------------------
EquationName = Literal["advection", "burgers", "gasdynamics", "kuramotosivashinsky", "reactiondiffusion", "wave"]
StructureName = Literal["grid"]
ResolutionName = Literal["low", "medium", "high"]


# ---------------------------
# Dataclasses
# ---------------------------
@dataclass(frozen=True)
class DataConfig:
    """
    Configuration for loading DynaBench data.

    Attributes
    ----------
    equation : str
        PDE name. Fixed to 'burgers'.
    structure : str
        Data structure type. Must be 'grid'.
    resolution : str
        Spatial resolution: 'low', 'medium', or 'high'.
    split : str
        Dataset split: 'train', 'val', or 'test'.
    n_trajectories : int
        Number of trajectories to load.
    trajectory_ids : sequence of int, optional
        Explicit trajectory indices to load.
    data_path : str
        Filesystem location where Dynabench data are stored.
    download_if_missing : bool
        If True, trigger Dynabench to download data to data_path when missing.
    lookback : int
        Number of past steps returned by DynabenchIterator.
    rollout : int
        Number of future steps returned by DynabenchIterator.
    squeeze_lookback_dim : bool
        Whether to remove lookback dimension when lookback == 1.
    time_stride : int
        Subsampling factor along the time axis.
    time_limit : int, optional
        Maximum number of time steps per trajectory.
    """

    equation: EquationName = "burgers"
    structure: StructureName = "grid"
    resolution: ResolutionName = "low"

    split: Literal["train", "val", "test"] = "train"

    n_trajectories: int = 1
    trajectory_ids: Optional[Sequence[int]] = None

    data_path: str | None = DATA_DIR.as_posix()
    download_if_missing: bool = False

    lookback: int = 1  # time dimension in X
    rollout: int = 1  # time dimension in y

    squeeze_lookback_dim: bool = True

    time_stride: int = 1
    time_limit: Optional[int] = None

    def __post_init__(self):

        if self.data_path is None:
            object.__setattr__(self, "data_path", DATA_DIR.as_posix())


@dataclass(frozen=True)
class PODConfig:
    """
    Configuration for POD (SVD) stage.

    Attributes
    ----------
    rank : int, optional
        Fixed POD rank.
    energy_fraction : float, optional
        Target cumulative energy fraction for rank selection.
    center : bool
        Whether to subtract the mean state before POD.
    """

    rank: Optional[int] = 20
    energy_fraction: Optional[float] = None

    center: bool = True


@dataclass(frozen=True)
class DerivConfig:
    """
    Configuration for time differentiation of POD coefficients.

    Attributes
    ----------
    method : str
        Differentiation method. Currently only 'finite_difference'.
    scheme : str
        Finite difference scheme. Currently 'central'.
    """

    method: Literal["finite_difference"] = "finite_difference"
    scheme: Literal["central"] = "central"


@dataclass(frozen=True)
class SINDyConfig:
    """
    Configuration for sparse system identification.

    Attributes
    ----------
    enabled : bool
        Whether to run the SINDy fitting and rollout stage.
    poly_order : int
        Polynomial library order.
    include_bias : bool
        Whether to include a constant term in the library.
    optimizer : str
        Sparse optimizer type: 'stlsq', 'sr3', 'ssr', or 'frols'.
    optimizer_params : dict
        Dictionary of optimizer-specific parameters.
    constrain_energy : bool
        Whether to enforce energy-preserving constraints.
    parallel : bool
        Whether to use parallel processing for fitting SINDy.
    parallel_procs : int, optional
        Number of processes for parallel fitting. If None, uses all cores.
    """

    enabled: bool = True
    poly_order: int = 2
    include_bias: bool = False

    optimizer: Literal["stlsq", "sr3", "ssr", "frols", "trappingsr3"] = "sr3"
    optimizer_params: Dict[str, Any] = field(default_factory=lambda: {"threshold": 0.1, "nu": 1.0})

    constrain_energy: bool = False  # NOTE UNIMPLEMENTED

    parallel: bool = False
    parallel_procs: Optional[int] = None


@dataclass(frozen=True)
class PlotConfig:
    """
    Configuration for generating figures and movies for a run.
    """
    enabled: bool = True

    figures_subdir: str = "figures"
    movies_subdir: str = "movies"

    dpi: int = 150

    # Basis visualization
    basis_n_modes: int = 8
    basis_cmap: str = "viridis"

    # POD decomposition matrix
    pod_decomposition_matrix: bool = True
    pod_decomposition_matrix_square_pixels: bool = True
    plot_true_state_matrix: bool = True

    # Rollout visualization
    rollout_cmap: str = "viridis"
    movie_fps: int = 15
    movie_every: int = 1
    movie_components: Optional[list[int]] = None  # None means all components

    # 3D surface movie
    movie_3d_surface: bool = False
    movie_3d_interp_factor: int = 3
    movie_3d_mode_contributions: bool = False
    movie_3d_parallel: bool = False
    movie_3d_parallel_procs: Optional[int] = None
    movie_3d_decomposition: bool = False
    movie_3d_decomposition_show_titles: bool = True
    movie_3d_reconstruction_comparison: bool = False
    movie_chronos_comparison: bool = False
    movie_3d_surface_with_state_matrix: bool = False
    movie_3d_surface_clean: bool = False

    # field error plots
    metrics_curves: bool = True
    dynabench_comparison_timesteps: int = 50

    # Coefficient diagnostics
    coeff_time_series: bool = True
    coeff_pair_phase: bool = True
    coeff_pair_max_pairs: int = 45  # caps number of (i,j) phase plots saved

    # Optional sympy labels
    sympy_labels: bool = False
    sympy_label_style: Literal["a_i", "x_i"] = "a_i"



@dataclass(frozen=True)
class RunConfig:
    """
    Top-level configuration for a single pipeline run.

    Attributes
    ----------
    name : str
        Human-readable run name.
    seed : int
        Random seed.
    data : DataConfig
        Data loading configuration.
    pod : PODConfig
        POD configuration.
    deriv : DerivConfig
        Time-derivative configuration.
    sindy : SINDyConfig
        Sparse regression configuration.
    outputs_dir : str | None
        Path to the outputs directory. Runs will be
        nested under this directory.
    """

    name: str = "burgers_rom"
    seed: int = 0  # NOTE: UNUSED SO FAR

    data: DataConfig = DataConfig()
    pod: PODConfig = PODConfig()
    deriv: DerivConfig = DerivConfig()
    sindy: SINDyConfig = SINDyConfig()
    # rollout: RolloutConfig = RolloutConfig()
    # !! above bandaidpatch
    plots: PlotConfig = PlotConfig()

    outputs_dir: str | None = OUTPUTS_DIR.as_posix()

    def __post_init__(self):
        if self.outputs_dir is None:
            object.__setattr__(self, "outputs_dir", OUTPUTS_DIR.as_posix())

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the configuration to a plain dictionary.

        Returns
        -------
        dict
            Nested dictionary representation of the configuration.
        """
        return asdict(self)
