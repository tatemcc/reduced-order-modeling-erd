"""Configuration models and YAML helpers for the toy ERD FiPy plant.

This module centralizes all tunable parameters for plant simulation, forcing,
metrics, and run-folder outputs. It also provides preset constructors and
load/save helpers for YAML-backed experiments.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Literal, Optional

import yaml


PresetName = Literal["smoke", "report"]



def _repo_root(start: Optional[Path] = None) -> Path:
    """Locate the ERD-capstone workspace root.

    Args:
        start: Optional path used as the search anchor.

    Returns:
        Absolute path to the workspace root containing ``erd_fipy``, ``model``,
        and ``control`` package directories.
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
class DomainConfig:
    """Spatial domain and grid resolution.

    Attributes:
        R_min: Inner radius of the annular domain.
        R_max: Outer radius of the annular domain.
        N_r: Number of radial cells.
        N_phi: Number of azimuthal cells.
    """

    R_min: float = 0.3
    R_max: float = 1.0
    N_r: int = 48
    N_phi: int = 96


@dataclass(frozen=True)
class TimeConfig:
    """Time discretization and snapshot cadence.

    Attributes:
        dt: Fixed plant time step.
        T_final: Final physical simulation time.
        snapshot_every: Save one field snapshot every this many steps.
    """

    dt: float = 2.5e-3
    T_final: float = 6.0
    snapshot_every: int = 1

    @property
    def n_steps(self) -> int:
        """Number of explicit stepping iterations derived from ``dt`` and ``T_final``.

        Args:
            None.

        Returns:
            Integer number of simulation steps.
        """

        return int(round(self.T_final / self.dt))


@dataclass(frozen=True)
class PDEConfig:
    """Physical coefficients for the toy PDE system.

    Attributes:
        nu: Diffusion coefficient in the vorticity equation.
        gamma: Linear damping coefficient in the vorticity equation.
        D_r: Radial diffusion coefficient in the density equation.
        D_phi: Azimuthal diffusion coefficient in the density equation.
        alpha: Relaxation rate toward the target ring profile.
        beta_instab: Coupled instability gain multiplying ``b0(r) * (n - n_eq(r))``
            in the vorticity equation. Set to ``0`` to recover the previous
            baseline model.
    """

    nu: float = 2.5e-4
    gamma: float = 0.004
    D_r: float = 4.0e-5
    D_phi: float = 2.0e-5
    alpha: float = 0.001
    beta_instab: float = 3.4


@dataclass(frozen=True)
class RingConfig:
    """Target axisymmetric ring definition for ``n_eq(r)``.

    Attributes:
        n_bg: Background density level.
        n_amp: Gaussian ring amplitude above background.
        r_star: Target ring center radius.
        sigma_star: Ring Gaussian width.
    """

    n_bg: float = 0.3
    n_amp: float = 1.0
    r_star: float = 0.72
    sigma_star: float = 0.08


@dataclass(frozen=True)
class WallConfig:
    """Wall-loss layer parameters.

    Attributes:
        kappa_0: Peak wall sink strength.
        delta_w: Width of the wall-loss Gaussian layer.
    """

    kappa_0: float = 0.13
    delta_w: float = 0.06


@dataclass(frozen=True)
class ForcingConfig:
    """Actuation envelope widths and control bounds.

    Attributes:
        sigma_0: Width for axisymmetric drive envelope ``b0``.
        sigma_1: Width for ``m=1`` control envelopes.
        sigma_2: Width for ``m=2`` control envelopes.
        drive_u0_base: Default open-loop value for channel ``u0``.
        u_bounds: Symmetric hard bounds for channels ``[u0, ..., u4]``.
    """

    sigma_0: float = 0.11
    sigma_1: float = 0.10
    sigma_2: float = 0.12
    drive_u0_base: float = 2.6
    u_bounds: tuple[float, float, float, float, float] = (3.0, 2.8, 2.8, 2.4, 2.4)


@dataclass(frozen=True)
class DisturbanceModeConfig:
    """Deterministic multi-sine schedule for one disturbance mode.

    Attributes:
        amplitudes: Per-tone amplitudes.
        frequencies: Per-tone frequencies in Hz-equivalent units.
        phases: Per-tone initial phase offsets.
        phase_rate: Linear drift rate for the carrier phase.
        phase_offset: Static phase offset.
        phase_mod_amp: Amplitude of slow sinusoidal phase modulation.
        phase_mod_freq: Frequency of slow sinusoidal phase modulation.
    """

    amplitudes: tuple[float, ...] = (0.30, 0.14)
    frequencies: tuple[float, ...] = (0.19, 0.43)
    phases: tuple[float, ...] = (0.0, 1.2)
    phase_rate: float = 0.4
    phase_offset: float = 0.0
    phase_mod_amp: float = 0.2
    phase_mod_freq: float = 0.13


@dataclass(frozen=True)
class DisturbanceConfig:
    """Complete deterministic disturbance specification.

    Attributes:
        mode1: Multi-sine schedule for the ``m=1`` disturbance component.
        mode2: Multi-sine schedule for the ``m=2`` disturbance component.
        warmup_fraction: Fraction of run time in low-disturbance warmup.
        excite_fraction: Fraction of run time in high-excitation stage.
        warmup_scale: Disturbance amplitude multiplier in warmup stage.
        excite_scale: Disturbance amplitude multiplier in excitation stage.
        hold_scale: Disturbance amplitude multiplier in final hold stage.
        mode_jitter_amp_frac: Per-trajectory deterministic amplitude jitter.
        mode_jitter_phase: Per-trajectory deterministic phase jitter (radians).
        mode_jitter_freq_frac: Per-trajectory deterministic frequency jitter.
        phase_c_boost_freq: Additional rapid phase modulation frequency during
            the final stage to create visually stronger late-time dynamics.
        multiscale_modes: Additional mode numbers to excite for small/medium scales.
        multiscale_amplitudes: Base amplitudes for each mode in ``multiscale_modes``.
        multiscale_frequencies: Time frequencies for each extra mode.
        multiscale_strength: Global scaling for multiscale deterministic forcing.
        noise_strength: Band-limited pseudo-random forcing amplitude.
        noise_m_min: Minimum azimuthal mode index in pseudo-random forcing band.
        noise_m_max: Maximum azimuthal mode index in pseudo-random forcing band.
        noise_freq_base: Base temporal frequency for pseudo-random forcing.
    """

    mode1: DisturbanceModeConfig = DisturbanceModeConfig(
        amplitudes=(3.2, 1.6),
        frequencies=(0.20, 0.47),
        phases=(0.2, 1.4),
        phase_rate=0.62,
        phase_offset=0.0,
        phase_mod_amp=0.30,
        phase_mod_freq=0.15,
    )
    mode2: DisturbanceModeConfig = DisturbanceModeConfig(
        amplitudes=(2.8, 1.4),
        frequencies=(0.16, 0.41),
        phases=(1.1, -0.5),
        phase_rate=-0.36,
        phase_offset=0.7,
        phase_mod_amp=0.24,
        phase_mod_freq=0.18,
    )
    warmup_fraction: float = 0.0
    excite_fraction: float = 0.70
    warmup_scale: float = 2.2
    excite_scale: float = 2.40
    hold_scale: float = 3.00
    mode_jitter_amp_frac: float = 0.20
    mode_jitter_phase: float = 0.45
    mode_jitter_freq_frac: float = 0.08
    phase_c_boost_freq: float = 4.0
    multiscale_modes: tuple[int, ...] = (5, 8, 12, 16, 20)
    multiscale_amplitudes: tuple[float, ...] = (1.2, 1.0, 0.8, 0.55, 0.4)
    multiscale_frequencies: tuple[float, ...] = (0.31, 0.43, 0.57, 0.71, 0.89)
    multiscale_strength: float = 1.4
    noise_strength: float = 0.4
    noise_m_min: int = 14
    noise_m_max: int = 30
    noise_freq_base: float = 0.65


@dataclass(frozen=True)
class CostWeights:
    """Stage-cost weights used by random-shooting MPC.

    The wobble-energy metric is numerically small (typically ``1e-6`` to
    ``1e-4``), so its default weight is intentionally large to balance against
    control-power and wall-loss terms.

    Attributes:
        w_wob_abs: Weight on absolute wobble-energy penalty.
        w_wob_growth: Weight on wobble-energy growth above reference.
        w_wall: Weight on wall-loss term.
        w_coh: Weight on coherence/thickness term.
        w_coh_growth: Weight on ring-thickness growth above reference.
        w_wall_slope: Weight on positive wall-loss slope penalty.
        w_pow: Weight on control-power penalty.
        w_rate: Weight on control-rate penalty.
    """

    w_wob_abs: float = 5.0e10
    w_wob_growth: float = 2.0e11
    w_wall: float = 0.4
    w_coh: float = 0.3
    w_coh_growth: float = 0.6
    w_wall_slope: float = 0.15
    w_pow: float = 0.02
    w_rate: float = 0.15

    @property
    def w_wob(self) -> float:
        """Backward-compatible alias for legacy ``w_wob`` consumers.

        Args:
            None.

        Returns:
            Absolute wobble penalty weight.
        """

        return self.w_wob_abs


@dataclass(frozen=True)
class ControlConfig:
    """Controller hyperparameters shared between fitting and control loops.

    Attributes:
        H: Prediction horizon length in steps.
        N_shoot: Number of random shooting candidates per control step.
        rate_penalty: Additional multiplier on the rate-penalty term.
        weights: Structured stage-cost weights.
    """

    H: int = 6
    N_shoot: int = 24
    rate_penalty: float = 1.0
    weights: CostWeights = CostWeights()


@dataclass(frozen=True)
class MetricsConfig:
    """Metric post-processing constants.

    Attributes:
        sigma_w: Width of radial weighting used in wobble-energy metric.
        lambda_1: Control-power scaling for ``m=1`` channels.
        lambda_2: Control-power scaling for ``m=2`` channels.
        eta_window_steps: Window length (in steps) for efficiency metric.
        eta_eps: Numerical floor in efficiency denominator.
        wE: Reporting weight for ``E_wob_excess`` in ``badness_score``.
        wL: Reporting weight for ``L_w_cum`` in ``badness_score``.
        wS: Reporting weight for ``sigma_r_excess`` in ``badness_score``.
    """

    sigma_w: float = 0.1
    lambda_1: float = 1.0
    lambda_2: float = 1.2
    eta_window_steps: int = 20
    eta_eps: float = 1e-6
    wE: float = 1.0
    wL: float = 1.0
    wS: float = 1.0


@dataclass(frozen=True)
class InitConfig:
    """Initial-condition perturbation settings.

    Attributes:
        eps_n: Relative random perturbation amplitude for density.
        eps_omega: Random perturbation amplitude for vorticity.
        mode1_amp: Coherent ``m=1`` initial density perturbation amplitude.
        mode2_amp: Coherent ``m=2`` initial density perturbation amplitude.
        mode5_amp: Coherent ``m=5`` initial density perturbation amplitude.
        seed: Random seed for reproducible initial fields.
    """

    eps_n: float = 0.2
    eps_omega: float = 1.0
    mode1_amp: float = 0.82
    mode2_amp: float = 0.62
    mode5_amp: float = 0.38
    seed: int = 0


@dataclass(frozen=True)
class OutputConfig:
    """Run-folder output controls.

    Attributes:
        outputs_root: Parent directory for timestamped run folders.
        tag: Suffix tag appended to the timestamped run folder name.
        preset: Preset name used to initialize defaults.
        run_dir: Optional explicit output directory override.
    """

    outputs_root: str = DEFAULT_OUTPUTS_ROOT
    tag: str = "erd"
    preset: PresetName = "smoke"
    run_dir: Optional[str] = None


@dataclass(frozen=True)
class RunConfig:
    """Top-level configuration for a plant run.

    Attributes:
        name: Human-readable run name.
        domain: Spatial domain and mesh settings.
        time: Time stepping and snapshot settings.
        pde: PDE coefficient set.
        ring: Equilibrium ring parameters.
        wall: Wall sink parameters.
        forcing: Actuation envelopes and bounds.
        disturbance: Deterministic disturbance schedules.
        control: MPC hyperparameters.
        metrics: Metric post-processing constants.
        init: Initial-condition noise parameters.
        output: Output-directory and naming settings.
    """

    name: str = "erd_toy"
    domain: DomainConfig = DomainConfig()
    time: TimeConfig = TimeConfig()
    pde: PDEConfig = PDEConfig()
    ring: RingConfig = RingConfig()
    wall: WallConfig = WallConfig()
    forcing: ForcingConfig = ForcingConfig()
    disturbance: DisturbanceConfig = DisturbanceConfig()
    control: ControlConfig = ControlConfig()
    metrics: MetricsConfig = MetricsConfig()
    init: InitConfig = InitConfig()
    output: OutputConfig = OutputConfig()

    def to_dict(self) -> Dict[str, Any]:
        """Convert the full configuration tree into plain serializable dictionaries.

        Args:
            None.

        Returns:
            Nested dictionary representation of this run configuration.
        """

        return asdict(self)


@dataclass(frozen=True)
class DatasetSplitConfig:
    """Trajectory counts for train/validation/test dataset generation.

    Attributes:
        n_train: Number of training trajectories.
        n_val: Number of validation trajectories.
        n_test: Number of test trajectories.
    """

    n_train: int
    n_val: int
    n_test: int


@dataclass(frozen=True)
class PresetConfig:
    """Preset bundle containing run defaults and dataset split sizes.

    Attributes:
        run: Base run configuration for this preset.
        split: Dataset split cardinalities for this preset.
        block_steps: Piecewise-constant control block length for dataset generation.
    """

    run: RunConfig
    split: DatasetSplitConfig
    block_steps: int



def _merge(base: Dict[str, Any], updates: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge ``updates`` into ``base``.

    Args:
        base: Existing dictionary values.
        updates: New values that should override ``base``.

    Returns:
        Merged dictionary with recursive handling for nested mappings.
    """

    out = dict(base)
    for k, v in updates.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out



def _construct_run_config(payload: Dict[str, Any]) -> RunConfig:
    """Build a :class:`RunConfig` from a nested dictionary payload.

    Args:
        payload: Dictionary produced from YAML or another serialized source.

    Returns:
        Fully typed run configuration.
    """

    d = payload.get("domain", {})
    t = payload.get("time", {})
    p = payload.get("pde", {})
    r = payload.get("ring", {})
    w = payload.get("wall", {})
    f = payload.get("forcing", {})
    dist = payload.get("disturbance", {})
    c = payload.get("control", {})
    cw = c.get("weights", {})
    m = payload.get("metrics", {})
    init = payload.get("init", {})
    out = payload.get("output", {})

    mode1 = dist.get("mode1", {})
    mode2 = dist.get("mode2", {})

    default_control = ControlConfig()
    default_weights = CostWeights()
    default_dist = DisturbanceConfig()

    if "w_wob" in cw and "w_wob_abs" not in cw:
        cw = dict(cw)
        cw["w_wob_abs"] = cw["w_wob"]

    merged_weights = {
        "w_wob_abs": cw.get("w_wob_abs", default_weights.w_wob_abs),
        "w_wob_growth": cw.get("w_wob_growth", default_weights.w_wob_growth),
        "w_wall": cw.get("w_wall", default_weights.w_wall),
        "w_coh": cw.get("w_coh", default_weights.w_coh),
        "w_coh_growth": cw.get("w_coh_growth", default_weights.w_coh_growth),
        "w_wall_slope": cw.get("w_wall_slope", default_weights.w_wall_slope),
        "w_pow": cw.get("w_pow", default_weights.w_pow),
        "w_rate": cw.get("w_rate", default_weights.w_rate),
    }

    return RunConfig(
        name=payload.get("name", "erd_toy"),
        domain=DomainConfig(**d),
        time=TimeConfig(**t),
        pde=PDEConfig(**p),
        ring=RingConfig(**r),
        wall=WallConfig(**w),
        forcing=ForcingConfig(**f),
        disturbance=DisturbanceConfig(
            mode1=DisturbanceModeConfig(**mode1),
            mode2=DisturbanceModeConfig(**mode2),
            warmup_fraction=dist.get("warmup_fraction", default_dist.warmup_fraction),
            excite_fraction=dist.get("excite_fraction", default_dist.excite_fraction),
            warmup_scale=dist.get("warmup_scale", default_dist.warmup_scale),
            excite_scale=dist.get("excite_scale", default_dist.excite_scale),
            hold_scale=dist.get("hold_scale", default_dist.hold_scale),
            mode_jitter_amp_frac=dist.get("mode_jitter_amp_frac", default_dist.mode_jitter_amp_frac),
            mode_jitter_phase=dist.get("mode_jitter_phase", default_dist.mode_jitter_phase),
            mode_jitter_freq_frac=dist.get("mode_jitter_freq_frac", default_dist.mode_jitter_freq_frac),
            phase_c_boost_freq=dist.get("phase_c_boost_freq", default_dist.phase_c_boost_freq),
            multiscale_modes=tuple(dist.get("multiscale_modes", default_dist.multiscale_modes)),
            multiscale_amplitudes=tuple(dist.get("multiscale_amplitudes", default_dist.multiscale_amplitudes)),
            multiscale_frequencies=tuple(dist.get("multiscale_frequencies", default_dist.multiscale_frequencies)),
            multiscale_strength=dist.get("multiscale_strength", default_dist.multiscale_strength),
            noise_strength=dist.get("noise_strength", default_dist.noise_strength),
            noise_m_min=dist.get("noise_m_min", default_dist.noise_m_min),
            noise_m_max=dist.get("noise_m_max", default_dist.noise_m_max),
            noise_freq_base=dist.get("noise_freq_base", default_dist.noise_freq_base),
        ),
        control=ControlConfig(
            H=c.get("H", default_control.H),
            N_shoot=c.get("N_shoot", default_control.N_shoot),
            rate_penalty=c.get("rate_penalty", default_control.rate_penalty),
            weights=CostWeights(**merged_weights),
        ),
        metrics=MetricsConfig(**m),
        init=InitConfig(**init),
        output=OutputConfig(**out),
    )



def run_config_from_dict(payload: Dict[str, Any]) -> RunConfig:
    """Build a :class:`RunConfig` from an in-memory dictionary.

    Args:
        payload: Nested dictionary in the same schema as ``RunConfig.to_dict()``.

    Returns:
        Typed run configuration object.
    """

    return _construct_run_config(payload)



def preset_config(name: PresetName) -> PresetConfig:
    """Return a named preset with run defaults and dataset split sizes.

    Args:
        name: Preset name (``"smoke"`` or ``"report"``).

    Returns:
        Preset bundle with run config, split sizes, and control block length.
    """

    if name == "smoke":
        base = RunConfig(
            time=TimeConfig(dt=2.0e-3, T_final=0.9, snapshot_every=2),
            domain=DomainConfig(N_r=36, N_phi=72),
            control=ControlConfig(H=8, N_shoot=64),
            output=OutputConfig(tag="smoke", preset="smoke"),
        )
        return PresetConfig(run=base, split=DatasetSplitConfig(4, 0, 2), block_steps=7)

    if name == "report":
        base = RunConfig(
            time=TimeConfig(dt=1.5e-3, T_final=1.5, snapshot_every=4),
            domain=DomainConfig(N_r=48, N_phi=96),
            control=ControlConfig(
                H=12,
                N_shoot=128,
                rate_penalty=1.0,
                weights=CostWeights(
                    w_wob_abs=5.0e10,
                    w_wob_growth=2.0e11,
                    w_wall=0.4,
                    w_coh=0.3,
                    w_coh_growth=0.6,
                    w_wall_slope=0.15,
                    w_pow=0.02,
                    w_rate=0.15,
                ),
            ),
            output=OutputConfig(tag="report", preset="report"),
        )
        return PresetConfig(run=base, split=DatasetSplitConfig(12, 0, 4), block_steps=8)

    raise ValueError(f"Unknown preset: {name}")



def load_run_config(path: Optional[str | Path], preset: Optional[PresetName] = None) -> RunConfig:
    """Load a run configuration from YAML layered on top of a preset.

    Args:
        path: Optional path to a YAML file with overrides.
        preset: Optional preset name; defaults to ``"smoke"``.

    Returns:
        Final merged and typed run configuration.
    """

    base = preset_config(preset or "smoke").run.to_dict()
    if path is None:
        return _construct_run_config(base)

    p = Path(path)
    with p.open("r", encoding="utf-8") as f:
        payload = yaml.safe_load(f) or {}

    merged = _merge(base, payload)
    return _construct_run_config(merged)



def save_run_config(path: str | Path, cfg: RunConfig) -> None:
    """Serialize a run configuration as YAML.

    Args:
        path: Destination path for the YAML document.
        cfg: Run configuration instance to serialize.

    Returns:
        None.
    """

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        yaml.safe_dump(cfg.to_dict(), f, sort_keys=False)
