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
DynamicsModelName = Literal["mhw", "wave_landau"]



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
    """Physical coefficients for the annular modified Hasegawa-Wakatani plant.

    Attributes:
        dynamics_model: Active plant closure. ``"mhw"`` keeps the current
            modified-Hasegawa-Wakatani-style branch. ``"wave_landau"`` enables
            a more direct near-threshold propagating-wave scaffold with
            thresholded vorticity growth/saturation.
        C: Adiabaticity coupling coefficient multiplying ``tilde(phi) - tilde(n)``.
        kappa_0: Peak amplitude of the drift-wave drive profile ``kappa(r)``.
        sigma_kappa: Radial width of ``kappa(r)`` around the ring center.
        gradient_drive_gain: Weight of the self-consistent zonal density-gradient
            contribution to the drift-wave drive. When the thresholded branch
            is disabled this acts as the legacy additive zonal-gradient drive.
            When the thresholded branch is active it scales the excess
            supercritical drift-wave drive above the baseline profile.
        critical_gradient_ratio: Dimensionless critical zonal-gradient ratio
            relative to ``n_amp / sigma_star``. Values above zero activate the
            near-threshold branch where stochastic perturbations only grow once
            the ring profile steepens past this threshold.
        shear_suppression_gain: Strength of zonal-shear suppression in the
            instability margin. Larger values keep the ring coherent for longer
            by letting ExB shear quench drift-wave growth.
        shear_ref: Reference ExB shear magnitude used to nondimensionalize the
            zonal-shear suppression term.
        threshold_width: Smoothness of the critical-gradient activation. Smaller
            values make the branch behave more like a hard threshold, while
            larger values produce a softer onset of breakup.
        shear_damping_gain: Additional damping applied to nonzonal density and
            vorticity fluctuations in proportion to zonal ExB shear. This is a
            simple quasilinear surrogate for shear decorrelation near the
            marginal operating point.
        curvature_omega_gain: Strength of the density-to-vorticity curvature /
            interchange-style coupling. This lets small azimuthal density
            perturbations seed local vorticity growth instead of relying only on
            external forcing.
        baroclinic_omega_gain: Strength of a nonlinear baroclinic-style
            density-to-vorticity coupling proportional to radial and azimuthal
            density gradients. This is intended to generate more localized,
            bidirectional vortical structures than the separable curvature term
            alone.
        flux_balance_omega_gain: Strength of the flux-balanced Hasegawa-
            Wakatani-style radial-flux feedback in the vorticity equation. This
            term couples the evolving zonal density gradient to fluctuating
            radial transport, making bursts emerge from the state rather than
            from an imposed scripted carrier.
        phase_advection_gain: Gain for a supercritical diamagnetic-advection
            term acting on nonzonal density and vorticity. When the ring moves
            above threshold, this increases the azimuthal phase speed of
            unstable structures using the signed zonal-gradient direction,
            rather than adding more externally scripted forcing.
        supercritical_coupling_gain: Gain for strengthening the mHW
            fluctuation coupling term once the instability margin becomes
            positive. This makes supercritical regions exchange density and
            vorticity faster through the existing two-field physics instead of
            only through added source terms.
        supercritical_feedback_gain: Gain for a local positive feedback loop
            between nonzonal potential and density once the instability margin
            becomes positive. This is the simplest way to let marginally stable
            turbulence stay coherent until stochastic perturbations cross
            threshold, after which fluctuation amplitudes can grow quickly.
        supercritical_transport_gain: Gain for amplifying nonlinear advection
            once the instability margin becomes positive. This is a surrogate
            for the jump from coherent drift-wave propagation to fast convective
            transport and filamentation after threshold crossing.
        wave_density_coupling_gain: Strength of the ``∂_phi omega`` wave term
            in the density equation on the ``wave_landau`` branch.
        wave_vorticity_coupling_gain: Strength of the ``∂_phi n`` wave term in
            the vorticity equation on the ``wave_landau`` branch.
        wave_packet_phi_gain: Strength of the self-generated azimuthal packet
            speed on the ``wave_landau`` branch. Positive values let local
            density/vorticity structure set the direction and magnitude of
            wave propagation around the ring.
        wave_packet_r_gain: Strength of the self-generated radial packet speed
            on the ``wave_landau`` branch. This is the main lever for letting
            unstable packets peel away from the ring instead of remaining as
            radially aligned stripes.
        wave_burst_speed_gain: Multiplier that increases packet speeds once the
            local instability activation becomes positive. This formalizes the
            desired behavior that unstable structures do not just grow in
            amplitude, but also move faster and transport more aggressively.
        packet_mode_cut: Highest azimuthal Fourier mode retained in the
            packet-coupling channels on the ``wave_landau`` branch. This lets
            the plant keep fine-scale stochastic residue in ``omega`` while
            driving density transport from a smoother, more coherent packet
            field that is closer to the desired turbulent-ring morphology.
        omega_high_damp: Additional damping applied only to the high-mode
            vorticity residue above ``packet_mode_cut``. This is a targeted
            way to suppress late diagonal stripe blow-up without flattening the
            lower-mode packet dynamics that SINDy needs to predict.
        omega_landau_gain: Linear supercritical growth gain for nonzonal
            vorticity on the ``wave_landau`` branch.
        omega_landau_sat: Cubic saturation strength paired with
            ``omega_landau_gain``.
        density_landau_gain: Optional weaker supercritical growth gain for
            nonzonal density on the ``wave_landau`` branch.
        density_landau_sat: Cubic saturation strength paired with
            ``density_landau_gain``.
        density_vortex_gain: Strength of a bounded vorticity-to-density
            trapping/corrugation term on the ``wave_landau`` branch. This is
            the main lever for making density visibly inherit the propagating
            vortical packets instead of immediately relaxing back to a smooth
            axisymmetric ring.
        density_packet_drive_gain: Strength of a zero-mean packet-locked
            density drive on the ``wave_landau`` branch. This keeps density
            blobs attached to the propagating vortical packets after the
            instability turns on, which is critical for avoiding the observed
            postcritical collapse back to a diffuse axisymmetric ring.
        density_retention_gain: Additional zero-mean density entrainment gain
            used when vorticity packets remain strong but density nonaxisymmetry
            starts collapsing. This targets the observed late relaminarization
            failure without forcing early-time symmetry breaking.
        collapse_guard_gain: Strength of a late-stage annular-support term
            that preserves each azimuthal column's mass while pulling
            over-diffused density back toward the ring band. This is intended
            to stop the pathological whole-rectangle spill without flattening
            the earlier burst morphology.
        collapse_edge_threshold: Edge-mass-fraction threshold above which the
            annular-collapse guard begins to activate.
        collapse_core_threshold: Core-mass-fraction threshold below which the
            annular-collapse guard begins to activate.
        collapse_omega_damp: Targeted vorticity damping gain used only when
            the collapse guard activates. This damps edge-localized
            high-mode vorticity that would otherwise drive the late radial
            blow-up and subsequent relaminarization.
        radial_confinement_gain: Strength of an inward density-confinement
            flux centered on the ring. This keeps the operating point
            marginally stable in radius instead of letting the ring hit the
            radial walls immediately.
        radial_release_gain: Weakening of the confinement flux once the local
            instability activation turns positive. This formalizes the desired
            transition from coherent ring transport to radial breakup.
        radial_damage_release_gain: Additional weakening of the confinement
            flux once the density profile is already visibly degraded. This
            prevents the active branch from instantly re-confining after one
            burst and helps produce a more gradual degradation process.
        zonal_profile_restore_gain: Strength of an axisymmetric profile-support
            term that nudges the zonal density back toward ``n_eq(r)``. This
            is the main operating-point stabilizer for the wave branch because
            it anchors the mean ring without directly suppressing nonzonal wave
            packets.
        zonal_profile_release_gain: Weakening of the zonal profile-support term
            once the local instability activation turns positive. This lets the
            mean ring degrade only after disturbances genuinely become
            supercritical instead of remaining servoed to the operating point.
        zonal_profile_damage_release_gain: Additional weakening of the zonal
            profile-support term once the zonal ring profile is already
            damaged. This prevents the active branch from immediately
            re-axisymmetrizing after a burst and helps sustain degradation into
            the later part of the run.
        landau_phi_gain: Coefficient of a large-scale potential feedback term in
            the vorticity equation. Positive values keep the ring lively and
            fast-moving near threshold instead of relaxing toward a static tube.
        S_0: Peak amplitude of the particle source profile ``S(r)``.
        sigma_S: Radial width of ``S(r)`` around the ring center.
        nu_n: Dissipation coefficient in the density equation.
        nu_omega: Dissipation coefficient in the vorticity equation.
        hyper_p: Dissipation order. ``1`` gives Laplacian diffusion and ``2``
            gives dissipative biharmonic regularization.
        gamma_omega: Linear drag coefficient in the vorticity equation.
        source_floor_scale: Minimum fueling fraction applied to ``S(r)``.
        source_balance_gain: Gain that matches fueling to instantaneous wall loss.
        source_mass_gain: Gain that restores total ring mass toward its target.
        source_scale_max: Upper clamp on the dynamic fueling multiplier.
    """

    dynamics_model: DynamicsModelName = "mhw"
    C: float = 1.4
    kappa_0: float = 3.2
    sigma_kappa: float = 0.11
    gradient_drive_gain: float = 0.0
    critical_gradient_ratio: float = 0.0
    shear_suppression_gain: float = 0.0
    shear_ref: float = 1.0
    threshold_width: float = 0.08
    shear_damping_gain: float = 0.0
    curvature_omega_gain: float = 0.0
    baroclinic_omega_gain: float = 0.0
    flux_balance_omega_gain: float = 0.0
    phase_advection_gain: float = 0.0
    supercritical_coupling_gain: float = 0.0
    supercritical_feedback_gain: float = 0.0
    supercritical_transport_gain: float = 0.0
    wave_density_coupling_gain: float = 0.0
    wave_vorticity_coupling_gain: float = 0.0
    wave_packet_phi_gain: float = 0.0
    wave_packet_r_gain: float = 0.0
    wave_burst_speed_gain: float = 0.0
    packet_mode_cut: int = 0
    omega_high_damp: float = 0.0
    omega_landau_gain: float = 0.0
    omega_landau_sat: float = 0.30
    density_landau_gain: float = 0.0
    density_landau_sat: float = 0.30
    density_vortex_gain: float = 0.0
    density_packet_drive_gain: float = 0.0
    density_retention_gain: float = 0.0
    collapse_guard_gain: float = 0.0
    collapse_edge_threshold: float = 0.18
    collapse_core_threshold: float = 0.42
    collapse_omega_damp: float = 0.0
    radial_confinement_gain: float = 0.0
    radial_release_gain: float = 0.0
    radial_damage_release_gain: float = 0.0
    zonal_profile_restore_gain: float = 0.0
    zonal_profile_release_gain: float = 0.0
    zonal_profile_damage_release_gain: float = 0.0
    landau_phi_gain: float = 0.0
    S_0: float = 0.14
    sigma_S: float = 0.08
    nu_n: float = 1.8e-4
    nu_omega: float = 1.6e-4
    hyper_p: int = 1
    gamma_omega: float = 0.02
    source_floor_scale: float = 0.15
    source_balance_gain: float = 1.0
    source_mass_gain: float = 0.8
    source_scale_max: float = 4.0


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
        noise_freq_base: Legacy base frequency scale retained for backwards
            compatibility with older configs. New stochastic noise uses the
            correlation and refresh settings below.
        noise_corr_time: Correlation time for stochastic modal amplitudes.
        noise_refresh_dt: Time spacing used to precompute stochastic modal
            amplitudes before interpolation in ``t``.
        noise_mode_decay: Power-law bias toward low azimuthal modes in the
            stochastic forcing band. Larger values favor large-scale motion.
        phase_drift_strength: RMS stochastic phase-speed scale for low-mode
            carrier drift. Larger values promote earlier direction changes.
        phase_drift_corr_time: Correlation time of stochastic phase drift.
        eddy_strength: Global amplitude of localized stochastic vorticity packets.
        eddy_count: Number of simultaneously active eddy packets.
        eddy_sigma_r: Radial width of each eddy packet.
        eddy_sigma_phi: Azimuthal width of each eddy packet in radians.
        eddy_speed_std: RMS angular speed of drifting eddy packets.
        eddy_corr_time: Correlation time for eddy speeds and amplitudes.
        eddy_radial_jitter: RMS radial wandering of eddy centers about ``r_star``.
        eddy_geometry: Geometry used for localized stochastic eddies. The
            default ``fixed_pair`` reproduces the stronger historical branch,
            while ``rotating_pair`` rotates the dipole orientation stochastically
            to reduce scripted-looking packets.
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
    noise_corr_time: float = 0.06
    noise_refresh_dt: float = 0.012
    noise_mode_decay: float = 0.8
    phase_drift_strength: float = 0.0
    phase_drift_corr_time: float = 0.08
    eddy_strength: float = 0.0
    eddy_count: int = 0
    eddy_sigma_r: float = 0.06
    eddy_sigma_phi: float = 0.35
    eddy_speed_std: float = 1.2
    eddy_corr_time: float = 0.10
    eddy_radial_jitter: float = 0.06
    eddy_geometry: str = "fixed_pair"


@dataclass(frozen=True)
class CostWeights:
    """Stage-cost weights used by random-shooting MPC.

    Attributes:
        w_j: Weight on absolute profile-deviation penalty ``J_prof``.
        w_j_growth: Weight on growth of ``J_prof`` above its reference value.
        w_e: Weight on absolute low-mode energy penalty ``E_low``.
        w_e_growth: Weight on growth of ``E_low`` above its reference value.
        w_l: Weight on wall-loss term ``L_w``.
        w_sigma: Weight on radial-thickness penalty ``sigma_r^2``.
        w_u: Weight on control-power penalty.
        w_delta_u: Weight on control-rate penalty.
    """

    w_j: float = 18.0
    w_j_growth: float = 36.0
    w_e: float = 4.0e4
    w_e_growth: float = 8.0e4
    w_l: float = 1.0
    w_sigma: float = 6.0
    w_u: float = 0.02
    w_delta_u: float = 0.15


@dataclass(frozen=True)
class ControlConfig:
    """Controller hyperparameters shared between fitting and control loops.

    Attributes:
        H: Prediction horizon length in steps.
        N_shoot: Number of random shooting candidates per control step.
        shoot_segments: Number of piecewise-constant segments per candidate.
        shoot_seg_len: Number of horizon steps in each piecewise-constant segment.
        rate_penalty: Additional multiplier on the rate-penalty term.
        weights: Structured stage-cost weights.
    """

    H: int = 6
    N_shoot: int = 24
    shoot_segments: int = 3
    shoot_seg_len: int = 2
    rate_penalty: float = 1.0
    weights: CostWeights = CostWeights()


@dataclass(frozen=True)
class MetricsConfig:
    """Metric post-processing constants.

    Attributes:
        sigma_w: Width of radial weighting used in the low-mode energy metric.
        lambda_1: Control-power scaling for ``m=1`` channels.
        lambda_2: Control-power scaling for ``m=2`` channels.
        eta_window_steps: Window length (in steps) for efficiency metric.
        eta_eps: Numerical floor in efficiency denominator.
        wJ: Reporting weight for ``J_prof_excess`` in ``badness_score``.
        wE: Reporting weight for ``E_low_excess`` in ``badness_score``.
        wL: Reporting weight for ``L_w_cum`` in ``badness_score``.
        wS: Reporting weight for ``sigma_r_excess`` in ``badness_score``.
        temporal_activity_baseline: Reference median frame-to-frame change used
            when no external baseline run is supplied for variation diagnostics.
    """

    sigma_w: float = 0.1
    lambda_1: float = 1.0
    lambda_2: float = 1.2
    eta_window_steps: int = 20
    eta_eps: float = 1e-6
    wJ: float = 1.0
    wE: float = 1.0
    wL: float = 1.0
    wS: float = 1.0
    temporal_activity_baseline: float = 7.5e-4


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

    if "w_wob" in cw and "w_e" not in cw:
        cw = dict(cw)
        cw["w_e"] = cw["w_wob"]
    if "w_wob_abs" in cw and "w_e" not in cw:
        cw = dict(cw)
        cw["w_e"] = cw["w_wob_abs"]
    if "w_wob_growth" in cw and "w_e_growth" not in cw:
        cw = dict(cw)
        cw["w_e_growth"] = cw["w_wob_growth"]
    if "w_wall" in cw and "w_l" not in cw:
        cw = dict(cw)
        cw["w_l"] = cw["w_wall"]
    if "w_coh" in cw and "w_sigma" not in cw:
        cw = dict(cw)
        cw["w_sigma"] = cw["w_coh"]
    if "w_pow" in cw and "w_u" not in cw:
        cw = dict(cw)
        cw["w_u"] = cw["w_pow"]
    if "w_rate" in cw and "w_delta_u" not in cw:
        cw = dict(cw)
        cw["w_delta_u"] = cw["w_rate"]

    merged_weights = {
        "w_j": cw.get("w_j", default_weights.w_j),
        "w_j_growth": cw.get("w_j_growth", default_weights.w_j_growth),
        "w_e": cw.get("w_e", default_weights.w_e),
        "w_e_growth": cw.get("w_e_growth", default_weights.w_e_growth),
        "w_l": cw.get("w_l", default_weights.w_l),
        "w_sigma": cw.get("w_sigma", default_weights.w_sigma),
        "w_u": cw.get("w_u", default_weights.w_u),
        "w_delta_u": cw.get("w_delta_u", default_weights.w_delta_u),
    }

    pde_keys = set(PDEConfig.__annotations__.keys())
    filtered_p = {k: v for k, v in p.items() if k in pde_keys}

    return RunConfig(
        name=payload.get("name", "erd_toy"),
        domain=DomainConfig(**d),
        time=TimeConfig(**t),
        pde=PDEConfig(**filtered_p),
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
            noise_corr_time=dist.get("noise_corr_time", default_dist.noise_corr_time),
            noise_refresh_dt=dist.get("noise_refresh_dt", default_dist.noise_refresh_dt),
            noise_mode_decay=dist.get("noise_mode_decay", default_dist.noise_mode_decay),
            phase_drift_strength=dist.get("phase_drift_strength", default_dist.phase_drift_strength),
            phase_drift_corr_time=dist.get("phase_drift_corr_time", default_dist.phase_drift_corr_time),
            eddy_strength=dist.get("eddy_strength", default_dist.eddy_strength),
            eddy_count=dist.get("eddy_count", default_dist.eddy_count),
            eddy_sigma_r=dist.get("eddy_sigma_r", default_dist.eddy_sigma_r),
            eddy_sigma_phi=dist.get("eddy_sigma_phi", default_dist.eddy_sigma_phi),
            eddy_speed_std=dist.get("eddy_speed_std", default_dist.eddy_speed_std),
            eddy_corr_time=dist.get("eddy_corr_time", default_dist.eddy_corr_time),
            eddy_radial_jitter=dist.get("eddy_radial_jitter", default_dist.eddy_radial_jitter),
            eddy_geometry=dist.get("eddy_geometry", default_dist.eddy_geometry),
        ),
        control=ControlConfig(
            H=c.get("H", default_control.H),
            N_shoot=c.get("N_shoot", default_control.N_shoot),
            shoot_segments=c.get("shoot_segments", default_control.shoot_segments),
            shoot_seg_len=c.get("shoot_seg_len", default_control.shoot_seg_len),
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
            time=TimeConfig(dt=1.5e-3, T_final=0.9, snapshot_every=2),
            domain=DomainConfig(N_r=36, N_phi=72),
            ring=RingConfig(n_bg=0.0, n_amp=1.0, r_star=0.72, sigma_star=0.08),
            wall=WallConfig(kappa_0=1.0, delta_w=0.06),
            pde=PDEConfig(
                C=1.1,
                kappa_0=4.8,
                sigma_kappa=0.11,
                S_0=0.10,
                sigma_S=0.16,
                nu_n=6.0e-6,
                nu_omega=6.0e-6,
                hyper_p=2,
                gamma_omega=0.02,
                source_floor_scale=0.20,
                source_balance_gain=1.0,
                source_mass_gain=1.1,
                source_scale_max=4.5,
            ),
            disturbance=DisturbanceConfig(
                mode1=DisturbanceModeConfig(
                    amplitudes=(4.4, 2.5),
                    frequencies=(0.20, 0.47),
                    phases=(0.2, 1.4),
                    phase_rate=0.62,
                    phase_offset=0.0,
                    phase_mod_amp=0.30,
                    phase_mod_freq=0.15,
                ),
                mode2=DisturbanceModeConfig(
                    amplitudes=(3.9, 2.3),
                    frequencies=(0.16, 0.41),
                    phases=(1.1, -0.5),
                    phase_rate=-0.36,
                    phase_offset=0.7,
                    phase_mod_amp=0.24,
                    phase_mod_freq=0.18,
                ),
                warmup_fraction=0.0,
                excite_fraction=0.70,
                warmup_scale=1.8,
                excite_scale=4.8,
                hold_scale=5.4,
                mode_jitter_amp_frac=0.20,
                mode_jitter_phase=0.45,
                mode_jitter_freq_frac=0.08,
                phase_c_boost_freq=4.0,
                multiscale_modes=(5, 8, 12, 16, 20),
                multiscale_amplitudes=(1.3, 1.1, 0.9, 0.65, 0.45),
                multiscale_frequencies=(0.31, 0.43, 0.57, 0.71, 0.89),
                multiscale_strength=1.75,
                noise_strength=0.45,
                noise_m_min=14,
                noise_m_max=30,
                noise_freq_base=0.65,
                phase_drift_strength=0.30,
                phase_drift_corr_time=0.09,
                eddy_strength=0.20,
                eddy_count=4,
                eddy_sigma_r=0.05,
                eddy_sigma_phi=0.25,
                eddy_speed_std=1.0,
                eddy_corr_time=0.10,
                eddy_radial_jitter=0.04,
            ),
            init=InitConfig(
                eps_n=0.05,
                eps_omega=0.9,
                mode1_amp=1.0,
                mode2_amp=0.70,
                mode5_amp=0.26,
                seed=0,
            ),
            control=ControlConfig(H=6, N_shoot=64, shoot_segments=3, shoot_seg_len=2),
            output=OutputConfig(tag="smoke", preset="smoke"),
        )
        return PresetConfig(run=base, split=DatasetSplitConfig(4, 0, 2), block_steps=8)

    if name == "report":
        base = RunConfig(
            time=TimeConfig(dt=1.5e-3, T_final=1.5, snapshot_every=4),
            domain=DomainConfig(N_r=48, N_phi=96),
            pde=PDEConfig(
                C=1.7,
                kappa_0=4.0,
                sigma_kappa=0.11,
                S_0=0.16,
                sigma_S=0.08,
                nu_n=1.1e-4,
                nu_omega=9.0e-5,
                hyper_p=1,
                gamma_omega=0.014,
                source_floor_scale=0.20,
                source_balance_gain=1.0,
                source_mass_gain=0.8,
                source_scale_max=4.0,
            ),
            control=ControlConfig(
                H=12,
                N_shoot=128,
                shoot_segments=4,
                shoot_seg_len=3,
                rate_penalty=1.0,
                weights=CostWeights(
                    w_j=18.0,
                    w_j_growth=36.0,
                    w_e=4.0e4,
                    w_e_growth=8.0e4,
                    w_l=1.0,
                    w_sigma=6.0,
                    w_u=0.02,
                    w_delta_u=0.15,
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
