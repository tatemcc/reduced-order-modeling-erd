"""Run a parameter scan for visually tuning the ERD plant.

The scan supports two modes:
- ``open_loop``: strict uncontrolled runs with a safer default variant family,
- ``random_piecewise``: training-like runs that can tolerate more aggressive
  parameter sweeps because the forcing schedule is richer.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from datetime import datetime
import json
import os
from pathlib import Path
import sys
from typing import Callable, Dict, Iterable, List, Sequence

import numpy as np
import yaml

if __package__ in (None, ""):
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)

from erd_fipy import ERDPlant, RunConfig, load_run_config


VariantBuilder = Callable[[RunConfig], tuple[RunConfig, Dict[str, object]]]


@dataclass(frozen=True)
class VariantSpec:
    """Named scan variant for one uncontrolled data run.

    Attributes:
        name: Stable identifier used in CLI filtering and run tags.
        description: Short explanation of what the variant is testing.
        build: Function that clones and modifies a base :class:`RunConfig`.
    """

    name: str
    description: str
    build: VariantBuilder


def build_parser(variant_names: Sequence[str]) -> argparse.ArgumentParser:
    """Build the CLI parser for the open-loop data scan.

    Args:
        variant_names: Available variant names for argument validation.

    Returns:
        Configured argument parser.
    """

    parser = argparse.ArgumentParser(
        description="Run a scan of aggressive uncontrolled ERD data variants with movies"
    )
    parser.add_argument("--config", type=str, default=None, help="Optional YAML override for the base plant config")
    parser.add_argument("--preset", type=str, default="report", choices=["smoke", "report"])
    parser.add_argument("--tag", type=str, default="single_data_scan", help="Human-friendly scan tag")
    parser.add_argument(
        "--variant",
        action="append",
        dest="variants",
        choices=list(variant_names),
        help="Run only the named variant; repeat to select several",
    )
    parser.add_argument("--limit", type=int, default=None, help="Optional limit on the number of variants to run")
    parser.add_argument("--manifest", type=str, default=None, help="Optional explicit scan manifest output path")
    parser.add_argument(
        "--outputs-root",
        type=str,
        default=None,
        help="Optional override for the run-folder root directory",
    )
    parser.add_argument(
        "--baseline-run",
        type=str,
        default=None,
        help="Optional baseline run directory used for variation diagnostics",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available scan variants and exit",
    )
    parser.add_argument(
        "--skip-media",
        action="store_true",
        help="Skip plots and GIFs for faster debugging of the scan script itself",
    )
    parser.add_argument(
        "--control-mode",
        choices=["open_loop", "random_piecewise"],
        default="open_loop",
        help="Use the default open-loop actuation or a training-style random piecewise schedule.",
    )
    parser.add_argument("--seed", type=int, default=0, help="Base RNG seed used for random-piecewise scans.")
    parser.add_argument(
        "--block-steps",
        type=int,
        default=8,
        help="Control block length for random-piecewise scans.",
    )
    return parser


def _mode_update(
    cfg_mode,
    *,
    amplitudes: Iterable[float] | None = None,
    frequencies: Iterable[float] | None = None,
    phases: Iterable[float] | None = None,
    phase_rate: float | None = None,
    phase_offset: float | None = None,
    phase_mod_amp: float | None = None,
    phase_mod_freq: float | None = None,
):
    """Return a modified disturbance-mode config from a base mode object."""

    updates = {}
    if amplitudes is not None:
        updates["amplitudes"] = tuple(float(v) for v in amplitudes)
    if frequencies is not None:
        updates["frequencies"] = tuple(float(v) for v in frequencies)
    if phases is not None:
        updates["phases"] = tuple(float(v) for v in phases)
    if phase_rate is not None:
        updates["phase_rate"] = float(phase_rate)
    if phase_offset is not None:
        updates["phase_offset"] = float(phase_offset)
    if phase_mod_amp is not None:
        updates["phase_mod_amp"] = float(phase_mod_amp)
    if phase_mod_freq is not None:
        updates["phase_mod_freq"] = float(phase_mod_freq)
    return replace(cfg_mode, **updates)


def _random_piecewise_schedule(
    rng: np.random.Generator,
    bounds: np.ndarray,
    base_u0: float,
    block_steps: int,
    t_final: float,
    warmup_fraction: float,
    excite_fraction: float,
):
    """Create a training-style piecewise-constant control schedule.

    The schedule intentionally keeps the axisymmetric ``u0`` channel active
    while asymmetry channels remain persistently exciting but somewhat gentler
    than their hard bounds, so turbulence and advection are not swamped by
    actuator noise.
    """

    current = np.array([base_u0, 0.0, 0.0, 0.0, 0.0], dtype=float)

    def schedule(k: int, t: float, _n: np.ndarray, _omega: np.ndarray) -> np.ndarray:
        nonlocal current
        if (k % max(1, block_steps)) == 0:
            sample = rng.uniform(-bounds, bounds)
            frac = float(np.clip(t / max(t_final, 1e-12), 0.0, 1.0))
            warm = float(np.clip(warmup_fraction, 0.0, 1.0))
            excite = float(np.clip(excite_fraction, 0.0, 1.0 - warm))

            if frac < warm:
                asym_scale = 0.18
                u0_scale = 0.55
            elif frac < warm + excite:
                asym_scale = 0.50
                u0_scale = 0.85
            else:
                asym_scale = 0.65
                u0_scale = 0.75

            sample[1:] = np.clip(sample[1:] * asym_scale, -bounds[1:], bounds[1:])
            sample[0] = np.clip(base_u0 + u0_scale * sample[0], -bounds[0], bounds[0])
            current = sample
        return current

    return schedule


def _deterministic_open_loop_schedule(
    *,
    base_u0: float,
    bound_u0: float,
    t_final: float,
    warmup_fraction: float,
    excite_fraction: float,
):
    """Create a deterministic open-loop ``u0(t)`` schedule with zero asymmetry channels.

    The scan uses this schedule instead of a constant ``u0`` so open-loop runs
    remain informative and numerically better behaved on the turbulent plant.

    Args:
        base_u0: Nominal axisymmetric drive level.
        bound_u0: Symmetric hard bound for the ``u0`` channel.
        t_final: Final simulation time.
        warmup_fraction: Fraction of run spent in the low-shear warmup.
        excite_fraction: Fraction of run spent in the stronger drive stage.

    Returns:
        Callback compatible with :meth:`ERDPlant.run`.
    """

    def schedule(_k: int, t: float, _n: np.ndarray, _omega: np.ndarray) -> np.ndarray:
        frac = float(np.clip(t / max(t_final, 1e-12), 0.0, 1.0))
        warm = float(np.clip(warmup_fraction, 0.0, 1.0))
        excite = float(np.clip(excite_fraction, 0.0, 1.0 - warm))
        low_base_reversal = base_u0 <= 0.40 * bound_u0
        if frac < warm:
            if low_base_reversal:
                amp = 0.30 * bound_u0
                carrier = 0.85 * np.sin(2.0 * np.pi * 0.20 * t)
            else:
                amp = 0.18 * bound_u0
                carrier = np.sin(2.0 * np.pi * 0.32 * t)
            u0 = base_u0 + amp * carrier
        elif frac < warm + excite:
            if low_base_reversal:
                amp = 0.62 * bound_u0
                carrier = (
                    0.80 * np.sin(2.0 * np.pi * 0.18 * t)
                    + 0.20 * np.sin(2.0 * np.pi * 0.49 * t + 0.6)
                )
            else:
                amp = 0.32 * bound_u0
                carrier = (
                    0.7 * np.sin(2.0 * np.pi * 0.41 * t)
                    + 0.3 * np.sin(2.0 * np.pi * 0.93 * t + 0.4)
                )
            u0 = base_u0 + amp * carrier
        else:
            if low_base_reversal:
                amp = 0.52 * bound_u0
                carrier = (
                    0.85 * np.sin(2.0 * np.pi * 0.16 * t + 0.3)
                    + 0.15 * np.sin(2.0 * np.pi * 0.37 * t)
                )
            else:
                amp = 0.24 * bound_u0
                carrier = (
                    0.6 * np.sin(2.0 * np.pi * 0.27 * t + 0.3)
                    + 0.4 * np.sin(2.0 * np.pi * 0.74 * t)
                )
            u0 = base_u0 + amp * carrier
        return np.array([np.clip(u0, -bound_u0, bound_u0), 0.0, 0.0, 0.0, 0.0], dtype=float)

    return schedule


def _stabilize_open_loop_base(cfg: RunConfig) -> RunConfig:
    """Retune the base config so strict open-loop scans stay finite.

    Open loop only applies the constant ``u0`` drive, so it cannot tolerate the
    same aggressive drift/noise settings as the training-style random-piecewise
    schedule. This helper moves the base into a smoke-like, hyperdiffusive
    regime that still develops visible transport and symmetry breaking.

    Args:
        cfg: Loaded base run configuration from the chosen preset.

    Returns:
        Safer base configuration for open-loop visual scans.
    """

    dist = replace(
        cfg.disturbance,
        mode1=_mode_update(
            cfg.disturbance.mode1,
            amplitudes=(0.50, 0.16),
            frequencies=(0.20, 0.46),
            phases=(0.15, 1.35),
            phase_rate=0.50,
            phase_mod_amp=0.18,
            phase_mod_freq=0.14,
        ),
        mode2=_mode_update(
            cfg.disturbance.mode2,
            amplitudes=(0.30, 0.10),
            frequencies=(0.16, 0.40),
            phases=(1.05, -0.45),
            phase_rate=-0.28,
            phase_offset=0.65,
            phase_mod_amp=0.16,
            phase_mod_freq=0.17,
        ),
        warmup_fraction=0.08,
        excite_fraction=0.60,
        warmup_scale=0.45,
        excite_scale=2.4,
        hold_scale=2.8,
        multiscale_modes=(6, 10, 15),
        multiscale_amplitudes=(0.35, 0.24, 0.16),
        multiscale_frequencies=(0.27, 0.43, 0.71),
        multiscale_strength=0.35,
        noise_strength=1.00,
        noise_m_min=1,
        noise_m_max=18,
        noise_freq_base=0.55,
        noise_corr_time=0.045,
        noise_refresh_dt=0.009,
        noise_mode_decay=0.95,
    )
    return replace(
        cfg,
        domain=replace(cfg.domain, N_r=min(int(cfg.domain.N_r), 36), N_phi=min(int(cfg.domain.N_phi), 72)),
        time=replace(
            cfg.time,
            dt=min(float(cfg.time.dt), 1.5e-3),
            T_final=min(float(cfg.time.T_final), 1.2),
            snapshot_every=max(4, int(cfg.time.snapshot_every)),
        ),
        ring=replace(cfg.ring, n_bg=0.0, n_amp=1.0, r_star=0.72, sigma_star=0.08),
        wall=replace(cfg.wall, kappa_0=1.0, delta_w=0.06),
        pde=replace(
            cfg.pde,
            C=1.0,
            kappa_0=4.2,
            sigma_kappa=0.11,
            S_0=0.09,
            sigma_S=0.16,
            nu_n=max(float(cfg.pde.nu_n), 2.0e-5),
            nu_omega=max(float(cfg.pde.nu_omega), 2.0e-5),
            hyper_p=2,
            gamma_omega=max(float(cfg.pde.gamma_omega), 0.028),
        ),
        forcing=replace(cfg.forcing, drive_u0_base=2.85),
        disturbance=dist,
        init=replace(cfg.init, eps_n=0.0, eps_omega=0.45, mode1_amp=0.0, mode2_amp=0.0, mode5_amp=0.0),
    )


def _with_symmetric_density_start(cfg: RunConfig, *, eps_omega: float) -> tuple[RunConfig, Dict[str, object]]:
    """Remove density asymmetry at ``t=0`` while keeping vorticity perturbations."""

    init_cfg = replace(
        cfg.init,
        eps_n=0.0,
        eps_omega=float(eps_omega),
        mode1_amp=0.0,
        mode2_amp=0.0,
        mode5_amp=0.0,
    )
    return replace(cfg, init=init_cfg), {
        "init.eps_n": 0.0,
        "init.eps_omega": float(eps_omega),
        "init.mode1_amp": 0.0,
        "init.mode2_amp": 0.0,
        "init.mode5_amp": 0.0,
    }


def _with_minimal_seed(cfg: RunConfig, *, eps_omega: float) -> tuple[RunConfig, Dict[str, object]]:
    """Use only weak random vorticity seeding and no coherent density seed.

    Args:
        cfg: Base run configuration to modify.
        eps_omega: Random vorticity seed amplitude.

    Returns:
        Modified configuration and a change log dictionary.
    """

    init_cfg = replace(
        cfg.init,
        eps_n=0.0,
        eps_omega=float(eps_omega),
        mode1_amp=0.0,
        mode2_amp=0.0,
        mode5_amp=0.0,
    )
    return replace(cfg, init=init_cfg), {
        "init.eps_n": 0.0,
        "init.eps_omega": float(eps_omega),
        "init.mode1_amp": 0.0,
        "init.mode2_amp": 0.0,
        "init.mode5_amp": 0.0,
    }


def _variant_base_reference(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Keep the loaded configuration unchanged as the comparison point."""

    return base, {}


def _variant_symmetric_emergence(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Delay visible symmetry breaking so the movie starts from a clean ring."""

    cfg, changes = _with_symmetric_density_start(base, eps_omega=1.15)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=1.8, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.90,
            kappa_0=4.0,
            S_0=0.08,
            nu_n=1.2e-5,
            nu_omega=1.2e-5,
            hyper_p=2,
            gamma_omega=0.03,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=2.9),
        disturbance=replace(
            cfg.disturbance,
            warmup_fraction=0.18,
            excite_fraction=0.55,
            warmup_scale=0.15,
            excite_scale=3.2,
            hold_scale=4.2,
            multiscale_modes=(5, 8, 12, 18, 24),
            multiscale_amplitudes=(1.5, 1.25, 1.0, 0.75, 0.55),
            multiscale_frequencies=(0.24, 0.36, 0.51, 0.69, 0.92),
            multiscale_strength=1.9,
            noise_strength=0.7,
            noise_m_min=14,
            noise_m_max=32,
            noise_freq_base=0.9,
        ),
    )
    changes.update(
        {
            "time.T_final": 1.8,
            "time.snapshot_every": cfg.time.snapshot_every,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.S_0": cfg.pde.S_0,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.hyper_p": cfg.pde.hyper_p,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "disturbance.warmup_fraction": cfg.disturbance.warmup_fraction,
            "disturbance.warmup_scale": cfg.disturbance.warmup_scale,
            "disturbance.excite_scale": cfg.disturbance.excite_scale,
            "disturbance.hold_scale": cfg.disturbance.hold_scale,
            "disturbance.multiscale_modes": list(cfg.disturbance.multiscale_modes),
            "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
        }
    )
    return cfg, changes


def _variant_stochastic_advection_breakup(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Favor stochastic symmetry breaking with stronger advection and minimal seeding."""

    cfg, changes = _with_minimal_seed(base, eps_omega=0.70)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=1.8, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.90,
            kappa_0=4.2,
            S_0=0.08,
            nu_n=1.2e-5,
            nu_omega=1.2e-5,
            hyper_p=2,
            gamma_omega=0.028,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=3.05),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(0.90, 0.30),
                frequencies=(0.18, 0.34),
                phases=(0.1, 1.2),
                phase_rate=0.22,
                phase_mod_amp=0.10,
                phase_mod_freq=0.09,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(0.55, 0.18),
                frequencies=(0.16, 0.30),
                phases=(0.8, -0.3),
                phase_rate=-0.14,
                phase_mod_amp=0.08,
                phase_mod_freq=0.08,
            ),
            warmup_fraction=0.06,
            excite_fraction=0.64,
            warmup_scale=0.30,
            excite_scale=3.6,
            hold_scale=4.4,
            multiscale_modes=(6, 10, 14, 20),
            multiscale_amplitudes=(0.35, 0.25, 0.18, 0.10),
            multiscale_frequencies=(0.23, 0.39, 0.66, 0.92),
            multiscale_strength=0.55,
            noise_strength=3.0,
            noise_m_min=1,
            noise_m_max=26,
            noise_freq_base=0.45,
            noise_corr_time=0.035,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.65,
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "disturbance.mode1.amplitudes": list(cfg.disturbance.mode1.amplitudes),
            "disturbance.mode2.amplitudes": list(cfg.disturbance.mode2.amplitudes),
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.noise_band": [cfg.disturbance.noise_m_min, cfg.disturbance.noise_m_max],
            "disturbance.noise_corr_time": cfg.disturbance.noise_corr_time,
            "disturbance.noise_mode_decay": cfg.disturbance.noise_mode_decay,
            "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
        }
    )
    return cfg, changes


def _variant_hybrid_stochastic_shear(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Blend weak coherent scaffold with stronger stochastic mid/high-band forcing.

    The intent is to preserve the strong progression speed of
    ``symmetric_emergence`` while making the disturbance look less phase-scripted.
    Large-scale structure is still supported by low-mode coherent forcing, but
    smaller scales are driven primarily by seeded stochastic content.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=1.10)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=1.8, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.89,
            kappa_0=4.15,
            S_0=0.08,
            nu_n=1.05e-5,
            nu_omega=1.05e-5,
            hyper_p=2,
            gamma_omega=0.027,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=3.10),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(1.30, 0.52),
                frequencies=(0.20, 0.41),
                phases=(0.18, 1.25),
                phase_rate=0.34,
                phase_mod_amp=0.10,
                phase_mod_freq=0.10,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(0.85, 0.34),
                frequencies=(0.16, 0.34),
                phases=(0.95, -0.35),
                phase_rate=-0.24,
                phase_offset=0.55,
                phase_mod_amp=0.08,
                phase_mod_freq=0.09,
            ),
            warmup_fraction=0.12,
            excite_fraction=0.58,
            warmup_scale=0.20,
            excite_scale=3.45,
            hold_scale=4.25,
            multiscale_modes=(6, 10, 14, 20),
            multiscale_amplitudes=(0.40, 0.28, 0.18, 0.12),
            multiscale_frequencies=(0.24, 0.36, 0.53, 0.78),
            multiscale_strength=0.45,
            noise_strength=2.60,
            noise_m_min=1,
            noise_m_max=22,
            noise_freq_base=0.70,
            noise_corr_time=0.055,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.18,
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "disturbance.mode1.amplitudes": list(cfg.disturbance.mode1.amplitudes),
            "disturbance.mode2.amplitudes": list(cfg.disturbance.mode2.amplitudes),
            "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.noise_band": [cfg.disturbance.noise_m_min, cfg.disturbance.noise_m_max],
            "disturbance.noise_corr_time": cfg.disturbance.noise_corr_time,
            "disturbance.noise_mode_decay": cfg.disturbance.noise_mode_decay,
        }
    )
    return cfg, changes


def _variant_textured_emergence(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Keep the old strong progression but add stochastic texture on top.

    This variant intentionally stays close to ``symmetric_emergence`` because
    that is still the strongest visual baseline. The only goal here is to make
    the forcing look less synthetic by layering seeded stochastic structure on
    top of the existing multiscale scaffold, not to replace the scaffold.
    """

    cfg, changes = _with_symmetric_density_start(base, eps_omega=1.10)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=1.8, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.90,
            kappa_0=4.0,
            S_0=0.08,
            nu_n=1.15e-5,
            nu_omega=1.15e-5,
            hyper_p=2,
            gamma_omega=0.029,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=2.95),
        disturbance=replace(
            cfg.disturbance,
            warmup_fraction=0.16,
            excite_fraction=0.56,
            warmup_scale=0.15,
            excite_scale=3.2,
            hold_scale=4.1,
            multiscale_modes=(5, 8, 12, 18, 24),
            multiscale_amplitudes=(1.25, 1.05, 0.85, 0.62, 0.44),
            multiscale_frequencies=(0.24, 0.36, 0.51, 0.69, 0.92),
            multiscale_strength=1.55,
            noise_strength=1.30,
            noise_m_min=4,
            noise_m_max=28,
            noise_freq_base=0.85,
            noise_corr_time=0.045,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.55,
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
            "disturbance.multiscale_amplitudes": list(cfg.disturbance.multiscale_amplitudes),
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.noise_band": [cfg.disturbance.noise_m_min, cfg.disturbance.noise_m_max],
            "disturbance.noise_corr_time": cfg.disturbance.noise_corr_time,
            "disturbance.noise_mode_decay": cfg.disturbance.noise_mode_decay,
        }
    )
    return cfg, changes


def _variant_azimuthal_shear_texture(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Bias the plant toward stronger azimuthal shear with stochastic texture.

    This variant is aimed at the current gap in the visuals: the ring broadens
    and breaks up, but too much of the visible transport still feels radial.
    The tuning keeps the symmetric density start, increases the axisymmetric
    drive that produces ``u_phi``, and shifts the stochastic band toward lower
    modes so the azimuthal transport looks less scripted than the older
    coherent-carrier variants.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=1.25)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=1.9, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.82,
            kappa_0=4.55,
            S_0=0.10,
            nu_n=9.0e-6,
            nu_omega=1.0e-5,
            hyper_p=2,
            gamma_omega=0.024,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=3.35, u_bounds=(3.8, 2.8, 2.8, 2.4, 2.4)),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(0.95, 0.34),
                frequencies=(0.16, 0.30),
                phases=(0.10, 1.10),
                phase_rate=0.14,
                phase_mod_amp=0.06,
                phase_mod_freq=0.08,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(0.55, 0.18),
                frequencies=(0.14, 0.28),
                phases=(0.72, -0.22),
                phase_rate=-0.10,
                phase_offset=0.40,
                phase_mod_amp=0.05,
                phase_mod_freq=0.07,
            ),
            warmup_fraction=0.04,
            excite_fraction=0.62,
            warmup_scale=1.00,
            excite_scale=3.4,
            hold_scale=4.0,
            multiscale_modes=(4, 7, 11, 16, 22),
            multiscale_amplitudes=(0.82, 0.64, 0.46, 0.28, 0.16),
            multiscale_frequencies=(0.20, 0.31, 0.45, 0.63, 0.84),
            multiscale_strength=1.15,
            noise_strength=2.8,
            noise_m_min=1,
            noise_m_max=20,
            noise_freq_base=0.55,
            noise_corr_time=0.075,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.10,
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.S_0": cfg.pde.S_0,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "forcing.u_bounds": list(cfg.forcing.u_bounds),
            "disturbance.mode1.amplitudes": list(cfg.disturbance.mode1.amplitudes),
            "disturbance.mode2.amplitudes": list(cfg.disturbance.mode2.amplitudes),
            "disturbance.multiscale_modes": list(cfg.disturbance.multiscale_modes),
            "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.noise_band": [cfg.disturbance.noise_m_min, cfg.disturbance.noise_m_max],
            "disturbance.noise_corr_time": cfg.disturbance.noise_corr_time,
            "disturbance.noise_mode_decay": cfg.disturbance.noise_mode_decay,
        }
    )
    return cfg, changes


def _variant_swirl_rich_texture(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Favor stronger circulating shear while keeping the disturbance stochastic.

    This branch starts from the azimuthal-shear variant and pushes the levers
    that most directly strengthen visible swirl:
    - stronger axisymmetric shear drive,
    - slightly lower damping and diffusion,
    - low-mode-biased stochastic forcing with longer correlation time.

    The intent is to make all fields rotate and shear more clearly without
    reverting to a heavily scripted multi-carrier disturbance.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=1.45)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=2.0, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.78,
            kappa_0=4.9,
            S_0=0.11,
            nu_n=8.5e-6,
            nu_omega=9.0e-6,
            hyper_p=2,
            gamma_omega=0.020,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=3.60, u_bounds=(4.4, 2.8, 2.8, 2.4, 2.4)),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(1.10, 0.42),
                frequencies=(0.15, 0.28),
                phases=(0.10, 1.00),
                phase_rate=0.18,
                phase_mod_amp=0.07,
                phase_mod_freq=0.08,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(0.68, 0.24),
                frequencies=(0.13, 0.25),
                phases=(0.65, -0.18),
                phase_rate=-0.12,
                phase_offset=0.35,
                phase_mod_amp=0.06,
                phase_mod_freq=0.07,
            ),
            warmup_fraction=0.02,
            excite_fraction=0.64,
            warmup_scale=1.20,
            excite_scale=3.8,
            hold_scale=4.4,
            multiscale_modes=(4, 6, 9, 13, 18),
            multiscale_amplitudes=(0.95, 0.76, 0.56, 0.34, 0.20),
            multiscale_frequencies=(0.18, 0.28, 0.40, 0.56, 0.76),
            multiscale_strength=1.35,
            noise_strength=3.0,
            noise_m_min=1,
            noise_m_max=16,
            noise_freq_base=0.50,
            noise_corr_time=0.095,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.45,
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.S_0": cfg.pde.S_0,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "forcing.u_bounds": list(cfg.forcing.u_bounds),
            "disturbance.mode1.amplitudes": list(cfg.disturbance.mode1.amplitudes),
            "disturbance.mode2.amplitudes": list(cfg.disturbance.mode2.amplitudes),
            "disturbance.warmup_scale": cfg.disturbance.warmup_scale,
            "disturbance.excite_scale": cfg.disturbance.excite_scale,
            "disturbance.hold_scale": cfg.disturbance.hold_scale,
            "disturbance.multiscale_modes": list(cfg.disturbance.multiscale_modes),
            "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.noise_band": [cfg.disturbance.noise_m_min, cfg.disturbance.noise_m_max],
            "disturbance.noise_corr_time": cfg.disturbance.noise_corr_time,
            "disturbance.noise_mode_decay": cfg.disturbance.noise_mode_decay,
        }
    )
    return cfg, changes


def _variant_swirl_drift_texture(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Increase visible rotation by combining stronger shear with drifting low modes.

    Compared with ``swirl_rich_texture``, this variant eases the radial drive a
    little and adds more explicit low-mode phase drift. The goal is to keep the
    stronger activity while making structures circulate around ``phi`` more
    clearly instead of primarily broadening in ``r``.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=1.40)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=2.0, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.76,
            kappa_0=4.65,
            S_0=0.10,
            nu_n=8.5e-6,
            nu_omega=8.5e-6,
            hyper_p=2,
            gamma_omega=0.018,
            source_floor_scale=0.05,
            source_balance_gain=0.90,
            source_mass_gain=2.4,
            source_scale_max=4.2,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=3.75, u_bounds=(4.6, 2.8, 2.8, 2.4, 2.4)),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(1.05, 0.40),
                frequencies=(0.15, 0.26),
                phases=(0.08, 0.92),
                phase_rate=0.30,
                phase_mod_amp=0.10,
                phase_mod_freq=0.10,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(0.62, 0.22),
                frequencies=(0.12, 0.24),
                phases=(0.58, -0.15),
                phase_rate=-0.22,
                phase_offset=0.45,
                phase_mod_amp=0.08,
                phase_mod_freq=0.09,
            ),
            warmup_fraction=0.03,
            excite_fraction=0.62,
            warmup_scale=1.10,
            excite_scale=3.7,
            hold_scale=4.3,
            multiscale_modes=(4, 6, 9, 13, 18),
            multiscale_amplitudes=(0.82, 0.64, 0.46, 0.28, 0.16),
            multiscale_frequencies=(0.19, 0.29, 0.42, 0.58, 0.80),
            multiscale_strength=1.15,
            noise_strength=2.8,
            noise_m_min=1,
            noise_m_max=14,
            noise_freq_base=0.50,
            noise_corr_time=0.080,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.35,
            phase_drift_strength=0.95,
            phase_drift_corr_time=0.11,
            eddy_strength=1.10,
            eddy_count=6,
            eddy_sigma_r=0.045,
            eddy_sigma_phi=0.22,
            eddy_speed_std=1.9,
            eddy_corr_time=0.12,
            eddy_radial_jitter=0.045,
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.S_0": cfg.pde.S_0,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "pde.source_mass_gain": cfg.pde.source_mass_gain,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "forcing.u_bounds": list(cfg.forcing.u_bounds),
            "disturbance.mode1.amplitudes": list(cfg.disturbance.mode1.amplitudes),
            "disturbance.mode1.phase_rate": cfg.disturbance.mode1.phase_rate,
            "disturbance.mode2.amplitudes": list(cfg.disturbance.mode2.amplitudes),
            "disturbance.mode2.phase_rate": cfg.disturbance.mode2.phase_rate,
            "disturbance.warmup_scale": cfg.disturbance.warmup_scale,
            "disturbance.excite_scale": cfg.disturbance.excite_scale,
            "disturbance.hold_scale": cfg.disturbance.hold_scale,
            "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.noise_band": [cfg.disturbance.noise_m_min, cfg.disturbance.noise_m_max],
            "disturbance.noise_corr_time": cfg.disturbance.noise_corr_time,
            "disturbance.noise_mode_decay": cfg.disturbance.noise_mode_decay,
            "disturbance.phase_drift_strength": cfg.disturbance.phase_drift_strength,
            "disturbance.eddy_count": cfg.disturbance.eddy_count,
            "disturbance.eddy_strength": cfg.disturbance.eddy_strength,
        }
    )
    return cfg, changes


def _variant_counter_swirl_texture(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Favor bidirectional swirl with weaker carriers and stronger eddy dipoles.

    This branch assumes the global one-way look is coming from two sources:
    a strong positive ``u0`` drive and low-mode coherent carriers that keep
    marching around the ring in one dominant direction. The retune keeps the
    plant active, but shifts more of the visible transport into stochastic
    dipolar eddies and lets the open-loop ``u0`` schedule cross zero.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=1.25)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=2.0, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.72,
            kappa_0=4.85,
            S_0=0.08,
            nu_n=7.8e-6,
            nu_omega=7.8e-6,
            hyper_p=2,
            gamma_omega=0.017,
            source_floor_scale=0.0,
            source_balance_gain=1.00,
            source_mass_gain=0.35,
            source_scale_max=3.2,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=1.30, u_bounds=(5.0, 2.8, 2.8, 2.4, 2.4)),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(0.44, 0.16),
                frequencies=(0.13, 0.23),
                phases=(0.05, 0.70),
                phase_rate=0.08,
                phase_mod_amp=0.06,
                phase_mod_freq=0.09,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(0.28, 0.10),
                frequencies=(0.11, 0.20),
                phases=(0.52, -0.10),
                phase_rate=-0.06,
                phase_offset=0.35,
                phase_mod_amp=0.05,
                phase_mod_freq=0.08,
            ),
            warmup_fraction=0.02,
            excite_fraction=0.62,
            warmup_scale=1.00,
            excite_scale=3.6,
            hold_scale=4.1,
            multiscale_modes=(3, 5, 8, 12, 17),
            multiscale_amplitudes=(0.42, 0.34, 0.24, 0.16, 0.10),
            multiscale_frequencies=(0.17, 0.27, 0.40, 0.57, 0.78),
            multiscale_strength=0.78,
            noise_strength=3.2,
            noise_m_min=1,
            noise_m_max=16,
            noise_freq_base=0.48,
            noise_corr_time=0.070,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.12,
            phase_drift_strength=1.45,
            phase_drift_corr_time=0.09,
            eddy_strength=1.85,
            eddy_count=10,
            eddy_sigma_r=0.040,
            eddy_sigma_phi=0.18,
            eddy_speed_std=2.3,
            eddy_corr_time=0.11,
            eddy_radial_jitter=0.035,
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.S_0": cfg.pde.S_0,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "forcing.u_bounds": list(cfg.forcing.u_bounds),
            "disturbance.mode1.amplitudes": list(cfg.disturbance.mode1.amplitudes),
            "disturbance.mode2.amplitudes": list(cfg.disturbance.mode2.amplitudes),
            "disturbance.multiscale_modes": list(cfg.disturbance.multiscale_modes),
            "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.phase_drift_strength": cfg.disturbance.phase_drift_strength,
            "disturbance.eddy_count": cfg.disturbance.eddy_count,
            "disturbance.eddy_strength": cfg.disturbance.eddy_strength,
        }
    )
    return cfg, changes


def _variant_alternating_swirl_texture(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Recover activity while keeping the lower-bias alternating swirl regime.

    This branch keeps the weaker coherent carriers and low-base open-loop
    drive from ``counter_swirl_texture``, but restores more motion through
    stronger drift-wave drive, stronger stochastic forcing, and more eddy
    packets. The intent is to keep the density visually alive without going
    back to a single dominant circulation direction.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=1.35)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=2.0, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.74,
            kappa_0=5.15,
            S_0=0.08,
            nu_n=7.2e-6,
            nu_omega=7.2e-6,
            hyper_p=2,
            gamma_omega=0.016,
            source_floor_scale=0.0,
            source_balance_gain=1.00,
            source_mass_gain=0.35,
            source_scale_max=3.4,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=1.35, u_bounds=(5.0, 2.8, 2.8, 2.4, 2.4)),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(0.52, 0.20),
                frequencies=(0.13, 0.24),
                phases=(0.08, 0.76),
                phase_rate=0.10,
                phase_mod_amp=0.08,
                phase_mod_freq=0.09,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(0.34, 0.12),
                frequencies=(0.11, 0.21),
                phases=(0.58, -0.08),
                phase_rate=-0.08,
                phase_offset=0.35,
                phase_mod_amp=0.06,
                phase_mod_freq=0.08,
            ),
            warmup_fraction=0.02,
            excite_fraction=0.62,
            warmup_scale=1.05,
            excite_scale=3.9,
            hold_scale=4.4,
            multiscale_modes=(3, 5, 8, 12, 17),
            multiscale_amplitudes=(0.50, 0.40, 0.28, 0.18, 0.11),
            multiscale_frequencies=(0.17, 0.27, 0.40, 0.57, 0.78),
            multiscale_strength=0.95,
            noise_strength=3.8,
            noise_m_min=1,
            noise_m_max=16,
            noise_freq_base=0.50,
            noise_corr_time=0.070,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.10,
            phase_drift_strength=1.60,
            phase_drift_corr_time=0.09,
            eddy_strength=2.35,
            eddy_count=12,
            eddy_sigma_r=0.040,
            eddy_sigma_phi=0.18,
            eddy_speed_std=2.45,
            eddy_corr_time=0.10,
            eddy_radial_jitter=0.035,
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.S_0": cfg.pde.S_0,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "forcing.u_bounds": list(cfg.forcing.u_bounds),
            "disturbance.mode1.amplitudes": list(cfg.disturbance.mode1.amplitudes),
            "disturbance.mode2.amplitudes": list(cfg.disturbance.mode2.amplitudes),
            "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.phase_drift_strength": cfg.disturbance.phase_drift_strength,
            "disturbance.eddy_count": cfg.disturbance.eddy_count,
            "disturbance.eddy_strength": cfg.disturbance.eddy_strength,
        }
    )
    return cfg, changes


def _variant_gradient_swirl_texture(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Use zonal density-gradient drive to reduce radial column coherence.

    The current alternating branch fixed mass drift and introduced sign changes
    in the azimuthal flow, but the density patterns remain too coherent in
    ``r``. This branch keeps the low-bias alternating swirl setup and adds a
    self-consistent drift-wave drive based on the evolving zonal density
    gradient. The inner and outer flanks of the ring can then propagate
    differently instead of sharing one nearly separable radial pattern.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=1.35)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=2.0, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.74,
            kappa_0=3.65,
            sigma_kappa=0.09,
            gradient_drive_gain=2.4,
            S_0=0.08,
            nu_n=7.0e-6,
            nu_omega=7.0e-6,
            hyper_p=2,
            gamma_omega=0.016,
            source_floor_scale=0.0,
            source_balance_gain=1.00,
            source_mass_gain=0.35,
            source_scale_max=3.4,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=1.30, u_bounds=(5.0, 2.8, 2.8, 2.4, 2.4)),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(0.48, 0.18),
                frequencies=(0.13, 0.24),
                phases=(0.08, 0.76),
                phase_rate=0.10,
                phase_mod_amp=0.08,
                phase_mod_freq=0.09,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(0.30, 0.11),
                frequencies=(0.11, 0.21),
                phases=(0.58, -0.08),
                phase_rate=-0.08,
                phase_offset=0.35,
                phase_mod_amp=0.06,
                phase_mod_freq=0.08,
            ),
            warmup_fraction=0.02,
            excite_fraction=0.62,
            warmup_scale=1.05,
            excite_scale=3.9,
            hold_scale=4.4,
            multiscale_modes=(3, 5, 8, 12, 17),
            multiscale_amplitudes=(0.54, 0.42, 0.30, 0.20, 0.12),
            multiscale_frequencies=(0.17, 0.27, 0.40, 0.57, 0.78),
            multiscale_strength=1.05,
            noise_strength=3.9,
            noise_m_min=1,
            noise_m_max=16,
            noise_freq_base=0.50,
            noise_corr_time=0.070,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.10,
            phase_drift_strength=1.65,
            phase_drift_corr_time=0.09,
            eddy_strength=2.45,
            eddy_count=12,
            eddy_sigma_r=0.034,
            eddy_sigma_phi=0.16,
            eddy_speed_std=2.55,
            eddy_corr_time=0.10,
            eddy_radial_jitter=0.050,
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.gradient_drive_gain": cfg.pde.gradient_drive_gain,
            "pde.S_0": cfg.pde.S_0,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "forcing.u_bounds": list(cfg.forcing.u_bounds),
            "disturbance.mode1.amplitudes": list(cfg.disturbance.mode1.amplitudes),
            "disturbance.mode2.amplitudes": list(cfg.disturbance.mode2.amplitudes),
            "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.phase_drift_strength": cfg.disturbance.phase_drift_strength,
            "disturbance.eddy_count": cfg.disturbance.eddy_count,
            "disturbance.eddy_strength": cfg.disturbance.eddy_strength,
        }
    )
    return cfg, changes


def _variant_curvature_swirl_texture(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Drive vorticity growth from density structure through curvature coupling.

    This is the main deeper-physics branch. Small density perturbations can now
    seed local vorticity through ``curvature_omega_gain * kappa(r) * dphi(n)``,
    so breakup can emerge from internal growth instead of relying as heavily on
    direct external vorticity forcing.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=1.20)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=2.0, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.76,
            kappa_0=4.95,
            sigma_kappa=0.10,
            gradient_drive_gain=0.0,
            curvature_omega_gain=0.72,
            S_0=0.08,
            nu_n=6.8e-6,
            nu_omega=6.8e-6,
            hyper_p=2,
            gamma_omega=0.015,
            source_floor_scale=0.0,
            source_balance_gain=1.00,
            source_mass_gain=0.35,
            source_scale_max=3.4,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=1.35, u_bounds=(5.0, 2.8, 2.8, 2.4, 2.4)),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(0.42, 0.16),
                frequencies=(0.13, 0.24),
                phases=(0.05, 0.70),
                phase_rate=0.08,
                phase_mod_amp=0.06,
                phase_mod_freq=0.09,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(0.28, 0.10),
                frequencies=(0.11, 0.21),
                phases=(0.52, -0.08),
                phase_rate=-0.06,
                phase_offset=0.35,
                phase_mod_amp=0.05,
                phase_mod_freq=0.08,
            ),
            warmup_fraction=0.03,
            excite_fraction=0.60,
            warmup_scale=0.90,
            excite_scale=3.8,
            hold_scale=4.3,
            multiscale_modes=(3, 5, 8, 12, 17),
            multiscale_amplitudes=(0.46, 0.36, 0.26, 0.17, 0.11),
            multiscale_frequencies=(0.17, 0.27, 0.40, 0.57, 0.78),
            multiscale_strength=0.90,
            noise_strength=3.4,
            noise_m_min=1,
            noise_m_max=16,
            noise_freq_base=0.50,
            noise_corr_time=0.070,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.10,
            phase_drift_strength=1.45,
            phase_drift_corr_time=0.09,
            eddy_strength=2.10,
            eddy_count=10,
            eddy_sigma_r=0.040,
            eddy_sigma_phi=0.18,
            eddy_speed_std=2.40,
            eddy_corr_time=0.10,
            eddy_radial_jitter=0.040,
            eddy_geometry="rotating_pair",
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.curvature_omega_gain": cfg.pde.curvature_omega_gain,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "forcing.u_bounds": list(cfg.forcing.u_bounds),
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.phase_drift_strength": cfg.disturbance.phase_drift_strength,
            "disturbance.eddy_count": cfg.disturbance.eddy_count,
            "disturbance.eddy_strength": cfg.disturbance.eddy_strength,
            "disturbance.eddy_geometry": cfg.disturbance.eddy_geometry,
        }
    )
    return cfg, changes


def _variant_interchange_swirl_texture(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Use nonlinear density-gradient coupling to create less columnar vortices.

    The curvature-only branch increased activity, but it also collapsed the
    swirl back toward a more one-way large-scale response. This branch keeps a
    weak curvature seed, then lets a stronger nonlinear baroclinic-style
    coupling generate local positive and negative vorticity from the evolving
    density gradients themselves.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=1.15)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=2.0, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.75,
            kappa_0=5.05,
            sigma_kappa=0.10,
            gradient_drive_gain=0.0,
            curvature_omega_gain=0.18,
            baroclinic_omega_gain=1.35,
            S_0=0.08,
            nu_n=6.9e-6,
            nu_omega=6.9e-6,
            hyper_p=2,
            gamma_omega=0.0155,
            source_floor_scale=0.0,
            source_balance_gain=1.00,
            source_mass_gain=0.35,
            source_scale_max=3.4,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=1.25, u_bounds=(5.0, 2.8, 2.8, 2.4, 2.4)),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(0.34, 0.13),
                frequencies=(0.12, 0.23),
                phases=(0.04, 0.64),
                phase_rate=0.07,
                phase_mod_amp=0.05,
                phase_mod_freq=0.08,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(0.22, 0.08),
                frequencies=(0.10, 0.19),
                phases=(0.46, -0.06),
                phase_rate=-0.05,
                phase_offset=0.28,
                phase_mod_amp=0.04,
                phase_mod_freq=0.07,
            ),
            warmup_fraction=0.02,
            excite_fraction=0.60,
            warmup_scale=0.95,
            excite_scale=4.0,
            hold_scale=4.5,
            multiscale_modes=(3, 5, 8, 12, 17),
            multiscale_amplitudes=(0.38, 0.30, 0.22, 0.14, 0.09),
            multiscale_frequencies=(0.17, 0.27, 0.40, 0.57, 0.78),
            multiscale_strength=0.72,
            noise_strength=3.6,
            noise_m_min=1,
            noise_m_max=18,
            noise_freq_base=0.50,
            noise_corr_time=0.072,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.08,
            phase_drift_strength=1.85,
            phase_drift_corr_time=0.09,
            eddy_strength=2.40,
            eddy_count=14,
            eddy_sigma_r=0.036,
            eddy_sigma_phi=0.17,
            eddy_speed_std=2.55,
            eddy_corr_time=0.10,
            eddy_radial_jitter=0.040,
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.curvature_omega_gain": cfg.pde.curvature_omega_gain,
            "pde.baroclinic_omega_gain": cfg.pde.baroclinic_omega_gain,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "forcing.u_bounds": list(cfg.forcing.u_bounds),
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.phase_drift_strength": cfg.disturbance.phase_drift_strength,
            "disturbance.eddy_count": cfg.disturbance.eddy_count,
            "disturbance.eddy_strength": cfg.disturbance.eddy_strength,
        }
    )
    return cfg, changes


def _variant_marginal_swirl_debug(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Dynamic operating-point branch that should propagate without breaking up.

    This is the requested debugging mode: the ring should remain coherent and
    clearly dynamic, but the stochastic perturbations should stay subcritical so
    the run does not fully break apart.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=0.45)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=2.0, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.82,
            kappa_0=3.90,
            sigma_kappa=0.10,
            gradient_drive_gain=0.0,
            curvature_omega_gain=0.12,
            baroclinic_omega_gain=0.18,
            S_0=0.08,
            nu_n=8.8e-6,
            nu_omega=8.8e-6,
            hyper_p=2,
            gamma_omega=0.020,
            source_floor_scale=0.0,
            source_balance_gain=1.00,
            source_mass_gain=0.30,
            source_scale_max=3.0,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=1.10, u_bounds=(4.5, 2.0, 2.0, 1.8, 1.8)),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(0.16, 0.06),
                frequencies=(0.12, 0.22),
                phases=(0.04, 0.50),
                phase_rate=0.05,
                phase_mod_amp=0.03,
                phase_mod_freq=0.07,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(0.10, 0.04),
                frequencies=(0.10, 0.18),
                phases=(0.40, -0.06),
                phase_rate=-0.04,
                phase_offset=0.25,
                phase_mod_amp=0.03,
                phase_mod_freq=0.07,
            ),
            warmup_fraction=0.10,
            excite_fraction=0.50,
            warmup_scale=0.40,
            excite_scale=1.15,
            hold_scale=1.30,
            multiscale_modes=(3, 5, 8),
            multiscale_amplitudes=(0.12, 0.09, 0.05),
            multiscale_frequencies=(0.16, 0.26, 0.40),
            multiscale_strength=0.30,
            noise_strength=0.80,
            noise_m_min=1,
            noise_m_max=10,
            noise_freq_base=0.45,
            noise_corr_time=0.080,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.25,
            phase_drift_strength=0.45,
            phase_drift_corr_time=0.10,
            eddy_strength=0.55,
            eddy_count=6,
            eddy_sigma_r=0.045,
            eddy_sigma_phi=0.20,
            eddy_speed_std=1.40,
            eddy_corr_time=0.11,
            eddy_radial_jitter=0.030,
            eddy_geometry="rotating_pair",
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.curvature_omega_gain": cfg.pde.curvature_omega_gain,
            "pde.baroclinic_omega_gain": cfg.pde.baroclinic_omega_gain,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "forcing.u_bounds": list(cfg.forcing.u_bounds),
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.phase_drift_strength": cfg.disturbance.phase_drift_strength,
            "disturbance.eddy_count": cfg.disturbance.eddy_count,
            "disturbance.eddy_strength": cfg.disturbance.eddy_strength,
            "disturbance.eddy_geometry": cfg.disturbance.eddy_geometry,
        }
    )
    return cfg, changes


def _variant_threshold_bhw_active(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Thresholded flux-driven branch with coherent motion and burst onset.

    The ring is driven by a narrower source that steepens the zonal density
    profile. Instability growth is then controlled by a critical-gradient
    activation opposed by zonal ExB shear. Disturbance forcing is present only
    to seed crossings of that threshold rather than to script the whole movie.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=0.70)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=2.1, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.88,
            kappa_0=4.85,
            sigma_kappa=0.10,
            gradient_drive_gain=3.2,
            critical_gradient_ratio=0.56,
            shear_suppression_gain=0.52,
            shear_ref=0.65,
            threshold_width=0.05,
            shear_damping_gain=0.22,
            curvature_omega_gain=0.14,
            baroclinic_omega_gain=0.26,
            flux_balance_omega_gain=1.45,
            phase_advection_gain=1.5,
            supercritical_coupling_gain=0.8,
            supercritical_feedback_gain=0.6,
            supercritical_transport_gain=1.0,
            landau_phi_gain=0.20,
            S_0=0.11,
            sigma_S=0.055,
            nu_n=6.5e-6,
            nu_omega=6.5e-6,
            hyper_p=2,
            gamma_omega=0.014,
            source_floor_scale=0.55,
            source_balance_gain=0.45,
            source_mass_gain=0.08,
            source_scale_max=3.8,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=1.45, u_bounds=(5.0, 2.8, 2.8, 2.4, 2.4)),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(0.10, 0.04),
                frequencies=(0.11, 0.20),
                phases=(0.02, 0.35),
                phase_rate=0.03,
                phase_mod_amp=0.04,
                phase_mod_freq=0.07,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(0.07, 0.03),
                frequencies=(0.09, 0.16),
                phases=(0.26, -0.05),
                phase_rate=-0.03,
                phase_offset=0.18,
                phase_mod_amp=0.03,
                phase_mod_freq=0.07,
            ),
            warmup_fraction=0.12,
            excite_fraction=0.54,
            warmup_scale=0.25,
            excite_scale=1.10,
            hold_scale=1.45,
            multiscale_modes=(3, 5, 8, 12),
            multiscale_amplitudes=(0.08, 0.06, 0.04, 0.02),
            multiscale_frequencies=(0.14, 0.24, 0.38, 0.56),
            multiscale_strength=0.16,
            noise_strength=1.45,
            noise_m_min=1,
            noise_m_max=14,
            noise_freq_base=0.45,
            noise_corr_time=0.070,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.04,
            phase_drift_strength=1.25,
            phase_drift_corr_time=0.09,
            eddy_strength=1.10,
            eddy_count=10,
            eddy_sigma_r=0.038,
            eddy_sigma_phi=0.18,
            eddy_speed_std=1.85,
            eddy_corr_time=0.10,
            eddy_radial_jitter=0.032,
            eddy_geometry="rotating_pair",
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.gradient_drive_gain": cfg.pde.gradient_drive_gain,
            "pde.critical_gradient_ratio": cfg.pde.critical_gradient_ratio,
            "pde.shear_suppression_gain": cfg.pde.shear_suppression_gain,
            "pde.shear_damping_gain": cfg.pde.shear_damping_gain,
            "pde.flux_balance_omega_gain": cfg.pde.flux_balance_omega_gain,
            "pde.phase_advection_gain": cfg.pde.phase_advection_gain,
            "pde.supercritical_coupling_gain": cfg.pde.supercritical_coupling_gain,
            "pde.supercritical_feedback_gain": cfg.pde.supercritical_feedback_gain,
            "pde.supercritical_transport_gain": cfg.pde.supercritical_transport_gain,
            "pde.landau_phi_gain": cfg.pde.landau_phi_gain,
            "pde.curvature_omega_gain": cfg.pde.curvature_omega_gain,
            "pde.baroclinic_omega_gain": cfg.pde.baroclinic_omega_gain,
            "pde.sigma_S": cfg.pde.sigma_S,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.phase_drift_strength": cfg.disturbance.phase_drift_strength,
            "disturbance.eddy_strength": cfg.disturbance.eddy_strength,
            "disturbance.eddy_count": cfg.disturbance.eddy_count,
            "disturbance.eddy_geometry": cfg.disturbance.eddy_geometry,
        }
    )
    return cfg, changes
def _variant_threshold_bhw_burst(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """More supercritical threshold branch intended to produce clear breakup.

    This keeps the same threshold/shear mechanism as
    :func:`_variant_threshold_bhw_active`, but makes the ring slightly more
    supercritical and reduces shear suppression so stochastic disturbances can
    accelerate into obvious bursts.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=0.58)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=2.1, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.86,
            kappa_0=5.05,
            sigma_kappa=0.10,
            gradient_drive_gain=3.8,
            critical_gradient_ratio=0.48,
            shear_suppression_gain=0.34,
            shear_ref=0.60,
            threshold_width=0.05,
            shear_damping_gain=0.14,
            curvature_omega_gain=0.18,
            baroclinic_omega_gain=0.32,
            flux_balance_omega_gain=2.10,
            phase_advection_gain=2.8,
            supercritical_coupling_gain=1.4,
            supercritical_feedback_gain=1.2,
            supercritical_transport_gain=2.0,
            landau_phi_gain=0.22,
            S_0=0.12,
            sigma_S=0.050,
            nu_n=5.8e-6,
            nu_omega=5.8e-6,
            hyper_p=2,
            gamma_omega=0.013,
            source_floor_scale=0.75,
            source_balance_gain=0.35,
            source_mass_gain=0.06,
            source_scale_max=3.9,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=1.50, u_bounds=(5.0, 2.8, 2.8, 2.4, 2.4)),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(0.12, 0.05),
                frequencies=(0.11, 0.20),
                phases=(0.02, 0.36),
                phase_rate=0.04,
                phase_mod_amp=0.04,
                phase_mod_freq=0.07,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(0.08, 0.03),
                frequencies=(0.09, 0.16),
                phases=(0.28, -0.05),
                phase_rate=-0.03,
                phase_offset=0.20,
                phase_mod_amp=0.03,
                phase_mod_freq=0.07,
            ),
            warmup_fraction=0.10,
            excite_fraction=0.56,
            warmup_scale=0.28,
            excite_scale=1.22,
            hold_scale=1.58,
            multiscale_modes=(3, 5, 8, 12, 17),
            multiscale_amplitudes=(0.09, 0.07, 0.05, 0.03, 0.02),
            multiscale_frequencies=(0.15, 0.25, 0.39, 0.56, 0.76),
            multiscale_strength=0.18,
            noise_strength=1.70,
            noise_m_min=1,
            noise_m_max=16,
            noise_freq_base=0.45,
            noise_corr_time=0.068,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.02,
            phase_drift_strength=1.45,
            phase_drift_corr_time=0.09,
            eddy_strength=1.35,
            eddy_count=12,
            eddy_sigma_r=0.038,
            eddy_sigma_phi=0.18,
            eddy_speed_std=1.95,
            eddy_corr_time=0.10,
            eddy_radial_jitter=0.034,
            eddy_geometry="rotating_pair",
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.gradient_drive_gain": cfg.pde.gradient_drive_gain,
            "pde.critical_gradient_ratio": cfg.pde.critical_gradient_ratio,
            "pde.shear_suppression_gain": cfg.pde.shear_suppression_gain,
            "pde.shear_damping_gain": cfg.pde.shear_damping_gain,
            "pde.flux_balance_omega_gain": cfg.pde.flux_balance_omega_gain,
            "pde.phase_advection_gain": cfg.pde.phase_advection_gain,
            "pde.supercritical_coupling_gain": cfg.pde.supercritical_coupling_gain,
            "pde.supercritical_feedback_gain": cfg.pde.supercritical_feedback_gain,
            "pde.supercritical_transport_gain": cfg.pde.supercritical_transport_gain,
            "pde.landau_phi_gain": cfg.pde.landau_phi_gain,
            "pde.curvature_omega_gain": cfg.pde.curvature_omega_gain,
            "pde.baroclinic_omega_gain": cfg.pde.baroclinic_omega_gain,
            "pde.sigma_S": cfg.pde.sigma_S,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.phase_drift_strength": cfg.disturbance.phase_drift_strength,
            "disturbance.eddy_strength": cfg.disturbance.eddy_strength,
            "disturbance.eddy_count": cfg.disturbance.eddy_count,
            "disturbance.eddy_geometry": cfg.disturbance.eddy_geometry,
        }
    )
    return cfg, changes
def _variant_threshold_bhw_marginal(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Thresholded marginal operating point with visible propagation only.

    This is the clean debug branch for the new paradigm: a narrow source keeps
    the profile close to critical, zonal shear suppresses the supercritical
    margin most of the time, and stochastic forcing remains small enough that
    the ring advects and fluctuates without fully breaking apart.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=0.35)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=2.1, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=0.90,
            kappa_0=4.40,
            sigma_kappa=0.10,
            gradient_drive_gain=2.7,
            critical_gradient_ratio=0.66,
            shear_suppression_gain=0.78,
            shear_ref=0.65,
            threshold_width=0.05,
            shear_damping_gain=0.34,
            curvature_omega_gain=0.10,
            baroclinic_omega_gain=0.18,
            flux_balance_omega_gain=0.95,
            phase_advection_gain=1.0,
            supercritical_coupling_gain=0.4,
            supercritical_feedback_gain=0.2,
            supercritical_transport_gain=0.5,
            landau_phi_gain=0.20,
            S_0=0.10,
            sigma_S=0.058,
            nu_n=7.8e-6,
            nu_omega=7.8e-6,
            hyper_p=2,
            gamma_omega=0.018,
            source_floor_scale=0.40,
            source_balance_gain=0.40,
            source_mass_gain=0.06,
            source_scale_max=3.4,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=1.18, u_bounds=(4.5, 2.0, 2.0, 1.8, 1.8)),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(0.06, 0.02),
                frequencies=(0.11, 0.19),
                phases=(0.02, 0.36),
                phase_rate=0.02,
                phase_mod_amp=0.03,
                phase_mod_freq=0.06,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(0.04, 0.02),
                frequencies=(0.09, 0.15),
                phases=(0.28, -0.04),
                phase_rate=-0.02,
                phase_offset=0.18,
                phase_mod_amp=0.02,
                phase_mod_freq=0.06,
            ),
            warmup_fraction=0.12,
            excite_fraction=0.50,
            warmup_scale=0.20,
            excite_scale=0.70,
            hold_scale=0.82,
            multiscale_modes=(3, 5, 8),
            multiscale_amplitudes=(0.06, 0.05, 0.03),
            multiscale_frequencies=(0.15, 0.24, 0.37),
            multiscale_strength=0.10,
            noise_strength=0.75,
            noise_m_min=1,
            noise_m_max=10,
            noise_freq_base=0.45,
            noise_corr_time=0.080,
            noise_refresh_dt=0.009,
            noise_mode_decay=0.25,
            phase_drift_strength=0.55,
            phase_drift_corr_time=0.10,
            eddy_strength=0.42,
            eddy_count=6,
            eddy_sigma_r=0.040,
            eddy_sigma_phi=0.20,
            eddy_speed_std=1.40,
            eddy_corr_time=0.11,
            eddy_radial_jitter=0.030,
            eddy_geometry="rotating_pair",
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.gradient_drive_gain": cfg.pde.gradient_drive_gain,
            "pde.critical_gradient_ratio": cfg.pde.critical_gradient_ratio,
            "pde.shear_suppression_gain": cfg.pde.shear_suppression_gain,
            "pde.shear_damping_gain": cfg.pde.shear_damping_gain,
            "pde.curvature_omega_gain": cfg.pde.curvature_omega_gain,
            "pde.baroclinic_omega_gain": cfg.pde.baroclinic_omega_gain,
            "pde.flux_balance_omega_gain": cfg.pde.flux_balance_omega_gain,
            "pde.phase_advection_gain": cfg.pde.phase_advection_gain,
            "pde.supercritical_coupling_gain": cfg.pde.supercritical_coupling_gain,
            "pde.supercritical_feedback_gain": cfg.pde.supercritical_feedback_gain,
            "pde.supercritical_transport_gain": cfg.pde.supercritical_transport_gain,
            "pde.landau_phi_gain": cfg.pde.landau_phi_gain,
            "pde.sigma_S": cfg.pde.sigma_S,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.phase_drift_strength": cfg.disturbance.phase_drift_strength,
            "disturbance.eddy_strength": cfg.disturbance.eddy_strength,
            "disturbance.eddy_count": cfg.disturbance.eddy_count,
            "disturbance.eddy_geometry": cfg.disturbance.eddy_geometry,
        }
    )
    return cfg, changes


def _variant_wave_landau_burst(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Near-threshold wave branch with strong stochastic breakup.

    This branch abandons the weak-response mHW threshold closure in favor of a
    direct propagating-wave density/vorticity coupling plus thresholded Landau
    growth and saturation in nonzonal vorticity. The intent is:

    - coherent fast ring motion in the warmup stage,
    - seed-dependent propagating disturbances around the ring,
    - rapid breakup once stochastic vorticity packets push the local state
      above threshold.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=0.78)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=1.0, snapshot_every=max(cfg.time.snapshot_every, 4)),
        pde=replace(
            cfg.pde,
            dynamics_model="wave_landau",
            C=0.34,
            kappa_0=4.2,
            sigma_kappa=0.10,
            critical_gradient_ratio=0.12,
            shear_suppression_gain=0.05,
            shear_ref=0.80,
            threshold_width=0.05,
            shear_damping_gain=0.0,
            phase_advection_gain=2.9,
            supercritical_coupling_gain=0.78,
            supercritical_feedback_gain=0.92,
            supercritical_transport_gain=1.95,
            wave_density_coupling_gain=3.6,
            wave_vorticity_coupling_gain=4.6,
            wave_packet_phi_gain=2.2,
            wave_packet_r_gain=0.14,
            wave_burst_speed_gain=1.05,
            packet_mode_cut=8,
            omega_high_damp=0.04,
            omega_landau_gain=1.72,
            omega_landau_sat=0.40,
            density_landau_gain=0.16,
            density_landau_sat=0.25,
            density_vortex_gain=0.18,
            density_packet_drive_gain=0.26,
            density_retention_gain=0.70,
            collapse_guard_gain=0.10,
            collapse_edge_threshold=0.30,
            collapse_core_threshold=0.0,
            collapse_omega_damp=0.08,
            radial_confinement_gain=0.04,
            radial_release_gain=0.80,
            radial_damage_release_gain=0.25,
            zonal_profile_restore_gain=0.0,
            zonal_profile_release_gain=0.0,
            zonal_profile_damage_release_gain=0.0,
            landau_phi_gain=0.08,
            S_0=0.12,
            sigma_S=0.050,
            nu_n=5.0e-6,
            nu_omega=5.0e-6,
            hyper_p=2,
            gamma_omega=0.018,
            source_floor_scale=0.42,
            source_balance_gain=0.09,
            source_mass_gain=0.04,
            source_scale_max=4.0,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=2.8, u_bounds=(5.0, 2.8, 2.8, 2.4, 2.4)),
        disturbance=replace(
            cfg.disturbance,
            warmup_fraction=0.22,
            excite_fraction=0.50,
            warmup_scale=0.24,
            excite_scale=1.28,
            hold_scale=1.78,
            multiscale_strength=0.04,
            noise_strength=2.55,
            noise_m_min=1,
            noise_m_max=18,
            noise_corr_time=0.07,
            noise_refresh_dt=0.006,
            phase_drift_strength=2.10,
            eddy_strength=1.90,
            eddy_count=14,
            eddy_speed_std=1.8,
            eddy_corr_time=0.07,
            eddy_geometry="rotating_pair",
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.dynamics_model": cfg.pde.dynamics_model,
            "pde.critical_gradient_ratio": cfg.pde.critical_gradient_ratio,
            "pde.shear_suppression_gain": cfg.pde.shear_suppression_gain,
            "pde.phase_advection_gain": cfg.pde.phase_advection_gain,
            "pde.supercritical_transport_gain": cfg.pde.supercritical_transport_gain,
            "pde.wave_density_coupling_gain": cfg.pde.wave_density_coupling_gain,
            "pde.wave_vorticity_coupling_gain": cfg.pde.wave_vorticity_coupling_gain,
            "pde.wave_packet_phi_gain": cfg.pde.wave_packet_phi_gain,
            "pde.wave_packet_r_gain": cfg.pde.wave_packet_r_gain,
            "pde.wave_burst_speed_gain": cfg.pde.wave_burst_speed_gain,
            "pde.packet_mode_cut": cfg.pde.packet_mode_cut,
            "pde.omega_high_damp": cfg.pde.omega_high_damp,
            "pde.omega_landau_gain": cfg.pde.omega_landau_gain,
            "pde.density_landau_gain": cfg.pde.density_landau_gain,
            "pde.density_vortex_gain": cfg.pde.density_vortex_gain,
            "pde.density_packet_drive_gain": cfg.pde.density_packet_drive_gain,
            "pde.density_retention_gain": cfg.pde.density_retention_gain,
            "pde.collapse_guard_gain": cfg.pde.collapse_guard_gain,
            "pde.collapse_edge_threshold": cfg.pde.collapse_edge_threshold,
            "pde.collapse_core_threshold": cfg.pde.collapse_core_threshold,
            "pde.collapse_omega_damp": cfg.pde.collapse_omega_damp,
            "pde.radial_confinement_gain": cfg.pde.radial_confinement_gain,
            "pde.radial_release_gain": cfg.pde.radial_release_gain,
            "pde.radial_damage_release_gain": cfg.pde.radial_damage_release_gain,
            "pde.zonal_profile_restore_gain": cfg.pde.zonal_profile_restore_gain,
            "pde.zonal_profile_release_gain": cfg.pde.zonal_profile_release_gain,
            "pde.zonal_profile_damage_release_gain": cfg.pde.zonal_profile_damage_release_gain,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.phase_drift_strength": cfg.disturbance.phase_drift_strength,
            "disturbance.eddy_strength": cfg.disturbance.eddy_strength,
            "disturbance.eddy_count": cfg.disturbance.eddy_count,
        }
    )
    return cfg, changes


def _variant_wave_landau_marginal(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Near-threshold wave branch below breakup.

    This is the operating-point debug branch for the new scaffold: clear
    azimuthal wave propagation and stochastic texture without full ring
    destruction.
    """

    cfg, changes = _with_minimal_seed(base, eps_omega=0.24)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=0.9, snapshot_every=max(cfg.time.snapshot_every, 4)),
        pde=replace(
            cfg.pde,
            dynamics_model="wave_landau",
            C=0.27,
            kappa_0=3.6,
            sigma_kappa=0.10,
            critical_gradient_ratio=0.24,
            shear_suppression_gain=0.28,
            shear_ref=0.80,
            threshold_width=0.05,
            shear_damping_gain=0.03,
            phase_advection_gain=2.15,
            supercritical_coupling_gain=0.55,
            supercritical_feedback_gain=0.50,
            supercritical_transport_gain=0.85,
            wave_density_coupling_gain=2.4,
            wave_vorticity_coupling_gain=3.1,
            wave_packet_phi_gain=1.55,
            wave_packet_r_gain=0.10,
            wave_burst_speed_gain=0.40,
            omega_landau_gain=0.90,
            omega_landau_sat=0.22,
            density_landau_gain=0.10,
            density_landau_sat=0.25,
            density_vortex_gain=0.10,
            radial_confinement_gain=0.10,
            radial_release_gain=0.70,
            radial_damage_release_gain=0.20,
            zonal_profile_restore_gain=1.65,
            zonal_profile_release_gain=1.35,
            zonal_profile_damage_release_gain=0.70,
            landau_phi_gain=0.05,
            S_0=0.11,
            sigma_S=0.052,
            nu_n=5.0e-6,
            nu_omega=5.0e-6,
            hyper_p=2,
            gamma_omega=0.012,
            source_floor_scale=0.20,
            source_balance_gain=0.08,
            source_mass_gain=0.04,
            source_scale_max=3.8,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=2.4, u_bounds=(4.5, 2.4, 2.4, 2.0, 2.0)),
        disturbance=replace(
            cfg.disturbance,
            warmup_fraction=0.24,
            excite_fraction=0.48,
            warmup_scale=0.18,
            excite_scale=0.88,
            hold_scale=0.98,
            multiscale_strength=0.03,
            noise_strength=1.15,
            noise_m_min=1,
            noise_m_max=14,
            noise_corr_time=0.06,
            noise_refresh_dt=0.006,
            phase_drift_strength=1.0,
            eddy_strength=0.80,
            eddy_count=8,
            eddy_speed_std=1.55,
            eddy_corr_time=0.08,
            eddy_geometry="rotating_pair",
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.dynamics_model": cfg.pde.dynamics_model,
            "pde.critical_gradient_ratio": cfg.pde.critical_gradient_ratio,
            "pde.shear_suppression_gain": cfg.pde.shear_suppression_gain,
            "pde.phase_advection_gain": cfg.pde.phase_advection_gain,
            "pde.supercritical_transport_gain": cfg.pde.supercritical_transport_gain,
            "pde.wave_density_coupling_gain": cfg.pde.wave_density_coupling_gain,
            "pde.wave_vorticity_coupling_gain": cfg.pde.wave_vorticity_coupling_gain,
            "pde.wave_packet_phi_gain": cfg.pde.wave_packet_phi_gain,
            "pde.wave_packet_r_gain": cfg.pde.wave_packet_r_gain,
            "pde.wave_burst_speed_gain": cfg.pde.wave_burst_speed_gain,
            "pde.omega_landau_gain": cfg.pde.omega_landau_gain,
            "pde.density_vortex_gain": cfg.pde.density_vortex_gain,
            "pde.radial_confinement_gain": cfg.pde.radial_confinement_gain,
            "pde.radial_release_gain": cfg.pde.radial_release_gain,
            "pde.radial_damage_release_gain": cfg.pde.radial_damage_release_gain,
            "pde.zonal_profile_restore_gain": cfg.pde.zonal_profile_restore_gain,
            "pde.zonal_profile_release_gain": cfg.pde.zonal_profile_release_gain,
            "pde.zonal_profile_damage_release_gain": cfg.pde.zonal_profile_damage_release_gain,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.phase_drift_strength": cfg.disturbance.phase_drift_strength,
            "disturbance.eddy_strength": cfg.disturbance.eddy_strength,
            "disturbance.eddy_count": cfg.disturbance.eddy_count,
        }
    )
    return cfg, changes


def _variant_shear_breakup(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Increase advective dominance so the ring shears and distorts sooner."""

    cfg, changes = _with_symmetric_density_start(base, eps_omega=1.55)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=1.9, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            C=1.05,
            kappa_0=5.4,
            S_0=0.10,
            nu_n=8.0e-6,
            nu_omega=8.0e-6,
            hyper_p=2,
            gamma_omega=0.018,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=3.2),
        disturbance=replace(
            cfg.disturbance,
            warmup_fraction=0.05,
            excite_fraction=0.65,
            warmup_scale=0.8,
            excite_scale=3.6,
            hold_scale=4.0,
            multiscale_modes=(4, 7, 10, 14, 20, 28),
            multiscale_amplitudes=(1.8, 1.6, 1.3, 1.0, 0.7, 0.5),
            multiscale_frequencies=(0.22, 0.31, 0.41, 0.53, 0.69, 0.87),
            multiscale_strength=2.6,
            noise_strength=0.95,
            noise_m_min=16,
            noise_m_max=34,
            noise_freq_base=1.05,
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.S_0": cfg.pde.S_0,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.hyper_p": cfg.pde.hyper_p,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "disturbance.multiscale_modes": list(cfg.disturbance.multiscale_modes),
            "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.noise_band": [cfg.disturbance.noise_m_min, cfg.disturbance.noise_m_max],
        }
    )
    return cfg, changes


def _variant_high_mode_filaments(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Push energy into higher azimuthal modes to test filament formation."""

    cfg, changes = _with_symmetric_density_start(base, eps_omega=1.7)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=2.2, snapshot_every=max(cfg.time.snapshot_every, 8)),
        pde=replace(
            cfg.pde,
            C=1.15,
            kappa_0=6.2,
            S_0=0.11,
            nu_n=6.0e-6,
            nu_omega=6.0e-6,
            hyper_p=2,
            gamma_omega=0.015,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=3.25),
        disturbance=replace(
            cfg.disturbance,
            warmup_fraction=0.0,
            excite_fraction=0.55,
            warmup_scale=1.2,
            excite_scale=4.0,
            hold_scale=5.0,
            multiscale_modes=(6, 9, 13, 17, 24, 32, 40),
            multiscale_amplitudes=(1.5, 1.45, 1.3, 1.0, 0.8, 0.6, 0.45),
            multiscale_frequencies=(0.25, 0.34, 0.47, 0.59, 0.73, 0.91, 1.07),
            multiscale_strength=3.3,
            noise_strength=1.45,
            noise_m_min=20,
            noise_m_max=44,
            noise_freq_base=1.2,
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "time.snapshot_every": cfg.time.snapshot_every,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.hyper_p": cfg.pde.hyper_p,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "disturbance.multiscale_modes": list(cfg.disturbance.multiscale_modes),
            "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.noise_band": [cfg.disturbance.noise_m_min, cfg.disturbance.noise_m_max],
        }
    )
    return cfg, changes


def _variant_long_horizon_sparse_frames(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Use a smaller internal time step and sparser saves to amplify visible motion."""

    cfg, changes = _with_symmetric_density_start(base, eps_omega=1.05)
    cfg = replace(
        cfg,
        time=replace(cfg.time, dt=1.0e-3, T_final=2.8, snapshot_every=10),
        pde=replace(
            cfg.pde,
            C=1.00,
            kappa_0=4.5,
            S_0=0.08,
            nu_n=1.1e-5,
            nu_omega=1.1e-5,
            hyper_p=2,
            gamma_omega=0.024,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=2.85),
        disturbance=replace(
            cfg.disturbance,
            warmup_fraction=0.12,
            excite_fraction=0.58,
            warmup_scale=0.35,
            excite_scale=3.0,
            hold_scale=3.8,
            multiscale_modes=(5, 8, 12, 16, 22, 30),
            multiscale_amplitudes=(1.4, 1.2, 1.0, 0.8, 0.6, 0.45),
            multiscale_frequencies=(0.18, 0.27, 0.39, 0.53, 0.71, 0.95),
            multiscale_strength=2.2,
            noise_strength=0.8,
            noise_m_min=14,
            noise_m_max=30,
            noise_freq_base=0.85,
        ),
    )
    changes.update(
        {
            "time.dt": cfg.time.dt,
            "time.T_final": cfg.time.T_final,
            "time.snapshot_every": cfg.time.snapshot_every,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.hyper_p": cfg.pde.hyper_p,
            "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
        }
    )
    return cfg, changes


def _variant_violent_burst(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Drive the strongest instability/noise combination in the scan."""

    cfg, changes = _with_symmetric_density_start(base, eps_omega=1.9)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=1.7, snapshot_every=max(cfg.time.snapshot_every, 5)),
        pde=replace(
            cfg.pde,
            C=1.20,
            kappa_0=6.6,
            S_0=0.13,
            nu_n=5.0e-6,
            nu_omega=5.0e-6,
            hyper_p=2,
            gamma_omega=0.012,
        ),
        forcing=replace(cfg.forcing, drive_u0_base=3.5),
        disturbance=replace(
            cfg.disturbance,
            mode1=_mode_update(
                cfg.disturbance.mode1,
                amplitudes=(4.8, 2.5),
                frequencies=(0.23, 0.58),
                phases=(0.1, 1.9),
                phase_rate=0.92,
                phase_mod_amp=0.45,
                phase_mod_freq=0.32,
            ),
            mode2=_mode_update(
                cfg.disturbance.mode2,
                amplitudes=(4.1, 2.1),
                frequencies=(0.18, 0.49),
                phases=(1.3, -0.2),
                phase_rate=-0.64,
                phase_mod_amp=0.38,
                phase_mod_freq=0.29,
            ),
            warmup_fraction=0.0,
            excite_fraction=0.50,
            warmup_scale=1.5,
            excite_scale=4.4,
            hold_scale=5.4,
            multiscale_modes=(5, 9, 14, 20, 27, 35, 44),
            multiscale_amplitudes=(2.0, 1.75, 1.5, 1.15, 0.85, 0.62, 0.45),
            multiscale_frequencies=(0.29, 0.41, 0.55, 0.73, 0.91, 1.08, 1.23),
            multiscale_strength=3.8,
            noise_strength=1.8,
            noise_m_min=18,
            noise_m_max=48,
            noise_freq_base=1.35,
        ),
    )
    changes.update(
        {
            "time.T_final": cfg.time.T_final,
            "pde.C": cfg.pde.C,
            "pde.kappa_0": cfg.pde.kappa_0,
            "pde.S_0": cfg.pde.S_0,
            "pde.nu_n": cfg.pde.nu_n,
            "pde.nu_omega": cfg.pde.nu_omega,
            "pde.hyper_p": cfg.pde.hyper_p,
            "pde.gamma_omega": cfg.pde.gamma_omega,
            "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
            "disturbance.mode1.amplitudes": list(cfg.disturbance.mode1.amplitudes),
            "disturbance.mode2.amplitudes": list(cfg.disturbance.mode2.amplitudes),
            "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
            "disturbance.noise_strength": cfg.disturbance.noise_strength,
            "disturbance.noise_band": [cfg.disturbance.noise_m_min, cfg.disturbance.noise_m_max],
        }
    )
    return cfg, changes


def _variant_tiny_seed_fast_break(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Reintroduce only a tiny density seed to compare against symmetric starts."""

    cfg = replace(
        base,
        init=replace(base.init, eps_n=0.008, eps_omega=0.9, mode1_amp=0.03, mode2_amp=0.015, mode5_amp=0.008),
        time=replace(base.time, T_final=1.8, snapshot_every=max(base.time.snapshot_every, 6)),
        pde=replace(
            base.pde,
            C=0.98,
            kappa_0=4.6,
            S_0=0.09,
            nu_n=1.3e-5,
            nu_omega=1.3e-5,
            hyper_p=2,
            gamma_omega=0.028,
        ),
        forcing=replace(base.forcing, drive_u0_base=2.85),
        disturbance=replace(
            base.disturbance,
            warmup_fraction=0.05,
            excite_fraction=0.60,
            warmup_scale=0.4,
            excite_scale=3.0,
            hold_scale=3.6,
            multiscale_modes=(5, 8, 11, 16, 22, 30),
            multiscale_amplitudes=(1.3, 1.12, 0.96, 0.74, 0.55, 0.4),
            multiscale_frequencies=(0.21, 0.33, 0.45, 0.61, 0.82, 1.01),
            multiscale_strength=1.8,
            noise_strength=0.55,
            noise_m_min=15,
            noise_m_max=34,
            noise_freq_base=0.9,
        ),
    )
    changes = {
        "init.eps_n": cfg.init.eps_n,
        "init.eps_omega": cfg.init.eps_omega,
        "init.mode1_amp": cfg.init.mode1_amp,
        "init.mode2_amp": cfg.init.mode2_amp,
        "init.mode5_amp": cfg.init.mode5_amp,
        "time.T_final": cfg.time.T_final,
        "pde.C": cfg.pde.C,
        "pde.kappa_0": cfg.pde.kappa_0,
        "pde.nu_n": cfg.pde.nu_n,
        "pde.nu_omega": cfg.pde.nu_omega,
        "pde.hyper_p": cfg.pde.hyper_p,
        "pde.gamma_omega": cfg.pde.gamma_omega,
        "forcing.drive_u0_base": cfg.forcing.drive_u0_base,
        "disturbance.multiscale_strength": cfg.disturbance.multiscale_strength,
        "disturbance.noise_strength": cfg.disturbance.noise_strength,
    }
    return cfg, changes


VARIANTS: tuple[VariantSpec, ...] = (
    VariantSpec(
        name="base_reference",
        description="Current loaded base config with no parameter changes.",
        build=_variant_base_reference,
    ),
    VariantSpec(
        name="symmetric_emergence",
        description="Exactly symmetric density at t=0, then delayed symmetry breaking and stronger late forcing.",
        build=_variant_symmetric_emergence,
    ),
    VariantSpec(
        name="stochastic_advection_breakup",
        description="Minimal seeding plus low-mode-biased stochastic forcing to favor advection-led breakup.",
        build=_variant_stochastic_advection_breakup,
    ),
    VariantSpec(
        name="hybrid_stochastic_shear",
        description="Weak low-mode scaffold plus stronger stochastic mid/high-band forcing for less scripted breakup.",
        build=_variant_hybrid_stochastic_shear,
    ),
    VariantSpec(
        name="textured_emergence",
        description="Keep the strong emergence scaffold but overlay stochastic texture to reduce scripted appearance.",
        build=_variant_textured_emergence,
    ),
    VariantSpec(
        name="azimuthal_shear_texture",
        description="Raise axisymmetric shear and use low/mid-band stochastic forcing to favor phi transport.",
        build=_variant_azimuthal_shear_texture,
    ),
    VariantSpec(
        name="swirl_rich_texture",
        description="Push the azimuthal-shear branch toward stronger circulating motion and longer-lived swirl.",
        build=_variant_swirl_rich_texture,
    ),
    VariantSpec(
        name="swirl_drift_texture",
        description="Favor stronger low-mode phase drift so structures rotate more clearly around phi.",
        build=_variant_swirl_drift_texture,
    ),
    VariantSpec(
        name="counter_swirl_texture",
        description="Reduce one-way carrier bias and rely more on stochastic dipole eddies for bidirectional swirl.",
        build=_variant_counter_swirl_texture,
    ),
    VariantSpec(
        name="alternating_swirl_texture",
        description="Keep the lower-bias swirl regime but restore stronger stochastic motion and eddy activity.",
        build=_variant_alternating_swirl_texture,
    ),
    VariantSpec(
        name="gradient_swirl_texture",
        description="Add self-consistent zonal density-gradient drive to break radial column coherence.",
        build=_variant_gradient_swirl_texture,
    ),
    VariantSpec(
        name="curvature_swirl_texture",
        description="Add density-to-vorticity curvature coupling so stochastic density structure can grow into swirl.",
        build=_variant_curvature_swirl_texture,
    ),
    VariantSpec(
        name="interchange_swirl_texture",
        description="Use nonlinear density-gradient coupling to create faster, less columnar stochastic swirl.",
        build=_variant_interchange_swirl_texture,
    ),
    VariantSpec(
        name="marginal_swirl_debug",
        description="Dynamic but subcritical operating point: propagating structures without full breakup.",
        build=_variant_marginal_swirl_debug,
    ),
    VariantSpec(
        name="threshold_bhw_active",
        description="Flux-driven threshold branch where radial-flux bursts can break a coherent turbulent ring.",
        build=_variant_threshold_bhw_active,
    ),
    VariantSpec(
        name="threshold_bhw_burst",
        description="More aggressive flux-driven threshold branch intended to trigger actual breakup bursts.",
        build=_variant_threshold_bhw_burst,
    ),
    VariantSpec(
        name="threshold_bhw_marginal",
        description="Flux-driven threshold branch below breakup: fast turbulent propagation without runaway bursts.",
        build=_variant_threshold_bhw_marginal,
    ),
    VariantSpec(
        name="wave_landau_burst",
        description="Near-threshold propagating-wave branch with thresholded vorticity growth and stochastic breakup.",
        build=_variant_wave_landau_burst,
    ),
    VariantSpec(
        name="wave_landau_marginal",
        description="Near-threshold propagating-wave branch below breakup: wave motion and texture without full destruction.",
        build=_variant_wave_landau_marginal,
    ),
    VariantSpec(
        name="shear_breakup",
        description="Lower damping and diffusion plus stronger multiscale forcing to trigger earlier breakup.",
        build=_variant_shear_breakup,
    ),
    VariantSpec(
        name="high_mode_filaments",
        description="Shift energy into higher azimuthal modes to test small-scale filament generation.",
        build=_variant_high_mode_filaments,
    ),
    VariantSpec(
        name="long_horizon_sparse_frames",
        description="Smaller internal dt, longer runtime, and sparser saves to amplify visible frame-to-frame change.",
        build=_variant_long_horizon_sparse_frames,
    ),
    VariantSpec(
        name="violent_burst",
        description="Most aggressive instability and forcing variant in the scan; tests the edge of acceptable dynamics.",
        build=_variant_violent_burst,
    ),
    VariantSpec(
        name="tiny_seed_fast_break",
        description="Near-symmetric density start with only tiny seed modes for comparison against fully symmetric starts.",
        build=_variant_tiny_seed_fast_break,
    ),
)

OPEN_LOOP_SAFE_VARIANTS: tuple[str, ...] = (
    "base_reference",
    "azimuthal_shear_texture",
    "counter_swirl_texture",
    "alternating_swirl_texture",
    "gradient_swirl_texture",
    "curvature_swirl_texture",
    "interchange_swirl_texture",
    "marginal_swirl_debug",
    "threshold_bhw_active",
    "threshold_bhw_burst",
    "threshold_bhw_marginal",
    "wave_landau_burst",
    "wave_landau_marginal",
    "hybrid_stochastic_shear",
    "symmetric_emergence",
    "textured_emergence",
    "long_horizon_sparse_frames",
    "tiny_seed_fast_break",
)


def _select_variants(args: argparse.Namespace) -> List[VariantSpec]:
    """Resolve the ordered list of variants requested on the CLI.

    Args:
        args: Parsed CLI arguments.

    Returns:
        Ordered list of scan variants to execute.
    """

    selected = list(VARIANTS)
    if args.variants:
        wanted = set(args.variants)
        selected = [variant for variant in VARIANTS if variant.name in wanted]
    elif args.control_mode == "open_loop":
        safe = set(OPEN_LOOP_SAFE_VARIANTS)
        selected = [variant for variant in VARIANTS if variant.name in safe]
    if args.limit is not None:
        selected = selected[: max(0, int(args.limit))]
    return selected


def _write_variant_note(run_dir: Path, variant: VariantSpec, changes: Dict[str, object]) -> None:
    """Write a per-run markdown note summarizing the scan variant.

    Args:
        run_dir: Completed run directory.
        variant: Variant metadata that produced the run.
        changes: Key parameter overrides applied to the base config.

    Returns:
        None.
    """

    lines = [
        "# Scan Variant",
        "",
        f"- name: `{variant.name}`",
        f"- description: {variant.description}",
        "",
        "## Key changes",
    ]
    if changes:
        for key, value in changes.items():
            lines.append(f"- `{key}`: `{value}`")
    else:
        lines.append("- none; this is the unchanged base reference")
    (run_dir / "notes" / "scan_variant.md").write_text("\n".join(lines), encoding="utf-8")


def _scan_payload_entry(
    variant: VariantSpec,
    changes: Dict[str, object],
    cfg: RunConfig,
    run_dir: str,
    summary: Dict[str, float],
    curves: Dict[str, List[float]],
) -> Dict[str, object]:
    """Build one manifest entry for a completed scan run."""

    all_finite_metrics = True
    for values in curves.values():
        arr = np.asarray(values, dtype=float)
        if arr.size and not np.isfinite(arr).all():
            all_finite_metrics = False
            break

    variation_path = Path(run_dir) / "metrics" / "variation_diagnostics.json"
    transport_path = Path(run_dir) / "metrics" / "transport_diagnostics.json"
    variation_summary: Dict[str, object] = {}
    transport_summary: Dict[str, object] = {}
    if variation_path.is_file():
        with variation_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        variation_summary = {
            "summary": payload.get("diagnostics", {}).get("summary", {}),
            "acceptance": payload.get("acceptance", {}),
        }
    if transport_path.is_file():
        with transport_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        transport_summary = payload.get("summary", {})

    return {
        "name": variant.name,
        "description": variant.description,
        "run_dir": run_dir,
        "config_path": str(Path(run_dir) / "config.yaml"),
        "movie_paths": {
            "n": str(Path(run_dir) / "movies" / "n.gif"),
            "omega": str(Path(run_dir) / "movies" / "omega.gif"),
            "u_mag": str(Path(run_dir) / "movies" / "u_mag.gif"),
            "u_r": str(Path(run_dir) / "movies" / "u_r.gif"),
            "u_phi": str(Path(run_dir) / "movies" / "u_phi.gif"),
        },
        "time": {
            "dt": cfg.time.dt,
            "T_final": cfg.time.T_final,
            "snapshot_every": cfg.time.snapshot_every,
            "n_steps": cfg.time.n_steps,
        },
        "changes": changes,
        "status": {
            "all_finite_metrics": bool(all_finite_metrics),
        },
        "summary": summary,
        "variation": variation_summary,
        "transport": transport_summary,
    }


def _write_scan_overview(manifest_path: Path, entries: Sequence[Dict[str, object]], scan_id: str) -> None:
    """Write a markdown overview next to the scan manifest.

    Args:
        manifest_path: Manifest file path for the full scan.
        entries: Completed run entries written into the manifest.
        scan_id: Human-readable scan identifier.

    Returns:
        None.
    """

    lines = [
        "# Data Scan Overview",
        "",
        f"- scan_id: `{scan_id}`",
        f"- manifest: `{manifest_path}`",
        "",
        "## Runs",
    ]
    for entry in entries:
        lines.append(f"- `{entry['name']}`: {entry['description']}")
        lines.append(f"  run_dir: `{entry['run_dir']}`")
        lines.append(f"  movie: `{entry['movie_paths']['n']}`")
    manifest_path.with_suffix(".md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Execute the uncontrolled data-variation scan and write a manifest.

    The function performs the following steps:
    1. Load one base plant configuration from a preset plus optional YAML file.
    2. Select the requested scan variants.
    3. Run one open-loop simulation per variant with full artifact writing.
    4. Save a scan manifest summarizing variant intent, run directories, and
       key parameter overrides.

    Args:
        None.

    Returns:
        None.
    """

    parser = build_parser([variant.name for variant in VARIANTS])
    args = parser.parse_args()

    if args.list:
        for variant in VARIANTS:
            print(f"{variant.name}: {variant.description}")
        return

    selected = _select_variants(args)
    if not selected:
        raise ValueError("No scan variants selected")

    base_cfg = load_run_config(args.config, preset=args.preset)
    if args.control_mode == "open_loop":
        base_cfg = _stabilize_open_loop_base(base_cfg)
    if args.outputs_root is not None:
        base_cfg = replace(base_cfg, output=replace(base_cfg.output, outputs_root=str(Path(args.outputs_root).resolve())))

    if args.baseline_run:
        os.environ["ERD_VARIATION_BASELINE_RUN"] = str(Path(args.baseline_run).resolve())
    else:
        os.environ.pop("ERD_VARIATION_BASELINE_RUN", None)

    if args.skip_media:
        os.environ["ERD_SKIP_MEDIA"] = "1"
    else:
        os.environ.pop("ERD_SKIP_MEDIA", None)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    scan_id = f"{stamp}_{args.tag}_{args.preset}"
    outputs_root = Path(base_cfg.output.outputs_root).resolve()
    outputs_root.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(args.manifest).resolve() if args.manifest else outputs_root / f"{scan_id}_manifest.yaml"

    entries: List[Dict[str, object]] = []

    for variant in selected:
        cfg_variant, changes = variant.build(base_cfg)
        cfg_variant = replace(
            cfg_variant,
            name=f"{base_cfg.name}_{variant.name}",
            output=replace(cfg_variant.output, tag=f"{scan_id}_{variant.name}", preset=args.preset),
        )
        print(f"[scan] {variant.name}: {variant.description}")
        plant = ERDPlant(cfg_variant)
        if args.control_mode == "random_piecewise":
            rng = np.random.default_rng(int(args.seed) + len(entries))
            schedule = _random_piecewise_schedule(
                rng=rng,
                bounds=np.asarray(cfg_variant.forcing.u_bounds, dtype=float),
                base_u0=cfg_variant.forcing.drive_u0_base,
                block_steps=max(1, int(args.block_steps)),
                t_final=cfg_variant.time.T_final,
                warmup_fraction=cfg_variant.disturbance.warmup_fraction,
                excite_fraction=cfg_variant.disturbance.excite_fraction,
            )
            result = plant.run(control_schedule=schedule, write_artifacts=True)
        else:
            schedule = _deterministic_open_loop_schedule(
                base_u0=cfg_variant.forcing.drive_u0_base,
                bound_u0=float(cfg_variant.forcing.u_bounds[0]),
                t_final=cfg_variant.time.T_final,
                warmup_fraction=cfg_variant.disturbance.warmup_fraction,
                excite_fraction=cfg_variant.disturbance.excite_fraction,
            )
            result = plant.run(control_schedule=schedule, write_artifacts=True)
        run_dir = Path(result.run_dir)
        _write_variant_note(run_dir, variant, changes)
        entry = _scan_payload_entry(variant, changes, cfg_variant, result.run_dir, result.summary, result.curves)
        entries.append(entry)
        print(f"  run_dir: {result.run_dir}")
        if not entry["status"]["all_finite_metrics"]:
            print("  warning: non-finite metrics detected; this variant may be numerically unstable")

    payload = {
        "scan_id": scan_id,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "preset": args.preset,
        "base_config_path": str(Path(args.config).resolve()) if args.config else None,
        "base_config": base_cfg.to_dict(),
        "variants": entries,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)

    _write_scan_overview(manifest_path, entries, scan_id)

    print(f"scan_id: {scan_id}")
    print(f"manifest: {manifest_path}")
    for entry in entries:
        print(f"{entry['name']}: {entry['run_dir']}")


if __name__ == "__main__":
    main()
