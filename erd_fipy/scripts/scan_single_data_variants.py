"""Run an open-loop parameter scan for visually tuning the ERD plant.

This script is intentionally focused on uncontrolled plant behavior. It loads a
single base plant configuration, applies a library of aggressive variant
overrides, and executes one open-loop run per variant with the standard ERD
artifacts, including GIF movies.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from datetime import datetime
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
            nu=1.8e-4,
            gamma=3.0e-3,
            D_r=3.0e-5,
            D_phi=1.2e-5,
            alpha=8.0e-4,
            beta_instab=4.2,
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
            "pde.nu": cfg.pde.nu,
            "pde.gamma": cfg.pde.gamma,
            "pde.D_r": cfg.pde.D_r,
            "pde.D_phi": cfg.pde.D_phi,
            "pde.alpha": cfg.pde.alpha,
            "pde.beta_instab": cfg.pde.beta_instab,
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


def _variant_shear_breakup(base: RunConfig) -> tuple[RunConfig, Dict[str, object]]:
    """Increase advective dominance so the ring shears and distorts sooner."""

    cfg, changes = _with_symmetric_density_start(base, eps_omega=1.55)
    cfg = replace(
        cfg,
        time=replace(cfg.time, T_final=1.9, snapshot_every=max(cfg.time.snapshot_every, 6)),
        pde=replace(
            cfg.pde,
            nu=1.1e-4,
            gamma=2.0e-3,
            D_r=1.8e-5,
            D_phi=6.0e-6,
            alpha=4.0e-4,
            beta_instab=5.4,
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
            "pde.nu": cfg.pde.nu,
            "pde.gamma": cfg.pde.gamma,
            "pde.D_r": cfg.pde.D_r,
            "pde.D_phi": cfg.pde.D_phi,
            "pde.alpha": cfg.pde.alpha,
            "pde.beta_instab": cfg.pde.beta_instab,
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
            nu=8.0e-5,
            gamma=2.0e-3,
            D_r=1.2e-5,
            D_phi=3.0e-6,
            alpha=2.0e-4,
            beta_instab=6.4,
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
            "pde.beta_instab": cfg.pde.beta_instab,
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
            nu=1.4e-4,
            gamma=2.5e-3,
            D_r=2.4e-5,
            D_phi=8.0e-6,
            alpha=6.0e-4,
            beta_instab=4.8,
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
            "pde.beta_instab": cfg.pde.beta_instab,
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
            nu=7.0e-5,
            gamma=1.5e-3,
            D_r=1.0e-5,
            D_phi=2.5e-6,
            alpha=1.0e-4,
            beta_instab=7.0,
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
            "pde.nu": cfg.pde.nu,
            "pde.gamma": cfg.pde.gamma,
            "pde.D_phi": cfg.pde.D_phi,
            "pde.beta_instab": cfg.pde.beta_instab,
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
        pde=replace(base.pde, nu=1.7e-4, gamma=3.2e-3, D_r=2.8e-5, D_phi=1.1e-5, alpha=7.0e-4, beta_instab=4.2),
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
        "pde.beta_instab": cfg.pde.beta_instab,
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

    return {
        "name": variant.name,
        "description": variant.description,
        "run_dir": run_dir,
        "config_path": str(Path(run_dir) / "config.yaml"),
        "movie_paths": {
            "n": str(Path(run_dir) / "movies" / "n.gif"),
            "omega": str(Path(run_dir) / "movies" / "omega.gif"),
            "u_mag": str(Path(run_dir) / "movies" / "u_mag.gif"),
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
        result = ERDPlant(cfg_variant).run(write_artifacts=True)
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
