"""End-to-end POD + controlled SINDy pipeline for ERD datasets.

This module owns model fitting orchestration and model-run artifact generation.
It consumes dataset manifests from ``erd_fipy``, trains a reduced model, runs
sanity rollouts, and writes lineage metadata so control runs can reference the
exact training source.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
import os
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List

import numpy as np

from .artifacts import save_json, save_npy, save_pickle, save_yaml
from .config import RunConfig
from .dataset_io import TrajectoryData, load_manifest, load_split
from .derivs import finite_difference_with_controls
from .metrics import compute_curves, summarize_aggregates
from .paths import make_run_dir
from .pod import PODResult, compute_pod
from .snapshot import SnapshotLayout, build_snapshot_matrix, fields_to_state_vec

if TYPE_CHECKING:
    from .rollout import RolloutResult
    from .sindy_model import SINDyControlFitResult


@dataclass(frozen=True)
class ModelRunResult:
    """Container returned by :func:`run`.

    Attributes:
        run_dir: Model run-folder path.
        pod: POD decomposition result.
        sindy: Fitted controlled SINDy model bundle.
        rollout: Representative held-out rollout result.
        aggregates: Aggregate validation metrics.
    """

    run_dir: str
    pod: PODResult
    sindy: "SINDyControlFitResult"
    rollout: "RolloutResult"
    aggregates: Dict[str, float]



def _setup_logger(run_dir: Path) -> logging.Logger:
    """Create file + stream logger for a model run.

    Args:
        run_dir: Run-folder path containing ``logs/``.

    Returns:
        Configured logger instance.
    """

    logger = logging.getLogger(f"erd_model_{run_dir.name}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    fh = logging.FileHandler(run_dir / "logs" / "run.log")
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(fh)

    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(sh)
    return logger



def _project_traj_to_coeffs(
    traj: TrajectoryData,
    layout: SnapshotLayout,
    basis: np.ndarray,
    mean_state: np.ndarray | None,
) -> np.ndarray:
    """Project one trajectory from full fields into reduced coefficients.

    Args:
        traj: Loaded trajectory data.
        layout: Snapshot stacking layout metadata.
        basis: POD basis matrix.
        mean_state: Optional centering vector used during POD.

    Returns:
        Coefficient trajectory array of shape ``(T, r)``.
    """

    T = traj.fields.shape[0]
    A = np.empty((T, basis.shape[1]), dtype=float)
    for k in range(T):
        x = fields_to_state_vec(traj.fields[k], layout)
        if mean_state is not None:
            x = x - mean_state
        A[k] = basis.T @ x
    return A



def _clip_rank(r: int) -> int:
    """Clamp configured POD rank to the accepted range ``[4, 8]``.

    Args:
        r: Requested POD rank.

    Returns:
        Clamped integer rank.
    """

    return int(min(8, max(4, r)))


def _coefficients_collapsed(Xi: np.ndarray, atol: float = 1e-12) -> bool:
    """Check whether fitted coefficients are effectively all zero/non-finite.

    Args:
        Xi: Coefficient matrix from sparse regression.
        atol: Absolute threshold for nonzero detection.

    Returns:
        ``True`` when coefficients are non-finite or all below ``atol``.
    """

    Xi_arr = np.asarray(Xi, dtype=float)
    if Xi_arr.size == 0:
        return True
    if not np.all(np.isfinite(Xi_arr)):
        return True
    return bool(np.count_nonzero(np.abs(Xi_arr) > atol) == 0)


def _quick_validation_mean_rel_l2(
    model: object,
    heldout_traj: TrajectoryData,
    layout: SnapshotLayout,
    pod: PODResult,
    mean_state: np.ndarray | None,
    dt: float,
    horizon_steps: int,
    state_clip: float,
) -> float:
    """Compute a quick held-out mean relative field error for fit triage.

    Args:
        model: Candidate controlled model with ``predict`` API.
        heldout_traj: One held-out trajectory.
        layout: Snapshot layout metadata.
        pod: POD decomposition output.
        mean_state: Optional POD centering vector.
        dt: Integration time step.
        horizon_steps: Maximum validation horizon.
        state_clip: Reduced-state clipping bound.

    Returns:
        Mean relative field L2 error for the candidate model.
    """

    from .rollout import rollout_one

    A_traj = _project_traj_to_coeffs(heldout_traj, layout, pod.U, mean_state)
    horizon = min(horizon_steps, A_traj.shape[0] - 1)
    ro = rollout_one(
        model=model,
        U_basis=pod.U,
        layout=layout,
        A_traj=A_traj,
        U_traj=heldout_traj.controls,
        dt=dt,
        horizon_steps=horizon,
        state_clip=state_clip,
        mean_state=mean_state,
        start_idx=0,
    )
    err, _ = compute_curves(ro.A_true, ro.A_pred, ro.fields_true, ro.fields_pred)
    return float(np.mean(err.field_rel_l2))



def run(cfg: RunConfig) -> ModelRunResult:
    """Execute the full ERD model-fitting workflow.

    The workflow is:
    1) load train/val/test trajectories from a manifest,
    2) build centered snapshot matrix and controls,
    3) fit POD and controlled SINDy,
    4) validate on held-out rollouts,
    5) write run-folder artifacts (arrays, metadata, plots, GIF).

    Args:
        cfg: Full model pipeline configuration.

    Returns:
        :class:`ModelRunResult` with key in-memory outputs and run path.
    """

    manifest = load_manifest(cfg.data.manifest_path)
    dataset_id = str(manifest.get("dataset_id", "dataset")).strip() or "dataset"

    run_tag = cfg.output.tag
    if run_tag.strip().lower() in {"", "auto"}:
        run_tag = f"model_{dataset_id}"

    run_dir = make_run_dir(cfg.output.outputs_root, run_tag)
    logger = _setup_logger(run_dir)

    logger.info("Loading manifest and trajectories")
    resolved_manifest = manifest.get("_resolved_manifest_path", cfg.data.manifest_path)
    logger.info(f"Using manifest: {resolved_manifest}")
    logger.info(f"Dataset ID: {dataset_id}")

    if cfg.deriv.order != 2:
        raise ValueError(f"Only central finite differences (order=2) are supported, got {cfg.deriv.order}")
    if cfg.sindy.optimizer.lower() not in {"stlsq", "trappingsr3", "trapping_sr3", "trapping-sr3"}:
        raise ValueError(f"Unsupported SINDy optimizer: {cfg.sindy.optimizer!r}")

    train = load_split(manifest, cfg.data.train_split, max_trajectories=cfg.data.max_train)
    val = load_split(manifest, cfg.data.val_split, max_trajectories=cfg.data.max_val)
    test = load_split(manifest, cfg.data.test_split, max_trajectories=cfg.data.max_test)

    if not train:
        raise ValueError("No training trajectories in manifest")
    heldout = test if test else val
    if not heldout:
        raise ValueError("Need at least one held-out trajectory (val or test)")

    dt = float(np.mean([tr.dt for tr in train]))
    logger.info(f"Train trajectories: {len(train)}, heldout: {len(heldout)}, dt={dt:.4e}")

    X_train, U_train, layout, mean_state, seg_lengths, _ = build_snapshot_matrix(
        trajectories=train,
        center=cfg.pod.center,
        time_stride=cfg.data.time_stride,
        time_limit=cfg.data.time_limit,
    )
    logger.info(
        "Snapshot matrix built: X=%s, U=%s, segments=%d",
        X_train.shape,
        U_train.shape,
        len(seg_lengths),
    )

    rank = _clip_rank(cfg.pod.rank)
    pod = compute_pod(X_train, rank=rank)
    logger.info(
        "POD complete: rank=%d, retained_energy=%.4f",
        pod.r,
        pod.energy_fraction,
    )
    if pod.r > int(cfg.acceptance.max_rank):
        raise ValueError(
            f"POD rank gate failed: rank={pod.r} > max_rank={cfg.acceptance.max_rank}"
        )
    if pod.energy_fraction < float(cfg.acceptance.min_pod_energy):
        raise ValueError(
            "POD retained-energy gate failed: "
            f"{pod.energy_fraction:.4f} < {cfg.acceptance.min_pod_energy:.4f}"
        )

    state_clip = float(max(25.0, 2.5 * np.percentile(np.abs(pod.A), 95)))
    logger.info("Using rollout state clip: %.3f", state_clip)

    deriv = finite_difference_with_controls(
        A=pod.A,
        U_controls=U_train,
        dt=dt,
        segment_lengths=seg_lengths,
    )
    logger.info(
        "Derivative samples prepared: A_used=%s, U_used=%s",
        deriv.A_used.shape,
        deriv.U_used.shape,
    )

    logger.info("Fitting controlled SINDy model")
    from .sindy_model import fit_sindy_control
    backoff_factor = float(cfg.sindy.backoff_factor)
    if not (0.0 < backoff_factor < 1.0):
        backoff_factor = 0.1
    retries = int(max(0, cfg.sindy.backoff_retries))
    spike_thresh = float(cfg.sindy.validation_spike_threshold)

    fit_attempts: List[Dict[str, object]] = []
    collapse_events = 0
    fallback_used = False
    fit_warning: str | None = None

    def _attempt_fit(optimizer_name: str) -> tuple[object | None, float]:
        """Try fitting with threshold backoff and quick held-out triage.

        Args:
            optimizer_name: Sparse optimizer mode for this attempt block.

        Returns:
            Tuple ``(result_or_none, used_threshold)``.
        """

        nonlocal collapse_events

        threshold_i = float(cfg.sindy.threshold)
        for retry_idx in range(retries + 1):
            try:
                candidate = fit_sindy_control(
                    A_used=deriv.A_used,
                    dA_dt=deriv.dA_dt,
                    U_used=deriv.U_used,
                    dt=dt,
                    threshold=threshold_i,
                    alpha=cfg.sindy.alpha,
                    max_iter=cfg.sindy.max_iter,
                    optimizer=optimizer_name,
                    trapping=cfg.sindy.trapping,
                )
            except Exception as exc:
                fit_attempts.append(
                    {
                        "optimizer": optimizer_name,
                        "retry": retry_idx,
                        "threshold": threshold_i,
                        "status": "fit_exception",
                        "message": str(exc),
                    }
                )
                threshold_i *= backoff_factor
                continue

            collapsed = _coefficients_collapsed(candidate.coefficient_matrix)
            if collapsed:
                collapse_events += 1

            try:
                quick_rel = _quick_validation_mean_rel_l2(
                    model=candidate.model,
                    heldout_traj=heldout[0],
                    layout=layout,
                    pod=pod,
                    mean_state=mean_state,
                    dt=dt,
                    horizon_steps=cfg.rollout.horizon_steps,
                    state_clip=state_clip,
                )
            except Exception:
                quick_rel = float("inf")

            spiking = (not np.isfinite(quick_rel)) or (quick_rel > spike_thresh)
            status = "accepted"
            if collapsed:
                status = "collapsed"
            elif spiking:
                status = "validation_spike"

            fit_attempts.append(
                {
                    "optimizer": optimizer_name,
                    "retry": retry_idx,
                    "threshold": threshold_i,
                    "status": status,
                    "quick_mean_field_rel_l2": quick_rel,
                    "collapsed": collapsed,
                }
            )

            if (not collapsed) and (not spiking):
                return candidate, threshold_i

            threshold_i *= backoff_factor

        return None, threshold_i

    used_optimizer = str(cfg.sindy.optimizer)
    sindy, used_threshold = _attempt_fit(used_optimizer)

    if sindy is None and used_optimizer.lower() != str(cfg.sindy.fallback_optimizer).lower():
        fallback_used = True
        used_optimizer = str(cfg.sindy.fallback_optimizer)
        fit_warning = (
            f"Primary optimizer {cfg.sindy.optimizer} failed after backoff; "
            f"falling back to {used_optimizer}."
        )
        logger.warning(fit_warning)
        sindy, used_threshold = _attempt_fit(used_optimizer)

    if sindy is None:
        save_yaml(run_dir / "config.yaml", cfg.to_dict())
        save_json(
            run_dir / "model" / "fit_diagnostics.json",
            {
                "optimizer_requested": cfg.sindy.optimizer,
                "optimizer_used": used_optimizer,
                "fallback_optimizer": cfg.sindy.fallback_optimizer,
                "fallback_used": fallback_used,
                "threshold_requested": cfg.sindy.threshold,
                "threshold_used": used_threshold,
                "backoff_factor": backoff_factor,
                "backoff_retries": retries,
                "validation_spike_threshold": spike_thresh,
                "coefficient_collapse_events": int(collapse_events),
                "fit_warning": fit_warning,
                "attempts": fit_attempts,
            },
        )
        raise RuntimeError(
            "SINDy fitting failed after backoff/fallback attempts; "
            f"see {run_dir / 'model' / 'fit_diagnostics.json'} for details."
        )

    logger.info("Running held-out validation rollouts")
    from .rollout import RolloutResult, rollout_one

    rollout_results: List[RolloutResult] = []
    agg_list: List[Dict[str, float]] = []

    for traj in heldout:
        A_traj = _project_traj_to_coeffs(traj, layout, pod.U, mean_state)
        horizon = min(cfg.rollout.horizon_steps, A_traj.shape[0] - 1)

        ro = rollout_one(
            model=sindy.model,
            U_basis=pod.U,
            layout=layout,
            A_traj=A_traj,
            U_traj=traj.controls,
            dt=dt,
            horizon_steps=horizon,
            state_clip=state_clip,
            mean_state=mean_state,
            start_idx=0,
        )
        err, energy = compute_curves(ro.A_true, ro.A_pred, ro.fields_true, ro.fields_pred)
        agg = summarize_aggregates(err, energy)

        rollout_results.append(ro)
        agg_list.append(agg)

    aggregates = {k: float(np.mean([a[k] for a in agg_list])) for k in agg_list[0].keys()}

    first_rollout = rollout_results[0]
    err0, energy0 = compute_curves(
        first_rollout.A_true,
        first_rollout.A_pred,
        first_rollout.fields_true,
        first_rollout.fields_pred,
    )

    split_eval_curves: Dict[str, Dict[str, List[float]]] = {}
    split_eval_summary: Dict[str, Dict[str, float]] = {}
    split_targets = {
        cfg.data.train_split: train,
        cfg.data.test_split: test if test else val,
    }
    for split_name, split_traj in split_targets.items():
        if not split_traj:
            continue
        traj = split_traj[0]
        A_traj = _project_traj_to_coeffs(traj, layout, pod.U, mean_state)
        horizon = min(cfg.rollout.horizon_steps, A_traj.shape[0] - 1)
        ro = rollout_one(
            model=sindy.model,
            U_basis=pod.U,
            layout=layout,
            A_traj=A_traj,
            U_traj=traj.controls,
            dt=dt,
            horizon_steps=horizon,
            state_clip=state_clip,
            mean_state=mean_state,
            start_idx=0,
        )
        err_split, energy_split = compute_curves(ro.A_true, ro.A_pred, ro.fields_true, ro.fields_pred)
        split_eval_curves[split_name] = {
            "coeff_mse": err_split.coeff_mse.tolist(),
            "field_l2": err_split.field_l2.tolist(),
            "field_rel_l2": err_split.field_rel_l2.tolist(),
            "energy_true": energy_split.energy_true.tolist(),
            "energy_pred": energy_split.energy_pred.tolist(),
        }
        split_eval_summary[split_name] = summarize_aggregates(err_split, energy_split)

    mean_heldout_rel = float(aggregates["mean_field_rel_l2"])
    max_rel_gate = float(cfg.acceptance.max_mean_field_rel_l2)
    acceptance = {
        "pod_rank": pod.r,
        "pod_rank_max": int(cfg.acceptance.max_rank),
        "pod_rank_ok": bool(pod.r <= int(cfg.acceptance.max_rank)),
        "pod_energy": float(pod.energy_fraction),
        "pod_energy_min": float(cfg.acceptance.min_pod_energy),
        "pod_energy_ok": bool(pod.energy_fraction >= float(cfg.acceptance.min_pod_energy)),
        "mean_field_rel_l2": mean_heldout_rel,
        "max_mean_field_rel_l2": max_rel_gate,
        "heldout_error_ok": bool(mean_heldout_rel <= max_rel_gate),
        "coefficient_collapse_events": int(collapse_events),
        "fallback_used": bool(fallback_used),
        "optimizer_used": used_optimizer,
    }
    acceptance["all_gates_passed"] = bool(
        acceptance["pod_rank_ok"] and acceptance["pod_energy_ok"] and acceptance["heldout_error_ok"]
    )

    logger.info("Writing artifacts")
    save_yaml(run_dir / "config.yaml", cfg.to_dict())
    manifest_to_save = dict(manifest)
    manifest_to_save.pop("_resolved_manifest_path", None)
    save_yaml(run_dir / "model" / "training_manifest.yaml", manifest_to_save)

    save_npy(run_dir / "model" / "basis_U.npy", pod.U)
    save_npy(run_dir / "model" / "singular_values.npy", pod.s)
    if mean_state is not None:
        save_npy(run_dir / "model" / "mean_state.npy", mean_state)
    save_npy(run_dir / "model" / "coeffs_A_train.npy", pod.A)

    save_pickle(run_dir / "model" / "sindy_model.pkl", sindy.model)
    save_npy(run_dir / "model" / "Xi.npy", sindy.coefficient_matrix)
    save_json(run_dir / "model" / "feature_names.json", {"features": sindy.feature_names})
    save_json(
        run_dir / "model" / "model_meta.json",
        {
            "dataset_id": dataset_id,
            "rank": pod.r,
            "dt": dt,
            "state_clip": state_clip,
            "optimizer_requested": cfg.sindy.optimizer,
            "optimizer_used": used_optimizer,
            "threshold_requested": cfg.sindy.threshold,
            "threshold_used": used_threshold,
            "alpha": cfg.sindy.alpha,
            "max_iter": cfg.sindy.max_iter,
            "fallback_optimizer": cfg.sindy.fallback_optimizer,
            "fallback_used": fallback_used,
            "validation_spike_threshold": cfg.sindy.validation_spike_threshold,
            "trapping": {
                "method": cfg.sindy.trapping.method,
                "reg_weight_lam": cfg.sindy.trapping.reg_weight_lam,
                "relax_coeff_nu": cfg.sindy.trapping.relax_coeff_nu,
                "tol": cfg.sindy.trapping.tol,
                "max_iter": cfg.sindy.trapping.max_iter,
                "eta": cfg.sindy.trapping.eta,
                "stability_alpha": cfg.sindy.trapping.stability_alpha,
                "stability_beta": cfg.sindy.trapping.stability_beta,
                "control_threshold": cfg.sindy.trapping.control_threshold,
                "control_alpha": cfg.sindy.trapping.control_alpha,
                "control_max_iter": cfg.sindy.trapping.control_max_iter,
            },
            "fit_info": sindy.fit_info,
            "n_train_traj": len(train),
            "n_heldout_traj": len(heldout),
            "layout": {"nr": layout.nr, "nphi": layout.nphi, "n_components": layout.n_components},
        },
    )
    save_json(
        run_dir / "model" / "fit_diagnostics.json",
        {
            "optimizer_requested": cfg.sindy.optimizer,
            "optimizer_used": used_optimizer,
            "fallback_optimizer": cfg.sindy.fallback_optimizer,
            "fallback_used": fallback_used,
            "threshold_requested": cfg.sindy.threshold,
            "threshold_used": used_threshold,
            "backoff_factor": backoff_factor,
            "backoff_retries": retries,
            "validation_spike_threshold": spike_thresh,
            "coefficient_collapse_events": int(collapse_events),
            "fit_warning": fit_warning,
            "attempts": fit_attempts,
            "fit_info": sindy.fit_info,
        },
    )
    save_json(
        run_dir / "model" / "lineage.json",
        {
            "dataset_id": dataset_id,
            "manifest_path": str(resolved_manifest),
            "train_split": cfg.data.train_split,
            "test_split": cfg.data.test_split,
            "train_run_dirs": [tr.run_dir for tr in train],
            "test_run_dirs": [tr.run_dir for tr in (test if test else val)],
        },
    )

    save_npy(run_dir / "fields" / "rollout_q_true.npy", first_rollout.q_true)
    save_npy(run_dir / "fields" / "rollout_q_pred.npy", first_rollout.q_pred)
    save_npy(run_dir / "controls" / "rollout_u.npy", first_rollout.U_used)

    save_json(
        run_dir / "metrics" / "curves.json",
        {
            "coeff_mse": err0.coeff_mse.tolist(),
            "field_l2": err0.field_l2.tolist(),
            "field_rel_l2": err0.field_rel_l2.tolist(),
            "energy_true": energy0.energy_true.tolist(),
            "energy_pred": energy0.energy_pred.tolist(),
        },
    )
    save_json(run_dir / "metrics" / "aggregates.json", aggregates)
    save_json(run_dir / "metrics" / "split_aggregates.json", split_eval_summary)
    save_json(run_dir / "metrics" / "acceptance.json", acceptance)
    for split_name, split_curves in split_eval_curves.items():
        save_json(run_dir / "metrics" / f"{split_name}_curves.json", split_curves)
    save_json(
        run_dir / "summary.json",
        {
            "run_dir": str(run_dir),
            "optimizer_used": used_optimizer,
            "fallback_used": fallback_used,
            "threshold_used": used_threshold,
            "coefficient_collapse_events": int(collapse_events),
            **aggregates,
            "acceptance_all_gates_passed": acceptance["all_gates_passed"],
        },
    )

    skip_media = os.environ.get("ERD_SKIP_MEDIA", "").strip() == "1"
    if not skip_media:
        from .plots_and_movies import generate_all_plots_and_movies

        generate_all_plots_and_movies(
            cfg=cfg.plots,
            run_dir=run_dir,
            pod=pod,
            sindy=sindy,
            rollout=first_rollout,
            err=err0,
            energy=energy0,
        )
    else:
        logger.info("Skipping model plot/movie rendering because ERD_SKIP_MEDIA=1")

    if not acceptance["all_gates_passed"]:
        raise ValueError(
            "Model acceptance gates failed: "
            f"pod_rank_ok={acceptance['pod_rank_ok']}, "
            f"pod_energy_ok={acceptance['pod_energy_ok']}, "
            f"heldout_error_ok={acceptance['heldout_error_ok']} "
            f"(mean_field_rel_l2={mean_heldout_rel:.4f}, max={max_rel_gate:.4f})."
        )

    logger.info(f"Model run complete: {run_dir}")
    return ModelRunResult(
        run_dir=str(run_dir),
        pod=pod,
        sindy=sindy,
        rollout=first_rollout,
        aggregates=aggregates,
    )
