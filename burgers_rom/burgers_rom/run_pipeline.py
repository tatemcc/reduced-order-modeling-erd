"""
End-to-end ROM + SINDy pipeline runner for DynaBench Burgers.

This module provides:
- A single run entry point that executes: load -> snapshots -> POD -> derivatives -> SINDy -> rollout -> metrics
- Artifact writing to a deterministic run directory
"""

from __future__ import annotations
import numpy as np

from dataclasses import dataclass
from typing import Dict, Optional, Sequence


from .artifacts import save_config_yaml, save_json, save_npy
from .config import RunConfig
from .dynabench_io import build_iterator, infer_grid_spacing, load_trajectories
from .metrics import compute_curves, summarize_aggregates
from .paths import config_hash, ensure_run_subdirs, run_dir
from .pod import PODResult, compute_pod, reconstruct_from_pod
from .rollout import RolloutResult, reshape_coeffs_by_trajectory, rollout_one, _coeff_segment
from .snapshot import build_snapshot_matrix
from .sindy_model import SINDyFitResult, fit_sindy_on_coeffs
from .derivs import finite_difference_coeff_derivative
from .plots_and_movies import generate_all_plots_and_movies


@dataclass(frozen=True)
class RunResult:
    """
    Container for a full pipeline run.

    Attributes
    ----------
    run_id : str
        Deterministic hash ID for the run configuration.
    rundir : str
        Path to the run output directory.
    pod : PODResult
        POD result object.
    sindy : SINDyFitResult
        Fitted SINDy model result.
    rollout : RolloutResult
        One rollout result from a chosen trajectory.
    aggregates : dict
        Scalar summary metrics.
    """

    run_id: str
    rundir: str
    pod: PODResult
    sindy: Optional[SINDyFitResult]
    rollout: RolloutResult
    aggregates: Dict[str, float]


def _select_trajectory_ids(
    cfg: RunConfig,
) -> Sequence[int]:
    """
    Select trajectory IDs according to config.

    Parameters
    ----------
    cfg : RunConfig
        Run configuration.

    Returns
    -------
    sequence of int
        Trajectory indices to load.
    """
    if cfg.data.trajectory_ids is not None:
        if len(cfg.data.trajectory_ids) < cfg.data.n_trajectories:
            raise ValueError("trajectory_ids shorter than n_trajectories")
        return list(cfg.data.trajectory_ids)[: cfg.data.n_trajectories]
    return list(range(cfg.data.n_trajectories))


def run(
    cfg: RunConfig,
    write_artifacts: bool = True,
) -> RunResult:
    """
    Execute the full pipeline for one configuration.

    Parameters
    ----------
    cfg : RunConfig
        Full run configuration.
    write_artifacts : bool
        If True, write all artifacts under outputs/ per the run ID.

    Returns
    -------
    RunResult
        Full run result including the fitted model and one rollout.
    """
    print("Building iterator...")
    iterator = build_iterator(
        split=cfg.data.split,
        equation=cfg.data.equation,
        structure=cfg.data.structure,
        resolution=cfg.data.resolution,
        base_path=cfg.data.data_path,
        lookback=cfg.data.lookback,
        # rollout=cfg.data.rollout
        rollout=1,
        # !! caution: bandaid patch
        squeeze_lookback_dim=cfg.data.squeeze_lookback_dim,
        download=cfg.data.download_if_missing,
    )

    print(f"Loading {cfg.data.n_trajectories} trajectories...")
    traj_ids = _select_trajectory_ids(cfg)
    trajectories = load_trajectories(
        iterator=iterator,
        n_trajectories=cfg.data.n_trajectories,
        trajectory_ids=traj_ids,
        time_stride=cfg.data.time_stride,
        time_limit=cfg.data.time_limit,
    )

    print("Building snapshot matrix...")
    dx, dy = infer_grid_spacing(cfg.data.resolution)
    dt = float(getattr(iterator, "dt", 1.0))

    X, layout, mean_state = build_snapshot_matrix(
        trajectories=trajectories,
        center=cfg.pod.center,
        time_stride=1,
        time_limit=None,
    )

    print(f"Computing POD (rank={cfg.pod.rank}, energy={cfg.pod.energy_fraction})...")
    pod = compute_pod(X, rank=cfg.pod.rank, energy_fraction=cfg.pod.energy_fraction)

    n_traj, T_used, _, _, _ = trajectories.shape
    A_by_traj = reshape_coeffs_by_trajectory(
        A=pod.A.T,
        n_trajectories=n_traj,
        T_per_traj=T_used,
    )

    print("Computing derivatives...")
    seg_lengths = [T_used for _ in range(n_traj)]
    deriv_res = finite_difference_coeff_derivative(
        A=pod.A,
        dt=dt,
        segment_lengths=seg_lengths,
        order=2,
    )
    
    sindy: Optional[SINDyFitResult] = None
    rollout: RolloutResult
    aggregates: Dict[str, float] = {}

    if cfg.sindy.enabled:
        print("Fitting SINDy model...")
        sindy = fit_sindy_on_coeffs(
            A_used=deriv_res.A_used,
            dA_dt=deriv_res.dA_dt,
            dt=dt,
            poly_order=cfg.sindy.poly_order,
            include_bias=cfg.sindy.include_bias,
            optimizer_name=cfg.sindy.optimizer,
            optimizer_params=cfg.sindy.optimizer_params,
            parallel=cfg.sindy.parallel,
            parallel_procs=cfg.sindy.parallel_procs,
            feature_names=None,
        )

        print("Running rollout...")
        rollout = rollout_one(
            model=sindy.model,
            U=pod.U,
            layout=layout,
            A_traj=A_by_traj[0],
            dt=dt,
            horizon_steps=cfg.data.rollout,
            mean_state=mean_state,
            start_idx=0,
        )

        print("Computing metrics...")
        err, energy = compute_curves(
            A_true=rollout.A_true,
            A_pred=rollout.A_pred,
            fields_true=rollout.fields_true,
            fields_pred=rollout.fields_pred,
            dx=dx,
            dy=dy,
            equation=cfg.data.equation,
        )
        aggregates = summarize_aggregates(err, energy)
    else:
        print("SINDy is disabled. Skipping fitting, rollout, and metrics.")
        from .snapshot import state_vec_to_fields

        A_true_segment = _coeff_segment(A_by_traj[0], start_idx=0, horizon_steps=cfg.data.rollout)
        q_true_mat = reconstruct_from_pod(pod.U, A_true_segment.T, mean_state=mean_state).T
        fields_true = np.stack(
            [state_vec_to_fields(q_true_mat[i], layout) for i in range(q_true_mat.shape[0])], axis=0
        )
        rollout = RolloutResult(
            A_true=A_true_segment, q_true=q_true_mat, fields_true=fields_true
        )

    
    cfg_dict = cfg.to_dict()
    run_id = config_hash(cfg_dict)

    rundir = run_dir(cfg.outputs_dir, cfg.data.equation, run_id)
    subdirs = ensure_run_subdirs(rundir)

    if write_artifacts:
        print(f"Saving artifacts to {rundir}...")
        save_config_yaml(rundir / "config.yaml", cfg)
        save_json(rundir / "summary.json", {"run_id": run_id, **aggregates})

        save_npy(subdirs["pod"] / "singular_values.npy", pod.s)
        save_npy(subdirs["pod"] / "basis_U.npy", pod.U)
        if mean_state is not None:
            save_npy(subdirs["pod"] / "mean_q.npy", mean_state)
        save_npy(subdirs["pod"] / "coeffs_A.npy", pod.A)

        if sindy is not None:
            save_npy(subdirs["sindy"] / "Xi.npy", sindy.coefficient_matrix)
            save_json(subdirs["sindy"] / "feature_names.json", {"features": sindy.feature_names})
            save_json(
                subdirs["sindy"] / "model_meta.json",
                {
                    "poly_order": cfg.sindy.poly_order,
                    "include_bias": cfg.sindy.include_bias,
                    "optimizer": cfg.sindy.optimizer,
                    "optimizer_params": cfg.sindy.optimizer_params,
                    "n_targets": sindy.n_targets,
                    "n_features": sindy.n_features,
                },
            )

        save_npy(subdirs["rollout"] / "A_true.npy", rollout.A_true)
        if rollout.A_pred is not None:
            save_npy(subdirs["rollout"] / "A_pred.npy", rollout.A_pred)
        save_npy(subdirs["rollout"] / "q_true.npy", rollout.q_true)
        if rollout.q_pred is not None:
            save_npy(subdirs["rollout"] / "q_pred.npy", rollout.q_pred)

        if cfg.sindy.enabled and sindy is not None:
            save_json(
                subdirs["metrics"] / "curves.json",
                {
                    "dt": dt,
                    "coeff_mse": err.coeff_mse.tolist(),
                    "field_mse": err.field_mse.tolist(),
                    "field_l2": err.field_l2.tolist(),
                    "field_rel_l2": err.field_rel_l2.tolist(),
                    "energy_true": energy.energy_true.tolist(),
                    "energy_pred": energy.energy_pred.tolist(),
                },
            )
            save_npy(subdirs["metrics"] / "field_mse.npy", err.field_mse)
            save_json(subdirs["metrics"] / "aggregates.json", aggregates)

        generate_all_plots_and_movies(
            cfg=cfg.plots,
            equation=cfg.data.equation,
            rundir=rundir,
            layout=layout,
            pod=pod,
            X=X,
            sindy=sindy,
            rollout=rollout,
            mean_state=mean_state,
            is_centered=cfg.pod.center,
        )
    
    return RunResult(
        run_id=run_id,
        rundir=str(rundir),
        pod=pod,
        sindy=sindy,
        rollout=rollout,
        aggregates=aggregates,
    )
