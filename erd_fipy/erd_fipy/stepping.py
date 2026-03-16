"""Time integration utilities for the ERD toy plant.

The :class:`ERDPlant` class orchestrates mesh/state initialization, one-step
updates through :class:`ToyERDOperators`, metric evaluation, and optional
artifact writing through :class:`RunWriter`.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable, Dict, List, Optional

import numpy as np

from .config import RunConfig
from .fields import assign_grid, make_initial_state, stack_state, variable_to_grid
from .io import RunWriter
from .mesh import make_mesh
from .metrics import add_degradation_metrics, add_window_efficiency, compute_metric_state
from .pdes import ToyERDOperators


ControlSchedule = Callable[[int, float, np.ndarray, np.ndarray], np.ndarray]
ControllerCallback = Callable[[np.ndarray, np.ndarray, float, int], np.ndarray]


@dataclass(frozen=True)
class PlantRunResult:
    """Summary object returned by :meth:`ERDPlant.run`.

    Attributes:
        run_dir: Output directory used for artifacts (empty when disabled).
        curves: Time-series metric curves accumulated during the run.
        summary: Small dictionary of aggregate scalar metrics.
    """

    run_dir: str
    curves: Dict[str, List[float]]
    summary: Dict[str, float]



def default_open_loop_control(cfg: RunConfig) -> np.ndarray:
    """Return the default open-loop control vector.

    Args:
        cfg: Plant run configuration containing forcing defaults.

    Returns:
        Array ``[u0, u1, u2, u3, u4]`` with only ``u0`` nonzero.
    """

    return np.array([cfg.forcing.drive_u0_base, 0.0, 0.0, 0.0, 0.0], dtype=float)


class ERDPlant:
    """Toy ERD plant wrapper exposing step-wise and full-run execution APIs."""

    def __init__(self, cfg: RunConfig):
        """Initialize mesh, state variables, and PDE operators.

        Args:
            cfg: Full plant run configuration.

        Returns:
            None.
        """

        self.cfg = cfg
        self.bundle = make_mesh(cfg)
        self.state = make_initial_state(self.bundle, cfg)
        self.ops = ToyERDOperators(cfg, self.bundle)

    def get_fields(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Read current ``(n, omega, psi)`` arrays from FiPy variables.

        Args:
            None.

        Returns:
            Tuple ``(n, omega, psi)`` with shape ``(N_r, N_phi)`` each.
        """

        n = variable_to_grid(self.state.n, self.bundle.shape)
        omega = variable_to_grid(self.state.omega, self.bundle.shape)
        psi = variable_to_grid(self.state.psi, self.bundle.shape)
        return n, omega, psi

    def set_fields(self, n: np.ndarray, omega: np.ndarray, psi: np.ndarray) -> None:
        """Write arrays back into FiPy state variables.

        Args:
            n: Density field array.
            omega: Vorticity field array.
            psi: Streamfunction field array.

        Returns:
            None.
        """

        assign_grid(self.state.n, n)
        assign_grid(self.state.omega, omega)
        assign_grid(self.state.psi, psi)

    def velocity_magnitude_from_omega(self, omega: np.ndarray) -> np.ndarray:
        """Compute ``|u|`` from a vorticity field using the diagnostic Poisson solve.

        Args:
            omega: Vorticity field ``omega(r, phi)``.

        Returns:
            Velocity magnitude field with shape ``(N_r, N_phi)``.
        """

        psi = self.ops.solve_poisson(omega)
        u_r, u_phi = self.ops.velocity(psi)
        return np.sqrt(u_r**2 + u_phi**2)

    def velocity_components_from_omega(self, omega: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Compute ``(u_r, u_phi)`` from a vorticity field.

        Args:
            omega: Vorticity field ``omega(r, phi)``.

        Returns:
            Tuple ``(u_r, u_phi)`` derived from the diagnostic Poisson solve.
        """

        psi = self.ops.solve_poisson(omega)
        return self.ops.velocity(psi)

    def state_vector(self) -> np.ndarray:
        """Return the stacked state vector ``[vec(n), vec(omega)]``.

        Args:
            None.

        Returns:
            One-dimensional state vector used by ROM/controller code.
        """

        n, omega, _ = self.get_fields()
        return stack_state(n, omega)

    def step(self, u: np.ndarray, t: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, float]]:
        """Advance the plant by one time step and compute instantaneous metrics.

        Args:
            u: Control vector at the current step.
            t: Current simulation time.

        Returns:
            Tuple ``(n_next, omega_next, psi_next, metric_dict)``.
        """

        n, omega, _ = self.get_fields()
        n_next, omega_next, psi_next, diag = self.ops.step(n, omega, u=u, t=t, dt=self.cfg.time.dt)
        self.set_fields(n_next, omega_next, psi_next)

        metrics = compute_metric_state(
            cfg=self.cfg,
            n=n_next,
            kappa_w_r=self.ops.kappa_w_r,
            area_weights=self.bundle.area_weights,
            u_r=diag.u_r,
            u_phi=diag.u_phi,
            u=u,
            r=self.bundle.r,
            dr=self.bundle.dr,
        )
        return n_next, omega_next, psi_next, {
            "J_prof": metrics.J_prof,
            "E_low": metrics.E_low,
            "L_w": metrics.L_w,
            "sigma_r": metrics.sigma_r,
            "r_mean": metrics.r_mean,
            "M": metrics.M,
            "E_u": metrics.E_u,
            "P_ctrl": metrics.P_ctrl,
        }

    def run(
        self,
        control_schedule: Optional[ControlSchedule] = None,
        controller_callback: Optional[ControllerCallback] = None,
        run_dir: Optional[str] = None,
        tag_override: Optional[str] = None,
        write_artifacts: bool = True,
    ) -> PlantRunResult:
        """Execute a full simulation with either open-loop or callback control.

        Args:
            control_schedule: Optional open-loop schedule ``u = f(k, t, n, omega)``.
            controller_callback: Optional closed-loop callback ``u = f(x, u_prev, t, k)``.
            run_dir: Optional explicit run directory override.
            tag_override: Optional run-folder tag override.
            write_artifacts: If ``True``, write HDF5/plots/GIF/log outputs.

        Returns:
            :class:`PlantRunResult` with the output path, curves, and aggregates.
        """

        if (control_schedule is not None) and (controller_callback is not None):
            raise ValueError("Provide either control_schedule or controller_callback, not both")

        cfg = self.cfg
        if tag_override is not None:
            cfg = replace(cfg, output=replace(cfg.output, tag=tag_override))

        writer = RunWriter(cfg, run_dir=run_dir) if write_artifacts else None
        if writer is not None:
            writer.log(f"Starting run: steps={cfg.time.n_steps}, dt={cfg.time.dt:.4e}")
            writer.write_grid(self.bundle.r, self.bundle.phi)

        curves = {
            "t": [],
            "J_prof": [],
            "E_low": [],
            "L_w": [],
            "sigma_r": [],
            "r_mean": [],
            "M": [],
            "E_u": [],
            "P_ctrl": [],
        }

        n, omega, _ = self.get_fields()
        if writer is not None:
            u_r0, u_phi0 = self.velocity_components_from_omega(omega)
            writer.append_snapshot(0.0, n, omega, u_mag=np.sqrt(u_r0**2 + u_phi0**2), u_r=u_r0, u_phi=u_phi0)

        u_prev = default_open_loop_control(cfg)

        for k in range(cfg.time.n_steps):
            t = k * cfg.time.dt
            if controller_callback is not None:
                xk = self.state_vector()
                u = np.asarray(controller_callback(xk, u_prev.copy(), t, k), dtype=float).reshape(5)
            elif control_schedule is not None:
                n_cur, omega_cur, _ = self.get_fields()
                u = np.asarray(control_schedule(k, t, n_cur, omega_cur), dtype=float).reshape(5)
            else:
                u = default_open_loop_control(cfg)

            bounds = np.asarray(cfg.forcing.u_bounds, dtype=float)
            u = np.clip(u, -bounds, bounds)

            if writer is not None:
                writer.append_control(t, u)

            n_next, omega_next, _, metric = self.step(u=u, t=t)
            u_prev = u

            t_next = (k + 1) * cfg.time.dt
            if writer is not None:
                writer.append_metrics(t_next, metric)
                if ((k + 1) % cfg.time.snapshot_every) == 0 or (k == cfg.time.n_steps - 1):
                    u_r_snap, u_phi_snap = self.velocity_components_from_omega(omega_next)
                    writer.append_snapshot(
                        t_next,
                        n_next,
                        omega_next,
                        u_mag=np.sqrt(u_r_snap**2 + u_phi_snap**2),
                        u_r=u_r_snap,
                        u_phi=u_phi_snap,
                    )

            curves["t"].append(float(t_next))
            for key in ("J_prof", "E_low", "L_w", "sigma_r", "r_mean", "M", "E_u", "P_ctrl"):
                curves[key].append(float(metric[key]))

        if writer is not None:
            # Persist a terminal control sample so control timestamps include
            # the final snapshot time.
            writer.append_control(cfg.time.T_final, u_prev)

        if writer is not None:
            curves = writer.curves

        add_window_efficiency(curves, cfg)
        add_degradation_metrics(curves, cfg)
        summary = {
            "final_J_prof": float(curves["J_prof"][-1]) if curves["J_prof"] else 0.0,
            "mean_J_prof": float(np.mean(curves["J_prof"])) if curves["J_prof"] else 0.0,
            "final_E_low": float(curves["E_low"][-1]) if curves["E_low"] else 0.0,
            "mean_E_low": float(np.mean(curves["E_low"])) if curves["E_low"] else 0.0,
            "final_L_w": float(curves["L_w"][-1]) if curves["L_w"] else 0.0,
            "final_sigma_r": float(curves["sigma_r"][-1]) if curves["sigma_r"] else 0.0,
            "final_L_w_cum": float(curves["L_w_cum"][-1]) if curves["L_w_cum"] else 0.0,
            "final_badness_score": float(curves["badness_score"][-1]) if curves["badness_score"] else 0.0,
            "mean_P_ctrl": float(np.mean(curves["P_ctrl"])) if curves["P_ctrl"] else 0.0,
        }

        if writer is not None:
            writer.finalize(curves, summary)
            return PlantRunResult(run_dir=str(writer.run_dir), curves=curves, summary=summary)

        return PlantRunResult(run_dir="", curves=curves, summary=summary)
