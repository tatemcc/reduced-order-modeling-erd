"""Random-shooting MPC controller over a POD + SINDy reduced model."""

from __future__ import annotations

from dataclasses import dataclass
import json
import pickle
from pathlib import Path
from typing import Tuple

import numpy as np

from erd_fipy.config import RunConfig


@dataclass(frozen=True)
class ModelBundle:
    """Serialized model artifacts required by the MPC controller.

    Attributes:
        basis: POD basis matrix ``Phi`` with shape ``(n_state, r)``.
        mean_state: POD centering vector ``x_bar``.
        model: Fitted controlled SINDy model object.
        layout: Spatial field layout ``(N_r, N_phi)``.
        dt: Reduced-model integration step.
        dataset_id: Dataset lineage identifier from model metadata.
        state_clip: Absolute clip bound for reduced-state predictions.
    """

    basis: np.ndarray
    mean_state: np.ndarray
    model: object
    layout: Tuple[int, int]
    dt: float
    dataset_id: str
    state_clip: float


class RandomShootingMPCController:
    """Finite-horizon random-shooting MPC controller using lifted field metrics."""

    def __init__(
        self,
        model_bundle: ModelBundle,
        plant_cfg: RunConfig,
        seed: int = 0,
    ):
        """Initialize controller state and precomputed metric weights.

        Args:
            model_bundle: Loaded POD/SINDy artifacts.
            plant_cfg: Plant configuration that defines bounds and metric constants.
            seed: Random seed for reproducible shooting samples.

        Returns:
            None.
        """

        self.model_bundle = model_bundle
        self.cfg = plant_cfg
        self.rng = np.random.default_rng(seed)

        self.nr, self.nphi = model_bundle.layout
        self.n_points = self.nr * self.nphi

        self.r = np.linspace(
            self.cfg.domain.R_min + 0.5 * (self.cfg.domain.R_max - self.cfg.domain.R_min) / self.cfg.domain.N_r,
            self.cfg.domain.R_max - 0.5 * (self.cfg.domain.R_max - self.cfg.domain.R_min) / self.cfg.domain.N_r,
            self.cfg.domain.N_r,
        )
        self.dr = (self.cfg.domain.R_max - self.cfg.domain.R_min) / self.cfg.domain.N_r

        self.kappa_w_r = self.cfg.wall.kappa_0 * np.exp(
            -((self.cfg.domain.R_max - self.r) ** 2) / (2.0 * self.cfg.wall.delta_w**2)
        )
        self.w_r = self.r * np.exp(
            -((self.r - self.cfg.ring.r_star) ** 2) / (2.0 * self.cfg.metrics.sigma_w**2)
        )

        self.bounds = np.asarray(self.cfg.forcing.u_bounds, dtype=float)

    @staticmethod
    def load_model_bundle(model_run_dir: str | Path) -> ModelBundle:
        """Load a model bundle from a completed ``erd_model`` run folder.

        Args:
            model_run_dir: Path to the model run directory.

        Returns:
            Populated :class:`ModelBundle` instance.
        """

        model_dir = Path(model_run_dir) / "model"

        basis = np.load(model_dir / "basis_U.npy")
        mean_path = model_dir / "mean_state.npy"
        mean_state = np.load(mean_path) if mean_path.is_file() else np.zeros(basis.shape[0])

        with (model_dir / "sindy_model.pkl").open("rb") as f:
            model = pickle.load(f)

        with (model_dir / "model_meta.json").open("r", encoding="utf-8") as f:
            meta = json.load(f)

        nr = int(meta["layout"]["nr"])
        nphi = int(meta["layout"]["nphi"])
        dt = float(meta["dt"])
        state_clip = float(meta.get("state_clip", 50.0))

        dataset_id = str(meta.get("dataset_id", "unknown_dataset"))
        return ModelBundle(
            basis=basis,
            mean_state=mean_state,
            model=model,
            layout=(nr, nphi),
            dt=dt,
            dataset_id=dataset_id,
            state_clip=state_clip,
        )

    def project(self, x: np.ndarray) -> np.ndarray:
        """Project a full state vector into reduced coordinates.

        Args:
            x: Full-order stacked state vector.

        Returns:
            Reduced state vector ``a``.
        """

        return self.model_bundle.basis.T @ (x - self.model_bundle.mean_state)

    def lift(self, a: np.ndarray) -> np.ndarray:
        """Lift reduced coordinates back to full-order stacked state space.

        Args:
            a: Reduced coordinates.

        Returns:
            Reconstructed full-order state vector.
        """

        return self.model_bundle.mean_state + self.model_bundle.basis @ a

    def _decode_n(self, x: np.ndarray) -> np.ndarray:
        """Extract ``n(r,phi)`` from stacked full-order state ``[n, omega]``."""

        return x[: self.n_points].reshape(self.nr, self.nphi)

    def _metrics_from_n(self, n: np.ndarray, u: np.ndarray) -> tuple[float, float, float, float, float]:
        """Evaluate control objective terms from a lifted ``n`` field and control.

        Args:
            n: Density field ``n(r, phi)``.
            u: Control vector for the current stage.

        Returns:
            Tuple ``(E_wob, L_w, sigma_r, sigma_r_sq, P_ctrl)``.
        """

        n_hat = np.fft.fft(n, axis=1) / n.shape[1]
        e1 = 0.5 * np.sum(self.w_r * np.abs(n_hat[:, 1]) ** 2) * self.dr
        e2 = 0.5 * np.sum(self.w_r * np.abs(n_hat[:, 2]) ** 2) * self.dr
        e_wob = float(e1 + e2)

        l_w = float(np.sum(self.kappa_w_r[:, None] * n * (self.r[:, None] * self.dr * (2.0 * np.pi / self.nphi))))

        nbar = np.mean(n, axis=1)
        mass = float(np.sum(nbar * self.r) * self.dr)
        mass = max(mass, 1e-12)
        r_mean = float(np.sum(self.r * nbar * self.r) * self.dr / mass)
        sigma_r2 = float(np.sum(((self.r - r_mean) ** 2) * nbar * self.r) * self.dr / mass)

        p_ctrl = float(
            u[0] ** 2
            + self.cfg.metrics.lambda_1 * (u[1] ** 2 + u[2] ** 2)
            + self.cfg.metrics.lambda_2 * (u[3] ** 2 + u[4] ** 2)
        )

        sigma_r2 = max(sigma_r2, 0.0)
        return e_wob, l_w, float(np.sqrt(sigma_r2)), sigma_r2, p_ctrl

    def _predict_step(self, a: np.ndarray, u: np.ndarray) -> np.ndarray:
        """Advance one reduced step with explicit Euler and finite-value guards.

        Args:
            a: Current reduced state.
            u: Control vector for this step.

        Returns:
            Next reduced state, or ``NaN`` vector when prediction is invalid.
        """

        try:
            adot = np.asarray(
                self.model_bundle.model.predict(a.reshape(1, -1), u=u.reshape(1, -1)),
                dtype=float,
            ).reshape(-1)
        except Exception:
            return np.full_like(a, np.nan, dtype=float)

        if not np.all(np.isfinite(adot)):
            return np.full_like(a, np.nan, dtype=float)

        a_next = a + self.model_bundle.dt * adot
        a_next = np.clip(a_next, -self.model_bundle.state_clip, self.model_bundle.state_clip)
        if not np.all(np.isfinite(a_next)):
            return np.full_like(a, np.nan, dtype=float)
        return a_next

    def select_u(self, a_k: np.ndarray, u_prev: np.ndarray, t_k: float) -> np.ndarray:
        """Choose the next control action via random-shooting MPC.

        Args:
            a_k: Current reduced state.
            u_prev: Previously applied control vector.
            t_k: Current simulation time (unused in base policy, kept for API parity).

        Returns:
            First control action from the minimum-cost candidate sequence.
        """

        _ = t_k
        H = self.cfg.control.H
        N_shoot = self.cfg.control.N_shoot
        w = self.cfg.control.weights
        rate_penalty = float(self.cfg.control.rate_penalty)

        x_ref = self.lift(np.asarray(a_k, dtype=float))
        n_ref = self._decode_n(x_ref)
        e_ref, l_ref, sigma_ref, _, _ = self._metrics_from_n(n_ref, np.asarray(u_prev, dtype=float))

        best_cost = np.inf
        best_u0 = np.asarray(u_prev, dtype=float)

        candidate0 = np.repeat(np.asarray(u_prev, dtype=float)[None, :], H, axis=0)
        sample_scale0 = np.asarray([0.05, 0.20, 0.20, 0.20, 0.20], dtype=float) * self.bounds
        sample_step = np.asarray([0.03, 0.12, 0.12, 0.12, 0.12], dtype=float) * self.bounds

        for n_cand in range(N_shoot + 1):
            if n_cand == 0:
                seq = candidate0.copy()
            else:
                # Sample smooth bounded sequences around the previous action
                # so the controller explores local corrections, not only
                # globally random bang-bang trajectories.
                seq = np.empty((H, 5), dtype=float)
                seq[0] = np.asarray(u_prev, dtype=float) + self.rng.normal(loc=0.0, scale=sample_scale0)
                for h in range(1, H):
                    seq[h] = seq[h - 1] + self.rng.normal(loc=0.0, scale=sample_step)
                seq = np.clip(seq, -self.bounds, self.bounds)

            a = np.asarray(a_k, dtype=float).copy()
            up = np.asarray(u_prev, dtype=float).copy()
            cost = 0.0
            l_prev = l_ref

            for h in range(H):
                u = np.clip(seq[h], -self.bounds, self.bounds)
                a = self._predict_step(a, u)
                if not np.all(np.isfinite(a)):
                    cost = np.inf
                    break
                x = self.lift(a)
                n_field = self._decode_n(x)

                e_wob, l_w, sigma_r, sigma_r2, p_ctrl = self._metrics_from_n(n_field, u)
                if not all(np.isfinite(v) for v in (e_wob, l_w, sigma_r, sigma_r2, p_ctrl)):
                    cost = np.inf
                    break
                rate = np.sum((u - up) ** 2)
                wob_growth = max(0.0, e_wob - e_ref)
                sig_growth = max(0.0, sigma_r - sigma_ref)
                lw_slope = max(0.0, l_w - l_prev)
                cost += (
                    w.w_wob_abs * (e_wob**2)
                    + w.w_wob_growth * (wob_growth**2)
                    + w.w_wall * l_w
                    + w.w_coh * sigma_r2
                    + w.w_coh_growth * (sig_growth**2)
                    + w.w_wall_slope * (lw_slope**2)
                    + w.w_pow * p_ctrl
                    + w.w_rate * rate_penalty * rate
                )
                up = u
                l_prev = l_w

            if cost < best_cost:
                best_cost = cost
                best_u0 = np.clip(seq[0], -self.bounds, self.bounds)

        return best_u0

    def select_u_from_state(self, x_k: np.ndarray, u_prev: np.ndarray, t_k: float) -> np.ndarray:
        """Project full state to reduced space and select control action.

        Args:
            x_k: Current full-order state vector.
            u_prev: Previously applied control vector.
            t_k: Current simulation time.

        Returns:
            Selected control vector for the next plant step.
        """

        a_k = self.project(np.nan_to_num(np.asarray(x_k, dtype=float), nan=0.0, posinf=1e6, neginf=-1e6))
        return self.select_u(a_k, u_prev, t_k)
