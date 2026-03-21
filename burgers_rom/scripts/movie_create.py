from __future__ import annotations

import numpy as np
import math
from datetime import datetime

from burgers_rom.snapshot import state_vec_to_fields, SnapshotLayout
from burgers_rom.config import OUTPUTS_DIR
from burgers_rom.plotting import (
    DEFAULT_FPS,
    animate_channels,
    animate_quiver,
    animate_vorticity,
    animate_streamlines,
    animate_moving_dashes,
    animate_amplitudes,
    plot_pod_basis,
)


def load_pod_data(run_id: str, equation: str = "burgers"):
    pod_dir = OUTPUTS_DIR / equation / run_id / "pod"

    if not pod_dir.exists():
        raise FileNotFoundError(f"POD directory not found: {pod_dir}")

    basis_U = np.load(pod_dir / "basis_U.npy", allow_pickle=True)
    coeffs_A = np.load(pod_dir / "coeffs_A.npy", allow_pickle=True)
    mean_q = np.load(pod_dir / "mean_q.npy", allow_pickle=True)
    singular_values = np.load(pod_dir / "singular_values.npy", allow_pickle=True)
    return basis_U, coeffs_A, mean_q, singular_values


def load_rollout_data(run_id: str, equation: str = "burgers"):
    rollout_dir = OUTPUTS_DIR / equation / run_id / "rollout"

    if not rollout_dir.exists():
        raise FileNotFoundError(f"Rollout directory not found: {rollout_dir}")

    q_true = np.load(rollout_dir / "q_true.npy")
    q_pred = np.load(rollout_dir / "q_pred.npy")
    A_true = np.load(rollout_dir / "A_true.npy")
    A_pred = np.load(rollout_dir / "A_pred.npy")
    return q_true, q_pred, A_true, A_pred


def reconstruct_fields(
    basis_U_flat: np.ndarray, mean_q_flat: np.ndarray, coeffs_A: np.ndarray, n_components: int = 2
) -> tuple[np.ndarray, np.ndarray, SnapshotLayout]:
    """
    Reinflates flat POD basis and mean state into spatial field tensors.

    Infers grid dimensions (ny, nx) from the state vector size, assuming a square grid.

    Parameters
    ----------
    basis_U_flat : np.ndarray
        Flat POD basis matrix of shape (n_state, r).
    mean_q_flat : np.ndarray
        Flat mean state vector of shape (n_state,).
    coeffs_A : np.ndarray
        POD coefficients of shape (r, nt).
    n_components : int
        Number of vector components (default 2 for Burgers).

    Returns
    -------
    basis_U_fields : np.ndarray
        Spatial basis modes of shape (r, C, ny, nx).
    mean_q_fields : np.ndarray
        Spatial mean field of shape (C, ny, nx).
    layout : SnapshotLayout
        Object containing grid dimensions (ny, nx) and component count.
    """
    # 1. Extract shapes
    n_state, n_modes = basis_U_flat.shape
    n_modes_coeffs, nt = coeffs_A.shape

    # 2. Validate consistency
    if n_modes != n_modes_coeffs:
        raise ValueError(
            f"Dimension mismatch: basis has {n_modes} modes, coeffs have {n_modes_coeffs}"
        )

    # 3. Infer spatial resolution (ny, nx)
    # n_state = C * ny * nx
    n_points = n_state // n_components
    ny = int(math.isqrt(n_points))
    nx = ny

    if n_components * ny * nx != n_state:
        raise ValueError(
            f"State vector length {n_state} is not compatible with C={n_components} on a square grid."
        )

    layout = SnapshotLayout(ny=ny, nx=nx, n_components=n_components)

    # 4. Reshape basis: (n_state, r) -> (r, C, ny, nx)
    # Transpose to (r, n_state) to iterate over modes
    basis_U_T = basis_U_flat.T
    basis_U_fields = np.array([state_vec_to_fields(basis_U_T[i], layout) for i in range(n_modes)])

    # 5. Reshape mean state: (n_state,) -> (C, ny, nx)
    mean_q_fields = state_vec_to_fields(mean_q_flat, layout)

    return basis_U_fields, mean_q_fields, layout


def reshape_trajectory(q: np.ndarray, n_components: int = 2):
    """
    Reshape trajectory (T, n_state) to (T, C, ny, nx).
    """
    T, n_state = q.shape
    n_points = n_state // n_components
    ny = int(math.isqrt(n_points))
    nx = ny

    if n_components * ny * nx != n_state:
        raise ValueError(
            f"State vector length {n_state} is not compatible with C={n_components} on a square grid."
        )

    layout = SnapshotLayout(ny=ny, nx=nx, n_components=n_components)
    fields = q.reshape((T, n_components, ny, nx), order="C")
    return fields, layout


def reconstruct_from_k_modes(
    U: np.ndarray,
    A: np.ndarray,
    k: int,
    layout: SnapshotLayout,
    mean_state: np.ndarray | None = None,
) -> np.ndarray:
    """
    Reconstruct trajectory fields using only the first k POD modes.

    Parameters
    ----------
    U : np.ndarray
        Flat POD basis of shape (n_state, r).
    A : np.ndarray
        POD coefficients of shape (r, T).
    k : int
        Number of modes to use.
    layout : SnapshotLayout
        Layout object for reshaping.
    mean_state : np.ndarray, optional
        Mean state vector of shape (n_state,).

    Returns
    -------
    np.ndarray
        Reconstructed fields of shape (T, C, ny, nx).
    """
    r, T = A.shape
    k = min(k, r)

    U_k = U[:, :k]
    A_k = A[:k, :]

    X_rec = U_k @ A_k
    if mean_state is not None:
        X_rec += mean_state[:, None]

    fields = np.empty((T, layout.n_components, layout.ny, layout.nx), dtype=X_rec.dtype)
    for t in range(T):
        fields[t] = state_vec_to_fields(X_rec[:, t], layout)

    return fields


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Generate animations for Burgers ROM rollout.")
    parser.add_argument("--run_id", type=str, required=True, help="Run ID hash.")
    parser.add_argument(
        "--equation", type=str, default="burgers", help="Equation name (default: burgers)."
    )
    parser.add_argument(
        "--style",
        type=str,
        default="channels",
        choices=["channels", "quiver", "vorticity", "streamlines", "dashes", "amplitudes", "all"],
        help="Animation style to generate.",
    )
    parser.add_argument("--fps", type=int, default=DEFAULT_FPS, help="Frames per second.")
    parser.add_argument(
        "--format", type=str, default="mp4", choices=["mp4", "gif"], help="Output format."
    )
    parser.add_argument("--n_bases", type=int, default=0, help="Number of POD bases to plot.")
    parser.add_argument(
        "--max_modes", type=int, default=0, help="Max modes for reconstruction animations."
    )
    args = parser.parse_args()

    print(f"Loading data for run: {args.run_id}")
    try:
        q_true_flat, q_pred_flat, A_true, A_pred = load_rollout_data(
            args.run_id, equation=args.equation
        )
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return

    # Reshape
    match args.equation:
        case "advection":
            n_components = 1
        case "burgers":
            n_components = 2
        case "gasdynamics":
            n_components = 4
        case "kuramotosivashinsky":
            n_components = 1
        case "reactiondiffusion":
            n_components = 2
        case "wave":
            n_components = 1
        case _:
            raise Exception("Unrecognized args.equation")
    q_true, layout = reshape_trajectory(q_true_flat, n_components=n_components)
    q_pred, _ = reshape_trajectory(q_pred_flat, n_components=n_components)

    output_dir = OUTPUTS_DIR / args.equation / args.run_id / "movies"
    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    GENERAL_STYLES = ["channels", "amplitudes"]
    BURGERS_STYLES = ["quiver", "vorticity", "streamlines", "dashes"]

    if args.style == "all":
        styles_to_run = GENERAL_STYLES.copy()
        if args.equation == "burgers":
            styles_to_run.extend(BURGERS_STYLES)
    else:
        if args.style in BURGERS_STYLES and args.equation != "burgers":
            print(f"Warning: Style '{args.style}' is specific to the Burgers equation. Skipping.")
            styles_to_run = []
        else:
            styles_to_run = [args.style]

    for style in styles_to_run:
        if style == "channels":
            print("Generating channel animation...")
            animate_channels(
                q_true,
                q_pred,
                output_dir / f"{timestamp}_{args.run_id}_truevpred_channels.{args.format}",
                fps=args.fps,
            )
        elif style == "quiver":
            print("Generating quiver animation...")
            animate_quiver(
                q_true,
                q_pred,
                output_dir / f"{timestamp}_{args.run_id}_truevpred_quiver.{args.format}",
                fps=args.fps,
            )
        elif style == "vorticity":
            print("Generating vorticity animation...")
            animate_vorticity(
                q_true,
                q_pred,
                output_dir / f"{timestamp}_{args.run_id}_truevpred_vorticity.{args.format}",
                fps=args.fps,
            )
        elif style == "streamlines":
            print("Generating streamline animation...")
            animate_streamlines(
                q_true,
                q_pred,
                output_dir / f"{timestamp}_{args.run_id}_truevpred_streamlines.{args.format}",
                fps=args.fps,
            )
        elif style == "amplitudes":
            print("Generating amplitude animation...")
            animate_amplitudes(
                A_true,
                A_pred,
                output_dir / f"{timestamp}_{args.run_id}_truevpred_amplitudes.{args.format}",
                fps=args.fps,
            )
        elif style == "dashes":
            print("Generating moving dashes animation...")
            animate_moving_dashes(
                q_true,
                q_pred,
                output_dir / f"{timestamp}_{args.run_id}_truevpred_dashes.{args.format}",
                fps=args.fps,
            )

    # 1) Plot POD Bases
    if args.n_bases > 0:
        print(f"Plotting first {args.n_bases} POD bases...")
        basis_U, _, _, _ = load_pod_data(args.run_id, equation=args.equation)
        bases_dir = OUTPUTS_DIR / args.equation / args.run_id / "figures" / "bases"
        bases_dir.mkdir(parents=True, exist_ok=True)

        for i in range(min(args.n_bases, basis_U.shape[1])):
            # basis_U is (n_state, r). Extract column i and reshape.
            mode_field = state_vec_to_fields(basis_U[:, i], layout)
            if args.equation == "burgers":
                plot_pod_basis(mode_field, bases_dir / f"mode_{i:03d}.png", title=f"Mode {i}")
            else:
                pass # plot_pod_basis currently only supports C=2
        print(f"Saved basis plots to {bases_dir}")

    # 2) Reconstruction Animations
    if args.max_modes > 0:
        print(f"Generating reconstruction animations up to {args.max_modes} modes...")
        basis_U, _, mean_q, _ = load_pod_data(args.run_id, equation=args.equation)
        recon_dir = OUTPUTS_DIR / args.equation / args.run_id / "movies" / "reconstruction"
        recon_dir.mkdir(parents=True, exist_ok=True)

        # Transpose A_pred to (r, T) for reconstruction utility
        A_pred_T = A_pred.T

        for k in range(1, min(args.max_modes, basis_U.shape[1]) + 1):
            # Reconstruct using k modes
            q_recon = reconstruct_from_k_modes(basis_U, A_pred_T, k, layout, mean_state=mean_q)

            if args.equation == "burgers":
                out_path = (
                    recon_dir / f"{timestamp}_{args.run_id}_recon_k{k:02d}_vorticity.{args.format}"
                )
                animate_vorticity(q_true, q_recon, out_path, fps=args.fps)
            else:
                out_path = recon_dir / f"{timestamp}_{args.run_id}_recon_k{k:02d}_channels.{args.format}"
                animate_channels(q_true, q_recon, out_path, fps=args.fps)
        print(f"Saved reconstruction animations to {recon_dir}")


if __name__ == "__main__":
    main()
