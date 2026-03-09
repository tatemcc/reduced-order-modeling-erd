import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from datetime import datetime

from burgers_rom.snapshot import state_vec_to_fields, SnapshotLayout
from burgers_rom.config import OUTPUTS_DIR

DEFAULT_FPS = 10


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
    return q_true, q_pred


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


def compute_vorticity(q: np.ndarray, L: float = 1.0) -> np.ndarray:
    """
    Compute scalar vorticity field from velocity field.

    omega = dv/dx - du/dy

    Parameters
    ----------
    q : np.ndarray
        Velocity field of shape (T, 2, ny, nx).
        q[:, 0] is u (x-velocity), q[:, 1] is v (y-velocity).
    L : float
        Domain size length (assumed square [0, L]x[0, L]).

    Returns
    -------
    omega : np.ndarray
        Vorticity field of shape (T, ny, nx).
    """
    T, C, ny, nx = q.shape
    dx = L / (nx - 1) if nx > 1 else 1.0
    dy = L / (ny - 1) if ny > 1 else 1.0

    u = q[:, 0]
    v = q[:, 1]

    # du/dy: gradient of u along axis 1 (y-axis)
    dudy = np.gradient(u, axis=1) / dy

    # dv/dx: gradient of v along axis 2 (x-axis)
    dvdx = np.gradient(v, axis=2) / dx

    omega = dvdx - dudy
    return omega


def animate_channels(q_true, q_pred, output_path, fps=DEFAULT_FPS):
    """
    Animate separate velocity components (channels) side-by-side.

    Creates a grid of subplots where rows correspond to channels (u, v)
    and columns correspond to True vs Predicted fields.

    Parameters
    ----------
    q_true : np.ndarray
        Ground truth trajectory of shape (T, C, ny, nx).
    q_pred : np.ndarray
        Predicted trajectory of shape (T, C, ny, nx).
    output_path : Path
        Destination path for the .mp4 file.
    fps : int
        Frames per second for the animation.
    """
    T, C, ny, nx = q_true.shape
    fig, axes = plt.subplots(C, 2, figsize=(10, 4 * C), squeeze=False)

    vmins = [min(q_true[:, c].min(), q_pred[:, c].min()) for c in range(C)]
    vmaxs = [max(q_true[:, c].max(), q_pred[:, c].max()) for c in range(C)]

    imgs = []
    for c in range(C):
        # True
        ax_true = axes[c, 0]
        ax_true.set_title(f"Channel {c} True")
        img_true = ax_true.imshow(
            q_true[0, c], origin="lower", vmin=vmins[c], vmax=vmaxs[c], cmap="viridis"
        )
        fig.colorbar(img_true, ax=ax_true)
        imgs.append(img_true)

        # Pred
        ax_pred = axes[c, 1]
        ax_pred.set_title(f"Channel {c} Pred")
        img_pred = ax_pred.imshow(
            q_pred[0, c], origin="lower", vmin=vmins[c], vmax=vmaxs[c], cmap="viridis"
        )
        fig.colorbar(img_pred, ax=ax_pred)
        imgs.append(img_pred)

    title = fig.suptitle("Time: 0")

    def update(frame):
        title.set_text(f"Time: {frame}")
        for c in range(C):
            imgs[2 * c].set_data(q_true[frame, c])
            imgs[2 * c + 1].set_data(q_pred[frame, c])
        return imgs + [title]

    ani = animation.FuncAnimation(fig, update, frames=T, blit=True, interval=1000 / fps)
    writer = "ffmpeg" if output_path.suffix == ".mp4" else "pillow"
    ani.save(output_path, writer=writer, fps=fps)
    plt.close(fig)
    print(f"Saved channel animation to {output_path}")


def animate_quiver(q_true, q_pred, output_path, fps=DEFAULT_FPS, stride=2):
    """
    Animate velocity vector fields using quiver plots side-by-side.

    Displays the vector field (u, v) as arrows on the grid.
    Left panel: True field. Right panel: Predicted field.

    Parameters
    ----------
    q_true : np.ndarray
        Ground truth trajectory of shape (T, C, ny, nx).
    q_pred : np.ndarray
        Predicted trajectory of shape (T, C, ny, nx).
    output_path : Path
        Destination path for the .mp4 file.
    fps : int
        Frames per second.
    stride : int
        Spatial subsampling factor to reduce arrow density.
    """
    T, C, ny, nx = q_true.shape
    if C < 2:
        print("Skipping quiver: need at least 2 components.")
        return

    Y, X = np.mgrid[0:ny, 0:nx]
    X_s = X[::stride, ::stride]
    Y_s = Y[::stride, ::stride]

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    # True
    ax_true = axes[0]
    ax_true.set_title("True Velocity Field")
    ax_true.set_xlim(0, nx - 1)
    ax_true.set_ylim(0, ny - 1)
    Q_true = ax_true.quiver(
        X_s, Y_s, q_true[0, 0, ::stride, ::stride], q_true[0, 1, ::stride, ::stride], scale=None
    )

    # Pred
    ax_pred = axes[1]
    ax_pred.set_title("Pred Velocity Field")
    ax_pred.set_xlim(0, nx - 1)
    ax_pred.set_ylim(0, ny - 1)
    Q_pred = ax_pred.quiver(
        X_s, Y_s, q_pred[0, 0, ::stride, ::stride], q_pred[0, 1, ::stride, ::stride], scale=None
    )

    title = fig.suptitle("Time: 0")

    def update(frame):
        title.set_text(f"Time: {frame}")
        Q_true.set_UVC(q_true[frame, 0, ::stride, ::stride], q_true[frame, 1, ::stride, ::stride])
        Q_pred.set_UVC(q_pred[frame, 0, ::stride, ::stride], q_pred[frame, 1, ::stride, ::stride])
        return Q_true, Q_pred

    ani = animation.FuncAnimation(fig, update, frames=T, blit=False, interval=1000 / fps)
    writer = "ffmpeg" if output_path.suffix == ".mp4" else "pillow"
    ani.save(output_path, writer=writer, fps=fps)
    plt.close(fig)
    print(f"Saved quiver animation to {output_path}")


def animate_vorticity(q_true, q_pred, output_path, fps=DEFAULT_FPS):
    """
    Animate scalar vorticity fields side-by-side.

    Computes the signed z-component of the curl (dv/dx - du/dy) and plots it
    using a diverging colormap centered at 0.

    Parameters
    ----------
    q_true : np.ndarray
        Ground truth trajectory of shape (T, C, ny, nx).
    q_pred : np.ndarray
        Predicted trajectory of shape (T, C, ny, nx).
    output_path : Path
        Destination path for the .mp4 file.
    fps : int
        Frames per second.
    """

    # Compute vorticity
    w_true = compute_vorticity(q_true)
    w_pred = compute_vorticity(q_pred)

    T, ny, nx = w_true.shape

    # Determine color limits centered on 0 for diverging map
    max_val = max(np.abs(w_true).max(), np.abs(w_pred).max())
    vmin, vmax = -max_val, max_val

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # True
    ax_true = axes[0]
    ax_true.set_title("True Vorticity")
    img_true = ax_true.imshow(w_true[0], origin="lower", vmin=vmin, vmax=vmax, cmap="RdBu_r")
    fig.colorbar(img_true, ax=ax_true)

    # Pred
    ax_pred = axes[1]
    ax_pred.set_title("Pred Vorticity")
    img_pred = ax_pred.imshow(w_pred[0], origin="lower", vmin=vmin, vmax=vmax, cmap="RdBu_r")
    fig.colorbar(img_pred, ax=ax_pred)

    title = fig.suptitle("Time: 0")

    def update(frame):
        title.set_text(f"Time: {frame}")
        img_true.set_data(w_true[frame])
        img_pred.set_data(w_pred[frame])
        return img_true, img_pred, title

    ani = animation.FuncAnimation(fig, update, frames=T, blit=True, interval=1000 / fps)
    writer = "ffmpeg" if output_path.suffix == ".mp4" else "pillow"
    ani.save(output_path, writer=writer, fps=fps)
    plt.close(fig)
    print(f"Saved vorticity animation to {output_path}")


def animate_streamlines(q_true, q_pred, output_path, fps=DEFAULT_FPS):
    """
    Animate streamlines side-by-side.

    Streamlines are colored and weighted by velocity magnitude.
    Note: Streamplot animations can be slow to render.

    Parameters
    ----------
    q_true : np.ndarray
        Ground truth trajectory of shape (T, C, ny, nx).
    q_pred : np.ndarray
        Predicted trajectory of shape (T, C, ny, nx).
    output_path : Path
        Destination path for the .mp4 file.
    fps : int
        Frames per second.
    """
    T, C, ny, nx = q_true.shape

    # Grid for streamplot
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    # Pre-compute magnitude limits for consistent coloring
    mag_true = np.sqrt(q_true[:, 0] ** 2 + q_true[:, 1] ** 2)
    mag_pred = np.sqrt(q_pred[:, 0] ** 2 + q_pred[:, 1] ** 2)
    vmax = max(mag_true.max(), mag_pred.max())
    if vmax < 1e-9:
        vmax = 1.0

    title = fig.suptitle("Time: 0")

    def update(frame):
        title.set_text(f"Time: {frame}")
        for ax in axes:
            ax.clear()

        # True
        ax_true = axes[0]
        ax_true.set_title("True Streamlines")
        ax_true.set_xlim(0, nx - 1)
        ax_true.set_ylim(0, ny - 1)

        u_t = q_true[frame, 0]
        v_t = q_true[frame, 1]
        speed_t = np.sqrt(u_t**2 + v_t**2)

        # Vary linewidth with speed to indicate magnitude
        lw_t = 2.5 * speed_t / vmax

        ax_true.streamplot(
            X,
            Y,
            u_t,
            v_t,
            color=speed_t,
            linewidth=lw_t,
            cmap="viridis",
            norm=plt.Normalize(0, vmax),
            density=1.5,
            arrowsize=1.0,
        )

        # Pred
        ax_pred = axes[1]
        ax_pred.set_title("Pred Streamlines")
        ax_pred.set_xlim(0, nx - 1)
        ax_pred.set_ylim(0, ny - 1)

        u_p = q_pred[frame, 0]
        v_p = q_pred[frame, 1]
        speed_p = np.sqrt(u_p**2 + v_p**2)

        lw_p = 2.5 * speed_p / vmax

        ax_pred.streamplot(
            X,
            Y,
            u_p,
            v_p,
            color=speed_p,
            linewidth=lw_p,
            cmap="viridis",
            norm=plt.Normalize(0, vmax),
            density=1.5,
            arrowsize=1.0,
        )

        return axes

    ani = animation.FuncAnimation(fig, update, frames=T, blit=False, interval=1000 / fps)
    writer = "ffmpeg" if output_path.suffix == ".mp4" else "pillow"
    ani.save(output_path, writer=writer, fps=fps)
    plt.close(fig)
    print(f"Saved streamline animation to {output_path}")


def animate_moving_dashes(q_true, q_pred, output_path, fps=DEFAULT_FPS):
    """
    Animate streamlines with moving dashes to indicate flow direction.

    Uses a shifting dash pattern on streamlines to simulate flow motion
    without using arrowheads.

    Parameters
    ----------
    q_true : np.ndarray
        Ground truth trajectory of shape (T, C, ny, nx).
    q_pred : np.ndarray
        Predicted trajectory of shape (T, C, ny, nx).
    output_path : Path
        Destination path for the .mp4 file.
    fps : int
        Frames per second.
    """
    T, C, ny, nx = q_true.shape

    # Grid for streamplot
    x = np.arange(nx)
    y = np.arange(ny)
    X, Y = np.meshgrid(x, y)

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    # Pre-compute magnitude limits
    mag_true = np.sqrt(q_true[:, 0] ** 2 + q_true[:, 1] ** 2)
    mag_pred = np.sqrt(q_pred[:, 0] ** 2 + q_pred[:, 1] ** 2)
    vmax = max(mag_true.max(), mag_pred.max())
    if vmax < 1e-9:
        vmax = 1.0

    title = fig.suptitle("Time: 0")

    def update(frame):
        title.set_text(f"Time: {frame}")
        for ax in axes:
            ax.clear()

        # Dash parameters
        dash_len = 10.0
        gap_len = 10.0
        total_len = dash_len + gap_len
        # Shift phase to animate (negative moves "downstream")
        phase = -(frame % 10) * (total_len / 10.0)

        # True
        ax_true = axes[0]
        ax_true.set_title("True Moving Dashes")
        ax_true.set_xlim(0, nx - 1)
        ax_true.set_ylim(0, ny - 1)

        u_t = q_true[frame, 0]
        v_t = q_true[frame, 1]
        speed_t = np.sqrt(u_t**2 + v_t**2)
        lw_t = 2.5 * speed_t / vmax

        st_true = ax_true.streamplot(
            X,
            Y,
            u_t,
            v_t,
            color=speed_t,
            linewidth=lw_t,
            cmap="viridis",
            norm=plt.Normalize(0, vmax),
            density=1.5,
            arrowsize=0.0,
        )
        st_true.lines.set_dashes([(phase, (dash_len, gap_len))])
        st_true.arrows.set_visible(False)

        # Pred
        ax_pred = axes[1]
        ax_pred.set_title("Pred Moving Dashes")
        ax_pred.set_xlim(0, nx - 1)
        ax_pred.set_ylim(0, ny - 1)

        u_p = q_pred[frame, 0]
        v_p = q_pred[frame, 1]
        speed_p = np.sqrt(u_p**2 + v_p**2)
        lw_p = 2.5 * speed_p / vmax

        st_pred = ax_pred.streamplot(
            X,
            Y,
            u_p,
            v_p,
            color=speed_p,
            linewidth=lw_p,
            cmap="viridis",
            norm=plt.Normalize(0, vmax),
            density=1.5,
            arrowsize=0.0,
        )
        st_pred.lines.set_dashes([(phase, (dash_len, gap_len))])
        st_pred.arrows.set_visible(False)

        return axes

    ani = animation.FuncAnimation(fig, update, frames=T, blit=False, interval=1000 / fps)
    writer = "ffmpeg" if output_path.suffix == ".mp4" else "pillow"
    ani.save(output_path, writer=writer, fps=fps)
    plt.close(fig)
    print(f"Saved moving dashes animation to {output_path}")


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
        default="all",
        choices=["channels", "quiver", "vorticity", "streamlines", "dashes", "all"],
        help="Animation style to generate.",
    )
    parser.add_argument("--fps", type=int, default=DEFAULT_FPS, help="Frames per second.")
    parser.add_argument(
        "--format", type=str, default="mp4", choices=["mp4", "gif"], help="Output format."
    )
    args = parser.parse_args()

    print(f"Loading data for run: {args.run_id}")
    try:
        q_true_flat, q_pred_flat = load_rollout_data(args.run_id, equation=args.equation)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return

    # Reshape
    q_true, layout = reshape_trajectory(q_true_flat, n_components=2)
    q_pred, _ = reshape_trajectory(q_pred_flat, n_components=2)

    output_dir = OUTPUTS_DIR / args.equation / args.run_id / "movies"
    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    if args.style in ["channels", "all"]:
        print("Generating channel animation...")
        animate_channels(
            q_true,
            q_pred,
            output_dir / f"{timestamp}_{args.run_id}_truevpred_channels.{args.format}",
            fps=args.fps,
        )

    if args.style in ["quiver", "all"]:
        print("Generating quiver animation...")
        animate_quiver(
            q_true,
            q_pred,
            output_dir / f"{timestamp}_{args.run_id}_truevpred_quiver.{args.format}",
            fps=args.fps,
        )

    if args.style in ["vorticity", "all"]:
        print("Generating vorticity animation...")
        animate_vorticity(
            q_true,
            q_pred,
            output_dir / f"{timestamp}_{args.run_id}_truevpred_vorticity.{args.format}",
            fps=args.fps,
        )

    if args.style in ["streamlines", "all"]:
        print("Generating streamline animation...")
        animate_streamlines(
            q_true,
            q_pred,
            output_dir / f"{timestamp}_{args.run_id}_truevpred_streamlines.{args.format}",
            fps=args.fps,
        )

    if args.style in ["dashes", "all"]:
        print("Generating moving dashes animation...")
        animate_moving_dashes(
            q_true,
            q_pred,
            output_dir / f"{timestamp}_{args.run_id}_truevpred_dashes.{args.format}",
            fps=args.fps,
        )


if __name__ == "__main__":
    main()
