"""
Plotting and animation utilities for the Burgers ROM pipeline.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

DEFAULT_FPS = 10


def animate_channels(
    q_true: np.ndarray,
    q_pred: np.ndarray,
    output_path: Path,
    fps: int = DEFAULT_FPS,
) -> None:
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


def animate_amplitudes(
    A_true: np.ndarray,
    A_pred: np.ndarray,
    output_path: Path,
    fps: int = DEFAULT_FPS,
) -> None:
    """
    Animate POD mode amplitudes over time using a grouped bar chart.

    Creates a single animation showing the temporal evolution of
    the true vs. predicted POD coefficients side-by-side for each mode.

    Parameters
    ----------
    A_true : np.ndarray
        Ground truth coefficients of shape (T, r).
    A_pred : np.ndarray
        Predicted coefficients of shape (T, r).
    output_path : Path
        Destination path for the .mp4 file.
    fps : int
        Frames per second for the animation.
    """
    T, r = A_true.shape
    if A_pred.shape != (T, r):
        raise ValueError("A_true and A_pred must have the same shape.")

    mode_indices = np.arange(r)
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    # Determine shared y-axis limits
    ymin = min(A_true.min(), A_pred.min())
    ymax = max(A_true.max(), A_pred.max())
    padding = (ymax - ymin) * 0.1
    ymin -= padding
    ymax += padding

    # Initial bar plots
    bar_true = ax.bar(
        mode_indices - width / 2, A_true[0], width, label="True", color="tab:blue"
    )
    bar_pred = ax.bar(
        mode_indices + width / 2, A_pred[0], width, label="Predicted", color="tab:orange"
    )

    # Add some text for labels, title and axes ticks
    ax.set_ylabel("Amplitude")
    ax.set_title("POD Mode Amplitudes")
    ax.set_xticks(mode_indices)
    ax.set_xlabel("Mode Index")
    ax.legend()
    ax.set_ylim(ymin, ymax)

    title = fig.suptitle("Time: 0")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    def update(frame):
        title.set_text(f"Time: {frame}")
        for i, b in enumerate(bar_true):
            b.set_height(A_true[frame, i])
        for i, b in enumerate(bar_pred):
            b.set_height(A_pred[frame, i])
        return list(bar_true) + list(bar_pred) + [title]

    ani = animation.FuncAnimation(fig, update, frames=T, blit=True, interval=1000 / fps)
    writer = "ffmpeg" if output_path.suffix == ".mp4" else "pillow"
    ani.save(output_path, writer=writer, fps=fps)
    plt.close(fig)
    print(f"Saved amplitude animation to {output_path}")


def plot_pod_basis(
    mode: np.ndarray,
    output_path: Path,
    title: str = "POD Mode",
) -> None:
    """
    Plot the spatial structure of a single POD basis mode (static).

    Parameters
    ----------
    mode : np.ndarray
        Field of shape (C, ny, nx).
    output_path : Path
        Destination path for the image.
    title : str
        Title prefix for the plot.
    """
    C, ny, nx = mode.shape
    if C != 2:
        print(f"Skipping basis plot: expected 2 components, got {C}")
        return

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), squeeze=True)

    # Center colormap at 0
    vmax = np.max(np.abs(mode))
    vmin = -vmax

    # u
    ax0 = axes[0]
    ax0.set_title(f"{title} - u")
    im0 = ax0.imshow(mode[0], origin="lower", cmap="RdBu_r", vmin=vmin, vmax=vmax)
    fig.colorbar(im0, ax=ax0)

    # v
    ax1 = axes[1]
    ax1.set_title(f"{title} - v")
    im1 = ax1.imshow(mode[1], origin="lower", cmap="RdBu_r", vmin=vmin, vmax=vmax)
    fig.colorbar(im1, ax=ax1)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close(fig)
