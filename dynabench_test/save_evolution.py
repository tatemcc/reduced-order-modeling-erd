import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np


def animate_video(u_true, u_pred, points, t, filename="sindy_forecast.gif"):
    print(f"      [Vis] Generating animation with {len(t)} frames...")

    # Setup Figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Flatten spatial dimensions to (Time, N_points)
    # This ensures u_t[frame] is Rank 1, matching the flattened points array
    u_t = u_true.reshape(u_true.shape[0], -1)
    u_p = u_pred.reshape(u_pred.shape[0], -1)

    # Determine global min/max for consistent colorbar (using nanmin/max for stability)
    vmin = np.nanmin([np.nanmin(u_t), np.nanmin(u_p)])
    vmax = np.nanmax([np.nanmax(u_t), np.nanmax(u_p)])

    # Initial Scatter Plots
    scat1 = ax1.scatter(points[:, 0], points[:, 1], c=u_t[0], cmap="viridis", vmin=vmin, vmax=vmax)
    ax1.set_title("Ground Truth")
    ax1.set_aspect("equal")
    fig.colorbar(scat1, ax=ax1)

    scat2 = ax2.scatter(points[:, 0], points[:, 1], c=u_p[0], cmap="viridis", vmin=vmin, vmax=vmax)
    ax2.set_title("SINDy Forecast")
    ax2.set_aspect("equal")
    fig.colorbar(scat2, ax=ax2)

    def update(frame):
        # Update colors for the current time step
        # set_array requires a 1D array
        scat1.set_array(u_t[frame])
        scat2.set_array(u_p[frame])

        fig.suptitle(f"Time: {t[frame]:.2f}s (Frame {frame}/{len(t)})")
        return scat1, scat2

    # Create Animation Object
    ani = FuncAnimation(fig, update, frames=len(t), blit=True)

    # Save to GIF using PillowWriter
    print(f"      [Vis] Saving to {filename}...")
    writer = PillowWriter(fps=15)

    try:
        ani.save(filename, writer=writer)
        print(f"      [Vis] Saved successfully: {filename}")
    except Exception as e:
        print(f"      [Vis] Error saving GIF: {e}")
        print("      Ensure 'pillow' is installed (pip install pillow).")

    plt.close(fig)
