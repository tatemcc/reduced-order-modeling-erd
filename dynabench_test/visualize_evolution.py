import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def animate_comparison(u_true, u_pred, points, t):
    plt.ion()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    u_t, u_p = u_true.squeeze(), u_pred.squeeze()
    vmin, vmax = u_t.min(), u_t.max()

    scat1 = ax1.scatter(points[:, 0], points[:, 1], c=u_t[0], cmap="viridis", vmin=vmin, vmax=vmax)
    ax1.set_title("Validation (Ground Truth)")
    plt.colorbar(scat1, ax=ax1)

    scat2 = ax2.scatter(points[:, 0], points[:, 1], c=u_p[0], cmap="viridis", vmin=vmin, vmax=vmax)
    ax2.set_title("SINDy Forecast")
    plt.colorbar(scat2, ax=ax2)

    def update(frame):
        scat1.set_array(u_t[frame])
        scat2.set_array(u_p[frame])
        fig.suptitle(f"Time Step: {frame}")
        return scat1, scat2

    _ani = FuncAnimation(fig, update, frames=len(t), interval=100, blit=True)
    plt.show(block=True)
