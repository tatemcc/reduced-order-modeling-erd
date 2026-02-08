import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


def animate_comparison(
    u_true,
    u_pred,
    t,
    optimizer_params=None,
    model_coeffs=None,
    filename="output/sindy_evolution.gif",
    save_anim=True,
):
    print(f"      [Vis] Generating animation with {len(t)} frames...")

    # Handle channels: (Time, X, Y, C) -> (Time, X, Y)
    if u_true.ndim == 4:
        u_true = u_true[..., 0]
    if u_pred.ndim == 4:
        u_pred = u_pred[..., 0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Adjust layout to make room for titles
    fig.subplots_adjust(top=0.80)

    vmin = min(u_true.min(), u_pred.min())
    vmax = max(u_true.max(), u_pred.max())
    extent = [0, 1, 0, 1]

    # Transpose (.T) to map (X, Y) data to (Row, Col) for imshow
    im1 = ax1.imshow(
        u_true[0].T, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax, extent=extent
    )
    ax1.set_title("Ground Truth")
    fig.colorbar(im1, ax=ax1)

    im2 = ax2.imshow(
        u_pred[0].T, origin="lower", cmap="viridis", vmin=vmin, vmax=vmax, extent=extent
    )
    ax2.set_title("SINDy Forecast")
    fig.colorbar(im2, ax=ax2)

    # Common text settings
    font_size = 11
    line_spacing = 0.04
    start_y = 0.96

    if optimizer_params:
        parts = [f"{k}: {v}" for k, v in optimizer_params.items()]
        # Use fig.text instead of suptitle so it aligns perfectly with the lines below
        fig.text(0.5, start_y, " | ".join(parts), ha="center", fontsize=font_size)

    if model_coeffs:
        coeff_str = ", ".join([f'"{k}": {v:.3f}' for k, v in model_coeffs.items()])
        fig.text(0.5, start_y - line_spacing, f"[{coeff_str}]", ha="center", fontsize=font_size)

    time_txt = fig.text(
        0.5, start_y - (line_spacing * 2), f"Time: {t[0]:.2f}s", ha="center", fontsize=font_size
    )

    def update(frame):
        im1.set_data(u_true[frame].T)
        im2.set_data(u_pred[frame].T)
        time_txt.set_text(f"Time: {t[frame]:.2f}s")
        return im1, im2, time_txt

    ani = FuncAnimation(fig, update, frames=len(t), blit=False)

    if save_anim:
        print(f"      [Vis] Saving to {filename}...")
        writer = PillowWriter(fps=15)
        ani.save(filename, writer=writer)

    plt.show()
