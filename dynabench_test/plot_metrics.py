import matplotlib.pyplot as plt
import numpy as np


def plot_mse_evolution(u_true, u_pred, t, filename="mse_vs_time.png"):
    print(f"      [Metrics] Calculating MSE over {len(t)} time steps...")

    # Ensure shapes match
    if u_true.shape != u_pred.shape:
        print(
            f"      [Warning] Shape mismatch: {u_true.shape} vs {u_pred.shape}. Truncating to shorter length."
        )
        min_len = min(u_true.shape[0], u_pred.shape[0])
        u_true = u_true[:min_len]
        u_pred = u_pred[:min_len]
        t = t[:min_len]

    # Calculate MSE at each time step (averaging over spatial dimensions)
    # Reshape to (Time, -1) to flatten all spatial dimensions
    diff = u_true.reshape(u_true.shape[0], -1) - u_pred.reshape(u_pred.shape[0], -1)
    mse_over_time = np.mean(diff**2, axis=1)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t, mse_over_time, label="MSE", linewidth=2, color="crimson")

    ax.set_yscale("log")  # <--- Log scale added

    ax.set_title("Forecast Error Evolution (MSE vs Time)")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Mean Squared Error (Log Scale)")
    ax.grid(True, linestyle="--", alpha=0.6, which="both")  # Grid for log scale
    ax.legend()

    # Save
    print(f"      [Metrics] Saving plot to {filename}...")
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)
