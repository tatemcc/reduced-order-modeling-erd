import os

# Force single-threaded execution to prevent deadlocks
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import numpy as np
import h5py
import glob

# Import your local modules
from sindy_dynamics import train_sindy_pde
from save_evolution import animate_video
from plot_metrics import plot_mse_evolution


def load_dynabench_data(
    equation="advection", structure="grid", resolution="low", split="train", base_path="data"
):
    print(f"      [Loader] Searching for {split} data ({structure}) in {base_path}...")
    search_path = os.path.join(base_path, equation, structure, resolution, f"*{split}*.h5")
    files = glob.glob(search_path)

    if not files:
        raise FileNotFoundError(f"No .h5 files found at: {search_path}")

    file_path = files[0]
    print(f"      [Loader] Loading file: {file_path}")

    with h5py.File(file_path, "r") as f:
        u_data = f["data"][0]

        if "points" in f:
            points = f["points"][0]
        else:
            # Generate implicit grid for 15x15
            nx, ny = u_data.shape[1], u_data.shape[2]
            x = np.linspace(0, 1, nx)
            y = np.linspace(0, 1, ny)
            X, Y = np.meshgrid(x, y, indexing="ij")
            points = np.stack([X, Y], axis=-1)

    return u_data, points


print("--- Starting SINDy Validation Script (Grid Mode) ---")

# 1. Load Training Data
print("[1/5] Loading Training Data...")
u_train, points = load_dynabench_data(split="train")

if points.ndim == 3:
    print(
        f"      [Pre-process] Flattening points from {points.shape} to ({points.shape[0]*points.shape[1]}, 2)"
    )
    points = points.reshape(-1, 2)

t_train = np.arange(u_train.shape[0]).astype(float)

# 2. Identify System Dynamics
print(f"[2/5] Training SINDy model on {len(t_train)} time steps...")
forecast_evolution, model = train_sindy_pde(u_train, t_train, points)

# 3. Load Validation Data
print("[3/5] Loading Validation Data...")
u_val_true, _ = load_dynabench_data(split="val")

# --- MODIFICATION: Limit to N steps ---
N_VAL_STEPS = 12  # <--- Set your desired number of steps here
if N_VAL_STEPS < u_val_true.shape[0]:
    print(f"      [Validation] Limiting validation to first {N_VAL_STEPS} steps.")
    u_val_true = u_val_true[:N_VAL_STEPS]

t_val = np.arange(u_val_true.shape[0]).astype(float)

# 4. Forecast
print("[4/5] Forecasting evolution on validation set...")
u_val_pred = forecast_evolution(u_val_true[0], t_val)

mse = np.mean((u_val_true.squeeze() - u_val_pred.squeeze()) ** 2)
print(f"\nValidation Forecast MSE: {mse:.6e}")

# 5. Visualize
print("[5/5] Generating visualization...")
animate_video(u_val_true, u_val_pred, points, t_val)

# 6. Plot Error Metrics
print("[6/6] Plotting error evolution...")
plot_mse_evolution(u_val_true, u_val_pred, t_val)
