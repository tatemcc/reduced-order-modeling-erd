import numpy as np

# from dynabench.dataset import DynabenchSimulationIterator

# local modules
from sindy_dynamics import train_sindy_pde
from save_evolution import animate_video
from plot_metrics import plot_mse_evolution
from dynabench_loader import load_dynabench_data

# limit number of timesteps
N_VAL_STEPS = 12

# 1. Load Training Data
print("[1/5] Loading Training Data...")

u_train, points = load_dynabench_data(split="train", sim_index=1)

t_train = np.arange(u_train.shape[0]).astype(float)

# Pre-process: Transpose (Time, Channels, X, Y) -> (X, Y, Time, Channels)
# This prepares the data for PySINDy's grid requirements
print(
    f"      [Pre-process] Reshaping training data from {u_train.shape} to (X, Y, Time, Channels)..."
)
u_train_sindy = np.transpose(u_train, (2, 3, 0, 1))

# 2. Identify System Dynamics
print(f"[2/5] Training SINDy model on {len(t_train)} time steps...")
forecast_evolution, model = train_sindy_pde(u_train_sindy, t_train, points)

# 3. Load Validation Data
print("[3/5] Loading Validation Data...")
# val_iterator = DynabenchSimulationIterator(
#     equation='advection',
#     structure='cloud',
#     resolution='low',
#     split='val'
# )
# val_batch = next(iter(val_iterator))

# u_val_true = val_batch.y[0]


u_val_true, _ = load_dynabench_data(split="val")

# limit to N timesteps
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
# Flatten points for scatter plot visualization: (X, Y, 2) -> (N_points, 2)
points_flat = points.reshape(-1, 2)

animate_video(u_val_true, u_val_pred, points_flat, t_val)

# 6. Plot Error Metrics
print("[6/6] Plotting error evolution...")
plot_mse_evolution(u_val_true, u_val_pred, t_val)
