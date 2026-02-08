import time
import itertools
import csv
import numpy as np
import os

# Local modules
from sindy_dynamics import train_sindy_pde
from dynabench_loader import load_dynabench_data

# Configuration Grid
PARAM_GRID = {
    "optimizer": ["STLSQ", "SR3"],
    "threshold": [0.01, 0.05, 0.1, 0.2],
    "alpha": [1e-5, 1e-3, 0.1],
    "poly_degree": [2],
    "derivative_order": [2],
    "ensemble": [False],  # Set to True if you have time/compute
}

OUTPUT_FILE = "output/benchmark_results.csv"
N_VAL_STEPS = 20  # Number of steps to forecast for validation


def run_benchmark():
    # 1. Load Data (Once)
    print("[Benchmark] Loading Data...")
    u_train, _ = load_dynabench_data(split="train", sim_index=1)
    u_val_true, _ = load_dynabench_data(split="val", sim_index=1)

    # Pre-process Training Data
    t_train = np.arange(u_train.shape[0]).astype(float)
    u_train_sindy = np.transpose(u_train, (2, 3, 0, 1))  # (X, Y, Time, Channels)

    # Pre-process Validation Data
    if N_VAL_STEPS < u_val_true.shape[0]:
        u_val_true = u_val_true[:N_VAL_STEPS]

    # Transpose Validation Data to (Time, X, Y, Channels)
    u_val_true = np.transpose(u_val_true, (0, 2, 3, 1))

    t_val = np.arange(u_val_true.shape[0]).astype(float)

    # 2. Generate Parameter Combinations
    keys, values = zip(*PARAM_GRID.items())
    combinations = [dict(zip(keys, v)) for v in itertools.product(*values)]

    print(f"[Benchmark] Starting sweep over {len(combinations)} configurations...")

    results = []

    # Ensure output directory exists
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

    # 3. Run Sweep
    for i, config in enumerate(combinations):
        print(f"\n--- Run {i+1}/{len(combinations)} ---")
        print(f"Config: {config}")

        start_time = time.time()
        run_data = config.copy()

        try:
            # Train
            forecast_evolution, model, _ = train_sindy_pde(u_train_sindy, t_train, config=config)
            train_duration = time.time() - start_time

            # Validate (Forecast)
            # u_val_true[0] is the initial condition for validation
            u_val_pred = forecast_evolution(u_val_true[0], t_val)

            # Calculate Metrics
            # MSE
            diff = u_val_true.squeeze() - u_val_pred.squeeze()
            mse = np.mean(diff**2)

            # Complexity (L0 Norm)
            coeffs = model.coefficients()
            l0_norm = np.count_nonzero(np.abs(coeffs) > 1e-6)

            # Check for explosion (NaNs)
            is_stable = not np.isnan(mse) and not np.isinf(mse)

            print(f"   -> MSE: {mse:.4e} | Terms: {l0_norm} | Time: {train_duration:.2f}s")

            run_data.update(
                {
                    "mse": mse if is_stable else 1e10,
                    "l0_norm": l0_norm,
                    "train_time": train_duration,
                    "status": "Success" if is_stable else "Unstable",
                }
            )

        except Exception as e:
            print(f"   -> Failed: {e}")
            run_data.update(
                {
                    "mse": None,
                    "l0_norm": None,
                    "train_time": time.time() - start_time,
                    "status": f"Error: {str(e)}",
                }
            )

        results.append(run_data)

    # 4. Save Results
    print(f"\n[Benchmark] Saving results to {OUTPUT_FILE}...")
    if results:
        keys = results[0].keys()
        with open(OUTPUT_FILE, "w", newline="") as f:
            dict_writer = csv.DictWriter(f, fieldnames=keys)
            dict_writer.writeheader()
            dict_writer.writerows(results)
        print("[Benchmark] Done.")
    else:
        print("[Benchmark] No results generated.")


if __name__ == "__main__":
    run_benchmark()
