import argparse
import numpy as np

# from dynabench.dataset import DynabenchSimulationIterator

# local modules
from sindy_dynamics import train_sindy_pde
from visualize_evolution import animate_comparison
from dynabench_loader import load_dynabench_data


def main():
    parser = argparse.ArgumentParser(description="Train and Validate SINDy Model")
    parser.add_argument(
        "--optimizer",
        type=str,
        default="STLSQ",
        choices=["STLSQ", "SR3", "SSR", "FROLS"],
        help="Optimizer type",
    )
    parser.add_argument("--threshold", type=float, default=0.1, help="Sparsity threshold")
    parser.add_argument("--alpha", type=float, default=1e-5, help="L2 regularization alpha")
    parser.add_argument("--poly_degree", type=int, default=2, help="Polynomial degree for library")
    parser.add_argument("--nu", type=float, default=1.0, help="Relaxation parameter for SR3")
    parser.add_argument("--kappa", type=int, default=None, help="Number of terms for FROLS")
    parser.add_argument(
        "--max_iter", type=int, default=2000, help="Maximum iterations for optimizer"
    )
    parser.add_argument(
        "--derivative_order", type=int, default=2, help="Derivative order for library"
    )
    parser.add_argument("--diff_order", type=int, default=4, help="Finite difference order")
    parser.add_argument("--ensemble", action="store_true", help="Enable ensemble SINDy")
    parser.add_argument("--n_models", type=int, default=10, help="Number of models for ensemble")
    parser.add_argument("--n_val_steps", type=int, default=12, help="Number of validation steps")
    parser.add_argument(
        "--save_animation", action="store_true", help="Save forecast animation as a .gif"
    )
    parser.add_argument("--index", type=int, default=0, help="Simulation index to load")

    args = parser.parse_args()

    # Create config dict
    config = vars(args)

    # 1. Load Training Data
    print("[1/5] Loading Training Data...")

    u_train, _ = load_dynabench_data(split="train", sim_index=args.index)

    t_train = np.arange(u_train.shape[0]).astype(float)

    # Pre-process: Transpose (Time, Channels, X, Y) -> (X, Y, Time, Channels)
    print(
        f"      [Pre-process] Reshaping training data from {u_train.shape} to (X, Y, Time, Channels)..."
    )
    u_train_sindy = np.transpose(u_train, (2, 3, 0, 1))

    # 2. Identify System Dynamics
    print(f"[2/5] Training SINDy model on {len(t_train)} time steps...")
    forecast_evolution, model, optimizer_params = train_sindy_pde(
        u_train_sindy, t_train, config=config
    )

    # 3. Load Validation Data
    print("[3/5] Loading Validation Data...")
    u_val_true, _ = load_dynabench_data(split="val")

    # Pre-process Validation Data: (Time, Channels, X, Y) -> (Time, X, Y, Channels)
    u_val_true = np.transpose(u_val_true, (0, 2, 3, 1))

    # limit to N timesteps
    if args.n_val_steps < u_val_true.shape[0]:
        print(f"      [Validation] Limiting validation to first {args.n_val_steps} steps.")
        u_val_true = u_val_true[: args.n_val_steps]

    t_val = np.arange(u_val_true.shape[0]).astype(float)

    # 4. Forecast
    print("[4/5] Forecasting evolution on validation set...")
    u_val_pred = forecast_evolution(u_val_true[0], t_val)

    mse = np.mean((u_val_true.squeeze() - u_val_pred.squeeze()) ** 2)
    print(f"\nValidation Forecast MSE: {mse:.6e}")

    # 5. Visualize
    print("[5/5] Generating visualization...")

    optimizer_params["Sim Index"] = args.index

    # Extract coefficients for visualization
    model_coeffs = {}
    if hasattr(model, "coefficients"):
        coeffs = model.coefficients()[0]
        features = model.get_feature_names()
        for name, val in zip(features, coeffs):
            if abs(val) > 0.001:
                model_coeffs[name] = val

    animate_comparison(
        u_val_true,
        u_val_pred,
        t_val,
        optimizer_params=optimizer_params,
        model_coeffs=model_coeffs,
        save_anim=args.save_animation,
    )

    # 6. Plot Error Metrics
    # print("[6/6] Plotting error evolution...")
    # plot_mse_evolution(u_val_true, u_val_pred, t_val)


if __name__ == "__main__":
    main()
