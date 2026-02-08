import numpy as np
import pysindy as ps
from scipy.integrate import solve_ivp

# Suppress warnings
# warnings.filterwarnings("ignore", category=UserWarning, module="pysindy")
# warnings.filterwarnings("ignore", category=RuntimeWarning)


def train_sindy_pde(u, t, config=None):
    print("      [SINDy] Starting GRID training pipeline (Native PDEFIND Mode)...")

    if config is None:
        config = {}

    # Extract config with defaults
    threshold = config.get("threshold", 0.1)
    alpha = config.get("alpha", 1e-05)
    ensemble = config.get("ensemble", False)
    n_models = config.get("n_models", 10)
    optimizer_name = config.get("optimizer", "STLSQ")
    poly_degree = config.get("poly_degree", 2)
    derivative_order = config.get("derivative_order", 2)
    diff_order = config.get("diff_order", 4)
    nu = config.get("nu", 1.0)
    kappa = config.get("kappa", None)
    max_iter = config.get("max_iter", 2000)

    # u is (X, Y, Time, Channels)
    n_channels = u.shape[-1]

    # 2. Define Grids
    nx, ny = u.shape[0], u.shape[1]
    # Use endpoint=False for periodic boundaries to match data spacing (dx = 1/N)
    grid_x = np.linspace(0, 1, nx, endpoint=False)
    grid_y = np.linspace(0, 1, ny, endpoint=False)
    X, Y = np.meshgrid(grid_x, grid_y, indexing="ij")
    spatial_grid = np.stack([X, Y], axis=-1)

    # 3. Define PDE Library
    # degree=2 generates [u, u^2]. Combined with derivatives, this allows terms like u*ux or u^2*ux.
    # Note: degree=1 is sufficient for u*ux, but degree=2 is needed for source terms like u^2.
    poly_library = ps.PolynomialLibrary(degree=poly_degree, include_bias=False)

    # Configure PDELibrary
    library = ps.PDELibrary(
        function_library=poly_library,
        derivative_order=derivative_order,
        spatial_grid=spatial_grid,
        include_bias=True,
        differentiation_method=ps.FiniteDifference,
        # order controls the finite difference stencil size (accuracy),
        # not the derivative order.
        diff_kwargs={"periodic": True, "order": diff_order, "is_uniform": True},
    )

    # Optimizer Configuration
    # - threshold: Minimum coefficient magnitude for sparsity (STLSQ, SR3).
    # - alpha: L2 regularization strength (STLSQ) or p-value threshold (SSR).
    optimizer_params = {"optimizer": optimizer_name}

    if optimizer_name == "STLSQ":
        # Standard SINDy optimizer: Sequential Thresholded Least Squares
        optimizer = ps.STLSQ(
            threshold=threshold, alpha=alpha, max_iter=max_iter, normalize_columns=True
        )
        optimizer_params["threshold"] = threshold
        optimizer_params["alpha"] = alpha
    elif optimizer_name == "SR3":
        # Sparse Relaxed Regularized Regression
        # - relax_coeff_nu: level of relaxation
        # - regularizer: regularization function. "L0" is the L0 norm
        # Note: Setting threshold via attribute to avoid __init__ issues in some versions
        optimizer = ps.SR3(
            relax_coeff_nu=nu, regularizer="L0", max_iter=max_iter, normalize_columns=True
        )
        optimizer.threshold = threshold
        optimizer_params["threshold"] = threshold
        optimizer_params["nu"] = nu
    elif optimizer_name == "SSR":
        # Stepwise Sparse Regression
        # - alpha: Used as the p-value significance level for term selection (default criteria='p_value').
        optimizer = ps.SSR(alpha=alpha, max_iter=max_iter, normalize_columns=True)
        optimizer_params["alpha"] = alpha
    elif optimizer_name == "FROLS":
        # Forward Regression Orthogonal Least Squares
        # - kappa: The exact number of non-zero terms to keep in the final model.
        optimizer = ps.FROLS(kappa=kappa, normalize_columns=True)
        optimizer_params["kappa"] = kappa
    else:
        print(f"      [Warning] Unknown optimizer '{optimizer_name}'. Defaulting to STLSQ.")
        optimizer = ps.STLSQ(threshold=threshold, alpha=alpha, normalize_columns=True)
        optimizer_params["optimizer"] = "STLSQ"
        optimizer_params["threshold"] = threshold
        optimizer_params["alpha"] = alpha

    feature_names = [f"u{i}" for i in range(n_channels)] if n_channels > 1 else ["u"]

    # 5. Fit Model
    print("      [SINDy] Fitting model...")
    model = ps.SINDy(feature_library=library, optimizer=optimizer)
    if ensemble:
        model.fit(
            u,
            t=t,
            feature_names=feature_names,
            ensemble=True,
            n_models=n_models,
            n_subset=int(len(t) * 0.6),
        )
    else:
        model.fit(u, t=t, feature_names=feature_names)

    # 6. Print Results
    print("\n" + "=" * 40)
    print("      IDENTIFIED DYNAMICS (Native PDEFIND)")
    print("=" * 40)
    model.print()

    print("\n--- Identified Parameters ---")
    if len(model.coefficients()) > 0:
        coeffs = model.coefficients()[0]
        features = model.get_feature_names()
        for name, val in zip(features, coeffs):
            if abs(val) > 0.001:
                print(f"  {name:<10} : {val:.5f}")

    # 7. Simulation Wrapper
    def forecast_evolution(u0, t_span):
        # u0 is expected to be (X, Y, Channels)
        u0_flat = u0.flatten()

        nx, ny = u0.shape[0], u0.shape[1]

        def rhs(t, u_flat):
            # Reshape to (X, Y, C)
            u_reshaped = u_flat.reshape(nx, ny, n_channels)

            # Add time axis for predict: (X, Y, 1, C)
            u_in = u_reshaped[..., np.newaxis, :]

            # Predict derivatives
            u_dot = model.predict(u_in)

            # Flatten
            return u_dot.flatten()

        print(f"      [SINDy] Integrating from t={t_span[0]} to {t_span[-1]}...")
        sol = solve_ivp(
            rhs,
            (t_span[0], t_span[-1]),
            u0_flat,
            t_eval=t_span,
            method="RK45",
            rtol=1e-5,
            atol=1e-8,
        )

        # Reshape: (Time, X*Y*C) -> (Time, X, Y, C)
        u_sim = sol.y.T.reshape(-1, nx, ny, n_channels)

        return u_sim

    return forecast_evolution, model, optimizer_params
