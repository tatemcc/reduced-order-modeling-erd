import numpy as np
import pysindy as ps
from scipy.integrate import solve_ivp

# Suppress warnings
# warnings.filterwarnings("ignore", category=UserWarning, module="pysindy")
# warnings.filterwarnings("ignore", category=RuntimeWarning)


def train_sindy_pde(u, t, spatial_grid, threshold=0.1, alpha=1e-05):
    print("      [SINDy] Starting GRID training pipeline (Native PDEFIND Mode)...")

    # u is (X, Y, Time, Channels)
    n_channels = u.shape[-1]

    # 2. Define Grids
    # Necessary to avoid singular matrix error in ps.FiniteDifference
    nx, ny = u.shape[0], u.shape[1]
    grid_x = np.linspace(0, 1, nx)
    grid_y = np.linspace(0, 1, ny)
    X, Y = np.meshgrid(grid_x, grid_y, indexing="ij")
    spatial_grid = np.stack([X, Y], axis=-1)

    # 3. Define PDE Library
    # degree=2 allows for quadratic terms (e.g. Burgers u*ux)
    poly_library = ps.PolynomialLibrary(degree=2, include_bias=False)

    # Configure PDELibrary
    # Pass differentiation_method as class and kwargs to allow PDELibrary to control 'd' (order)
    library = ps.PDELibrary(
        function_library=poly_library,
        derivative_order=2,
        spatial_grid=spatial_grid,
        include_bias=True,
        differentiation_method=ps.FiniteDifference,
        diff_kwargs={"periodic": True, "is_uniform": True},
    )

    # STLSQ optimizer (Page 8, 18)
    optimizer = ps.STLSQ(threshold=threshold, alpha=alpha, normalize_columns=True)

    feature_names = [f"u{i}" for i in range(n_channels)] if n_channels > 1 else ["u"]

    # 5. Fit Model
    print("      [SINDy] Fitting model...")
    model = ps.SINDy(feature_library=library, optimizer=optimizer)
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
        # u0 shape: (Channels, X, Y) -> Transpose to (X, Y, Channels)
        u0_native = np.transpose(u0, (1, 2, 0))
        u0_flat = u0_native.flatten()

        nx, ny = u0_native.shape[0], u0_native.shape[1]

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
        sol = solve_ivp(rhs, (t_span[0], t_span[-1]), u0_flat, t_eval=t_span, method="RK45")

        # Reshape: (Time, X*Y*C) -> (Time, X, Y, C)
        u_sim_native = sol.y.T.reshape(-1, nx, ny, n_channels)
        # Transpose back to (Time, Channels, X, Y)
        u_sim = np.transpose(u_sim_native, (0, 3, 1, 2))

        return u_sim

    return forecast_evolution, model
