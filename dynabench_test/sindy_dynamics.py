import numpy as np
import pysindy as ps
from scipy.integrate import solve_ivp

# Suppress warnings
# warnings.filterwarnings("ignore", category=UserWarning, module="pysindy")
# warnings.filterwarnings("ignore", category=RuntimeWarning)


def train_sindy_pde(u, t, points):
    print("      [SINDy] Starting GRID training pipeline (Native PDEFIND Mode)...")

    # 1. Transpose Data to (X, Y, Time, Channels)
    if u.ndim == 4 and u.shape[1] == 1:
        u_native = np.transpose(u, (2, 3, 0, 1))
    elif u.ndim == 4 and u.shape[3] == 1:
        u_native = np.transpose(u, (1, 2, 0, 3))
    elif u.ndim == 3:
        u_expanded = u[..., np.newaxis]
        u_native = np.transpose(u_expanded, (1, 2, 0, 3))
    else:
        u_native = u

    print(f"      [Pre-process] Transposed data to {u_native.shape} (X, Y, T, C)")

    # 2. Define Grids
    nx, ny = u_native.shape[0], u_native.shape[1]
    grid_x = np.linspace(0, 1, nx)
    grid_y = np.linspace(0, 1, ny)
    X, Y = np.meshgrid(grid_x, grid_y, indexing="ij")
    spatial_grid = np.stack([X, Y], axis=-1)

    # 3. Define PDE Library
    library_functions = [lambda x: x]
    function_names = [lambda x: "u"]

    custom_lib = ps.CustomLibrary(
        library_functions=library_functions, function_names=function_names
    )

    # diff_method = ps.SmoothedFiniteDifference(smoother_kws={'window_length': 5})

    library = ps.PDELibrary(
        function_library=custom_lib,
        derivative_order=2,
        spatial_grid=spatial_grid,
        include_bias=True,
        is_uniform=True,
        periodic=True,
    )

    # STLSQ is often more stable for physical systems
    optimizer = ps.STLSQ(
        threshold=0.1,  # Increase this to prune weak/unstable terms
        alpha=1e-05,  # L2 regularization to keep coefficients small
        normalize_columns=True,
    )

    # 5. Fit Model
    print("      [SINDy] Fitting model using SR3...")
    dt = t[1] - t[0]
    model = ps.SINDy(feature_library=library, optimizer=optimizer)
    model.fit(u_native, t=dt)

    # 6. Print Results
    print("\n" + "=" * 40)
    print("      IDENTIFIED DYNAMICS (Native PDEFIND)")
    print("=" * 40)
    model.print()

    print("\n--- Advection Parameters ---")
    if len(model.coefficients()) > 0:
        coeffs = model.coefficients()[0]
        features = model.get_feature_names()
        for name, val in zip(features, coeffs):
            if abs(val) > 0.001:
                print(f"  {name:<10} : {val:.5f}")
                if "x0_1" in name and "x0_11" not in name:
                    print(f"             ^ c_x ≈ {-val:.4f}")
                if "x0_2" in name and "x0_22" not in name:
                    print(f"             ^ c_y ≈ {-val:.4f}")

    # 7. Simulation Wrapper (Crash-Proof)
    def forecast_evolution(u0, t_span):
        u0_flat = u0.flatten()
        lib_grid = model.feature_library.spatial_grid
        nx, ny = lib_grid.shape[0], lib_grid.shape[1]

        def rhs(t, u_flat):
            # 1. Stability Check: If values explode, return 0 to stop integration safely
            if not np.all(np.isfinite(u_flat)) or np.max(np.abs(u_flat)) > 1e3:
                return np.zeros_like(u_flat)

            u_reshaped = u_flat.reshape(nx, ny, 1)

            # 2. Predict derivative
            # Note: We must suppress errors here if the library generates NaNs
            try:
                u_dot = model.predict(u_reshaped)
            except ValueError:
                return np.zeros_like(u_flat)

            return u_dot.flatten()

        # Wrap integration in try-except to catch Sklearn/LSODA crashes
        try:
            # Switched back to RK45 (less strict about stiffness warnings than LSODA)
            sol = solve_ivp(rhs, (t_span[0], t_span[-1]), u0_flat, t_eval=t_span, method="RK45")

            u_sim = sol.y.T.reshape(sol.y.shape[1], nx, ny, 1)
            steps_computed = u_sim.shape[0]

        except Exception as e:
            print(f"      [Warning] Integration crashed: {e}")
            steps_computed = 0
            u_sim = np.zeros((0, nx, ny, 1))

        # Pad with NaNs if incomplete
        if steps_computed < len(t_span):
            print(
                f"      [Warning] Simulation unstable! Stopped at step {steps_computed}/{len(t_span)}."
            )
            pad_len = len(t_span) - steps_computed
            if steps_computed == 0:
                # If it crashed immediately, return all NaNs
                u_sim = np.full((len(t_span), nx, ny, 1), np.nan)
            else:
                padding = np.full((pad_len, nx, ny, 1), np.nan)
                u_sim = np.vstack([u_sim, padding])

        return u_sim

    return forecast_evolution, model
