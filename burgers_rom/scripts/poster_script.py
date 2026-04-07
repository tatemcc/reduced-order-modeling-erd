#!/usr/bin/env python3
"""
Standalone script to load and plot a specific POD basis mode for a specific run ID.
Enhanced with customizable visualization parameters.
"""

import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, Normalize
import numpy as np
from scipy.interpolate import RectBivariateSpline


def main():
    # --- DATA PARAMETERS ---
    run_id = "5b007be0862a"
    equation = "burgers"
    modes_to_plot = [0, 1, 2]  # 0-indexed modes to plot
    mode_colors = ["#97279A", "#C7507B", "#F49C58"]  # Hex codes for each mode's chronos plot
    timestep = 20  # Time step of the true field to plot
    out_filename_true = f"true_field_t{timestep}_custom.png"

    # --- EDITABLE VISUALIZATION PARAMETERS ---
    FIG_SIZE = (2.5, 2)  # Width, Height in inches for 3D plots
    CHRONOS_FIG_SIZE = (5, 2)  # Width, Height in inches for chronos plots
    ENERGY_FIG_SIZE = (4, 3)  # Width, Height in inches for energy scatter plot
    CHRONOS_LINE_WIDTH = 3.0  # Line thickness for chronos plots
    DPI = 150  # Resolution

    # Perspective
    VIEW_ELEV = 20  # Elevation angle in degrees
    VIEW_AZIM = 285  # Azimuthal angle in degrees

    # Geometry & Scaling
    # set_box_aspect (x, y, z) where z controls the relative "height" of the plot area
    BOX_ASPECT = (1, 1, 0.4)
    INTERP_FACTOR = 6  # Higher = smoother surface, more memory usage

    # Appearance
    SHOW_AXES = False  # Set to False to remove grid, ticks, and labels
    CMAP_NAME = "plasma"  # Colormap (e.g., 'viridis', 'magma', 'inferno')
    VERT_EXAG = 0.2  # Vertical exaggeration for the lighting/shading effect
    LIGHT_AZIM = 225  # Angle of the light source
    LIGHT_ELEV = 45  # Elevation of the light source
    BG_COLOR = "white"  # Background color of the figure

    # -----------------------------------------

    # Resolve repository root assuming this script is in `scripts/`
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent
    outputs_dir = repo_root / "outputs"

    # Construct path to the POD basis
    pod_dir = outputs_dir / equation / run_id / "pod"
    basis_path = pod_dir / "basis_U.npy"

    if not basis_path.exists():
        print(f"Error: Could not find POD basis at {basis_path}")
        return

    print(f"Loading POD basis from {basis_path}...")
    U = np.load(basis_path)

    # Infer grid dimensions from the state vector
    n_state = U.shape[0]
    n_components = 2

    # Infer grid dimensions
    n_points = n_state // n_components
    ny = int(math.isqrt(n_points))
    nx = ny

    if n_components * ny * nx != n_state:
        print(f"Error: State vector size ({n_state}) incompatible with {n_components} components.")
        return

    def plot_vorticity_surface(flat_vector, title, filename):
        # Reshape the flat vector
        fields = flat_vector.reshape((n_components, ny, nx), order="C")
        u, v = fields[0], fields[1]

        # Calculate vorticity for height
        L = 1.0
        dx, dy = L / (nx - 1), L / (ny - 1)
        dudy = np.gradient(u, axis=0) / dy
        dvdx = np.gradient(v, axis=1) / dx
        z = dvdx - dudy

        # Calculate color data (velocity direction/angle)
        color_data = np.arctan2(v, u)

        # 3D Interpolation for smoothness
        x = np.arange(nx)
        y = np.arange(ny)
        x_fine = np.linspace(0, nx - 1, nx * INTERP_FACTOR)
        y_fine = np.linspace(0, ny - 1, ny * INTERP_FACTOR)
        xx_fine, yy_fine = np.meshgrid(x_fine, y_fine)

        spline_z = RectBivariateSpline(y, x, z, kx=3, ky=3)
        z_fine = spline_z(y_fine, x_fine)

        spline_color = RectBivariateSpline(y, x, color_data, kx=3, ky=3)
        color_fine = spline_color(y_fine, x_fine)

        # Apply lighting and shading
        light = LightSource(azdeg=LIGHT_AZIM, altdeg=LIGHT_ELEV)
        illumination = light.hillshade(z_fine, vert_exag=VERT_EXAG)
        rescaled_illumination = 0.5 + 0.5 * illumination

        color_min, color_max = color_data.min(), color_data.max()
        if color_max - color_min < 1e-9:
            color_max = color_min + 1.0
        color_norm = Normalize(vmin=color_min, vmax=color_max)
        cmap_obj = plt.get_cmap(CMAP_NAME)

        colors_from_map = cmap_obj(color_norm(color_fine))[:, :, :3]
        rgb_vertex = colors_from_map * rescaled_illumination[..., np.newaxis]

        # Calculate face colors (averaging corners)
        rgb_face = (
            rgb_vertex[:-1, :-1, :]
            + rgb_vertex[1:, :-1, :]
            + rgb_vertex[:-1, 1:, :]
            + rgb_vertex[1:, 1:, :]
        ) / 4.0

        # Initialize Plot
        fig = plt.figure(figsize=FIG_SIZE, dpi=DPI, facecolor=BG_COLOR)
        ax = fig.add_subplot(1, 1, 1, projection="3d", facecolor=BG_COLOR)

        # Plot Surface
        _ = ax.plot_surface(
            xx_fine,
            yy_fine,
            z_fine,
            facecolors=rgb_face,
            rstride=1,
            cstride=1,
            antialiased=True,
            shade=False,
        )

        # Apply View Controls
        ax.view_init(elev=VIEW_ELEV, azim=VIEW_AZIM)
        ax.set_box_aspect(BOX_ASPECT)  # Controls the height relative to width

        ax.set_xlim(0, nx - 1)
        ax.set_ylim(0, ny - 1)
        ax.dist = 8
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

        # Visibility Controls
        if not SHOW_AXES:
            ax.set_axis_off()
        else:
            ax.set_title(title)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("Vorticity")

        plt.tight_layout()
        plt.savefig(filename, dpi=DPI, bbox_inches="tight", pad_inches=0 if not SHOW_AXES else 0.1)
        plt.close(fig)
        print(f"Saved plot to {Path(filename).resolve()}")

    # 1. Plot the true field at the requested timestep
    q_true_path = outputs_dir / equation / run_id / "rollout" / "q_true.npy"
    if not q_true_path.exists():
        print(f"Error: Could not find true trajectory at {q_true_path}")
        return

    print(f"Loading true trajectory from {q_true_path}...")
    q_true = np.load(q_true_path)

    if timestep >= q_true.shape[0]:
        print(f"Error: Timestep {timestep} exceeds trajectory length {q_true.shape[0]}")
        return

    true_flat = q_true[timestep]
    plot_vorticity_surface(
        true_flat, f"True Field (t={timestep}) - Vorticity (Run: {run_id[:8]})", out_filename_true
    )

    # 1.5 Plot Mode Energy vs. Index
    s_path = pod_dir / "singular_values.npy"
    if s_path.exists():
        print(f"Loading singular values from {s_path}...")
        s = np.load(s_path)
        energy = s**2
        max_energy = energy[0] if energy.size > 0 else 0.0
        normalized_energy = energy / max_energy if max_energy > 0 else energy
        ranks = np.arange(1, len(s) + 1)

        cum_energy = np.cumsum(energy)
        total_energy = cum_energy[-1] if cum_energy.size > 0 else 1.0
        # The rank where cumulative energy reaches or exceeds 80%
        r_80 = int(np.searchsorted(cum_energy / total_energy, 0.80, side="left")) + 1

        fig, ax = plt.subplots(figsize=ENERGY_FIG_SIZE, dpi=DPI, facecolor=BG_COLOR)
        ax.scatter(ranks, normalized_energy, color="#97279A", s=20)

        if r_80 <= len(s):
            ax.axvline(x=r_80 + 0.5, color="#C7507B", linestyle="--", linewidth=3)

        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(2.0)

        out_filename_energy = "pod_energy_scatter_custom.png"
        fig.tight_layout()
        plt.savefig(
            out_filename_energy, dpi=DPI, bbox_inches="tight", pad_inches=0, facecolor=BG_COLOR
        )
        plt.close(fig)
        print(f"Saved energy scatter plot to {Path(out_filename_energy).resolve()}")
    else:
        print(f"Warning: singular_values.npy not found at {s_path}. Skipping energy plot.")

    # Load true and predicted chronos once
    a_true_path = outputs_dir / equation / run_id / "rollout" / "A_true.npy"
    a_pred_path = outputs_dir / equation / run_id / "rollout" / "A_pred.npy"

    A_true = None
    A_pred = None
    if a_true_path.exists() and a_pred_path.exists():
        print(f"Loading true and predicted coefficients from {a_true_path.parent}...")
        A_true = np.load(a_true_path)
        A_pred = np.load(a_pred_path)
    else:
        print(
            f"Warning: A_true.npy or A_pred.npy not found in {a_true_path.parent}. Skipping chronos plots."
        )

    # 2. Loop over modes to plot POD modes and chronos
    for mode_idx, line_color in zip(modes_to_plot, mode_colors):
        if U.shape[1] <= mode_idx:
            print(f"Warning: POD basis contains {U.shape[1]} modes, cannot plot mode {mode_idx}.")
            continue

        out_filename_mode = f"pod_mode_{mode_idx}_custom.png"
        out_filename_chronos = f"chronos_mode_{mode_idx}_custom.png"

        # Plot the requested POD mode
        mode_flat = U[:, mode_idx]
        plot_vorticity_surface(
            mode_flat, f"POD Mode {mode_idx} - Vorticity (Run: {run_id[:8]})", out_filename_mode
        )

        # Plot the true vs predicted chronos
        if A_true is not None and A_pred is not None:
            if mode_idx < A_true.shape[1] and mode_idx < A_pred.shape[1]:
                T = A_true.shape[0]
                time_axis = np.arange(T)

                fig, ax = plt.subplots(figsize=CHRONOS_FIG_SIZE, dpi=DPI, facecolor=BG_COLOR)

                ax.plot(
                    time_axis,
                    A_true[:, mode_idx],
                    color=line_color,
                    linestyle="-",
                    linewidth=CHRONOS_LINE_WIDTH,
                )
                ax.plot(
                    time_axis,
                    A_pred[:, mode_idx],
                    color=line_color,
                    linestyle="--",
                    linewidth=CHRONOS_LINE_WIDTH,
                )

                # ax.set_axis_off()
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_visible(False)
                    # spine.set_linewidth(2.0)

                fig.tight_layout()
                plt.savefig(out_filename_chronos, dpi=DPI, bbox_inches="tight", pad_inches=0)
                plt.close(fig)
                print(f"Saved chronos plot to {Path(out_filename_chronos).resolve()}")
            else:
                print(f"Error: Mode {mode_idx} out of bounds for coefficient matrices.")


if __name__ == "__main__":
    main()
