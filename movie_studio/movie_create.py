import numpy as np
from pathlib import Path
import sys, os
import math
import plotly.graph_objects as go
import subprocess
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'burgers_rom', 'burgers_rom'))
from snapshot import state_vec_to_fields, SnapshotLayout

# Settings: 
# ONLY_TOP: Boolean, just the top mode (True) or all included modes (False)
ONLY_TOP = False
# VORTICITY: Boolean, plots velocity vorticity on z axis (True), plots velocity magnitude on z axis (False)
VORTICITY = True

# Choose data files
# Note: by default, stuff happens relative to where you run python, not where the script lives
data_dir_path = "./3bdba31b7f76/pod/"
basis_U_path = data_dir_path + "basis_U.npy"
coeffs_A_path = data_dir_path + "coeffs_A.npy"
mean_q_path = data_dir_path + "mean_q.npy"
singular_values_path = data_dir_path + "singular_values.npy"
data_filepaths = [basis_U_path, coeffs_A_path, mean_q_path, singular_values_path]

# Check these all exist
assert(Path(data_dir_path).is_dir())
for path in data_filepaths:
    assert(Path(path).is_file())

# Unpack the datas into numpy arrays
def npy2arr(fp):
    return np.load(fp, allow_pickle=True)
basis_U_arr_flat = npy2arr(basis_U_path)
coeffs_A_arr = npy2arr(coeffs_A_path)
mean_q_arr_flat = npy2arr(mean_q_path)
singular_values_arr = npy2arr(singular_values_path)

# Reinflate the spatial data (1D back to 3D)
# assert() statements are used to check dimensionality so pay attention to those
space_product = (basis_U_arr_flat.shape)[0]   # C*ny*nx
nt = (coeffs_A_arr.shape)[1]
total_num_modes = (singular_values_arr.shape)[0]
assert(total_num_modes == min(space_product, nt))
nsm = (coeffs_A_arr.shape)[0]   # number of selected modes (top nsm ranked by singular values)
assert(nsm == (basis_U_arr_flat.shape)[1])
# !!Caution: making assumptions here related to space_product, C, ny, and nx (hard-coding)
C = 2
ny = 15
nx = 15
assert(space_product == C*ny*nx)
layout = SnapshotLayout(ny=ny, nx=nx, n_components=C)
# basis_U_arr_flat has dimensions (space_product, nsm)
# for the space_product to be unflattened, the raw array can be transposed first and then each row can be reinflated
# (I think reinflated sounds cooler than unflattened)
basis_U_afT = np.transpose(basis_U_arr_flat)
basis_U_arr = np.array([state_vec_to_fields(basis_U_afT[mdx,:], layout) for mdx in range(0, nsm)])
assert(basis_U_arr.shape == (nsm, C, ny, nx))
# mean_q_path_flat has dimensions (space_product,)
mean_q_arr = state_vec_to_fields(mean_q_arr_flat, layout)
assert(mean_q_arr.shape == (C, ny, nx))

# Add stuff for what want movie out of
data_xy = np.zeros((nt, C, ny, nx), float)
for tdx in range(nt):
    if ONLY_TOP:
        data_xy[tdx, :, :, :] = coeffs_A_arr[0, tdx] * basis_U_arr[0, :, :, :]
    else:
        for mdx in range(nsm):
            data_xy[tdx, :, :, :] += coeffs_A_arr[mdx, tdx] * basis_U_arr[mdx, :, :, :]
        data_xy[tdx, :, :, :] += mean_q_arr

# Convert spatial (x, y) into (magnitude, direction) or (vorticity, direction)
data = np.zeros((nt, C, ny, nx), float)
for tdx in range(nt):
    for ydx in range(ny):
        for xdx in range(nx):
            x_comp = data_xy[tdx, 0, ydx, xdx]
            y_comp = data_xy[tdx, 1, ydx, xdx]
            if VORTICITY:
                # Vorticity
                # Described by the 2D curl (this is a 2D vector field)
                # d(vy)/d(x) - d(dx)/d(y)
                # derivative given by formula (x_(n+1) - x_(n-1)) / (2h)
                # For edge cases, given by the single-sided version of that
                h = 1
                if (xdx > 0) and (xdx < nx-1):
                    vy_nm1 = data_xy[tdx, 1, ydx, xdx-1]
                    vy_np1 = data_xy[tdx, 1, ydx, xdx+1]
                    dvydx = (vy_np1 - vy_nm1) / (2*h)
                elif (xdx == 0):
                    vy_0 = data_xy[tdx, 1, ydx, 0]
                    vy_1 = data_xy[tdx, 1, ydx, 1]
                    dvydx = (vy_1 - vy_0) / h
                elif (xdx == nx-1):
                    vy_n2 = data_xy[tdx, 1, ydx, -2]
                    vy_n1 = data_xy[tdx, 1, ydx, -1]
                    dvydx = (vy_n1 - vy_n2) / h
                if (ydx > 0) and (ydx < ny-1):
                    vx_nm1 = data_xy[tdx, 0, ydx-1, xdx]
                    vx_np1 = data_xy[tdx, 0, ydx+1, xdx]
                    dvxdy = (vx_np1 - vx_nm1) / (2*h)
                elif (ydx == 0):
                    vx_0 = data_xy[tdx, 0, 0, xdx]
                    vx_1 = data_xy[tdx, 0, 1, xdx]
                    dvxdy = (vx_1 - vx_0) / h
                elif (ydx == ny-1):
                    vx_n2 = data_xy[tdx, 0, -2, xdx]
                    vx_n1 = data_xy[tdx, 0, -1, xdx]
                    dvxdy = (vx_n1 - vx_n2) / h
                data[tdx, 0, ydx, xdx] = dvydx - dvxdy
            else:
                # Magnitude
                data[tdx, 0, ydx, xdx] = math.sqrt(math.pow(x_comp,2) + math.pow(y_comp,2))
            # Direction
            data[tdx, 1, ydx, xdx] = math.atan2(y_comp, x_comp)


# Got lots of help from Copilot on this part:
FPS = 30
FRAME_DURATION_MS = int(1000 / FPS)
OUTDIR = Path("frames")
MP4_NAME = "plasma-movie-"
if ONLY_TOP:
    MP4_NAME += "top-mode-"
if VORTICITY:
    MP4_NAME += "vorticity"
else:
    MP4_NAME += "magnitude"
MP4_NAME += ".mp4"


# Choose your output resolution:
W, H = 1280, 720
SCALE = 2  # multiplies pixel density (higher => sharper, larger files)

# ----------------------------
# CYCLIC "PHASE" COLOR SCALE
# (-pi and +pi are same color)
# ----------------------------
phase_colorscale = [
    [0.0,   "rgb(255,0,0)"],     # red
    [1/6,   "rgb(255,255,0)"],   # yellow
    [2/6,   "rgb(0,255,0)"],     # green
    [3/6,   "rgb(0,255,255)"],   # cyan
    [4/6,   "rgb(0,0,255)"],     # blue
    [5/6,   "rgb(255,0,255)"],   # magenta
    [1.0,   "rgb(255,0,0)"],     # back to red
]


def make_figure(data, x=None, y=None):
    nt, _, ny, nx = data.shape
    if x is None:
        x = np.arange(nx)
    if y is None:
        y = np.arange(ny)

    # Heights and colours for initial timestep
    z0 = data[0, 0]
    c0 = data[0, 1]

    surf = go.Surface(
            x=x, y=y, z=z0,
            surfacecolor=c0,
            colorscale=phase_colorscale,
            cmin=-math.pi, cmax=math.pi,
            colorbar=dict(title="dir (rad)"),
            showscale=True,
            # Optional: a bit of lighting can improve "continuous" look
            lighting=dict(ambient=0.35, diffuse=0.7, specular=0.2, roughness=0.8),
    )

    fig = go.Figure(data=[surf])

    # Lock z scale
    if VORTICITY:
        absmaxvort = max(abs(np.min(data[:, 0, :, :])), abs(np.max(data[:, 0, :, :])))
        zmin = -absmaxvort
        zmax = +absmaxvort
    else:
        zmin = 0
        zmax = np.max(data[:, 0, :, :])
    fig.update_layout(
            scene=dict(
                zaxis=dict(
                    range=[zmin, zmax]
                )
            )
    )

    # Build frames (only update z & surfacecolor per frame)
    frames = []
    for tdx in range(nt):
        frames.append(
                go.Frame(
                    data=[go.Surface(z=data[tdx, 0], surfacecolor=data[tdx, 1])],
                    name=str(tdx)
                )
        )
    fig.frames = frames

    return fig



def export_mp4(fig, data):
    """
    Writes PNG frames using Kaleido, then stitches them into MP4 with ffmpeg.
    Plotly supports write_image via Kaleido. [3](https://plotly.com/python-api-reference/generated/plotly.io.write_image.html)[4](https://plotly.com/python/static-image-export/)
    Plotly animations aren't directly saved to mp4, so we do frame export. [2](https://community.plotly.com/t/how-to-export-animation-and-save-it-in-a-video-format-like-mp4-mpeg-or/64621)
    """
    nt = data.shape[0]
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # Render frames
    for t in range(nt):
        fig.update_traces(
            z=data[t, 0],
            surfacecolor=data[t, 1],
            selector=dict(type="surface"),
        )

        fig.write_image(str(OUTDIR / f"frame_{t:05d}.png"),
                        width=W, height=H, scale=SCALE)

    # Stitch with ffmpeg (must be installed on your system)
    # On macOS: brew install ffmpeg
    # On Ubuntu: sudo apt-get install ffmpeg
    # On Windows: install ffmpeg and ensure it's on PATH
    cmd = [
        "ffmpeg", "-y",
        "-framerate", str(FPS),
        "-i", str(OUTDIR / "frame_%05d.png"),
        "-c:v", "libx264",
        "-pix_fmt", "yuv420p",
        MP4_NAME,
    ]
    subprocess.run(cmd, check=True)
    print(f"Saved {MP4_NAME}")


fig = make_figure(data)
fig.show()           # interactive preview
export_mp4(fig, data)


print("so far so good")
