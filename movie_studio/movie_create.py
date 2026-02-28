import numpy as np
from pathlib import Path
import sys, os
import math
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'burgers_rom', 'burgers_rom'))
from snapshot import state_vec_to_fields, SnapshotLayout

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

# Convert spatial (x, y) into (magnitude, direction)
basis_U_md = np.zeros((nsm, C, ny, nx), float)
mean_q_md = np.zeros((C, ny, nx), float)
for mdx in range(nsm):
    for ydx in range(ny):
        for xdx in range(nx):
            x_comp = basis_U_arr[mdx, 0, ydx, xdx]
            y_comp = basis_U_arr[mdx, 1, ydx, xdx]
            # Magnitude
            basis_U_md[mdx, 0, ydx, xdx] = math.sqrt(math.pow(x_comp,2) + math.pow(y_comp,2))
            # Direction
            basis_U_md[mdx, 1, ydx, xdx] = math.atan2(y_comp, x_comp)
for ydx in range(ny):
    for xdx in range(nx):
        x_comp = mean_q_arr[0, ydx, xdx]
        y_comp = mean_q_arr[1, ydx, xdx]
        # Magnitude
        mean_q_md[0, ydx, xdx] = math.sqrt(math.pow(x_comp,2) + math.pow(y_comp,2))
        # Direction
        mean_q_md[1, ydx, xdx] = math.atan2(y_comp, x_comp)

print("so far so good")
