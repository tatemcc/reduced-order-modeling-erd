import numpy as np
from pathlib import Path
import sys, os
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

