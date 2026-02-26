import numpy as np
from pathlib import Path

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
mean_q_arr = npy2arr(mean_q_path)
singular_values_arr_flat = npy2arr(singular_values_path)

print("kk")
