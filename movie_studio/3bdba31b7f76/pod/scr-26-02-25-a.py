import numpy as np

arr = np.load("coeffs_A.npy", allow_pickle=True)
print(arr.shape, arr.dtype)
print(arr[:5])
