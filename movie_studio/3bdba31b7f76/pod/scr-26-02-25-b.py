import numpy as np
from pathlib import Path

files = [p.name for p in Path.cwd().glob("*.npy")]
for fname in files:
    print("\n" + fname)
    arr = np.load(fname, allow_pickle=True)
    print(arr.shape, arr.dtype)
