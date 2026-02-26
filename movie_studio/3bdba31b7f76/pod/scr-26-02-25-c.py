import numpy as np
import sys

# first arg: filename
# second arg: print the first n elements (-1 for all)

arr = np.load(sys.argv[1], allow_pickle=True)
print(arr.shape, arr.dtype)
print(arr[:int(sys.argv[2])])
