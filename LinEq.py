from scipy.io import mmread
import numpy as np

#A = mmread('small_matrix.mtx')
A = mmread('matrices/big_matrix.mtx')

A_arr = A.toarray().astype('float32')
size = A_arr.shape[0]
b = np.random.randn(size)
# Numpy method 
x = np.linalg.solve(A_arr,b)

# Write to /tmp/ans
with open("answer", "wb") as f:
    f.write(x.data)

print(x[:10])