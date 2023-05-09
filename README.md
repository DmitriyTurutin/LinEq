<p align="center">
    <img src="https://github.com/DmitriyTurutin/LinEq/blob/master/docs/logo.png?raw=true">
</p>

------------------------

### About library

LinEq is a C Hader library to solve systems of linear equtaions


## About CG(conjugate gradient)

The conjugate gradient (CG) algorithm is a popular iterative method for solving linear systems of equations of the form Ax = b, where A is a symmetric, positive-definite matrix, and b is a vector. The algorithm seeks to find a solution x that minimizes the quadratic form 1/2 *x^T* A *x - b^T* x.

When A is a sparse matrix stored in the compressed sparse row (CSR) format, the CG algorithm can be adapted to take advantage of the sparsity structure of the matrix. Here's a rough outline of how the algorithm works:

Initialize the solution vector x to zero, and the residual vector r = b - Ax to b.
Initialize the search direction d = r.
Repeat until convergence:

1. Compute the product Ad = A *d using sparse matrix-vector multiplication.
2. Compute the step size alpha = (r^T* r) / (d^T *Ad).
3. Update the solution x = x + alpha* d.
4. Update the residual r = r - alpha *Ad.
5. Compute the squared residual norm rho = r^T* r.
6. If rho is small enough, terminate the iteration.
7. Compute the new search direction beta = (rho_new / rho_old).
8. Update the search direction d = r + beta * d.

In step 3a, the matrix-vector multiplication can be efficiently computed using the CSR format, which stores the non-zero entries of the matrix in three arrays: the values array, the column indices array, and the row pointer array. The algorithm only needs to access the non-zero entries of A and can skip over the zero entries, which can greatly reduce the computational cost.

The CG algorithm is guaranteed to converge to the exact solution in at most n iterations, where n is the dimension of the matrix A. However, the rate of convergence can depend on the condition number of the matrix, which is the ratio of the largest to the smallest eigenvalue. If the condition number is large, the CG algorithm may converge slowly, and preconditioning techniques can be used to improve the convergence rate.
