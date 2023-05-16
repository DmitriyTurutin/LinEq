#include <stdio.h>
#include <stdlib.h>

#define N 3
typedef struct csr_matrix {
    int nrows;
    int ncols;
    int nnz;
    float *values;
    int *col_idx;
    int *row_ptr;

} csr_matrix;
// {30, 70, 180, 80}
float *matvec_product(csr_matrix *mtx, float *vec) {
    // Allocate result vector
    float *result = malloc(mtx->ncols * sizeof(float));
    // Max value to iterate through array of row pointers
    size_t size = mtx->nrows + 1;
    // Iterate through matrix and multily it on vector b
    for (size_t i = 0; i < size; ++i) {
        for (int j = mtx->row_ptr[i]; j < mtx->row_ptr[i + 1]; ++j) {
            // Use _m256_fmadd_pd instead  
            result[i] += mtx->values[j] * vec[mtx->col_idx[j]];
        }
    }
    return result;
}

void print_vec(float *vec, size_t size) {
    for (int i = 0; i < size; ++i) {
        printf("%f ", vec[i]);
    }
    printf("\n");
}

// TODO: research how to implement CG
// Ax=b
// P = {p_1, ..., p_n} -- set of mutualy conjugate vectors with respect to A
// p_i^TAp_j = 0 ∀ i != j
// x* = Σa_n p_i
// a_k = (p_k, b)/(p_k, p_k)_A
float *cg(csr_matrix *A, float *b) {
    float *x0 = malloc(N * sizeof(float));
    float *r0 = matvec_product(A, x0);
    return x0;
}

int main() {
    csr_matrix print_mat = {.nrows = 4,
                            .ncols = 6,
                            .nnz = 8,
                            .values = (float[]){10, 20, 30, 40, 50, 60, 70, 80},
                            .col_idx = (int[]){0, 1, 1, 3, 2, 3, 4, 5},
                            .row_ptr = (int[]){0, 2, 4, 7, 8}};

    float b[] = {1, 2, 3, 4, 5, 1};
    printf("\n\n<====================PRINT PRODUCT====================>\n");
    float *product = matvec_product(&print_mat, b);
    print_vec(product, 4);
    printf("<====================PRINT PRODUCT====================>\n\n");

    free(product);

    for (int i = 0; i < print_mat.nrows + 1; ++i) {
        for (int j = print_mat.row_ptr[i]; j < print_mat.row_ptr[i + 1]; ++j) {
            printf("%f ", print_mat.values[j]);
        }
        printf("\n");
    }

    return 0;

    // MATRIX
    csr_matrix test_matrix = {
        .nrows = N,
        .ncols = N,
        .nnz = 5,
        .values = (float[]){2., 1., 1., 1., 3.}, // is of size nnz
        .col_idx = (int[]){0, 1, 1, 0, 2, 2},    // is of size nnz
        .row_ptr = (int[]){0, 2, 4, 6}};         // is of size n + 1

    float test_b[] = {3, 7, 2.5};

    float *test_x = cg(&test_matrix, test_b);
    for (int i = 0; i < N; ++i) {
        printf("%e ", test_x[i]);
    }
    return 0;
}
