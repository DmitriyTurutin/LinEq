#ifndef PARSE_MATRIX_H
#define PARSE_MATRIX_H
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef struct csr_matrix {
    int nrows;
    int ncols;
    int nnz;
    float *values;
    int *col_idx;
    int *row_ptr;
} csr_matrix;

/// @brief Concatenates current working directory with filename
/// @param filename
/// @return
char *get_path(char *filename) {
    char current_dir[PATH_MAX];
    if (getcwd(current_dir, sizeof(current_dir)) == NULL) {
        exit(-1);
    }

    char *full_path = (char *)malloc(PATH_MAX);
    if (full_path == NULL) {
        exit(-1);
    }

    sprintf(full_path, "%s/%s", current_dir, filename);
    return full_path;
}

/// @brief Read .mtx file
/// @param filename
/// @return file in char* format without comment
char *read_mtx_file(char *filename) {
    FILE *file;
    char *buffer;
    long file_len;

    // Open the file in binary-read mode
    file = fopen(filename, "rb");
    if (!file) {
        fprintf(stderr, "Failed to open file\n");
        return NULL;
    }

    char line[256];

    // Get file size
    fseek(file, 0, SEEK_END);
    file_len = ftell(file);
    fseek(file, 0, SEEK_SET);

    // Allocate memory for the buffer
    buffer = (char *)malloc(file_len + 1);
    if (!buffer) {
        fprintf(stderr, "Failed to allocate memory\n");
        fclose(file);
        return NULL;
    }

    // Skip first line with comments
    fgets(line, sizeof(line), file);

    // Read the rest of the file into the buffer
    while (fgets(line, sizeof(line), file)) {
        strcat(buffer, line);
    }

    // Close the file
    fclose(file);

    return buffer;
}

csr_matrix new_csr(char *filename) {
    // Get info
    char *input = read_mtx_file(filename);
    char *info = strtok(input, "\n");

    // Parse the info string to extract the matrix dimensions and number of
    // non-zero elements
    int nrows, ncols, nnz;
    sscanf(info, "%d %d %d", &nrows, &ncols, &nnz);

    csr_matrix csr;

    csr.nrows = nrows;
    csr.ncols = ncols;
    csr.nnz = nnz;

    csr.row_ptr = (int *)malloc(sizeof(int) * (nrows + 1));
    csr.col_idx = (int *)malloc(sizeof(int) * nnz);
    csr.values = (float *)malloc(sizeof(float) * nnz);

    char *line;
    char *tok;
    int row = 0;
    int col = 0;
    float val = 0.0;
    int idx = 0;

    for (int i = 0; i < nnz; i++) {
        line = strtok(NULL, "\n");
        int row, col;
        float val;
        sscanf(line, "%d %d %f", &row, &col, &val);
        // Convert into format starting from 0
        row--;
        col--;
        // Add the value to the values array
        csr.values[idx] = val;

        // Add the column index to the col_idx array
        csr.col_idx[idx] = col;

        // If this is the first non-zero element in the current row, update
        // the corresponding entry in the row_ptr array
        if (i == 0 || row > (idx - 1) / (nrows)) {
            csr.row_ptr[row] = idx;
        }

        idx++;
    }

    // Set the final entry in the row_ptr array to nnz
    csr.row_ptr[nrows] = nnz;

    free(input);

    return csr;
}

void print_csr_matrix(csr_matrix matrix) {

    printf("Number of rows: %d\n", matrix.nrows);
    printf("Number of cols: %d\n", matrix.ncols);
    printf("Number of non zero elements: %d\n", matrix.nnz);
    printf("=======================Values:========================\n");
    const int size = matrix.nnz;
    for (int i = 0; i < size; i++) {
        printf("%f\n", matrix.values[i]);
    }
    printf("=======================END OF VALUES========================\n");
}

float *conjugate_gradient(csr_matrix *matrix, float *b) {
    float *x, *r, *p, *Ap;
    float alpha, beta, rsold, rsnew;
    int n = matrix->ncols;
    int max_iter = 1000;
    float tol = 1e-6;

    x = (float *)malloc(n * sizeof(float));
    r = (float *)malloc(n * sizeof(float));
    p = (float *)malloc(n * sizeof(float));
    Ap = (float *)malloc(n * sizeof(float));

    // Set initial guess x to zero
    for (int i = 0; i < n; i++) {
        x[i] = 0;
    }

    // Calculate initial residual r = b - Ax
    for (int i = 0; i < n; i++) {
        r[i] = b[i];
        for (int j = matrix->row_ptr[i]; j < matrix->row_ptr[i + 1]; j++) {
            r[i] -= matrix->values[j] * x[matrix->col_idx[j]];
        }
        p[i] = r[i];
    }

    rsold = 0;
    for (int i = 0; i < n; i++) {
        rsold += r[i] * r[i];
    }

    for (int k = 0; k < max_iter; k++) {
        // Calculate Ap = A * p using sparse matrix-vector multiplication
        for (int i = 0; i < n; i++) {
            Ap[i] = 0;
            for (int j = matrix->row_ptr[i]; j < matrix->row_ptr[i + 1]; j++) {
                Ap[i] += matrix->values[j] * p[matrix->col_idx[j]];
            }
        }

        // Calculate alpha = rsold / p' * Ap
        alpha = 0;
        for (int i = 0; i < n; i++) {
            alpha += p[i] * Ap[i];
        }
        alpha = rsold / alpha;

        // Update solution x = x + alpha * p
        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
        }

        // Update residual r = r - alpha * Ap
        for (int i = 0; i < n; i++) {
            r[i] -= alpha * Ap[i];
        }

        rsnew = 0;
        for (int i = 0; i < n; i++) {
            rsnew += r[i] * r[i];
        }

        // Check for convergence
        if (sqrt(rsnew) < tol) {
            break;
        }

        // Calculate beta = rsnew / rsold
        beta = rsnew / rsold;

        // Update direction p = r + beta * p
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * p[i];
        }

        rsold = rsnew;
    }

    free(r);
    free(p);
    free(Ap);

    return x;
}

float *generate_random_b(size_t size) {
    float *b = malloc(size * sizeof(float));
    if (b == NULL) {
        return NULL;
    }

    srand(time(NULL));

    for (size_t i = 0; i < size; i += 2) {
        float u1 = ((float)rand() + 1) / ((float)RAND_MAX + 2);
        float u2 = ((float)rand() + 1) / ((float)RAND_MAX + 2);

        float z1 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        float z2 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);

        b[i] = z1;
        if (i + 1 < size) {
            b[i + 1] = z2;
        }
    }

    return b;
}

#endif /* PARSE_MATRIX_H */
