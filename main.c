#include "parse_mtx.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

uint64_t nanos() {
  struct timespec start;
  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  return (uint64_t)start.tv_sec*1000000000 + (uint64_t)start.tv_nsec;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Provide only 1 argument!");
        exit(-1);
    }

    char *filename = argv[1];

    // Parse file and get matrix in CSR format
    csr_matrix matrix = new_csr(filename);
    // print_csr_matrix(matrix);
    float *b = generate_random_b(matrix.nrows);
    uint64_t start = nanos();
    float *x = conjugate_gradient(&matrix, b);
    uint64_t end = nanos();

    printf("\n\n%ld\n\n", end-start);

    for (int i = 0; i < 10; ++i) {
        printf("%e ", x[i]);
    }
    printf("\n");

    /// TEST IF RIGHT
    csr_matrix test_matrix = {.nrows = 3,
                              .ncols = 3,
                              .nnz = 5,
                              .values = (float[]){2.0, 1.0, 1.0, 1.0, 3.0},
                              .col_idx = (int[]){0, 1, 1, 0, 2, 2},
                              .row_ptr = (int[]){0, 2, 4, 6}};

    float test_b[] = {3, 7, 2.5};
    float *test_x = conjugate_gradient(&test_matrix, test_b);
    for (int i = 0; i < 3; ++i) {
        printf("%e ", test_x[i]);
    }
    // Free memory
    free(x);
    free(b);

    return 0;
}
