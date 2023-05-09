#include "parse_mtx.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Provide only 1 argument!");
        exit(-1);
    }

    char *filename = argv[1];

    char *path = get_path(filename);
    char *file = read_mtx_file(filename);
    csr_matrix matrix = new_csr(filename);
    print_csr_matrix(matrix);
    // printf("%s", file);
    float *b = generate_random_b(matrix.nrows);



    free(x);
    free(b);
    free(file);

    return 0;
}
