#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "parasail.h"

int main(void)
{
    const parasail_matrix_t *internal_matrix = NULL;
    parasail_matrix_t *user_matrix = NULL;
    
    internal_matrix = parasail_matrix_lookup("blosum62");
    if (NULL == internal_matrix) {
        fprintf(stderr, "matrix lookup failed");
        exit(EXIT_FAILURE);
    }

    user_matrix = parasail_matrix_create("ACGT", 2, -1);
    if (NULL == user_matrix) {
        fprintf(stderr, "matrix create failed");
        exit(EXIT_FAILURE);
    }

    parasail_matrix_free(user_matrix);

    user_matrix = parasail_matrix_copy(internal_matrix);
    if (NULL == user_matrix) {
        fprintf(stderr, "matrix copy failed");
        exit(EXIT_FAILURE);
    }

    parasail_matrix_set_value(user_matrix, 10, 10, 100);
    if (100 != user_matrix->max) {
        fprintf(stderr, "matrix set value failed");
        exit(EXIT_FAILURE);
    }

    parasail_matrix_free(user_matrix);

    return 0;
}

