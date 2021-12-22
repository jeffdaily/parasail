#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "parasail.h"
#include "parasail/io.h"
#include "parasail/matrices/blosum62.h"

static void print_matrix(const parasail_matrix_t *matrix) {
    int i = 0;
    int j = 0;
    const char *alphabet = NULL;
    alphabet = matrix->alphabet;

    printf("matrix '%s'\n", matrix->name);
    printf("  ");
    for (i=0; i<matrix->size; ++i) {
        printf("%4c", alphabet[i]);
    }
    printf("\n");
    for (j=0; j<matrix->length; ++j) {
        if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
            printf("%c ", alphabet[j]);
        }
        else {
            if (matrix->query != NULL) {
                printf("%c ", matrix->query[j]);
            }
            else {
                printf("* ");
            }
        }
        for (i=0; i<matrix->size; ++i) {
            printf("%4d", matrix->matrix[j*matrix->size + i]);
        }
        printf("\n");
    }
    printf("max = %d\n", matrix->max);
    printf("min = %d\n", matrix->min);
}

static void print_mapper(const parasail_matrix_t *matrix) {
    int i = 0;
    const int *mapper = matrix->mapper;

    for (i=0; i<256; ++i) {
        if (i%16 == 0) printf("\n");
        printf("%5d", mapper[i]);
    }
    printf("\n");
}

int main(int argc, char **argv)
{
    parasail_matrix_t *matrix = NULL;
    const parasail_matrix_t *internal_matrix = NULL;
    parasail_matrix_t *user_matrix = NULL;
    parasail_matrix_t *user_matrix_copy = NULL;
    const char *pssm_alphabet = "abcdef";
    const int pssm_values[] = {
        0,  1,  2,  3,   4,  5,
        6,  7,  8,  9,  10, 11,
        12, 13, 14, 15, 16, 17,
        18, 19, 20, 21, 22, 23
    };
    const int pssm_length = 4;

    if (argc == 2) {
        matrix = parasail_matrix_from_file(argv[1]);
        print_matrix(matrix);
        print_mapper(matrix);
        parasail_matrix_free(matrix);
    }
    else {
        printf("missing matrix file argument, skipping that part of the test\n");
    }

    print_matrix(&parasail_blosum62);
    print_mapper(&parasail_blosum62);

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
    print_matrix(user_matrix);
    print_mapper(user_matrix);
    parasail_matrix_free(user_matrix);

    user_matrix = parasail_matrix_create_case_sensitive("AcgT", 2, -1);
    if (NULL == user_matrix) {
        fprintf(stderr, "matrix create failed");
        exit(EXIT_FAILURE);
    }
    print_matrix(user_matrix);
    print_mapper(user_matrix);
    parasail_matrix_free(user_matrix);

    user_matrix = parasail_matrix_create_case_sensitive("ACGTacgt", 2, -1);
    if (NULL == user_matrix) {
        fprintf(stderr, "matrix create failed");
        exit(EXIT_FAILURE);
    }
    print_matrix(user_matrix);
    print_mapper(user_matrix);
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

    user_matrix = parasail_matrix_pssm_create(
            pssm_alphabet,
            pssm_values,
            pssm_length);
    if (NULL == user_matrix) {
        fprintf(stderr, "pssm matrix create failed");
        exit(EXIT_FAILURE);
    }
    print_matrix(user_matrix);
    print_mapper(user_matrix);

    user_matrix_copy = parasail_matrix_copy(user_matrix);
    if (NULL == user_matrix_copy) {
        fprintf(stderr, "pssm matrix copy failed");
        exit(EXIT_FAILURE);
    }
    parasail_matrix_free(user_matrix);

    print_matrix(user_matrix_copy);
    print_mapper(user_matrix_copy);
    parasail_matrix_free(user_matrix_copy);

    user_matrix = parasail_matrix_convert_square_to_pssm(internal_matrix, "asdf", 4);
    print_matrix(user_matrix);
    print_mapper(user_matrix);
    parasail_matrix_free(user_matrix);

    print_matrix(internal_matrix);

    return 0;
}

