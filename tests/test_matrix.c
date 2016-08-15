#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "parasail.h"
#include "parasail/io.h"
#include "parasail/matrices/blosum62.h"

static char* get_alphabet(const parasail_matrix_t *matrix) {
    int i = 0;
    char *alphabet = NULL;

    alphabet = (char*)malloc(sizeof(char)*(matrix->size+1));
    for (i=0; i<matrix->size; ++i) {
        alphabet[i] = '*';
    }
    alphabet[matrix->size+1] = '\0';

    for (i=65; i<91; ++i) {
        if (matrix->mapper[i] < matrix->size) {
            alphabet[matrix->mapper[i]] = i;
        }
    }

    return alphabet;
}

static void print_matrix(const parasail_matrix_t *matrix) {
    int i = 0;
    int j = 0;
    char *alphabet = NULL;

    alphabet = get_alphabet(matrix);

    printf("matrix '%s'\n", matrix->name);
    printf("  ");
    for (i=0; i<matrix->size; ++i) {
        printf("%4c", alphabet[i]);
    }
    printf("\n");
    for (j=0; j<matrix->size; ++j) {
        printf("%c ", alphabet[j]);
        for (i=0; i<matrix->size; ++i) {
            printf("%4d", matrix->matrix[j*matrix->size + i]);
        }
        printf("\n");
    }
    printf("max = %d\n", matrix->max);
    printf("min = %d\n", matrix->min);

    free(alphabet);
}

int main(int argc, char **argv)
{
    parasail_matrix_t *matrix = NULL;

    if (argc == 2) {
        matrix = parasail_matrix_from_file(argv[1]);
    }
    else {
        printf("missing matrix file argument\n");
        return -1;
    }

    print_matrix(matrix);
    parasail_matrix_free(matrix);

    print_matrix(&parasail_blosum62);

    return 0;
}

