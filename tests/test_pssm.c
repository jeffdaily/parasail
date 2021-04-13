#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <parasail.h>
#include <parasail/matrices/blosum62.h>

int main(int argc, char ** argv)
{
    int width = 10;
    int height = 20;
    int r = 0;
    int c = 0;
    int i = 0;
    const char *alphabet = "abcdefghij";
    const char *seq = "deadbeef";
    int len = (int)strlen(seq);
    int *values = malloc(sizeof(int)*width*height);
    parasail_matrix_t *matrix = NULL;

    assert(strlen(alphabet) == (size_t)width);
    assert(values);

    printf("values = {\n");
    for (r=0; r<height; ++r) {
        for (c=0; c<width; ++c) {
            values[i++] = (rand()%30)-10;
            printf("%4d ", values[i-1]);
        }
        printf("\n");
    }
    printf("}\n");

    matrix = parasail_matrix_pssm_create(alphabet, values, height);

    /* done with values; copied into matrix */
    free(values);

    printf("matrix = {\n");
    for (r=0; r<matrix->length; ++r) {
        for (c=0; c<matrix->size; ++c) {
            printf("%4d ", matrix->matrix[r*matrix->size + c]);
        }
        printf("\n");
    }
    printf("}\n");
    printf("matrix min %d\n", matrix->min);
    printf("matrix max %d\n", matrix->max);
    printf("matrix size %d\n", matrix->size);
    printf("matrix length %d\n", matrix->length);

    parasail_matrix_free(matrix);

    matrix = parasail_matrix_convert_square_to_pssm(&parasail_blosum62, seq, len);

    printf("pssm_from_blosum62 = {\n");
    for (r=0; r<matrix->length; ++r) {
        for (c=0; c<matrix->size; ++c) {
            printf("%4d ", matrix->matrix[r*matrix->size + c]);
        }
        printf("\n");
    }
    printf("}\n");
    printf("matrix min %d\n", matrix->min);
    printf("matrix max %d\n", matrix->max);
    printf("matrix size %d\n", matrix->size);
    printf("matrix length %d\n", matrix->length);

    parasail_matrix_free(matrix);
    return 0;
}
