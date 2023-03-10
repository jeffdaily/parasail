#include <stdio.h>
#include <stdlib.h>

#include <parasail.h>
#include <parasail/matrices/blosum62.h>
#include <parasail/matrices/pam50.h>

int main(int argc, char ** argv)
{
    int major = 0;
    int minor = 0;
    int patch = 0;
    parasail_result_t *result = NULL;
    parasail_matrix_t *matrix = NULL;
    parasail_profile_t *profile = NULL;

    parasail_version(&major, &minor, &patch);
    printf("parasail is using C lib version %d.%d.%d\n", major, minor, patch);

    printf("\ntest1\n");
    printf("result = parasail_sw(\"asdf\", 4, \"asdf\", 4, 10, 1, &parasail_blosum62)\n");
    result = parasail_sw("asdf", 4, "asdf", 4, 10, 1, &parasail_blosum62);
    if (result->score != 20) {
        printf("failed\n");
        return EXIT_FAILURE;
    }
    else {
        printf("pass\n");
    }
    parasail_result_free(result);

    printf("\ntest2\n");
    printf("result = parasail_sw(\"asdf\", 4, \"asdf\", 4, 10, 1, &parasail_pam50)\n");
    result = parasail_sw("asdf", 4, "asdf", 4, 10, 1, &parasail_pam50);
    if (result->score != 27) {
        printf("failed\n");
        return EXIT_FAILURE;
    }
    else {
        printf("pass\n");
    }
    parasail_result_free(result);

    printf("\ntest3\n");
    printf("matrix = parasail_matrix_create(\"acgt\", 1, -1)\n");
    matrix = parasail_matrix_create("acgt", 1, -1);
    printf("result = parasail_sw(\"acgt\", 4, \"acgt\", 4, 10, 1, matrix)\n");
    result = parasail_sw("acgt", 4, "acgt", 4, 10, 1, matrix);
    if (result->score != 4) {
        printf("failed\n");
        return EXIT_FAILURE;
    }
    else {
        printf("pass\n");
    }
    parasail_result_free(result);
    parasail_matrix_free(matrix);

    /* if no vector ISA detected at build-time and/or run-time, this will fail */
    printf("\ntest4\n");
    printf("profile = parasail_profile_create_8(\"asdf\", 4, &parasail_blosum62)\n");
    profile = parasail_profile_create_8("asdf", 4, &parasail_blosum62);
    if (!profile) {
        printf("NULL profile returned; vector ISA not supported\n");
        printf("pass\n");
    }
    else {
        printf("result = parasail_sw_striped_profile_8(profile, \"asdf\", 4, 10, 1)\n");
        result = parasail_sw_striped_profile_8(profile, "asdf", 4, 10, 1);
        if (!result) {
            printf("NULL result returned; vector ISA not supported\n");
            printf("pass\n");
        }
        else {
            if (result->score != 20) {
                printf("failed\n");
                return EXIT_FAILURE;
            }
            else {
                printf("pass\n");
            }
            parasail_result_free(result);
        }
    }
    parasail_profile_free(profile);

    return 0;

}
