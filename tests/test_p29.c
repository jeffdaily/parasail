#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parasail.h"

int main()
{
    parasail_result_t* result = NULL;
    parasail_matrix_t* user_matrix = NULL;
    parasail_cigar_t* cigar = NULL;
    int match = 0;
    int mismatch = 0;
    int gap_opening = 0;
    int gap_extension = 0;
    const char *s1 = NULL;
    const char *s2 = NULL;
    int s1_len = 0;
    int s2_len = 0;
    char * decoded = NULL;

    //match = 1;
    //mismatch = -5;
    match = 5;
    mismatch = -4;
    gap_opening = -5;
    gap_extension = -5;
    user_matrix = parasail_matrix_create("ACGT", match, mismatch);

    s1 ="CAGAACAGGACCACGACACTC";
    s2 ="ACACAGAAAAGGAACCCGCACTCCATC";
    s1_len = (int)strlen(s1);
    s2_len = (int)strlen(s2);

    /* semi-global */

    result = parasail_sg_trace_scan_sat(s1, s1_len, s2, s2_len, -gap_opening, -gap_extension, user_matrix);
    if (parasail_result_is_saturated(result)) {
        printf("the result saturated\n");
        parasail_result_free(result);
        return EXIT_FAILURE;
    }

    cigar = parasail_result_get_cigar(result, s1, s1_len, s2, s2_len, user_matrix);
    decoded = parasail_cigar_decode(cigar);
    printf("('%s', %d)\n", decoded, parasail_result_get_score(result));
    parasail_traceback_generic(s1, s1_len, s2, s2_len, "query", "subject",
            user_matrix, result, '|', ':', '.', 60, 8, 0);

    free(decoded);
    parasail_cigar_free(cigar);
    parasail_result_free(result);

    /* global */

    result = parasail_nw_trace_scan_sat(s1, s1_len, s2, s2_len, -gap_opening, -gap_extension, user_matrix);
    if (parasail_result_is_saturated(result)) {
        printf("the result saturated\n");
        parasail_result_free(result);
        return EXIT_FAILURE;
    }

    cigar = parasail_result_get_cigar(result, s1, s1_len, s2, s2_len, user_matrix);
    decoded = parasail_cigar_decode(cigar);
    printf("('%s', %d)\n", decoded, parasail_result_get_score(result));
    parasail_traceback_generic(s1, s1_len, s2, s2_len, "query", "subject",
            user_matrix, result, '|', ':', '.', 60, 8, 0);

    free(decoded);
    parasail_cigar_free(cigar);
    parasail_result_free(result);

    parasail_matrix_free(user_matrix);

    return EXIT_SUCCESS;
}

