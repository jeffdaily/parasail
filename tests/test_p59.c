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
    int64_t s1Limit = 0;
    int64_t s2Limit = 0;
    char * decoded = NULL;

    //match = 1;
    //mismatch = -5;
    match = 1;
    mismatch = -3;
    gap_opening = -4;
    gap_extension = -1;
    user_matrix = parasail_matrix_create("ACGT", match, mismatch);

    // load sequence data, assigns s1 and s2
#include "test_p59.h"

    s1_len = (int)strlen(s1);
    s2_len = (int)strlen(s2);
    s1Limit = ((int64_t)s1_len)*gap_extension + gap_opening;
    s2Limit = ((int64_t)s2_len)*gap_extension + gap_opening;
    printf("s1_len=%d limit=%ld\n", s1_len, s1Limit);
    printf("s2_len=%d limit=%ld\n", s2_len, s2Limit);
    printf("INT16_MIN=%d\n", (int)INT16_MIN);

    /* semi-global striped */

    result = parasail_sg_trace_striped_16(s1, s1_len, s2, s2_len, -gap_opening, -gap_extension, user_matrix);
    //result = parasail_sg_trace_striped_32(s1, s1_len, s2, s2_len, -gap_opening, -gap_extension, user_matrix);
    if (parasail_result_is_saturated(result)) {
        printf("the result saturated\n");
        parasail_result_free(result);
        //return EXIT_FAILURE;
    }
    else {
        cigar = parasail_result_get_cigar(result, s1, s1_len, s2, s2_len, user_matrix);
        decoded = parasail_cigar_decode(cigar);
        printf("('%s', %d)\n", decoded, parasail_result_get_score(result));
        parasail_traceback_generic(s1, s1_len, s2, s2_len, "query", "subject",
                user_matrix, result, '|', ':', '.', 60, 8, 0);
        free(decoded);
        parasail_cigar_free(cigar);
        parasail_result_free(result);
    }

    /* semi-global scan */

#if 1
    result = parasail_sg_trace_scan_16(s1, s1_len, s2, s2_len, -gap_opening, -gap_extension, user_matrix);
    //result = parasail_sg_trace_scan_32(s1, s1_len, s2, s2_len, -gap_opening, -gap_extension, user_matrix);
    if (parasail_result_is_saturated(result)) {
        printf("the result saturated\n");
        parasail_result_free(result);
        //return EXIT_FAILURE;
    }
    else {
        cigar = parasail_result_get_cigar(result, s1, s1_len, s2, s2_len, user_matrix);
        decoded = parasail_cigar_decode(cigar);
        printf("('%s', %d)\n", decoded, parasail_result_get_score(result));
        parasail_traceback_generic(s1, s1_len, s2, s2_len, "query", "subject",
                user_matrix, result, '|', ':', '.', 60, 8, 0);
        free(decoded);
        parasail_cigar_free(cigar);
        parasail_result_free(result);
    }
#endif

    parasail_matrix_free(user_matrix);

    return EXIT_SUCCESS;
}

