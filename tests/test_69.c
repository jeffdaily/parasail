#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parasail.h"
#include "parasail/matrices/nuc44.h"

int main()
{
    parasail_result_t* result;
    
    const char* s =     "CCGTACGTACGTCCCCC"; // recombinant gene query
    const char* t = "GGGGCCGTACGTACGTGGGGG"; // reference database

    printf("\nsg alignment\n");
    result = parasail_sg_trace(s, strlen(s), t, strlen(t), 2, 2, &parasail_nuc44);
    parasail_traceback_generic(
        s, strlen(s), t, strlen(t),
        "Query:", "Target:", &parasail_nuc44,
        result, '|', '*', '*', 60, 7, 1);
	parasail_result_free(result);

    printf("\nsw alignment\n");
    result = parasail_sw_trace(s, strlen(s), t, strlen(t), 2, 2, &parasail_nuc44);
    parasail_traceback_generic(
        s, strlen(s), t, strlen(t),
        "Query:", "Target:", &parasail_nuc44,
        result, '|', '*', '*', 60, 7, 1);
	parasail_result_free(result);

    printf("\nsg_qe_db alignment\n");
    result = parasail_sg_qe_db_trace(s, strlen(s), t, strlen(t), 2, 2, &parasail_nuc44);
    parasail_traceback_generic(
        s, strlen(s), t, strlen(t),
        "Query:", "Target:", &parasail_nuc44,
        result, '|', '*', '*', 60, 7, 1);
	parasail_result_free(result);

    printf("\nsg_qe_dx alignment\n");
    result = parasail_sg_qe_dx_trace(s, strlen(s), t, strlen(t), 2, 2, &parasail_nuc44);
    parasail_traceback_generic(
        s, strlen(s), t, strlen(t),
        "Query:", "Target:", &parasail_nuc44,
        result, '|', '*', '*', 60, 7, 1);
	parasail_result_free(result);

    return 0;
}

