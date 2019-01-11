#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parasail.h"
#include "parasail/matrices/nuc44.h"

int main()
{
    parasail_result_t* result;
    
    const char *s = "AAAAAAAAAATTTAAAAAAAAAA";
    const char *t = "AAAAAAAAAAAAAAAAAAAA";

    /* SSE4 alignment */
    result = parasail_sg_trace_striped_sse41_128_16(
        s, strlen(s), t, strlen(t), 1, 4, &parasail_nuc44);

    parasail_traceback_generic(
        s, strlen(s), t, strlen(t),
        "A", "B", &parasail_nuc44,
        result, '|', '+', '-', 79, 10, 1);

    /* AVX2 alignment */
    result = parasail_sg_trace_striped_avx2_256_16(
        s, strlen(s), t, strlen(t), 1, 4, &parasail_nuc44);

    parasail_traceback_generic(
        s, strlen(s), t, strlen(t),
        "A", "B", &parasail_nuc44,
        result, '|', '+', '-', 79, 10, 1);
    
    return 0;
}

