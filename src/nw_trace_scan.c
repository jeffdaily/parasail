/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include "parasail.h"
#include "parasail/memory.h"

#define NEG_INF_32 (INT32_MIN/2)
#define MAX(a,b) ((a)>(b)?(a):(b))

#define ENAME parasail_nw_trace_scan

parasail_result_t* ENAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_result_t *result = parasail_result_new_trace(s1Len, s2Len, sizeof(int));
    int * const restrict s1 = parasail_memalign_int(16, s1Len);
    int * const restrict s2 = parasail_memalign_int(16, s2Len);
    int * const restrict HB = parasail_memalign_int(16, s1Len+1);
    int * const restrict H  = HB+1;
    int * const restrict E  = parasail_memalign_int(16, s1Len);
    int * const restrict HtB= parasail_memalign_int(16, s1Len+1);
    int * const restrict Ht = HtB+1;
    int * const restrict FtB= parasail_memalign_int(16, s1Len+1);
    int * const restrict Ft = FtB+1;
    int i = 0;
    int j = 0;

    for (i=0; i<s1Len; ++i) {
        s1[i] = matrix->mapper[(unsigned char)_s1[i]];
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = matrix->mapper[(unsigned char)_s2[j]];
    }

    /* initialize H */
    H[-1] = 0;
    Ht[-1] = 0;
    for (i=0; i<s1Len; ++i) {
        H[i] = -open - i*gap;
    }

    /* initialize E */
    for (i=0; i<s1Len; ++i) {
        E[i] = NEG_INF_32;
    }

#if 1
    /* iterate over database */
    for (j=0; j<s2Len; ++j) {
        const int * const restrict matcol = &matrix->matrix[matrix->size*s2[j]];
        /* calculate E */
        for (i=0; i<s1Len; ++i) {
            E[i] = MAX(E[i]-gap, H[i]-open);
        }
        /* calculate Ht */
        for (i=0; i<s1Len; ++i) {
            Ht[i] = MAX(H[i-1]+matcol[s1[i]], E[i]);
        }
        Ht[-1] = -open -j*gap;
        Ft[-1] = NEG_INF_32;
        /* calculate Ft */
        for (i=0; i<s1Len; ++i) {
            Ft[i] = MAX(Ft[i-1]-gap, Ht[i-1]);
        }
        /* calculate H */
        for (i=0; i<s1Len; ++i) {
            H[i] = MAX(Ht[i], Ft[i]-open);
        }
        H[-1] = -open - j*gap;
    }
#else
    /* iterate over database */
    Ft[-1] = NEG_INF_32;
    for (j=0; j<s2Len; ++j) {
        const int * const restrict matcol = matrix->matrix_[s2[j]];
        int Hp = H[-1];
        for (i=0; i<s1Len; ++i) {
            E[i] = MAX(E[i]-gap, H[i]-open);
            Ht[i] = MAX(Hp+matcol[s1[i]], E[i]);
            Ft[i] = MAX(Ft[i-1]-gap, Ht[i-1]);
            Hp = H[i];
            H[i] = MAX(Ht[i], Ft[i]-open);
        }
        H[-1] = -open - j*gap;
        Ht[-1] = -open -j*gap;
    }
#endif

    result->score = H[s1Len-1];
    result->end_query = s1Len-1;
    result->end_ref = s2Len-1;
    result->flag = PARASAIL_FLAG_NW | PARASAIL_FLAG_NOVEC_SCAN
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_INT | PARASAIL_FLAG_LANES_1;

    parasail_free(FtB);
    parasail_free(HtB);
    parasail_free(E);
    parasail_free(HB);
    parasail_free(s2);
    parasail_free(s1);

    return result;
}

