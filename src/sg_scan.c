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

#ifdef PARASAIL_TABLE
#define ENAME parasail_sg_table_scan
#else
#ifdef PARASAIL_ROWCOL
#define ENAME parasail_sg_rowcol_scan
#else
#define ENAME parasail_sg_scan
#endif
#endif

parasail_result_t* ENAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif
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
    int score = NEG_INF_32;
    int end_query = s1Len;
    int end_ref = s2Len;

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
        H[i] = 0;
    }
    Ft[-1] = NEG_INF_32;

    /* initialize E */
    for (i=0; i<s1Len; ++i) {
        E[i] = NEG_INF_32;
    }

    /* iterate over database */
    for (j=0; j<s2Len-1; ++j) {
        const int * const restrict matcol = &matrix->matrix[matrix->size*s2[j]];
        /* calculate E */
        for (i=0; i<s1Len; ++i) {
            E[i] = MAX(E[i]-gap, H[i]-open);
        }
        /* calculate Ht */
        for (i=0; i<s1Len; ++i) {
            Ht[i] = MAX(H[i-1]+matcol[s1[i]], E[i]);
        }
        /* calculate Ft */
        for (i=0; i<s1Len; ++i) {
            Ft[i] = MAX(Ft[i-1]-gap, Ht[i-1]);
        }
        /* calculate H */
        for (i=0; i<s1Len; ++i) {
            H[i] = MAX(Ht[i], Ft[i]-open);
#ifdef PARASAIL_TABLE
            result->score_table[i*s2Len + j] = H[i];
#endif
        }
#ifdef PARASAIL_ROWCOL
        result->score_row[j] = H[s1Len-1];
#endif
        if (H[s1Len-1] > score) {
            score = H[s1Len-1];
            end_query = s1Len-1;
            end_ref = j;
        }
    }
    j = s2Len - 1;
    {
        const int * const restrict matcol = &matrix->matrix[matrix->size*s2[j]];
        /* calculate E */
        for (i=0; i<s1Len; ++i) {
            E[i] = MAX(E[i]-gap, H[i]-open);
        }
        /* calculate Ht */
        for (i=0; i<s1Len; ++i) {
            Ht[i] = MAX(H[i-1]+matcol[s1[i]], E[i]);
        }
        /* calculate Ft */
        for (i=0; i<s1Len; ++i) {
            Ft[i] = MAX(Ft[i-1]-gap, Ht[i-1]);
        }
        /* calculate H */
        for (i=0; i<s1Len; ++i) {
            H[i] = MAX(Ht[i], Ft[i]-open);
#ifdef PARASAIL_TABLE
            result->score_table[i*s2Len + j] = H[i];
#endif
#ifdef PARASAIL_ROWCOL
            result->score_col[i] = H[i];
#endif
            if (H[i] > score) {
                score = H[i];
                end_query = i;
                end_ref = j;
            }
        }
#ifdef PARASAIL_ROWCOL
        result->score_row[j] = H[s1Len-1];
#endif
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(FtB);
    parasail_free(HtB);
    parasail_free(E);
    parasail_free(HB);
    parasail_free(s2);
    parasail_free(s1);

    return result;
}
