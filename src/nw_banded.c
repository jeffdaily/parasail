/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
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
#define ENAME parasail_nw_banded_table
#else
#ifdef PARASAIL_ROWCOL
#define ENAME parasail_nw_banded_rowcol
#else
#define ENAME parasail_nw_banded
#endif
#endif

parasail_result_t* ENAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const int k,
        const parasail_matrix_t *matrix)
{
    /* declare local variables */
    parasail_result_t *result = NULL;
    int * restrict s1 = NULL;
    int * restrict s2 = NULL;
    int * HB = NULL;
    int * H = NULL;
    int * EB = NULL;
    int * E = NULL;
    int i = 0;
    int j = 0;
    int low = 0;
    int colLen = 0;
    int colOff = 0;

    /* validate inputs */
    PARASAIL_CHECK_NULL(_s1);
    PARASAIL_CHECK_GT0(s1Len);
    PARASAIL_CHECK_NULL(_s2);
    PARASAIL_CHECK_GT0(s2Len);
    PARASAIL_CHECK_GE0(open);
    PARASAIL_CHECK_GE0(gap);
    PARASAIL_CHECK_GT0(k);
    PARASAIL_CHECK_NULL(matrix);

    /* initialize stack variables */

    /* initialize result */
#ifdef PARASAIL_TABLE
    result = parasail_result_new_table1(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    result = parasail_result_new_rowcol1(s1Len, s2Len);
#else
    result = parasail_result_new();
#endif
#endif
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_NOVEC
        | PARASAIL_FLAG_BANDED
        | PARASAIL_FLAG_BITS_INT | PARASAIL_FLAG_LANES_1;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    /* initalize heap variables */
    s1 = parasail_memalign_int(16, s1Len);
    s2 = parasail_memalign_int(16, s2Len);

    /* validate heap variables */
    if (!s1) return NULL;
    if (!s2) return NULL;

    low = s1Len - s2Len - k;
    if (s1Len > s2Len) {
        low = s2Len - s1Len - k;
    }
    colLen = k-low+1;
    colOff = s1Len > s2Len ? -k : low;

    HB = parasail_memalign_int(16, colLen+2);
    if (!HB) return NULL;
    H = HB+1;
    EB = parasail_memalign_int(16, colLen+2);
    if (!EB) return NULL;
    E = EB+1;
    parasail_memset_int(HB, 0, colLen+2);
    parasail_memset_int(EB, 0, colLen+2);

    for (i=0; i<s1Len; ++i) {
        s1[i] = matrix->mapper[(unsigned char)_s1[i]];
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = matrix->mapper[(unsigned char)_s2[j]];
    }

#ifdef PARASAIL_TABLE
    /* fill the table with zeros */
    {
        int limit = s1Len*s2Len;
        for (i=0; i<limit; ++i) {
            result->tables->score_table[i] = 0;
        }
    }
#endif

    /* initialize H */
    /* initialize E */
    for (i=-colOff+1,j=0; i<colLen; ++i, ++j) {
        H[i] = -open - j*gap;
        E[i] = NEG_INF_32;
    }
    H[-colOff-1] = -open;
    E[-colOff-1] = NEG_INF_32;
    H[-colOff] = 0;
    E[-colOff] = NEG_INF_32;
    H[-1] = NEG_INF_32;
    E[-1] = NEG_INF_32;
    H[colLen] = NEG_INF_32;
    E[colLen] = NEG_INF_32;

    /* iter over db */
    for (j=0; j<s2Len; ++j) {
        /* substitution matrix row (really we wanted a column) */
        int F = NEG_INF_32;
        if (colOff < 0) {
            H[-colOff-1] = -open - j*gap;
        }
        /* iter over query */
        for (i=0; i<colLen; ++i) {
            int pos = 0;
            int H_dag = 0;
            int E_opn = 0;
            int E_ext = 0;
            int F_opn = 0;
            int F_ext = 0;
            int matval = 0;
            pos = colOff+i;
            if (pos < 0 || pos >= s1Len) {
                continue;
            }
            matval = matrix->type == PARASAIL_MATRIX_TYPE_SQUARE ?
                     matrix->matrix[matrix->size*s1[pos]+s2[j]] :
                     matrix->matrix[matrix->size*pos+s2[j]];
            H_dag = H[i] + matval;
            E_opn = H[i+1] - open;
            E_ext = E[i+1] - gap;
            F_opn = H[i-1] - open;
            F_ext = F - gap;
            E[i] = MAX(E_opn, E_ext);
            F = MAX(F_opn, F_ext);
            H[i] = MAX(E[i], F);
            H[i] = MAX(H[i], H_dag);
#ifdef PARASAIL_TABLE
            result->tables->score_table[pos*s2Len + j] = H[i];
#endif
        }
        colOff += 1;
    }

    result->score = s1Len > s2Len ? H[-low] : H[k];
    result->end_query = s1Len-1;
    result->end_ref = s2Len-1;

    parasail_free(EB);
    parasail_free(HB);
    parasail_free(s2);
    parasail_free(s1);

    return result;
}

