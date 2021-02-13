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
#define ENAME parasail_sw_table_scan
#else
#ifdef PARASAIL_ROWCOL
#define ENAME parasail_sw_rowcol_scan
#else
#define ENAME parasail_sw_scan
#endif
#endif

parasail_result_t* ENAME(
        const char * const restrict _s1, const int _s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    /* declare local variables */
    parasail_result_t *result = NULL;
    int * restrict s1 = NULL;
    int * restrict s2 = NULL;
    int * restrict HB = NULL;
    int * restrict H  = NULL;
    int * restrict E  = NULL;
    int * restrict HtB= NULL;
    int * restrict Ht = NULL;
    int * restrict FtB= NULL;
    int * restrict Ft = NULL;
    int s1Len = 0;
    int i = 0;
    int j = 0;
    int score = 0;
    int end_query = 0;
    int end_ref = 0;

    /* validate inputs */
    PARASAIL_CHECK_NULL(_s2);
    PARASAIL_CHECK_GT0(s2Len);
    PARASAIL_CHECK_GE0(open);
    PARASAIL_CHECK_GE0(gap);
    PARASAIL_CHECK_NULL(matrix);
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        PARASAIL_CHECK_NULL(_s1);
        PARASAIL_CHECK_GT0(_s1Len);
    }

    /* initialize stack variables */
    s1Len = matrix->type == PARASAIL_MATRIX_TYPE_SQUARE ? _s1Len : matrix->length;
    i = 0;
    j = 0;
    score = NEG_INF_32;
    end_query = s1Len;
    end_ref = s2Len;

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
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_NOVEC_SCAN
        | PARASAIL_FLAG_BITS_INT | PARASAIL_FLAG_LANES_1;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    /* initialize heap variables */
    s2 = parasail_memalign_int(16, s2Len);
    HB = parasail_memalign_int(16, s1Len+1);
    H  = HB+1;
    E  = parasail_memalign_int(16, s1Len);
    HtB= parasail_memalign_int(16, s1Len+1);
    Ht = HtB+1;
    FtB= parasail_memalign_int(16, s1Len+1);
    Ft = FtB+1;

    /* validate heap variables */
    if (!s2) return NULL;
    if (!HB) return NULL;
    if (!E) return NULL;
    if (!HtB) return NULL;
    if (!FtB) return NULL;

    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        s1 = parasail_memalign_int(16, s1Len);
        if (!s1) return NULL;
        for (i=0; i<s1Len; ++i) {
            s1[i] = matrix->mapper[(unsigned char)_s1[i]];
        }
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
    for (j=0; j<s2Len; ++j) {
        /* calculate E */
        for (i=0; i<s1Len; ++i) {
            E[i] = MAX(E[i]-gap, H[i]-open);
        }
        /* calculate Ht */
        for (i=0; i<s1Len; ++i) {
            int matval = matrix->type == PARASAIL_MATRIX_TYPE_SQUARE ?
                         matrix->matrix[matrix->size*s1[i]+s2[j]] :
                         matrix->matrix[matrix->size*i+s2[j]];
            Ht[i] = MAX(H[i-1]+matval, E[i]);
        }
        /* calculate Ft */
        for (i=0; i<s1Len; ++i) {
            Ft[i] = MAX(Ft[i-1]-gap, Ht[i-1]);
        }
        /* calculate H */
        for (i=0; i<s1Len; ++i) {
            H[i] = MAX(Ht[i], Ft[i]-open);
            H[i] = MAX(H[i], 0);
#ifdef PARASAIL_TABLE
            result->tables->score_table[i*s2Len + j] = H[i];
#endif
            if (H[i] > score) {
                score = H[i];
                end_query = i;
                end_ref = j;
            }
        }
#ifdef PARASAIL_ROWCOL
        if (j == s2Len-1) {
            for (i=0; i<s1Len; ++i) {
                result->rowcols->score_col[i] = H[i];
            }
        }
        result->rowcols->score_row[j] = H[s1Len-1];
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
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        parasail_free(s1);
    }

    return result;
}

