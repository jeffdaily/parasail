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

#define ENAME parasail_nw_trace_scan

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
    int * restrict H = NULL;
    int * restrict E = NULL;
    int * restrict HtB = NULL;
    int * restrict Ht = NULL;
    int8_t * restrict HT = NULL;
    int * restrict Ex = NULL;
    int i = 0;
    int j = 0;
    int s1Len = 0;

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

    /* initialize result */
    result = parasail_result_new_trace(s1Len, s2Len, 16, sizeof(int8_t));
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_NOVEC_SCAN
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_INT | PARASAIL_FLAG_LANES_1;

    /* initialize heap variables */
    s2 = parasail_memalign_int(16, s2Len);
    HB = parasail_memalign_int(16, s1Len+1);
    H  = HB+1;
    E  = parasail_memalign_int(16, s1Len);
    HtB= parasail_memalign_int(16, s1Len+1);
    Ht = HtB+1;
    HT = (int8_t* const restrict)result->trace->trace_table;
    Ex = parasail_memalign_int(16, s1Len);

    /* validate heap variables */
    if (!s2) return NULL;
    if (!HB) return NULL;
    if (!E) return NULL;
    if (!HtB) return NULL;
    if (!HT) return NULL;
    if (!Ex) return NULL;

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
        H[i] = -open - i*gap;
    }

    /* initialize E */
    for (i=0; i<s1Len; ++i) {
        E[i] = NEG_INF_32;
    }

    /* iterate over database */
    for (j=0; j<s2Len; ++j) {
        int Ft = NEG_INF_32;
        /* calculate E */
        for (i=0; i<s1Len; ++i) {
            int E_opn = H[i]-open;
            int E_ext = E[i]-gap;
            E[i] = MAX(E_ext, E_opn);
            HT[1LL*i*s2Len + j] = (E_opn > E_ext) ? PARASAIL_DIAG_E
                                                  : PARASAIL_INS_E;
        }
        /* calculate Ht */
        for (i=0; i<s1Len; ++i) {
            int matval = matrix->type == PARASAIL_MATRIX_TYPE_SQUARE ?
                         matrix->matrix[matrix->size*s1[i]+s2[j]] :
                         matrix->matrix[matrix->size*i+s2[j]];
            int H_dag = H[i-1]+matval;
            Ht[i] = MAX(H_dag, E[i]);
            Ex[i] = (E[i] > H_dag);
        }
        Ht[-1] = -open -j*gap;
        /* calculate H */
        for (i=0; i<s1Len; ++i) {
            int Ft_opn;
            int Ht_pre = Ht[i-1];
            int Ft_ext = Ft-gap;
            if (Ht_pre >= Ft_ext) {
                Ft = Ht_pre;
            }
            else {
                Ft = Ft_ext;
            }
            Ft_opn = Ft-open;
            if (H[i-1] > Ft_ext) {
                HT[1LL*i*s2Len + j] |= PARASAIL_DIAG_F;
            }
            else {
                HT[1LL*i*s2Len + j] |= PARASAIL_DEL_F;
            }
            if (Ht[i] > Ft_opn) {
                H[i] = Ht[i];
                HT[1LL*i*s2Len + j] |= Ex[i] ? PARASAIL_INS : PARASAIL_DIAG;
            }
            else {
                H[i] = Ft_opn;
                if (Ht[i] == Ft_opn) {
                    if (Ex[i]) {
                        HT[1LL*i*s2Len + j] |= PARASAIL_DEL;
                    }
                    else {
                        HT[1LL*i*s2Len + j] |= PARASAIL_DIAG;
                    }
                }
                else {
                    HT[1LL*i*s2Len + j] |= PARASAIL_DEL;
                }
            }
        }
        H[-1] = -open - j*gap;
    }

    result->score = H[s1Len-1];
    result->end_query = s1Len-1;
    result->end_ref = s2Len-1;

    parasail_free(Ex);
    parasail_free(HtB);
    parasail_free(E);
    parasail_free(HB);
    parasail_free(s2);
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        parasail_free(s1);
    }

    return result;
}

