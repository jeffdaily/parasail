/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>



#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_altivec.h"



#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        vec128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[1LL*(i+0)*s2Len + (j-0)] = (int32_t)_mm_extract_epi32(vWH, 3);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[1LL*(i+1)*s2Len + (j-1)] = (int32_t)_mm_extract_epi32(vWH, 2);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[1LL*(i+2)*s2Len + (j-2)] = (int32_t)_mm_extract_epi32(vWH, 1);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[1LL*(i+3)*s2Len + (j-3)] = (int32_t)_mm_extract_epi32(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_rowcol(
        int *row,
        int *col,
        vec128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (i+0 == s1Len-1 && 0 <= j-0 && j-0 < s2Len) {
        row[j-0] = (int32_t)_mm_extract_epi32(vWH, 3);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int32_t)_mm_extract_epi32(vWH, 3);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int32_t)_mm_extract_epi32(vWH, 2);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int32_t)_mm_extract_epi32(vWH, 2);
    }
    if (i+2 == s1Len-1 && 0 <= j-2 && j-2 < s2Len) {
        row[j-2] = (int32_t)_mm_extract_epi32(vWH, 1);
    }
    if (j-2 == s2Len-1 && 0 <= i+2 && i+2 < s1Len) {
        col[(i+2)] = (int32_t)_mm_extract_epi32(vWH, 1);
    }
    if (i+3 == s1Len-1 && 0 <= j-3 && j-3 < s2Len) {
        row[j-3] = (int32_t)_mm_extract_epi32(vWH, 0);
    }
    if (j-3 == s2Len-1 && 0 <= i+3 && i+3 < s1Len) {
        col[(i+3)] = (int32_t)_mm_extract_epi32(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_table_diag_altivec_128_32
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_rowcol_diag_altivec_128_32
#else
#define FNAME parasail_sw_diag_altivec_128_32
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int _s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    /* declare local variables */
    int32_t N = 0;
    int32_t PAD = 0;
    int32_t PAD2 = 0;
    int32_t s1Len = 0;
    int32_t s1Len_PAD = 0;
    int32_t s2Len_PAD = 0;
    int32_t * restrict s1 = NULL;
    int32_t * restrict s2B = NULL;
    int32_t * restrict _H_pr = NULL;
    int32_t * restrict _F_pr = NULL;
    int32_t * restrict s2 = NULL;
    int32_t * restrict H_pr = NULL;
    int32_t * restrict F_pr = NULL;
    parasail_result_t *result = NULL;
    int32_t i = 0;
    int32_t j = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int32_t NEG_LIMIT = 0;
    int32_t POS_LIMIT = 0;
    int32_t score = 0;
    vec128i vNegLimit;
    vec128i vPosLimit;
    vec128i vSaturationCheckMin;
    vec128i vSaturationCheckMax;
    vec128i vNegInf;
    vec128i vNegInf0;
    vec128i vOpen;
    vec128i vGap;
    vec128i vZero;
    vec128i vOne;
    vec128i vN;
    vec128i vNegOne;
    vec128i vI;
    vec128i vJreset;
    vec128i vMaxH;
    vec128i vEndI;
    vec128i vEndJ;
    vec128i vILimit;
    vec128i vJLimit;

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
    N = 4; /* number of values in vector */
    PAD = N-1;
    PAD2 = PAD*2;
    s1Len = matrix->type == PARASAIL_MATRIX_TYPE_SQUARE ? _s1Len : matrix->length;
    s1Len_PAD = s1Len+PAD;
    s2Len_PAD = s2Len+PAD;
    i = 0;
    j = 0;
    end_query = 0;
    end_ref = 0;
    NEG_LIMIT = (-open < matrix->min ? INT32_MIN + open : INT32_MIN - matrix->min) + 1;
    POS_LIMIT = INT32_MAX - matrix->max - 1;
    score = NEG_LIMIT;
    vNegLimit = _mm_set1_epi32(NEG_LIMIT);
    vPosLimit = _mm_set1_epi32(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;
    vNegInf = _mm_set1_epi32(NEG_LIMIT);
    vNegInf0 = _mm_srli_si128(vNegInf, 4); /* shift in a 0 */
    vOpen = _mm_set1_epi32(open);
    vGap  = _mm_set1_epi32(gap);
    vZero = _mm_set1_epi32(0);
    vOne = _mm_set1_epi32(1);
    vN = _mm_set1_epi32(N);
    vNegOne = _mm_set1_epi32(-1);
    vI = _mm_set_epi32(0,1,2,3);
    vJreset = _mm_set_epi32(0,-1,-2,-3);
    vMaxH = vNegInf;
    vEndI = vNegInf;
    vEndJ = vNegInf;
    vILimit = _mm_set1_epi32(s1Len);
    vJLimit = _mm_set1_epi32(s2Len);

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
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_BITS_32 | PARASAIL_FLAG_LANES_4;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    /* initialize heap variables */
    s2B= parasail_memalign_int32_t(16, s2Len+PAD2);
    _H_pr = parasail_memalign_int32_t(16, s2Len+PAD2);
    _F_pr = parasail_memalign_int32_t(16, s2Len+PAD2);
    s2 = s2B+PAD; /* will allow later for negative indices */
    H_pr = _H_pr+PAD;
    F_pr = _F_pr+PAD;

    /* validate heap variables */
    if (!s2B) return NULL;
    if (!_H_pr) return NULL;
    if (!_F_pr) return NULL;

    /* convert _s1 from char to int in range 0-23 */
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        s1 = parasail_memalign_int32_t(16, s1Len+PAD);
        if (!s1) return NULL;
        for (i=0; i<s1Len; ++i) {
            s1[i] = matrix->mapper[(unsigned char)_s1[i]];
        }
        /* pad back of s1 with dummy values */
        for (i=s1Len; i<s1Len_PAD; ++i) {
            s1[i] = 0; /* point to first matrix row because we don't care */
        }
    }

    /* convert _s2 from char to int in range 0-23 */
    for (j=0; j<s2Len; ++j) {
        s2[j] = matrix->mapper[(unsigned char)_s2[j]];
    }
    /* pad front of s2 with dummy values */
    for (j=-PAD; j<0; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }
    /* pad back of s2 with dummy values */
    for (j=s2Len; j<s2Len_PAD; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }

    /* set initial values for stored row */
    for (j=0; j<s2Len; ++j) {
        H_pr[j] = 0;
        F_pr[j] = NEG_LIMIT;
    }
    /* pad front of stored row values */
    for (j=-PAD; j<0; ++j) {
        H_pr[j] = NEG_LIMIT;
        F_pr[j] = NEG_LIMIT;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        H_pr[j] = NEG_LIMIT;
        F_pr[j] = NEG_LIMIT;
    }

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        vec128i vNH = vNegInf0;
        vec128i vWH = vNegInf0;
        vec128i vE = vNegInf;
        vec128i vF = vNegInf;
        vec128i vJ = vJreset;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+0] : ((i+0 >= s1Len) ? s1Len-1 : i+0))];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+1] : ((i+1 >= s1Len) ? s1Len-1 : i+1))];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+2] : ((i+2 >= s1Len) ? s1Len-1 : i+2))];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+3] : ((i+3 >= s1Len) ? s1Len-1 : i+3))];
        vec128i vIltLimit = _mm_cmplt_epi32(vI, vILimit);
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            vec128i vMat;
            vec128i vNWH = vNH;
            vNH = _mm_srli_si128(vWH, 4);
            vNH = _mm_insert_epi32(vNH, H_pr[j], 3);
            vF = _mm_srli_si128(vF, 4);
            vF = _mm_insert_epi32(vF, F_pr[j], 3);
            vF = _mm_max_epi32(
                    _mm_sub_epi32(vNH, vOpen),
                    _mm_sub_epi32(vF, vGap));
            vE = _mm_max_epi32(
                    _mm_sub_epi32(vWH, vOpen),
                    _mm_sub_epi32(vE, vGap));
            vMat = _mm_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]]
                    );
            vNWH = _mm_add_epi32(vNWH, vMat);
            vWH = _mm_max_epi32(vNWH, vE);
            vWH = _mm_max_epi32(vWH, vF);
            vWH = _mm_max_epi32(vWH, vZero);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                vec128i cond = _mm_cmpeq_epi32(vJ,vNegOne);
                vWH = _mm_andnot_si128(cond, vWH);
                vF = _mm_blendv_epi8(vF, vNegInf, cond);
                vE = _mm_blendv_epi8(vE, vNegInf, cond);
            }
            /* cannot start checking sat until after J clears boundary */
            if (j > PAD) {
                vSaturationCheckMin = _mm_min_epi32(vSaturationCheckMin, vWH);
                vSaturationCheckMax = _mm_max_epi32(vSaturationCheckMax, vWH);
            }
#ifdef PARASAIL_TABLE
            arr_store_si128(result->tables->score_table, vWH, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->rowcols->score_row, result->rowcols->score_col, vWH, i, s1Len, j, s2Len);
#endif
            H_pr[j-3] = (int32_t)_mm_extract_epi32(vWH,0);
            F_pr[j-3] = (int32_t)_mm_extract_epi32(vF,0);
            /* as minor diagonal vector passes across table, extract
             * max values within the i,j bounds */
            {
                vec128i cond_valid_J = _mm_and_si128(
                        _mm_cmpgt_epi32(vJ, vNegOne),
                        _mm_cmplt_epi32(vJ, vJLimit));
                vec128i cond_valid_IJ = _mm_and_si128(cond_valid_J, vIltLimit);
                vec128i cond_eq = _mm_cmpeq_epi32(vWH, vMaxH);
                vec128i cond_max = _mm_cmpgt_epi32(vWH, vMaxH);
                vec128i cond_all = _mm_and_si128(cond_max, cond_valid_IJ);
                vec128i cond_Jlt = _mm_cmplt_epi32(vJ, vEndJ);
                vMaxH = _mm_blendv_epi8(vMaxH, vWH, cond_all);
                vEndI = _mm_blendv_epi8(vEndI, vI, cond_all);
                vEndJ = _mm_blendv_epi8(vEndJ, vJ, cond_all);
                cond_all = _mm_and_si128(cond_Jlt, cond_eq);
                cond_all = _mm_and_si128(cond_all, cond_valid_IJ);
                vEndI = _mm_blendv_epi8(vEndI, vI, cond_all);
                vEndJ = _mm_blendv_epi8(vEndJ, vJ, cond_all);
            }
            vJ = _mm_add_epi32(vJ, vOne);
        }
        vI = _mm_add_epi32(vI, vN);
    }

    /* alignment ending position */
    {
        int32_t *t = (int32_t*)&vMaxH;
        int32_t *i = (int32_t*)&vEndI;
        int32_t *j = (int32_t*)&vEndJ;
        int32_t k;
        for (k=0; k<N; ++k, ++t, ++i, ++j) {
            if (*t > score) {
                score = *t;
                end_query = *i;
                end_ref = *j;
            }
            else if (*t == score) {
                if (*j < end_ref) {
                    end_query = *i;
                    end_ref = *j;
                }
                else if (*j == end_ref && *i < end_query) {
                    end_query = *i;
                    end_ref = *j;
                }
            }
        }
    }

    if (_mm_movemask_epi8(_mm_or_si128(
            _mm_cmplt_epi32(vSaturationCheckMin, vNegLimit),
            _mm_cmpgt_epi32(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(_F_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        parasail_free(s1);
    }

    return result;
}


