/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>

#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <emmintrin.h>
#include <smmintrin.h>
#endif

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_sse.h"

#define SG_TRACE
#define SG_SUFFIX _diag_sse41_128_16
#include "sg_helper.h"



static inline void arr_store_si128(
        int8_t *array,
        __m128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[1LL*(i+0)*s2Len + (j-0)] = (int8_t)_mm_extract_epi16(vWH, 7);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[1LL*(i+1)*s2Len + (j-1)] = (int8_t)_mm_extract_epi16(vWH, 6);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[1LL*(i+2)*s2Len + (j-2)] = (int8_t)_mm_extract_epi16(vWH, 5);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[1LL*(i+3)*s2Len + (j-3)] = (int8_t)_mm_extract_epi16(vWH, 4);
    }
    if (0 <= i+4 && i+4 < s1Len && 0 <= j-4 && j-4 < s2Len) {
        array[1LL*(i+4)*s2Len + (j-4)] = (int8_t)_mm_extract_epi16(vWH, 3);
    }
    if (0 <= i+5 && i+5 < s1Len && 0 <= j-5 && j-5 < s2Len) {
        array[1LL*(i+5)*s2Len + (j-5)] = (int8_t)_mm_extract_epi16(vWH, 2);
    }
    if (0 <= i+6 && i+6 < s1Len && 0 <= j-6 && j-6 < s2Len) {
        array[1LL*(i+6)*s2Len + (j-6)] = (int8_t)_mm_extract_epi16(vWH, 1);
    }
    if (0 <= i+7 && i+7 < s1Len && 0 <= j-7 && j-7 < s2Len) {
        array[1LL*(i+7)*s2Len + (j-7)] = (int8_t)_mm_extract_epi16(vWH, 0);
    }
}

#define FNAME parasail_sg_flags_trace_diag_sse41_128_16

parasail_result_t* FNAME(
        const char * const restrict _s1, const int _s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix,
        int s1_beg, int s1_end, int s2_beg, int s2_end)
{
    /* declare local variables */
    int32_t N = 0;
    int32_t PAD = 0;
    int32_t PAD2 = 0;
    int32_t s1Len_PAD = 0;
    int32_t s2Len_PAD = 0;
    int16_t * restrict s1 = NULL;
    int16_t * restrict s2B = NULL;
    int16_t * restrict _H_pr = NULL;
    int16_t * restrict _F_pr = NULL;
    int16_t * restrict s2 = NULL;
    int16_t * restrict H_pr = NULL;
    int16_t * restrict F_pr = NULL;
    parasail_result_t *result = NULL;
    int32_t i = 0;
    int32_t j = 0;
    int32_t s1Len = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int16_t NEG_LIMIT = 0;
    int16_t POS_LIMIT = 0;
    int16_t score = 0;
    __m128i vNegLimit;
    __m128i vPosLimit;
    __m128i vSaturationCheckMin;
    __m128i vSaturationCheckMax;
    __m128i vNegInf;
    __m128i vOpen;
    __m128i vGap;
    __m128i vOne;
    __m128i vN;
    __m128i vGapN;
    __m128i vNegOne;
    __m128i vI;
    __m128i vJreset;
    __m128i vMaxHRow;
    __m128i vMaxHCol;
    __m128i vLastVal;
    __m128i vEndI;
    __m128i vEndJ;
    __m128i vILimit;
    __m128i vILimit1;
    __m128i vJLimit;
    __m128i vJLimit1;
    __m128i vIBoundary;
    __m128i vTDiag;
    __m128i vTIns;
    __m128i vTDel;
    __m128i vTDiagE;
    __m128i vTInsE;
    __m128i vTDiagF;
    __m128i vTDelF;

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
    N = 8; /* number of values in vector */
    PAD = N-1;
    PAD2 = PAD*2;
    s1Len = matrix->type == PARASAIL_MATRIX_TYPE_SQUARE ? _s1Len : matrix->length;
    s1Len_PAD = s1Len+PAD;
    s2Len_PAD = s2Len+PAD;
    i = 0;
    j = 0;
    end_query = s1Len-1;
    end_ref = s2Len-1;
    NEG_LIMIT = (-open < matrix->min ? INT16_MIN + open : INT16_MIN - matrix->min) + 1;
    POS_LIMIT = INT16_MAX - matrix->max - 1;
    score = NEG_LIMIT;
    vNegLimit = _mm_set1_epi16(NEG_LIMIT);
    vPosLimit = _mm_set1_epi16(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;
    vNegInf = _mm_set1_epi16(NEG_LIMIT);
    vOpen = _mm_set1_epi16(open);
    vGap  = _mm_set1_epi16(gap);
    vOne = _mm_set1_epi16(1);
    vN = _mm_set1_epi16(N);
    vGapN = s1_beg ? _mm_set1_epi16(0) : _mm_set1_epi16(gap*N);
    vNegOne = _mm_set1_epi16(-1);
    vI = _mm_set_epi16(0,1,2,3,4,5,6,7);
    vJreset = _mm_set_epi16(0,-1,-2,-3,-4,-5,-6,-7);
    vMaxHRow = vNegInf;
    vMaxHCol = vNegInf;
    vLastVal = vNegInf;
    vEndI = vNegInf;
    vEndJ = vNegInf;
    vILimit = _mm_set1_epi16(s1Len);
    vILimit1 = _mm_subs_epi16(vILimit, vOne);
    vJLimit = _mm_set1_epi16(s2Len);
    vJLimit1 = _mm_subs_epi16(vJLimit, vOne);
    vIBoundary = s1_beg ? _mm_set1_epi16(0) : _mm_set_epi16(
            -open-0*gap,
            -open-1*gap,
            -open-2*gap,
            -open-3*gap,
            -open-4*gap,
            -open-5*gap,
            -open-6*gap,
            -open-7*gap);
    vTDiag = _mm_set1_epi16(PARASAIL_DIAG);
    vTIns = _mm_set1_epi16(PARASAIL_INS);
    vTDel = _mm_set1_epi16(PARASAIL_DEL);
    vTDiagE = _mm_set1_epi16(PARASAIL_DIAG_E);
    vTInsE = _mm_set1_epi16(PARASAIL_INS_E);
    vTDiagF = _mm_set1_epi16(PARASAIL_DIAG_F);
    vTDelF = _mm_set1_epi16(PARASAIL_DEL_F);

    /* initialize result */
    result = parasail_result_new_trace(s1Len, s2Len, 16, sizeof(int8_t));
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_16 | PARASAIL_FLAG_LANES_8;
    result->flag |= s1_beg ? PARASAIL_FLAG_SG_S1_BEG : 0;
    result->flag |= s1_end ? PARASAIL_FLAG_SG_S1_END : 0;
    result->flag |= s2_beg ? PARASAIL_FLAG_SG_S2_BEG : 0;
    result->flag |= s2_end ? PARASAIL_FLAG_SG_S2_END : 0;

    /* initialize heap variables */
    s2B= parasail_memalign_int16_t(16, s2Len+PAD2);
    _H_pr = parasail_memalign_int16_t(16, s2Len+PAD2);
    _F_pr = parasail_memalign_int16_t(16, s2Len+PAD2);
    s2 = s2B+PAD; /* will allow later for negative indices */
    H_pr = _H_pr+PAD;
    F_pr = _F_pr+PAD;

    /* validate heap variables */
    if (!s2B) return NULL;
    if (!_H_pr) return NULL;
    if (!_F_pr) return NULL;

    /* convert _s1 from char to int in range 0-23 */
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        s1 = parasail_memalign_int16_t(16, s1Len+PAD);
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
    if (s2_beg) {
        for (j=0; j<s2Len; ++j) {
            H_pr[j] = 0;
            F_pr[j] = NEG_LIMIT;
        }
    }
    else {
        for (j=0; j<s2Len; ++j) {
            H_pr[j] = -open - j*gap;
            F_pr[j] = NEG_LIMIT;
        }
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
    H_pr[-1] = 0; /* upper left corner */

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        __m128i vNH = vNegInf;
        __m128i vWH = vNegInf;
        __m128i vE = vNegInf;
        __m128i vE_opn = vNegInf;
        __m128i vE_ext = vNegInf;
        __m128i vF = vNegInf;
        __m128i vF_opn = vNegInf;
        __m128i vF_ext = vNegInf;
        __m128i vJ = vJreset;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+0] : ((i+0 >= s1Len) ? s1Len-1 : i+0))];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+1] : ((i+1 >= s1Len) ? s1Len-1 : i+1))];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+2] : ((i+2 >= s1Len) ? s1Len-1 : i+2))];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+3] : ((i+3 >= s1Len) ? s1Len-1 : i+3))];
        const int * const restrict matrow4 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+4] : ((i+4 >= s1Len) ? s1Len-1 : i+4))];
        const int * const restrict matrow5 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+5] : ((i+5 >= s1Len) ? s1Len-1 : i+5))];
        const int * const restrict matrow6 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+6] : ((i+6 >= s1Len) ? s1Len-1 : i+6))];
        const int * const restrict matrow7 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+7] : ((i+7 >= s1Len) ? s1Len-1 : i+7))];
        __m128i vIltLimit = _mm_cmplt_epi16(vI, vILimit);
        __m128i vIeqLimit1 = _mm_cmpeq_epi16(vI, vILimit1);
        vNH = _mm_srli_si128(vNH, 2);
        vNH = _mm_insert_epi16(vNH, H_pr[-1], 7);
        vWH = _mm_srli_si128(vWH, 2);
        vWH = _mm_insert_epi16(vWH, s1_beg ? 0 : (-open - i*gap), 7);
        H_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m128i vMat;
            __m128i vNWH = vNH;
            vNH = _mm_srli_si128(vWH, 2);
            vNH = _mm_insert_epi16(vNH, H_pr[j], 7);
            vF = _mm_srli_si128(vF, 2);
            vF = _mm_insert_epi16(vF, F_pr[j], 7);
            vF_opn = _mm_subs_epi16(vNH, vOpen);
            vF_ext = _mm_subs_epi16(vF, vGap);
            vF = _mm_max_epi16(vF_opn, vF_ext);
            vE_opn = _mm_subs_epi16(vWH, vOpen);
            vE_ext = _mm_subs_epi16(vE, vGap);
            vE = _mm_max_epi16(vE_opn, vE_ext);
            vMat = _mm_set_epi16(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]],
                    matrow4[s2[j-4]],
                    matrow5[s2[j-5]],
                    matrow6[s2[j-6]],
                    matrow7[s2[j-7]]
                    );
            vNWH = _mm_adds_epi16(vNWH, vMat);
            vWH = _mm_max_epi16(vNWH, vE);
            vWH = _mm_max_epi16(vWH, vF);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m128i cond = _mm_cmpeq_epi16(vJ,vNegOne);
                vWH = _mm_blendv_epi8(vWH, vIBoundary, cond);
                vF = _mm_blendv_epi8(vF, vNegInf, cond);
                vE = _mm_blendv_epi8(vE, vNegInf, cond);
            }
            /* cannot start checking sat until after J clears boundary */
            if (j > PAD) {
                vSaturationCheckMin = _mm_min_epi16(vSaturationCheckMin, vWH);
                vSaturationCheckMax = _mm_max_epi16(vSaturationCheckMax, vWH);
            }
            /* trace table */
            {
                __m128i case1 = _mm_cmpeq_epi16(vWH, vNWH);
                __m128i case2 = _mm_cmpeq_epi16(vWH, vF);
                __m128i vT = _mm_blendv_epi8(
                        _mm_blendv_epi8(vTIns, vTDel, case2),
                        vTDiag,
                        case1);
                __m128i condE = _mm_cmpgt_epi16(vE_opn, vE_ext);
                __m128i condF = _mm_cmpgt_epi16(vF_opn, vF_ext);
                __m128i vET = _mm_blendv_epi8(vTInsE, vTDiagE, condE);
                __m128i vFT = _mm_blendv_epi8(vTDelF, vTDiagF, condF);
                vT = _mm_or_si128(vT, vET);
                vT = _mm_or_si128(vT, vFT);
                arr_store_si128(result->trace->trace_table, vT, i, s1Len, j, s2Len);
            }
            H_pr[j-7] = (int16_t)_mm_extract_epi16(vWH,0);
            F_pr[j-7] = (int16_t)_mm_extract_epi16(vF,0);
            /* as minor diagonal vector passes across the i or j limit
             * boundary, extract the last value of the column or row */
            {
                __m128i vJeqLimit1 = _mm_cmpeq_epi16(vJ, vJLimit1);
                __m128i vJgtNegOne = _mm_cmpgt_epi16(vJ, vNegOne);
                __m128i vJltLimit = _mm_cmplt_epi16(vJ, vJLimit);
                __m128i cond_j = _mm_and_si128(vIltLimit, vJeqLimit1);
                __m128i cond_i = _mm_and_si128(vIeqLimit1,
                        _mm_and_si128(vJgtNegOne, vJltLimit));
                __m128i cond_max_row = _mm_cmpgt_epi16(vWH, vMaxHRow);
                __m128i cond_max_col = _mm_cmpgt_epi16(vWH, vMaxHCol);
                __m128i cond_last_val = _mm_and_si128(vIeqLimit1, vJeqLimit1);
                __m128i cond_all_row = _mm_and_si128(cond_max_row, cond_i);
                __m128i cond_all_col = _mm_and_si128(cond_max_col, cond_j);
                vMaxHRow = _mm_blendv_epi8(vMaxHRow, vWH, cond_all_row);
                vMaxHCol = _mm_blendv_epi8(vMaxHCol, vWH, cond_all_col);
                vLastVal = _mm_blendv_epi8(vLastVal, vWH, cond_last_val);
                vEndI = _mm_blendv_epi8(vEndI, vI, cond_all_col);
                vEndJ = _mm_blendv_epi8(vEndJ, vJ, cond_all_row);
            }
            vJ = _mm_adds_epi16(vJ, vOne);
        }
        vI = _mm_adds_epi16(vI, vN);
        vIBoundary = _mm_subs_epi16(vIBoundary, vGapN);
    }

    /* alignment ending position */
    {
        int16_t max_row = NEG_LIMIT;
        int16_t max_col = NEG_LIMIT;
        int16_t last_val = NEG_LIMIT;
        int16_t *s = (int16_t*)&vMaxHRow;
        int16_t *t = (int16_t*)&vMaxHCol;
        int16_t *u = (int16_t*)&vLastVal;
        int16_t *i = (int16_t*)&vEndI;
        int16_t *j = (int16_t*)&vEndJ;
        int32_t k;
        for (k=0; k<N; ++k, ++s, ++t, ++u, ++i, ++j) {
            if (*t > max_col || (*t == max_col && *i < end_query)) {
                max_col = *t;
                end_query = *i;
            }
            if (*s > max_row) {
                max_row = *s;
                end_ref = *j;
            }
            if (*u > last_val) {
                last_val = *u;
            }
        }
        if (s1_end && s2_end) {
            if (max_col > max_row || (max_col == max_row && end_ref == s2Len-1)) {
                score = max_col;
                end_ref = s2Len-1;
            }
            else {
                score = max_row;
                end_query = s1Len-1;
            }
        }
        else if (s1_end) {
            score = max_col;
            end_ref = s2Len-1;
        }
        else if (s2_end) {
            score = max_row;
            end_query = s1Len-1;
        }
        else {
            score = last_val;
            end_query = s1Len-1;
            end_ref = s2Len-1;
        }
    }

    if (_mm_movemask_epi8(_mm_or_si128(
            _mm_cmplt_epi16(vSaturationCheckMin, vNegLimit),
            _mm_cmpgt_epi16(vSaturationCheckMax, vPosLimit)))) {
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

SG_IMPL_ALL


