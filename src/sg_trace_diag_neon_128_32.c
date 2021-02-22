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
#include "parasail/internal_neon.h"

#define SG_TRACE
#define SG_SUFFIX _diag_neon_128_32
#include "sg_helper.h"



static inline void arr_store_si128(
        int8_t *array,
        simde__m128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[1LL*(i+0)*s2Len + (j-0)] = (int8_t)simde_mm_extract_epi32(vWH, 3);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[1LL*(i+1)*s2Len + (j-1)] = (int8_t)simde_mm_extract_epi32(vWH, 2);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[1LL*(i+2)*s2Len + (j-2)] = (int8_t)simde_mm_extract_epi32(vWH, 1);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[1LL*(i+3)*s2Len + (j-3)] = (int8_t)simde_mm_extract_epi32(vWH, 0);
    }
}

#define FNAME parasail_sg_flags_trace_diag_neon_128_32

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
    int32_t s1Len = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int32_t NEG_LIMIT = 0;
    int32_t POS_LIMIT = 0;
    int32_t score = 0;
    simde__m128i vNegLimit;
    simde__m128i vPosLimit;
    simde__m128i vSaturationCheckMin;
    simde__m128i vSaturationCheckMax;
    simde__m128i vNegInf;
    simde__m128i vOpen;
    simde__m128i vGap;
    simde__m128i vOne;
    simde__m128i vN;
    simde__m128i vGapN;
    simde__m128i vNegOne;
    simde__m128i vI;
    simde__m128i vJreset;
    simde__m128i vMaxHRow;
    simde__m128i vMaxHCol;
    simde__m128i vLastVal;
    simde__m128i vEndI;
    simde__m128i vEndJ;
    simde__m128i vILimit;
    simde__m128i vILimit1;
    simde__m128i vJLimit;
    simde__m128i vJLimit1;
    simde__m128i vIBoundary;
    simde__m128i vTDiag;
    simde__m128i vTIns;
    simde__m128i vTDel;
    simde__m128i vTDiagE;
    simde__m128i vTInsE;
    simde__m128i vTDiagF;
    simde__m128i vTDelF;

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
    end_query = s1Len-1;
    end_ref = s2Len-1;
    NEG_LIMIT = (-open < matrix->min ? INT32_MIN + open : INT32_MIN - matrix->min) + 1;
    POS_LIMIT = INT32_MAX - matrix->max - 1;
    score = NEG_LIMIT;
    vNegLimit = simde_mm_set1_epi32(NEG_LIMIT);
    vPosLimit = simde_mm_set1_epi32(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;
    vNegInf = simde_mm_set1_epi32(NEG_LIMIT);
    vOpen = simde_mm_set1_epi32(open);
    vGap  = simde_mm_set1_epi32(gap);
    vOne = simde_mm_set1_epi32(1);
    vN = simde_mm_set1_epi32(N);
    vGapN = s1_beg ? simde_mm_set1_epi32(0) : simde_mm_set1_epi32(gap*N);
    vNegOne = simde_mm_set1_epi32(-1);
    vI = simde_mm_set_epi32(0,1,2,3);
    vJreset = simde_mm_set_epi32(0,-1,-2,-3);
    vMaxHRow = vNegInf;
    vMaxHCol = vNegInf;
    vLastVal = vNegInf;
    vEndI = vNegInf;
    vEndJ = vNegInf;
    vILimit = simde_mm_set1_epi32(s1Len);
    vILimit1 = simde_mm_sub_epi32(vILimit, vOne);
    vJLimit = simde_mm_set1_epi32(s2Len);
    vJLimit1 = simde_mm_sub_epi32(vJLimit, vOne);
    vIBoundary = s1_beg ? simde_mm_set1_epi32(0) : simde_mm_set_epi32(
            -open-0*gap,
            -open-1*gap,
            -open-2*gap,
            -open-3*gap);
    vTDiag = simde_mm_set1_epi32(PARASAIL_DIAG);
    vTIns = simde_mm_set1_epi32(PARASAIL_INS);
    vTDel = simde_mm_set1_epi32(PARASAIL_DEL);
    vTDiagE = simde_mm_set1_epi32(PARASAIL_DIAG_E);
    vTInsE = simde_mm_set1_epi32(PARASAIL_INS_E);
    vTDiagF = simde_mm_set1_epi32(PARASAIL_DIAG_F);
    vTDelF = simde_mm_set1_epi32(PARASAIL_DEL_F);

    /* initialize result */
    result = parasail_result_new_trace(s1Len, s2Len, 16, sizeof(int8_t));
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_32 | PARASAIL_FLAG_LANES_4;
    result->flag |= s1_beg ? PARASAIL_FLAG_SG_S1_BEG : 0;
    result->flag |= s1_end ? PARASAIL_FLAG_SG_S1_END : 0;
    result->flag |= s2_beg ? PARASAIL_FLAG_SG_S2_BEG : 0;
    result->flag |= s2_end ? PARASAIL_FLAG_SG_S2_END : 0;

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
        simde__m128i vNH = vNegInf;
        simde__m128i vWH = vNegInf;
        simde__m128i vE = vNegInf;
        simde__m128i vE_opn = vNegInf;
        simde__m128i vE_ext = vNegInf;
        simde__m128i vF = vNegInf;
        simde__m128i vF_opn = vNegInf;
        simde__m128i vF_ext = vNegInf;
        simde__m128i vJ = vJreset;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+0] : ((i+0 >= s1Len) ? s1Len-1 : i+0))];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+1] : ((i+1 >= s1Len) ? s1Len-1 : i+1))];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+2] : ((i+2 >= s1Len) ? s1Len-1 : i+2))];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+3] : ((i+3 >= s1Len) ? s1Len-1 : i+3))];
        simde__m128i vIltLimit = simde_mm_cmplt_epi32(vI, vILimit);
        simde__m128i vIeqLimit1 = simde_mm_cmpeq_epi32(vI, vILimit1);
        vNH = simde_mm_srli_si128(vNH, 4);
        vNH = simde_mm_insert_epi32(vNH, H_pr[-1], 3);
        vWH = simde_mm_srli_si128(vWH, 4);
        vWH = simde_mm_insert_epi32(vWH, s1_beg ? 0 : (-open - i*gap), 3);
        H_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            simde__m128i vMat;
            simde__m128i vNWH = vNH;
            vNH = simde_mm_srli_si128(vWH, 4);
            vNH = simde_mm_insert_epi32(vNH, H_pr[j], 3);
            vF = simde_mm_srli_si128(vF, 4);
            vF = simde_mm_insert_epi32(vF, F_pr[j], 3);
            vF_opn = simde_mm_sub_epi32(vNH, vOpen);
            vF_ext = simde_mm_sub_epi32(vF, vGap);
            vF = simde_mm_max_epi32(vF_opn, vF_ext);
            vE_opn = simde_mm_sub_epi32(vWH, vOpen);
            vE_ext = simde_mm_sub_epi32(vE, vGap);
            vE = simde_mm_max_epi32(vE_opn, vE_ext);
            vMat = simde_mm_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]]
                    );
            vNWH = simde_mm_add_epi32(vNWH, vMat);
            vWH = simde_mm_max_epi32(vNWH, vE);
            vWH = simde_mm_max_epi32(vWH, vF);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                simde__m128i cond = simde_mm_cmpeq_epi32(vJ,vNegOne);
                vWH = simde_mm_blendv_epi8(vWH, vIBoundary, cond);
                vF = simde_mm_blendv_epi8(vF, vNegInf, cond);
                vE = simde_mm_blendv_epi8(vE, vNegInf, cond);
            }
            /* cannot start checking sat until after J clears boundary */
            if (j > PAD) {
                vSaturationCheckMin = simde_mm_min_epi32(vSaturationCheckMin, vWH);
                vSaturationCheckMax = simde_mm_max_epi32(vSaturationCheckMax, vWH);
            }
            /* trace table */
            {
                simde__m128i case1 = simde_mm_cmpeq_epi32(vWH, vNWH);
                simde__m128i case2 = simde_mm_cmpeq_epi32(vWH, vF);
                simde__m128i vT = simde_mm_blendv_epi8(
                        simde_mm_blendv_epi8(vTIns, vTDel, case2),
                        vTDiag,
                        case1);
                simde__m128i condE = simde_mm_cmpgt_epi32(vE_opn, vE_ext);
                simde__m128i condF = simde_mm_cmpgt_epi32(vF_opn, vF_ext);
                simde__m128i vET = simde_mm_blendv_epi8(vTInsE, vTDiagE, condE);
                simde__m128i vFT = simde_mm_blendv_epi8(vTDelF, vTDiagF, condF);
                vT = simde_mm_or_si128(vT, vET);
                vT = simde_mm_or_si128(vT, vFT);
                arr_store_si128(result->trace->trace_table, vT, i, s1Len, j, s2Len);
            }
            H_pr[j-3] = (int32_t)simde_mm_extract_epi32(vWH,0);
            F_pr[j-3] = (int32_t)simde_mm_extract_epi32(vF,0);
            /* as minor diagonal vector passes across the i or j limit
             * boundary, extract the last value of the column or row */
            {
                simde__m128i vJeqLimit1 = simde_mm_cmpeq_epi32(vJ, vJLimit1);
                simde__m128i vJgtNegOne = simde_mm_cmpgt_epi32(vJ, vNegOne);
                simde__m128i vJltLimit = simde_mm_cmplt_epi32(vJ, vJLimit);
                simde__m128i cond_j = simde_mm_and_si128(vIltLimit, vJeqLimit1);
                simde__m128i cond_i = simde_mm_and_si128(vIeqLimit1,
                        simde_mm_and_si128(vJgtNegOne, vJltLimit));
                simde__m128i cond_max_row = simde_mm_cmpgt_epi32(vWH, vMaxHRow);
                simde__m128i cond_max_col = simde_mm_cmpgt_epi32(vWH, vMaxHCol);
                simde__m128i cond_last_val = simde_mm_and_si128(vIeqLimit1, vJeqLimit1);
                simde__m128i cond_all_row = simde_mm_and_si128(cond_max_row, cond_i);
                simde__m128i cond_all_col = simde_mm_and_si128(cond_max_col, cond_j);
                vMaxHRow = simde_mm_blendv_epi8(vMaxHRow, vWH, cond_all_row);
                vMaxHCol = simde_mm_blendv_epi8(vMaxHCol, vWH, cond_all_col);
                vLastVal = simde_mm_blendv_epi8(vLastVal, vWH, cond_last_val);
                vEndI = simde_mm_blendv_epi8(vEndI, vI, cond_all_col);
                vEndJ = simde_mm_blendv_epi8(vEndJ, vJ, cond_all_row);
            }
            vJ = simde_mm_add_epi32(vJ, vOne);
        }
        vI = simde_mm_add_epi32(vI, vN);
        vIBoundary = simde_mm_sub_epi32(vIBoundary, vGapN);
    }

    /* alignment ending position */
    {
        int32_t max_row = NEG_LIMIT;
        int32_t max_col = NEG_LIMIT;
        int32_t last_val = NEG_LIMIT;
        int32_t *s = (int32_t*)&vMaxHRow;
        int32_t *t = (int32_t*)&vMaxHCol;
        int32_t *u = (int32_t*)&vLastVal;
        int32_t *i = (int32_t*)&vEndI;
        int32_t *j = (int32_t*)&vEndJ;
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

    if (simde_mm_movemask_epi8(simde_mm_or_si128(
            simde_mm_cmplt_epi32(vSaturationCheckMin, vNegLimit),
            simde_mm_cmpgt_epi32(vSaturationCheckMax, vPosLimit)))) {
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


