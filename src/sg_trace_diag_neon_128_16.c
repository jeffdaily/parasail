/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>



#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_neon.h"

#define SG_TRACE
#define SG_SUFFIX _diag_neon_128_16
#include "sg_helper.h"

#define NEG_INF (INT16_MIN/(int16_t)(2))


static inline void arr_store_si128(
        int8_t *array,
        simde__m128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[1LL*(i+0)*s2Len + (j-0)] = (int8_t)simde_mm_extract_epi16(vWH, 7);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[1LL*(i+1)*s2Len + (j-1)] = (int8_t)simde_mm_extract_epi16(vWH, 6);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[1LL*(i+2)*s2Len + (j-2)] = (int8_t)simde_mm_extract_epi16(vWH, 5);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[1LL*(i+3)*s2Len + (j-3)] = (int8_t)simde_mm_extract_epi16(vWH, 4);
    }
    if (0 <= i+4 && i+4 < s1Len && 0 <= j-4 && j-4 < s2Len) {
        array[1LL*(i+4)*s2Len + (j-4)] = (int8_t)simde_mm_extract_epi16(vWH, 3);
    }
    if (0 <= i+5 && i+5 < s1Len && 0 <= j-5 && j-5 < s2Len) {
        array[1LL*(i+5)*s2Len + (j-5)] = (int8_t)simde_mm_extract_epi16(vWH, 2);
    }
    if (0 <= i+6 && i+6 < s1Len && 0 <= j-6 && j-6 < s2Len) {
        array[1LL*(i+6)*s2Len + (j-6)] = (int8_t)simde_mm_extract_epi16(vWH, 1);
    }
    if (0 <= i+7 && i+7 < s1Len && 0 <= j-7 && j-7 < s2Len) {
        array[1LL*(i+7)*s2Len + (j-7)] = (int8_t)simde_mm_extract_epi16(vWH, 0);
    }
}

#define FNAME parasail_sg_flags_trace_diag_neon_128_16

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix,
        int s1_beg, int s1_end, int s2_beg, int s2_end)
{
    const int32_t N = 8; /* number of values in vector */
    const int32_t PAD = N-1;
    const int32_t PAD2 = PAD*2;
    const int32_t s1Len_PAD = s1Len+PAD;
    const int32_t s2Len_PAD = s2Len+PAD;
    int16_t * const restrict s1 = parasail_memalign_int16_t(16, s1Len+PAD);
    int16_t * const restrict s2B= parasail_memalign_int16_t(16, s2Len+PAD2);
    int16_t * const restrict _H_pr = parasail_memalign_int16_t(16, s2Len+PAD2);
    int16_t * const restrict _F_pr = parasail_memalign_int16_t(16, s2Len+PAD2);
    int16_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int16_t * const restrict H_pr = _H_pr+PAD;
    int16_t * const restrict F_pr = _F_pr+PAD;
    parasail_result_t *result = parasail_result_new_trace(s1Len, s2Len, 16, sizeof(int8_t));
    int32_t i = 0;
    int32_t j = 0;
    int32_t end_query = s1Len-1;
    int32_t end_ref = s2Len-1;
    int16_t score = NEG_INF;
    simde__m128i vNegInf = simde_mm_set1_epi16(NEG_INF);
    simde__m128i vOpen = simde_mm_set1_epi16(open);
    simde__m128i vGap  = simde_mm_set1_epi16(gap);
    simde__m128i vOne = simde_mm_set1_epi16(1);
    simde__m128i vN = simde_mm_set1_epi16(N);
    simde__m128i vGapN = s1_beg ? simde_mm_set1_epi16(0) : simde_mm_set1_epi16(gap*N);
    simde__m128i vNegOne = simde_mm_set1_epi16(-1);
    simde__m128i vI = simde_mm_set_epi16(0,1,2,3,4,5,6,7);
    simde__m128i vJreset = simde_mm_set_epi16(0,-1,-2,-3,-4,-5,-6,-7);
    simde__m128i vMaxHRow = vNegInf;
    simde__m128i vMaxHCol = vNegInf;
    simde__m128i vLastVal = vNegInf;
    simde__m128i vEndI = vNegInf;
    simde__m128i vEndJ = vNegInf;
    simde__m128i vILimit = simde_mm_set1_epi16(s1Len);
    simde__m128i vILimit1 = simde_mm_sub_epi16(vILimit, vOne);
    simde__m128i vJLimit = simde_mm_set1_epi16(s2Len);
    simde__m128i vJLimit1 = simde_mm_sub_epi16(vJLimit, vOne);
    simde__m128i vIBoundary = s1_beg ? simde_mm_set1_epi16(0) : simde_mm_set_epi16(
            -open-0*gap,
            -open-1*gap,
            -open-2*gap,
            -open-3*gap,
            -open-4*gap,
            -open-5*gap,
            -open-6*gap,
            -open-7*gap
            );
    simde__m128i vTDiag = simde_mm_set1_epi16(PARASAIL_DIAG);
    simde__m128i vTIns = simde_mm_set1_epi16(PARASAIL_INS);
    simde__m128i vTDel = simde_mm_set1_epi16(PARASAIL_DEL);
    simde__m128i vTDiagE = simde_mm_set1_epi16(PARASAIL_DIAG_E);
    simde__m128i vTInsE = simde_mm_set1_epi16(PARASAIL_INS_E);
    simde__m128i vTDiagF = simde_mm_set1_epi16(PARASAIL_DIAG_F);
    simde__m128i vTDelF = simde_mm_set1_epi16(PARASAIL_DEL_F);
    

    /* convert _s1 from char to int in range 0-23 */
    for (i=0; i<s1Len; ++i) {
        s1[i] = matrix->mapper[(unsigned char)_s1[i]];
    }
    /* pad back of s1 with dummy values */
    for (i=s1Len; i<s1Len_PAD; ++i) {
        s1[i] = 0; /* point to first matrix row because we don't care */
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
            F_pr[j] = NEG_INF;
        }
    }
    else {
        for (j=0; j<s2Len; ++j) {
            H_pr[j] = -open - j*gap;
            F_pr[j] = NEG_INF;
        }
    }
    /* pad front of stored row values */
    for (j=-PAD; j<0; ++j) {
        H_pr[j] = NEG_INF;
        F_pr[j] = NEG_INF;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        H_pr[j] = NEG_INF;
        F_pr[j] = NEG_INF;
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
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size*s1[i+2]];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size*s1[i+3]];
        const int * const restrict matrow4 = &matrix->matrix[matrix->size*s1[i+4]];
        const int * const restrict matrow5 = &matrix->matrix[matrix->size*s1[i+5]];
        const int * const restrict matrow6 = &matrix->matrix[matrix->size*s1[i+6]];
        const int * const restrict matrow7 = &matrix->matrix[matrix->size*s1[i+7]];
        simde__m128i vIltLimit = simde_mm_cmplt_epi16(vI, vILimit);
        simde__m128i vIeqLimit1 = simde_mm_cmpeq_epi16(vI, vILimit1);
        vNH = simde_mm_srli_si128(vNH, 2);
        vNH = simde_mm_insert_epi16(vNH, H_pr[-1], 7);
        vWH = simde_mm_srli_si128(vWH, 2);
        vWH = simde_mm_insert_epi16(vWH, s1_beg ? 0 : (-open - i*gap), 7);
        H_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            simde__m128i vMat;
            simde__m128i vNWH = vNH;
            vNH = simde_mm_srli_si128(vWH, 2);
            vNH = simde_mm_insert_epi16(vNH, H_pr[j], 7);
            vF = simde_mm_srli_si128(vF, 2);
            vF = simde_mm_insert_epi16(vF, F_pr[j], 7);
            vF_opn = simde_mm_sub_epi16(vNH, vOpen);
            vF_ext = simde_mm_sub_epi16(vF, vGap);
            vF = simde_mm_max_epi16(vF_opn, vF_ext);
            vE_opn = simde_mm_sub_epi16(vWH, vOpen);
            vE_ext = simde_mm_sub_epi16(vE, vGap);
            vE = simde_mm_max_epi16(vE_opn, vE_ext);
            vMat = simde_mm_set_epi16(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]],
                    matrow4[s2[j-4]],
                    matrow5[s2[j-5]],
                    matrow6[s2[j-6]],
                    matrow7[s2[j-7]]
                    );
            vNWH = simde_mm_add_epi16(vNWH, vMat);
            vWH = simde_mm_max_epi16(vNWH, vE);
            vWH = simde_mm_max_epi16(vWH, vF);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                simde__m128i cond = simde_mm_cmpeq_epi16(vJ,vNegOne);
                vWH = simde_mm_blendv_epi8(vWH, vIBoundary, cond);
                vF = simde_mm_blendv_epi8(vF, vNegInf, cond);
                vE = simde_mm_blendv_epi8(vE, vNegInf, cond);
            }
            
            /* trace table */
            {
                simde__m128i case1 = simde_mm_cmpeq_epi16(vWH, vNWH);
                simde__m128i case2 = simde_mm_cmpeq_epi16(vWH, vF);
                simde__m128i vT = simde_mm_blendv_epi8(
                        simde_mm_blendv_epi8(vTIns, vTDel, case2),
                        vTDiag,
                        case1);
                simde__m128i condE = simde_mm_cmpgt_epi16(vE_opn, vE_ext);
                simde__m128i condF = simde_mm_cmpgt_epi16(vF_opn, vF_ext);
                simde__m128i vET = simde_mm_blendv_epi8(vTInsE, vTDiagE, condE);
                simde__m128i vFT = simde_mm_blendv_epi8(vTDelF, vTDiagF, condF);
                vT = simde_mm_or_si128(vT, vET);
                vT = simde_mm_or_si128(vT, vFT);
                arr_store_si128(result->trace->trace_table, vT, i, s1Len, j, s2Len);
            }
            H_pr[j-7] = (int16_t)simde_mm_extract_epi16(vWH,0);
            F_pr[j-7] = (int16_t)simde_mm_extract_epi16(vF,0);
            /* as minor diagonal vector passes across the i or j limit
             * boundary, extract the last value of the column or row */
            {
                simde__m128i vJeqLimit1 = simde_mm_cmpeq_epi16(vJ, vJLimit1);
                simde__m128i vJgtNegOne = simde_mm_cmpgt_epi16(vJ, vNegOne);
                simde__m128i vJltLimit = simde_mm_cmplt_epi16(vJ, vJLimit);
                simde__m128i cond_j = simde_mm_and_si128(vIltLimit, vJeqLimit1);
                simde__m128i cond_i = simde_mm_and_si128(vIeqLimit1,
                        simde_mm_and_si128(vJgtNegOne, vJltLimit));
                simde__m128i cond_max_row = simde_mm_cmpgt_epi16(vWH, vMaxHRow);
                simde__m128i cond_max_col = simde_mm_cmpgt_epi16(vWH, vMaxHCol);
                simde__m128i cond_last_val = simde_mm_and_si128(vIeqLimit1, vJeqLimit1);
                simde__m128i cond_all_row = simde_mm_and_si128(cond_max_row, cond_i);
                simde__m128i cond_all_col = simde_mm_and_si128(cond_max_col, cond_j);
                vMaxHRow = simde_mm_blendv_epi8(vMaxHRow, vWH, cond_all_row);
                vMaxHCol = simde_mm_blendv_epi8(vMaxHCol, vWH, cond_all_col);
                vLastVal = simde_mm_blendv_epi8(vLastVal, vWH, cond_last_val);
                vEndI = simde_mm_blendv_epi8(vEndI, vI, cond_all_col);
                vEndJ = simde_mm_blendv_epi8(vEndJ, vJ, cond_all_row);
            }
            vJ = simde_mm_add_epi16(vJ, vOne);
        }
        vI = simde_mm_add_epi16(vI, vN);
        vIBoundary = simde_mm_sub_epi16(vIBoundary, vGapN);
    }

    /* alignment ending position */
    {
        int16_t max_row = NEG_INF;
        int16_t max_col = NEG_INF;
        int16_t last_val = NEG_INF;
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

    

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_16 | PARASAIL_FLAG_LANES_8;
    result->flag |= s1_beg ? PARASAIL_FLAG_SG_S1_BEG : 0;
    result->flag |= s1_end ? PARASAIL_FLAG_SG_S1_END : 0;
    result->flag |= s2_beg ? PARASAIL_FLAG_SG_S2_BEG : 0;
    result->flag |= s2_end ? PARASAIL_FLAG_SG_S2_END : 0;

    parasail_free(_F_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}

SG_IMPL_ALL


