/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
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

#define SG_SUFFIX _diag_sse41_128_64
#include "sg_helper.h"

#define NEG_INF (INT64_MIN/(int64_t)(2))

static inline __m128i _mm_cmpgt_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? 0xFFFFFFFFFFFFFFFF : 0;
    A.v[1] = (A.v[1]>B.v[1]) ? 0xFFFFFFFFFFFFFFFF : 0;
    return A.m;
}

#if HAVE_SSE41_MM_INSERT_EPI64
#define _mm_insert_epi64_rpl _mm_insert_epi64
#else
static inline __m128i _mm_insert_epi64_rpl(__m128i a, int64_t i, int imm) {
    __m128i_64_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif

static inline __m128i _mm_max_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]>B.v[1]) ? A.v[1] : B.v[1];
    return A.m;
}

#if HAVE_SSE2_MM_SET_EPI64X
#define _mm_set_epi64x_rpl _mm_set_epi64x
#else
static inline __m128i _mm_set_epi64x_rpl(int64_t e1, int64_t e0) {
    __m128i_64_t A;
    A.v[0] = e0;
    A.v[1] = e1;
    return A.m;
}
#endif

#if HAVE_SSE41_MM_EXTRACT_EPI64
#define _mm_extract_epi64_rpl _mm_extract_epi64
#else
static inline int64_t _mm_extract_epi64_rpl(__m128i a, int imm) {
    __m128i_64_t A;
    A.m = a;
    return A.v[imm];
}
#endif

static inline __m128i _mm_cmplt_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]<B.v[0]) ? 0xFFFFFFFFFFFFFFFF : 0;
    A.v[1] = (A.v[1]<B.v[1]) ? 0xFFFFFFFFFFFFFFFF : 0;
    return A.m;
}

#if HAVE_SSE2_MM_SET1_EPI64X
#define _mm_set1_epi64x_rpl _mm_set1_epi64x
#else
static inline __m128i _mm_set1_epi64x_rpl(int64_t i) {
    __m128i_64_t A;
    A.v[0] = i;
    A.v[1] = i;
    return A.m;
}
#endif


#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        __m128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[1LL*(i+0)*s2Len + (j-0)] = (int64_t)_mm_extract_epi64_rpl(vWH, 1);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[1LL*(i+1)*s2Len + (j-1)] = (int64_t)_mm_extract_epi64_rpl(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_rowcol(
        int *row,
        int *col,
        __m128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (i+0 == s1Len-1 && 0 <= j-0 && j-0 < s2Len) {
        row[j-0] = (int64_t)_mm_extract_epi64_rpl(vWH, 1);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int64_t)_mm_extract_epi64_rpl(vWH, 1);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int64_t)_mm_extract_epi64_rpl(vWH, 0);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int64_t)_mm_extract_epi64_rpl(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sg_flags_table_diag_sse41_128_64
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_flags_rowcol_diag_sse41_128_64
#else
#define FNAME parasail_sg_flags_diag_sse41_128_64
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix,
        int s1_beg, int s1_end, int s2_beg, int s2_end)
{
    const int32_t N = 2; /* number of values in vector */
    const int32_t PAD = N-1;
    const int32_t PAD2 = PAD*2;
    const int32_t s1Len_PAD = s1Len+PAD;
    const int32_t s2Len_PAD = s2Len+PAD;
    int64_t * const restrict s1 = parasail_memalign_int64_t(16, s1Len+PAD);
    int64_t * const restrict s2B= parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict _H_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict _F_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int64_t * const restrict H_pr = _H_pr+PAD;
    int64_t * const restrict F_pr = _F_pr+PAD;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif
    int32_t i = 0;
    int32_t j = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int64_t score = NEG_INF;
    __m128i vNegInf = _mm_set1_epi64x_rpl(NEG_INF);
    __m128i vOpen = _mm_set1_epi64x_rpl(open);
    __m128i vGap  = _mm_set1_epi64x_rpl(gap);
    __m128i vOne = _mm_set1_epi64x_rpl(1);
    __m128i vN = _mm_set1_epi64x_rpl(N);
    __m128i vGapN = s1_beg ? _mm_set1_epi64x_rpl(0) : _mm_set1_epi64x_rpl(gap*N);
    __m128i vNegOne = _mm_set1_epi64x_rpl(-1);
    __m128i vI = _mm_set_epi64x_rpl(0,1);
    __m128i vJreset = _mm_set_epi64x_rpl(0,-1);
    __m128i vMaxHRow = vNegInf;
    __m128i vMaxHCol = vNegInf;
    __m128i vLastVal = vNegInf;
    __m128i vEndI = vNegInf;
    __m128i vEndJ = vNegInf;
    __m128i vILimit = _mm_set1_epi64x_rpl(s1Len);
    __m128i vILimit1 = _mm_sub_epi64(vILimit, vOne);
    __m128i vJLimit = _mm_set1_epi64x_rpl(s2Len);
    __m128i vJLimit1 = _mm_sub_epi64(vJLimit, vOne);
    __m128i vIBoundary = s1_beg ? _mm_set1_epi64x_rpl(0) : _mm_set_epi64x_rpl(
            -open-0*gap,
            -open-1*gap
            );
    

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
        __m128i vNH = vNegInf;
        __m128i vWH = vNegInf;
        __m128i vE = vNegInf;
        __m128i vF = vNegInf;
        __m128i vJ = vJreset;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        __m128i vIltLimit = _mm_cmplt_epi64_rpl(vI, vILimit);
        __m128i vIeqLimit1 = _mm_cmpeq_epi64(vI, vILimit1);
        vNH = _mm_srli_si128(vNH, 8);
        vNH = _mm_insert_epi64_rpl(vNH, H_pr[-1], 1);
        vWH = _mm_srli_si128(vWH, 8);
        vWH = _mm_insert_epi64_rpl(vWH, s1_beg ? 0 : (-open - i*gap), 1);
        H_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m128i vMat;
            __m128i vNWH = vNH;
            vNH = _mm_srli_si128(vWH, 8);
            vNH = _mm_insert_epi64_rpl(vNH, H_pr[j], 1);
            vF = _mm_srli_si128(vF, 8);
            vF = _mm_insert_epi64_rpl(vF, F_pr[j], 1);
            vF = _mm_max_epi64_rpl(
                    _mm_sub_epi64(vNH, vOpen),
                    _mm_sub_epi64(vF, vGap));
            vE = _mm_max_epi64_rpl(
                    _mm_sub_epi64(vWH, vOpen),
                    _mm_sub_epi64(vE, vGap));
            vMat = _mm_set_epi64x_rpl(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]]
                    );
            vNWH = _mm_add_epi64(vNWH, vMat);
            vWH = _mm_max_epi64_rpl(vNWH, vE);
            vWH = _mm_max_epi64_rpl(vWH, vF);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m128i cond = _mm_cmpeq_epi64(vJ,vNegOne);
                vWH = _mm_blendv_epi8(vWH, vIBoundary, cond);
                vF = _mm_blendv_epi8(vF, vNegInf, cond);
                vE = _mm_blendv_epi8(vE, vNegInf, cond);
            }
            
#ifdef PARASAIL_TABLE
            arr_store_si128(result->tables->score_table, vWH, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->rowcols->score_row, result->rowcols->score_col, vWH, i, s1Len, j, s2Len);
#endif
            H_pr[j-1] = (int64_t)_mm_extract_epi64_rpl(vWH,0);
            F_pr[j-1] = (int64_t)_mm_extract_epi64_rpl(vF,0);
            /* as minor diagonal vector passes across the i or j limit
             * boundary, extract the last value of the column or row */
            {
                __m128i vJeqLimit1 = _mm_cmpeq_epi64(vJ, vJLimit1);
                __m128i vJgtNegOne = _mm_cmpgt_epi64_rpl(vJ, vNegOne);
                __m128i vJltLimit = _mm_cmplt_epi64_rpl(vJ, vJLimit);
                __m128i cond_j = _mm_and_si128(vIltLimit, vJeqLimit1);
                __m128i cond_i = _mm_and_si128(vIeqLimit1,
                        _mm_and_si128(vJgtNegOne, vJltLimit));
                __m128i cond_max_row = _mm_cmpgt_epi64_rpl(vWH, vMaxHRow);
                __m128i cond_max_col = _mm_cmpgt_epi64_rpl(vWH, vMaxHCol);
                __m128i cond_last_val = _mm_and_si128(vIeqLimit1, vJeqLimit1);
                __m128i cond_all_row = _mm_and_si128(cond_max_row, cond_i);
                __m128i cond_all_col = _mm_and_si128(cond_max_col, cond_j);
                vMaxHRow = _mm_blendv_epi8(vMaxHRow, vWH, cond_all_row);
                vMaxHCol = _mm_blendv_epi8(vMaxHCol, vWH, cond_all_col);
                vLastVal = _mm_blendv_epi8(vLastVal, vWH, cond_last_val);
                vEndI = _mm_blendv_epi8(vEndI, vI, cond_all_col);
                vEndJ = _mm_blendv_epi8(vEndJ, vJ, cond_all_row);
            }
            vJ = _mm_add_epi64(vJ, vOne);
        }
        vI = _mm_add_epi64(vI, vN);
        vIBoundary = _mm_sub_epi64(vIBoundary, vGapN);
    }

    /* alignment ending position */
    {
        int64_t max_row = NEG_INF;
        int64_t max_col = NEG_INF;
        int64_t last_val = NEG_INF;
        int64_t *s = (int64_t*)&vMaxHRow;
        int64_t *t = (int64_t*)&vMaxHCol;
        int64_t *u = (int64_t*)&vLastVal;
        int64_t *i = (int64_t*)&vEndI;
        int64_t *j = (int64_t*)&vEndJ;
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
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_2;
    result->flag |= s1_beg ? PARASAIL_FLAG_SG_S1_BEG : 0;
    result->flag |= s1_end ? PARASAIL_FLAG_SG_S1_END : 0;
    result->flag |= s2_beg ? PARASAIL_FLAG_SG_S2_BEG : 0;
    result->flag |= s2_end ? PARASAIL_FLAG_SG_S2_END : 0;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(_F_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}

SG_IMPL_ALL


