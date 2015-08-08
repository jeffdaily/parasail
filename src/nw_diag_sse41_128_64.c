/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>

#include <emmintrin.h>
#include <smmintrin.h>

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_sse.h"

#define NEG_INF (INT64_MIN/(int64_t)(2))

static inline __m128i _mm_max_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]>B.v[1]) ? A.v[1] : B.v[1];
    return A.m;
}


#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        __m128i vWscore,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[(i+0)*s2Len + (j-0)] = (int64_t)_mm_extract_epi64(vWscore, 1);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[(i+1)*s2Len + (j-1)] = (int64_t)_mm_extract_epi64(vWscore, 0);
    }
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_rowcol(
        int *row,
        int *col,
        __m128i vWscore,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (i+0 == s1Len-1 && 0 <= j-0 && j-0 < s2Len) {
        row[j-0] = (int64_t)_mm_extract_epi64(vWscore, 1);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int64_t)_mm_extract_epi64(vWscore, 1);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int64_t)_mm_extract_epi64(vWscore, 0);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int64_t)_mm_extract_epi64(vWscore, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_nw_table_diag_sse41_128_64
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_rowcol_diag_sse41_128_64
#else
#define FNAME parasail_nw_diag_sse41_128_64
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    const int32_t N = 2; /* number of values in vector */
    const int32_t PAD = N-1;
    const int32_t PAD2 = PAD*2;
    const int32_t s1Len_PAD = s1Len+PAD;
    const int32_t s2Len_PAD = s2Len+PAD;
    int64_t * const restrict s1 = parasail_memalign_int64_t(16, s1Len+PAD);
    int64_t * const restrict s2B= parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict _tbl_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict _del_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int64_t * const restrict tbl_pr = _tbl_pr+PAD;
    int64_t * const restrict del_pr = _del_pr+PAD;
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
    __m128i vNegInf = _mm_set1_epi64x(NEG_INF);
    __m128i vOpen = _mm_set1_epi64x(open);
    __m128i vGap  = _mm_set1_epi64x(gap);
    __m128i vOne = _mm_set1_epi64x(1);
    __m128i vN = _mm_set1_epi64x(N);
    __m128i vGapN = _mm_set1_epi64x(gap*N);
    __m128i vNegOne = _mm_set1_epi64x(-1);
    __m128i vI = _mm_set_epi64x(0,1);
    __m128i vJreset = _mm_set_epi64x(0,-1);
    __m128i vMax = vNegInf;
    __m128i vILimit = _mm_set1_epi64x(s1Len);
    __m128i vILimit1 = _mm_sub_epi64(vILimit, vOne);
    __m128i vJLimit = _mm_set1_epi64x(s2Len);
    __m128i vJLimit1 = _mm_sub_epi64(vJLimit, vOne);
    __m128i vIBoundary = _mm_set_epi64x(
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
    for (j=0; j<s2Len; ++j) {
        tbl_pr[j] = -open - j*gap;
        del_pr[j] = NEG_INF;
    }
    /* pad front of stored row values */
    for (j=-PAD; j<0; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
    }
    tbl_pr[-1] = 0; /* upper left corner */

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        __m128i vNscore = vNegInf;
        __m128i vWscore = vNegInf;
        __m128i vIns = vNegInf;
        __m128i vDel = vNegInf;
        __m128i vJ = vJreset;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        vNscore = _mm_srli_si128(vNscore, 8);
        vNscore = _mm_insert_epi64(vNscore, tbl_pr[-1], 1);
        vWscore = _mm_srli_si128(vWscore, 8);
        vWscore = _mm_insert_epi64(vWscore, -open - i*gap, 1);
        tbl_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            vNscore = _mm_srli_si128(vWscore, 8);
            vNscore = _mm_insert_epi64(vNscore, tbl_pr[j], 1);
            vDel = _mm_srli_si128(vDel, 8);
            vDel = _mm_insert_epi64(vDel, del_pr[j], 1);
            vDel = _mm_max_epi64_rpl(
                    _mm_sub_epi64(vNscore, vOpen),
                    _mm_sub_epi64(vDel, vGap));
            vIns = _mm_max_epi64_rpl(
                    _mm_sub_epi64(vWscore, vOpen),
                    _mm_sub_epi64(vIns, vGap));
            vMat = _mm_set_epi64x(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]]
                    );
            vNWscore = _mm_add_epi64(vNWscore, vMat);
            vWscore = _mm_max_epi64_rpl(vNWscore, vIns);
            vWscore = _mm_max_epi64_rpl(vWscore, vDel);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m128i cond = _mm_cmpeq_epi64(vJ,vNegOne);
                vWscore = _mm_blendv_epi8(vWscore, vIBoundary, cond);
                vDel = _mm_blendv_epi8(vDel, vNegInf, cond);
                vIns = _mm_blendv_epi8(vIns, vNegInf, cond);
            }
            
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->score_row, result->score_col, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-1] = (int64_t)_mm_extract_epi64(vWscore,0);
            del_pr[j-1] = (int64_t)_mm_extract_epi64(vDel,0);
            /* as minor diagonal vector passes across table, extract
               last table value at the i,j bound */
            {
                __m128i cond_valid_I = _mm_cmpeq_epi64(vI, vILimit1);
                __m128i cond_valid_J = _mm_cmpeq_epi64(vJ, vJLimit1);
                __m128i cond_all = _mm_and_si128(cond_valid_I, cond_valid_J);
                vMax = _mm_blendv_epi8(vMax, vWscore, cond_all);
            }
            vJ = _mm_add_epi64(vJ, vOne);
        }
        vI = _mm_add_epi64(vI, vN);
        vIBoundary = _mm_sub_epi64(vIBoundary, vGapN);
    }

    /* max in vMax */
    for (i=0; i<N; ++i) {
        int64_t value;
        value = (int64_t) _mm_extract_epi64(vMax, 1);
        if (value > score) {
            score = value;
        }
        vMax = _mm_slli_si128(vMax, 8);
    }

    

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(_del_pr);
    parasail_free(_tbl_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}


