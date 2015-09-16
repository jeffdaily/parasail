/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>

#include <immintrin.h>

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_avx.h"

#define NEG_INF (INT64_MIN/(int64_t)(2))

#if HAVE_AVX2_MM256_INSERT_EPI64
#define _mm256_insert_epi64_rpl _mm256_insert_epi64
#else
static inline __m256i _mm256_insert_epi64_rpl(__m256i a, int64_t i, int imm) {
    __m256i_64_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif

static inline __m256i _mm256_max_epi64_rpl(__m256i a, __m256i b) {
    __m256i_64_t A;
    __m256i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]>B.v[1]) ? A.v[1] : B.v[1];
    A.v[2] = (A.v[2]>B.v[2]) ? A.v[2] : B.v[2];
    A.v[3] = (A.v[3]>B.v[3]) ? A.v[3] : B.v[3];
    return A.m;
}

#if HAVE_AVX2_MM256_EXTRACT_EPI64
#define _mm256_extract_epi64_rpl _mm256_extract_epi64
#else
static inline int64_t _mm256_extract_epi64_rpl(__m256i a, int imm) {
    __m256i_64_t A;
    A.m = a;
    return A.v[imm];
}
#endif

#define _mm256_srli_si256_rpl(a,imm) _mm256_or_si256(_mm256_slli_si256(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(3,0,0,1)), 16-imm), _mm256_srli_si256(a, imm))

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)


#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vWscore,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[(i+0)*s2Len + (j-0)] = (int64_t)_mm256_extract_epi64_rpl(vWscore, 3);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[(i+1)*s2Len + (j-1)] = (int64_t)_mm256_extract_epi64_rpl(vWscore, 2);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[(i+2)*s2Len + (j-2)] = (int64_t)_mm256_extract_epi64_rpl(vWscore, 1);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[(i+3)*s2Len + (j-3)] = (int64_t)_mm256_extract_epi64_rpl(vWscore, 0);
    }
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_rowcol(
        int *row,
        int *col,
        __m256i vWscore,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (i+0 == s1Len-1 && 0 <= j-0 && j-0 < s2Len) {
        row[j-0] = (int64_t)_mm256_extract_epi64_rpl(vWscore, 3);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int64_t)_mm256_extract_epi64_rpl(vWscore, 3);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int64_t)_mm256_extract_epi64_rpl(vWscore, 2);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int64_t)_mm256_extract_epi64_rpl(vWscore, 2);
    }
    if (i+2 == s1Len-1 && 0 <= j-2 && j-2 < s2Len) {
        row[j-2] = (int64_t)_mm256_extract_epi64_rpl(vWscore, 1);
    }
    if (j-2 == s2Len-1 && 0 <= i+2 && i+2 < s1Len) {
        col[(i+2)] = (int64_t)_mm256_extract_epi64_rpl(vWscore, 1);
    }
    if (i+3 == s1Len-1 && 0 <= j-3 && j-3 < s2Len) {
        row[j-3] = (int64_t)_mm256_extract_epi64_rpl(vWscore, 0);
    }
    if (j-3 == s2Len-1 && 0 <= i+3 && i+3 < s1Len) {
        col[(i+3)] = (int64_t)_mm256_extract_epi64_rpl(vWscore, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_nw_table_diag_avx2_256_64
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_rowcol_diag_avx2_256_64
#else
#define FNAME parasail_nw_diag_avx2_256_64
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    const int32_t N = 4; /* number of values in vector */
    const int32_t PAD = N-1;
    const int32_t PAD2 = PAD*2;
    const int32_t s1Len_PAD = s1Len+PAD;
    const int32_t s2Len_PAD = s2Len+PAD;
    int64_t * const restrict s1 = parasail_memalign_int64_t(32, s1Len+PAD);
    int64_t * const restrict s2B= parasail_memalign_int64_t(32, s2Len+PAD2);
    int64_t * const restrict _tbl_pr = parasail_memalign_int64_t(32, s2Len+PAD2);
    int64_t * const restrict _del_pr = parasail_memalign_int64_t(32, s2Len+PAD2);
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
    __m256i vNegInf = _mm256_set1_epi64x(NEG_INF);
    __m256i vOpen = _mm256_set1_epi64x(open);
    __m256i vGap  = _mm256_set1_epi64x(gap);
    __m256i vOne = _mm256_set1_epi64x(1);
    __m256i vN = _mm256_set1_epi64x(N);
    __m256i vGapN = _mm256_set1_epi64x(gap*N);
    __m256i vNegOne = _mm256_set1_epi64x(-1);
    __m256i vI = _mm256_set_epi64x(0,1,2,3);
    __m256i vJreset = _mm256_set_epi64x(0,-1,-2,-3);
    __m256i vMax = vNegInf;
    __m256i vILimit = _mm256_set1_epi64x(s1Len);
    __m256i vILimit1 = _mm256_sub_epi64(vILimit, vOne);
    __m256i vJLimit = _mm256_set1_epi64x(s2Len);
    __m256i vJLimit1 = _mm256_sub_epi64(vJLimit, vOne);
    __m256i vIBoundary = _mm256_set_epi64x(
            -open-0*gap,
            -open-1*gap,
            -open-2*gap,
            -open-3*gap
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
        __m256i vNscore = vNegInf;
        __m256i vWscore = vNegInf;
        __m256i vIns = vNegInf;
        __m256i vDel = vNegInf;
        __m256i vJ = vJreset;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size*s1[i+2]];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size*s1[i+3]];
        vNscore = _mm256_srli_si256_rpl(vNscore, 8);
        vNscore = _mm256_insert_epi64_rpl(vNscore, tbl_pr[-1], 3);
        vWscore = _mm256_srli_si256_rpl(vWscore, 8);
        vWscore = _mm256_insert_epi64_rpl(vWscore, -open - i*gap, 3);
        tbl_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m256i vMat;
            __m256i vNWscore = vNscore;
            vNscore = _mm256_srli_si256_rpl(vWscore, 8);
            vNscore = _mm256_insert_epi64_rpl(vNscore, tbl_pr[j], 3);
            vDel = _mm256_srli_si256_rpl(vDel, 8);
            vDel = _mm256_insert_epi64_rpl(vDel, del_pr[j], 3);
            vDel = _mm256_max_epi64_rpl(
                    _mm256_sub_epi64(vNscore, vOpen),
                    _mm256_sub_epi64(vDel, vGap));
            vIns = _mm256_max_epi64_rpl(
                    _mm256_sub_epi64(vWscore, vOpen),
                    _mm256_sub_epi64(vIns, vGap));
            vMat = _mm256_set_epi64x(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]]
                    );
            vNWscore = _mm256_add_epi64(vNWscore, vMat);
            vWscore = _mm256_max_epi64_rpl(vNWscore, vIns);
            vWscore = _mm256_max_epi64_rpl(vWscore, vDel);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m256i cond = _mm256_cmpeq_epi64(vJ,vNegOne);
                vWscore = _mm256_blendv_epi8(vWscore, vIBoundary, cond);
                vDel = _mm256_blendv_epi8(vDel, vNegInf, cond);
                vIns = _mm256_blendv_epi8(vIns, vNegInf, cond);
            }
            
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->score_row, result->score_col, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-3] = (int64_t)_mm256_extract_epi64_rpl(vWscore,0);
            del_pr[j-3] = (int64_t)_mm256_extract_epi64_rpl(vDel,0);
            /* as minor diagonal vector passes across table, extract
               last table value at the i,j bound */
            {
                __m256i cond_valid_I = _mm256_cmpeq_epi64(vI, vILimit1);
                __m256i cond_valid_J = _mm256_cmpeq_epi64(vJ, vJLimit1);
                __m256i cond_all = _mm256_and_si256(cond_valid_I, cond_valid_J);
                vMax = _mm256_blendv_epi8(vMax, vWscore, cond_all);
            }
            vJ = _mm256_add_epi64(vJ, vOne);
        }
        vI = _mm256_add_epi64(vI, vN);
        vIBoundary = _mm256_sub_epi64(vIBoundary, vGapN);
    }

    /* max in vMax */
    for (i=0; i<N; ++i) {
        int64_t value;
        value = (int64_t) _mm256_extract_epi64_rpl(vMax, 3);
        if (value > score) {
            score = value;
        }
        vMax = _mm256_slli_si256_rpl(vMax, 8);
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


