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

#define NEG_INF (INT16_MIN/(int16_t)(2))

#define _mm256_cmplt_epi16_rpl(a,b) _mm256_cmpgt_epi16(b,a)

#if HAVE_AVX2_MM256_INSERT_EPI16
#define _mm256_insert_epi16_rpl _mm256_insert_epi16
#else
static inline __m256i _mm256_insert_epi16_rpl(__m256i a, int16_t i, int imm) {
    __m256i_16_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI16
#define _mm256_extract_epi16_rpl _mm256_extract_epi16
#else
static inline int16_t _mm256_extract_epi16_rpl(__m256i a, int imm) {
    __m256i_16_t A;
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
        array[(i+0)*s2Len + (j-0)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 15);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[(i+1)*s2Len + (j-1)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 14);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[(i+2)*s2Len + (j-2)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 13);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[(i+3)*s2Len + (j-3)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 12);
    }
    if (0 <= i+4 && i+4 < s1Len && 0 <= j-4 && j-4 < s2Len) {
        array[(i+4)*s2Len + (j-4)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 11);
    }
    if (0 <= i+5 && i+5 < s1Len && 0 <= j-5 && j-5 < s2Len) {
        array[(i+5)*s2Len + (j-5)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 10);
    }
    if (0 <= i+6 && i+6 < s1Len && 0 <= j-6 && j-6 < s2Len) {
        array[(i+6)*s2Len + (j-6)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 9);
    }
    if (0 <= i+7 && i+7 < s1Len && 0 <= j-7 && j-7 < s2Len) {
        array[(i+7)*s2Len + (j-7)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 8);
    }
    if (0 <= i+8 && i+8 < s1Len && 0 <= j-8 && j-8 < s2Len) {
        array[(i+8)*s2Len + (j-8)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 7);
    }
    if (0 <= i+9 && i+9 < s1Len && 0 <= j-9 && j-9 < s2Len) {
        array[(i+9)*s2Len + (j-9)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 6);
    }
    if (0 <= i+10 && i+10 < s1Len && 0 <= j-10 && j-10 < s2Len) {
        array[(i+10)*s2Len + (j-10)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 5);
    }
    if (0 <= i+11 && i+11 < s1Len && 0 <= j-11 && j-11 < s2Len) {
        array[(i+11)*s2Len + (j-11)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 4);
    }
    if (0 <= i+12 && i+12 < s1Len && 0 <= j-12 && j-12 < s2Len) {
        array[(i+12)*s2Len + (j-12)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 3);
    }
    if (0 <= i+13 && i+13 < s1Len && 0 <= j-13 && j-13 < s2Len) {
        array[(i+13)*s2Len + (j-13)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 2);
    }
    if (0 <= i+14 && i+14 < s1Len && 0 <= j-14 && j-14 < s2Len) {
        array[(i+14)*s2Len + (j-14)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 1);
    }
    if (0 <= i+15 && i+15 < s1Len && 0 <= j-15 && j-15 < s2Len) {
        array[(i+15)*s2Len + (j-15)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 0);
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
        row[j-0] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 15);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 15);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 14);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 14);
    }
    if (i+2 == s1Len-1 && 0 <= j-2 && j-2 < s2Len) {
        row[j-2] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 13);
    }
    if (j-2 == s2Len-1 && 0 <= i+2 && i+2 < s1Len) {
        col[(i+2)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 13);
    }
    if (i+3 == s1Len-1 && 0 <= j-3 && j-3 < s2Len) {
        row[j-3] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 12);
    }
    if (j-3 == s2Len-1 && 0 <= i+3 && i+3 < s1Len) {
        col[(i+3)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 12);
    }
    if (i+4 == s1Len-1 && 0 <= j-4 && j-4 < s2Len) {
        row[j-4] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 11);
    }
    if (j-4 == s2Len-1 && 0 <= i+4 && i+4 < s1Len) {
        col[(i+4)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 11);
    }
    if (i+5 == s1Len-1 && 0 <= j-5 && j-5 < s2Len) {
        row[j-5] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 10);
    }
    if (j-5 == s2Len-1 && 0 <= i+5 && i+5 < s1Len) {
        col[(i+5)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 10);
    }
    if (i+6 == s1Len-1 && 0 <= j-6 && j-6 < s2Len) {
        row[j-6] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 9);
    }
    if (j-6 == s2Len-1 && 0 <= i+6 && i+6 < s1Len) {
        col[(i+6)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 9);
    }
    if (i+7 == s1Len-1 && 0 <= j-7 && j-7 < s2Len) {
        row[j-7] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 8);
    }
    if (j-7 == s2Len-1 && 0 <= i+7 && i+7 < s1Len) {
        col[(i+7)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 8);
    }
    if (i+8 == s1Len-1 && 0 <= j-8 && j-8 < s2Len) {
        row[j-8] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 7);
    }
    if (j-8 == s2Len-1 && 0 <= i+8 && i+8 < s1Len) {
        col[(i+8)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 7);
    }
    if (i+9 == s1Len-1 && 0 <= j-9 && j-9 < s2Len) {
        row[j-9] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 6);
    }
    if (j-9 == s2Len-1 && 0 <= i+9 && i+9 < s1Len) {
        col[(i+9)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 6);
    }
    if (i+10 == s1Len-1 && 0 <= j-10 && j-10 < s2Len) {
        row[j-10] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 5);
    }
    if (j-10 == s2Len-1 && 0 <= i+10 && i+10 < s1Len) {
        col[(i+10)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 5);
    }
    if (i+11 == s1Len-1 && 0 <= j-11 && j-11 < s2Len) {
        row[j-11] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 4);
    }
    if (j-11 == s2Len-1 && 0 <= i+11 && i+11 < s1Len) {
        col[(i+11)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 4);
    }
    if (i+12 == s1Len-1 && 0 <= j-12 && j-12 < s2Len) {
        row[j-12] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 3);
    }
    if (j-12 == s2Len-1 && 0 <= i+12 && i+12 < s1Len) {
        col[(i+12)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 3);
    }
    if (i+13 == s1Len-1 && 0 <= j-13 && j-13 < s2Len) {
        row[j-13] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 2);
    }
    if (j-13 == s2Len-1 && 0 <= i+13 && i+13 < s1Len) {
        col[(i+13)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 2);
    }
    if (i+14 == s1Len-1 && 0 <= j-14 && j-14 < s2Len) {
        row[j-14] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 1);
    }
    if (j-14 == s2Len-1 && 0 <= i+14 && i+14 < s1Len) {
        col[(i+14)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 1);
    }
    if (i+15 == s1Len-1 && 0 <= j-15 && j-15 < s2Len) {
        row[j-15] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 0);
    }
    if (j-15 == s2Len-1 && 0 <= i+15 && i+15 < s1Len) {
        col[(i+15)] = (int16_t)_mm256_extract_epi16_rpl(vWscore, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_nw_stats_table_diag_avx2_256_16
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_stats_rowcol_diag_avx2_256_16
#else
#define FNAME parasail_nw_stats_diag_avx2_256_16
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    const int32_t N = 16; /* number of values in vector */
    const int32_t PAD = N-1;
    const int32_t PAD2 = PAD*2;
    const int32_t s1Len_PAD = s1Len+PAD;
    const int32_t s2Len_PAD = s2Len+PAD;
    int16_t * const restrict s1      = parasail_memalign_int16_t(32, s1Len+PAD);
    int16_t * const restrict s2B     = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict _tbl_pr = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict _del_pr = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict _mch_pr = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict _sim_pr = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict _len_pr = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int16_t * const restrict tbl_pr = _tbl_pr+PAD;
    int16_t * const restrict del_pr = _del_pr+PAD;
    int16_t * const restrict mch_pr = _mch_pr+PAD;
    int16_t * const restrict sim_pr = _sim_pr+PAD;
    int16_t * const restrict len_pr = _len_pr+PAD;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif
    int32_t i = 0;
    int32_t j = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int16_t score = NEG_INF;
    int16_t matches = NEG_INF;
    int16_t similar = NEG_INF;
    int16_t length = NEG_INF;
    
    __m256i vNegInf = _mm256_set1_epi16(NEG_INF);
    __m256i vOpen = _mm256_set1_epi16(open);
    __m256i vGap  = _mm256_set1_epi16(gap);
    __m256i vZero = _mm256_set1_epi16(0);
    __m256i vOne = _mm256_set1_epi16(1);
    __m256i vN = _mm256_set1_epi16(N);
    __m256i vGapN = _mm256_set1_epi16(gap*N);
    __m256i vNegOne = _mm256_set1_epi16(-1);
    __m256i vI = _mm256_set_epi16(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
    __m256i vJreset = _mm256_set_epi16(0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15);
    __m256i vMaxScore = vNegInf;
    __m256i vMaxMatch = vNegInf;
    __m256i vMaxSimilar = vNegInf;
    __m256i vMaxLength = vNegInf;
    __m256i vILimit = _mm256_set1_epi16(s1Len);
    __m256i vILimit1 = _mm256_sub_epi16(vILimit, vOne);
    __m256i vJLimit = _mm256_set1_epi16(s2Len);
    __m256i vJLimit1 = _mm256_sub_epi16(vJLimit, vOne);
    __m256i vIBoundary = _mm256_set_epi16(
            -open-0*gap,
            -open-1*gap,
            -open-2*gap,
            -open-3*gap,
            -open-4*gap,
            -open-5*gap,
            -open-6*gap,
            -open-7*gap,
            -open-8*gap,
            -open-9*gap,
            -open-10*gap,
            -open-11*gap,
            -open-12*gap,
            -open-13*gap,
            -open-14*gap,
            -open-15*gap);

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
        mch_pr[j] = 0;
        sim_pr[j] = 0;
        len_pr[j] = 0;
    }
    /* pad front of stored row values */
    for (j=-PAD; j<0; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        sim_pr[j] = 0;
        len_pr[j] = 0;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        sim_pr[j] = 0;
        len_pr[j] = 0;
    }
    tbl_pr[-1] = 0; /* upper left corner */

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        __m256i vNscore = vNegInf;
        __m256i vNmatch = vZero;
        __m256i vNsimilar = vZero;
        __m256i vNlength = vZero;
        __m256i vWscore = vNegInf;
        __m256i vWmatch = vZero;
        __m256i vWsimilar = vZero;
        __m256i vWlength = vZero;
        __m256i vIns = vNegInf;
        __m256i vDel = vNegInf;
        __m256i vJ = vJreset;
        __m256i vs1 = _mm256_set_epi16(
                s1[i+0],
                s1[i+1],
                s1[i+2],
                s1[i+3],
                s1[i+4],
                s1[i+5],
                s1[i+6],
                s1[i+7],
                s1[i+8],
                s1[i+9],
                s1[i+10],
                s1[i+11],
                s1[i+12],
                s1[i+13],
                s1[i+14],
                s1[i+15]);
        __m256i vs2 = vNegInf;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size*s1[i+2]];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size*s1[i+3]];
        const int * const restrict matrow4 = &matrix->matrix[matrix->size*s1[i+4]];
        const int * const restrict matrow5 = &matrix->matrix[matrix->size*s1[i+5]];
        const int * const restrict matrow6 = &matrix->matrix[matrix->size*s1[i+6]];
        const int * const restrict matrow7 = &matrix->matrix[matrix->size*s1[i+7]];
        const int * const restrict matrow8 = &matrix->matrix[matrix->size*s1[i+8]];
        const int * const restrict matrow9 = &matrix->matrix[matrix->size*s1[i+9]];
        const int * const restrict matrow10 = &matrix->matrix[matrix->size*s1[i+10]];
        const int * const restrict matrow11 = &matrix->matrix[matrix->size*s1[i+11]];
        const int * const restrict matrow12 = &matrix->matrix[matrix->size*s1[i+12]];
        const int * const restrict matrow13 = &matrix->matrix[matrix->size*s1[i+13]];
        const int * const restrict matrow14 = &matrix->matrix[matrix->size*s1[i+14]];
        const int * const restrict matrow15 = &matrix->matrix[matrix->size*s1[i+15]];
        vNscore = _mm256_srli_si256_rpl(vNscore, 2);
        vNscore = _mm256_insert_epi16_rpl(vNscore, tbl_pr[-1], 15);
        vNmatch = _mm256_srli_si256_rpl(vNmatch, 2);
        vNmatch = _mm256_insert_epi16_rpl(vNmatch, 0, 15);
        vNsimilar = _mm256_srli_si256_rpl(vNsimilar, 2);
        vNsimilar = _mm256_insert_epi16_rpl(vNsimilar, 0, 15);
        vNlength = _mm256_srli_si256_rpl(vNlength, 2);
        vNlength = _mm256_insert_epi16_rpl(vNlength, 0, 15);
        vWscore = _mm256_srli_si256_rpl(vWscore, 2);
        vWscore = _mm256_insert_epi16_rpl(vWscore, -open - i*gap, 15);
        vWmatch = _mm256_srli_si256_rpl(vWmatch, 2);
        vWmatch = _mm256_insert_epi16_rpl(vWmatch, 0, 15);
        vWsimilar = _mm256_srli_si256_rpl(vWsimilar, 2);
        vWsimilar = _mm256_insert_epi16_rpl(vWsimilar, 0, 15);
        vWlength = _mm256_srli_si256_rpl(vWlength, 2);
        vWlength = _mm256_insert_epi16_rpl(vWlength, 0, 15);
        tbl_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m256i vMat;
            __m256i vNWscore = vNscore;
            __m256i vNWmatch = vNmatch;
            __m256i vNWsimilar = vNsimilar;
            __m256i vNWlength = vNlength;
            vNscore = _mm256_srli_si256_rpl(vWscore, 2);
            vNscore = _mm256_insert_epi16_rpl(vNscore, tbl_pr[j], 15);
            vNmatch = _mm256_srli_si256_rpl(vWmatch, 2);
            vNmatch = _mm256_insert_epi16_rpl(vNmatch, mch_pr[j], 15);
            vNsimilar = _mm256_srli_si256_rpl(vWsimilar, 2);
            vNsimilar = _mm256_insert_epi16_rpl(vNsimilar, sim_pr[j], 15);
            vNlength = _mm256_srli_si256_rpl(vWlength, 2);
            vNlength = _mm256_insert_epi16_rpl(vNlength, len_pr[j], 15);
            vDel = _mm256_srli_si256_rpl(vDel, 2);
            vDel = _mm256_insert_epi16_rpl(vDel, del_pr[j], 15);
            vDel = _mm256_max_epi16(
                    _mm256_sub_epi16(vNscore, vOpen),
                    _mm256_sub_epi16(vDel, vGap));
            vIns = _mm256_max_epi16(
                    _mm256_sub_epi16(vWscore, vOpen),
                    _mm256_sub_epi16(vIns, vGap));
            vs2 = _mm256_srli_si256_rpl(vs2, 2);
            vs2 = _mm256_insert_epi16_rpl(vs2, s2[j], 15);
            vMat = _mm256_set_epi16(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]],
                    matrow4[s2[j-4]],
                    matrow5[s2[j-5]],
                    matrow6[s2[j-6]],
                    matrow7[s2[j-7]],
                    matrow8[s2[j-8]],
                    matrow9[s2[j-9]],
                    matrow10[s2[j-10]],
                    matrow11[s2[j-11]],
                    matrow12[s2[j-12]],
                    matrow13[s2[j-13]],
                    matrow14[s2[j-14]],
                    matrow15[s2[j-15]]
                    );
            vNWscore = _mm256_add_epi16(vNWscore, vMat);
            vWscore = _mm256_max_epi16(vNWscore, vIns);
            vWscore = _mm256_max_epi16(vWscore, vDel);
            /* conditional block */
            {
                __m256i case1not;
                __m256i case2not;
                __m256i case2;
                __m256i case3;
                __m256i vCmatch;
                __m256i vCsimilar;
                __m256i vClength;
                case1not = _mm256_or_si256(
                        _mm256_cmplt_epi16_rpl(vNWscore,vDel),
                        _mm256_cmplt_epi16_rpl(vNWscore,vIns));
                case2not = _mm256_cmplt_epi16_rpl(vDel,vIns);
                case2 = _mm256_andnot_si256(case2not,case1not);
                case3 = _mm256_and_si256(case1not,case2not);
                vCmatch = _mm256_andnot_si256(case1not,
                        _mm256_add_epi16(vNWmatch, _mm256_and_si256(
                                _mm256_cmpeq_epi16(vs1,vs2),vOne)));
                vCmatch = _mm256_or_si256(vCmatch, _mm256_and_si256(case2, vNmatch));
                vCmatch = _mm256_or_si256(vCmatch, _mm256_and_si256(case3, vWmatch));
                vCsimilar = _mm256_andnot_si256(case1not,
                        _mm256_add_epi16(vNWsimilar, _mm256_and_si256(
                                _mm256_cmpgt_epi16(vMat,vZero),vOne)));
                vCsimilar = _mm256_or_si256(vCsimilar, _mm256_and_si256(case2, vNsimilar));
                vCsimilar = _mm256_or_si256(vCsimilar, _mm256_and_si256(case3, vWsimilar));
                vClength= _mm256_andnot_si256(case1not,
                        _mm256_add_epi16(vNWlength, vOne));
                vClength= _mm256_or_si256(vClength,_mm256_and_si256(case2,
                            _mm256_add_epi16(vNlength, vOne)));
                vClength= _mm256_or_si256(vClength,_mm256_and_si256(case3,
                            _mm256_add_epi16(vWlength, vOne)));
                vWmatch = vCmatch;
                vWsimilar = vCsimilar;
                vWlength = vClength;
            }

            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m256i cond = _mm256_cmpeq_epi16(vJ,vNegOne);
                vWscore = _mm256_blendv_epi8(vWscore, vIBoundary, cond);
                vWmatch = _mm256_andnot_si256(cond, vWmatch);
                vWsimilar = _mm256_andnot_si256(cond, vWsimilar);
                vWlength = _mm256_andnot_si256(cond, vWlength);
                vDel = _mm256_blendv_epi8(vDel, vNegInf, cond);
                vIns = _mm256_blendv_epi8(vIns, vNegInf, cond);
            }
            
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vWscore, i, s1Len, j, s2Len);
            arr_store_si256(result->matches_table, vWmatch, i, s1Len, j, s2Len);
            arr_store_si256(result->similar_table, vWsimilar, i, s1Len, j, s2Len);
            arr_store_si256(result->length_table, vWlength, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->score_row, result->score_col, vWscore, i, s1Len, j, s2Len);
            arr_store_rowcol(result->matches_row, result->matches_col, vWmatch, i, s1Len, j, s2Len);
            arr_store_rowcol(result->similar_row, result->similar_col, vWsimilar, i, s1Len, j, s2Len);
            arr_store_rowcol(result->length_row, result->length_col, vWlength, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-15] = (int16_t)_mm256_extract_epi16_rpl(vWscore,0);
            mch_pr[j-15] = (int16_t)_mm256_extract_epi16_rpl(vWmatch,0);
            sim_pr[j-15] = (int16_t)_mm256_extract_epi16_rpl(vWsimilar,0);
            len_pr[j-15] = (int16_t)_mm256_extract_epi16_rpl(vWlength,0);
            del_pr[j-15] = (int16_t)_mm256_extract_epi16_rpl(vDel,0);
            /* as minor diagonal vector passes across table, extract
               last table value at the i,j bound */
            {
                __m256i cond_valid_I = _mm256_cmpeq_epi16(vI, vILimit1);
                __m256i cond_valid_J = _mm256_cmpeq_epi16(vJ, vJLimit1);
                __m256i cond_all = _mm256_and_si256(cond_valid_I, cond_valid_J);
                vMaxScore = _mm256_blendv_epi8(vMaxScore, vWscore, cond_all);
                vMaxMatch = _mm256_blendv_epi8(vMaxMatch, vWmatch, cond_all);
                vMaxSimilar = _mm256_blendv_epi8(vMaxSimilar, vWsimilar, cond_all);
                vMaxLength = _mm256_blendv_epi8(vMaxLength, vWlength, cond_all);
            }
            vJ = _mm256_add_epi16(vJ, vOne);
        }
        vI = _mm256_add_epi16(vI, vN);
        vIBoundary = _mm256_sub_epi16(vIBoundary, vGapN);
    }

    /* max in vMaxScore */
    for (i=0; i<N; ++i) {
        int16_t value;
        value = (int16_t) _mm256_extract_epi16_rpl(vMaxScore, 15);
        if (value > score) {
            score = value;
            matches = (int16_t) _mm256_extract_epi16_rpl(vMaxMatch, 15);
            similar = (int16_t) _mm256_extract_epi16_rpl(vMaxSimilar, 15);
            length= (int16_t) _mm256_extract_epi16_rpl(vMaxLength, 15);
        }
        vMaxScore = _mm256_slli_si256_rpl(vMaxScore, 2);
        vMaxMatch = _mm256_slli_si256_rpl(vMaxMatch, 2);
        vMaxSimilar = _mm256_slli_si256_rpl(vMaxSimilar, 2);
        vMaxLength = _mm256_slli_si256_rpl(vMaxLength, 2);
    }

    

    result->score = score;
    result->matches = matches;
    result->similar = similar;
    result->length = length;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(_len_pr);
    parasail_free(_sim_pr);
    parasail_free(_mch_pr);
    parasail_free(_del_pr);
    parasail_free(_tbl_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}


