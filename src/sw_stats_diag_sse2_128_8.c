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

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_sse.h"

#define NEG_INF INT8_MIN

static inline __m128i _mm_blendv_epi8_rpl(__m128i a, __m128i b, __m128i mask) {
    a = _mm_andnot_si128(mask, a);
    a = _mm_or_si128(a, _mm_and_si128(mask, b));
    return a;
}

static inline __m128i _mm_insert_epi8_rpl(__m128i a, int8_t i, const int imm) {
    __m128i_8_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}

static inline __m128i _mm_max_epi8_rpl(__m128i a, __m128i b) {
    __m128i mask = _mm_cmpgt_epi8(a, b);
    a = _mm_and_si128(a, mask);
    b = _mm_andnot_si128(mask, b);
    return _mm_or_si128(a, b);
}

static inline int8_t _mm_extract_epi8_rpl(__m128i a, const int imm) {
    __m128i_8_t A;
    A.m = a;
    return A.v[imm];
}

static inline __m128i _mm_min_epi8_rpl(__m128i a, __m128i b) {
    __m128i mask = _mm_cmpgt_epi8(b, a);
    a = _mm_and_si128(a, mask);
    b = _mm_andnot_si128(mask, b);
    return _mm_or_si128(a, b);
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
        array[(i+0)*s2Len + (j-0)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 15);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[(i+1)*s2Len + (j-1)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 14);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[(i+2)*s2Len + (j-2)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 13);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[(i+3)*s2Len + (j-3)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 12);
    }
    if (0 <= i+4 && i+4 < s1Len && 0 <= j-4 && j-4 < s2Len) {
        array[(i+4)*s2Len + (j-4)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 11);
    }
    if (0 <= i+5 && i+5 < s1Len && 0 <= j-5 && j-5 < s2Len) {
        array[(i+5)*s2Len + (j-5)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 10);
    }
    if (0 <= i+6 && i+6 < s1Len && 0 <= j-6 && j-6 < s2Len) {
        array[(i+6)*s2Len + (j-6)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 9);
    }
    if (0 <= i+7 && i+7 < s1Len && 0 <= j-7 && j-7 < s2Len) {
        array[(i+7)*s2Len + (j-7)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 8);
    }
    if (0 <= i+8 && i+8 < s1Len && 0 <= j-8 && j-8 < s2Len) {
        array[(i+8)*s2Len + (j-8)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 7);
    }
    if (0 <= i+9 && i+9 < s1Len && 0 <= j-9 && j-9 < s2Len) {
        array[(i+9)*s2Len + (j-9)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 6);
    }
    if (0 <= i+10 && i+10 < s1Len && 0 <= j-10 && j-10 < s2Len) {
        array[(i+10)*s2Len + (j-10)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 5);
    }
    if (0 <= i+11 && i+11 < s1Len && 0 <= j-11 && j-11 < s2Len) {
        array[(i+11)*s2Len + (j-11)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 4);
    }
    if (0 <= i+12 && i+12 < s1Len && 0 <= j-12 && j-12 < s2Len) {
        array[(i+12)*s2Len + (j-12)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 3);
    }
    if (0 <= i+13 && i+13 < s1Len && 0 <= j-13 && j-13 < s2Len) {
        array[(i+13)*s2Len + (j-13)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 2);
    }
    if (0 <= i+14 && i+14 < s1Len && 0 <= j-14 && j-14 < s2Len) {
        array[(i+14)*s2Len + (j-14)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 1);
    }
    if (0 <= i+15 && i+15 < s1Len && 0 <= j-15 && j-15 < s2Len) {
        array[(i+15)*s2Len + (j-15)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 0);
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
        row[j-0] = (int8_t)_mm_extract_epi8_rpl(vWscore, 15);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 15);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int8_t)_mm_extract_epi8_rpl(vWscore, 14);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 14);
    }
    if (i+2 == s1Len-1 && 0 <= j-2 && j-2 < s2Len) {
        row[j-2] = (int8_t)_mm_extract_epi8_rpl(vWscore, 13);
    }
    if (j-2 == s2Len-1 && 0 <= i+2 && i+2 < s1Len) {
        col[(i+2)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 13);
    }
    if (i+3 == s1Len-1 && 0 <= j-3 && j-3 < s2Len) {
        row[j-3] = (int8_t)_mm_extract_epi8_rpl(vWscore, 12);
    }
    if (j-3 == s2Len-1 && 0 <= i+3 && i+3 < s1Len) {
        col[(i+3)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 12);
    }
    if (i+4 == s1Len-1 && 0 <= j-4 && j-4 < s2Len) {
        row[j-4] = (int8_t)_mm_extract_epi8_rpl(vWscore, 11);
    }
    if (j-4 == s2Len-1 && 0 <= i+4 && i+4 < s1Len) {
        col[(i+4)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 11);
    }
    if (i+5 == s1Len-1 && 0 <= j-5 && j-5 < s2Len) {
        row[j-5] = (int8_t)_mm_extract_epi8_rpl(vWscore, 10);
    }
    if (j-5 == s2Len-1 && 0 <= i+5 && i+5 < s1Len) {
        col[(i+5)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 10);
    }
    if (i+6 == s1Len-1 && 0 <= j-6 && j-6 < s2Len) {
        row[j-6] = (int8_t)_mm_extract_epi8_rpl(vWscore, 9);
    }
    if (j-6 == s2Len-1 && 0 <= i+6 && i+6 < s1Len) {
        col[(i+6)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 9);
    }
    if (i+7 == s1Len-1 && 0 <= j-7 && j-7 < s2Len) {
        row[j-7] = (int8_t)_mm_extract_epi8_rpl(vWscore, 8);
    }
    if (j-7 == s2Len-1 && 0 <= i+7 && i+7 < s1Len) {
        col[(i+7)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 8);
    }
    if (i+8 == s1Len-1 && 0 <= j-8 && j-8 < s2Len) {
        row[j-8] = (int8_t)_mm_extract_epi8_rpl(vWscore, 7);
    }
    if (j-8 == s2Len-1 && 0 <= i+8 && i+8 < s1Len) {
        col[(i+8)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 7);
    }
    if (i+9 == s1Len-1 && 0 <= j-9 && j-9 < s2Len) {
        row[j-9] = (int8_t)_mm_extract_epi8_rpl(vWscore, 6);
    }
    if (j-9 == s2Len-1 && 0 <= i+9 && i+9 < s1Len) {
        col[(i+9)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 6);
    }
    if (i+10 == s1Len-1 && 0 <= j-10 && j-10 < s2Len) {
        row[j-10] = (int8_t)_mm_extract_epi8_rpl(vWscore, 5);
    }
    if (j-10 == s2Len-1 && 0 <= i+10 && i+10 < s1Len) {
        col[(i+10)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 5);
    }
    if (i+11 == s1Len-1 && 0 <= j-11 && j-11 < s2Len) {
        row[j-11] = (int8_t)_mm_extract_epi8_rpl(vWscore, 4);
    }
    if (j-11 == s2Len-1 && 0 <= i+11 && i+11 < s1Len) {
        col[(i+11)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 4);
    }
    if (i+12 == s1Len-1 && 0 <= j-12 && j-12 < s2Len) {
        row[j-12] = (int8_t)_mm_extract_epi8_rpl(vWscore, 3);
    }
    if (j-12 == s2Len-1 && 0 <= i+12 && i+12 < s1Len) {
        col[(i+12)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 3);
    }
    if (i+13 == s1Len-1 && 0 <= j-13 && j-13 < s2Len) {
        row[j-13] = (int8_t)_mm_extract_epi8_rpl(vWscore, 2);
    }
    if (j-13 == s2Len-1 && 0 <= i+13 && i+13 < s1Len) {
        col[(i+13)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 2);
    }
    if (i+14 == s1Len-1 && 0 <= j-14 && j-14 < s2Len) {
        row[j-14] = (int8_t)_mm_extract_epi8_rpl(vWscore, 1);
    }
    if (j-14 == s2Len-1 && 0 <= i+14 && i+14 < s1Len) {
        col[(i+14)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 1);
    }
    if (i+15 == s1Len-1 && 0 <= j-15 && j-15 < s2Len) {
        row[j-15] = (int8_t)_mm_extract_epi8_rpl(vWscore, 0);
    }
    if (j-15 == s2Len-1 && 0 <= i+15 && i+15 < s1Len) {
        col[(i+15)] = (int8_t)_mm_extract_epi8_rpl(vWscore, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_stats_table_diag_sse2_128_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_stats_rowcol_diag_sse2_128_8
#else
#define FNAME parasail_sw_stats_diag_sse2_128_8
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
    int8_t * const restrict s1      = parasail_memalign_int8_t(16, s1Len+PAD);
    int8_t * const restrict s2B     = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _tbl_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _del_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _mch_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _sim_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict _len_pr = parasail_memalign_int8_t(16, s2Len+PAD2);
    int8_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int8_t * const restrict tbl_pr = _tbl_pr+PAD;
    int8_t * const restrict del_pr = _del_pr+PAD;
    int8_t * const restrict mch_pr = _mch_pr+PAD;
    int8_t * const restrict sim_pr = _sim_pr+PAD;
    int8_t * const restrict len_pr = _len_pr+PAD;
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
    int8_t score = NEG_INF;
    int8_t matches = NEG_INF;
    int8_t similar = NEG_INF;
    int8_t length = NEG_INF;
    __m128i vNegLimit = _mm_set1_epi8(INT8_MIN);
    __m128i vPosLimit = _mm_set1_epi8(INT8_MAX);
    __m128i vSaturationCheckMin = vPosLimit;
    __m128i vSaturationCheckMax = vNegLimit;
    __m128i vNegInf = _mm_set1_epi8(NEG_INF);
    __m128i vNegInf0 = _mm_srli_si128(vNegInf, 1); /* shift in a 0 */
    __m128i vOpen = _mm_set1_epi8(open);
    __m128i vGap  = _mm_set1_epi8(gap);
    __m128i vZero = _mm_set1_epi8(0);
    __m128i vOne = _mm_set1_epi8(1);
    __m128i vOne16 = _mm_set1_epi16(1);
    __m128i vNegOne16 = _mm_set1_epi16(-1);
    __m128i vN16 = _mm_set1_epi16(N);
    __m128i vILo16 = _mm_set_epi16(8,9,10,11,12,13,14,15);
    __m128i vIHi16 = _mm_set_epi16(0,1,2,3,4,5,6,7);
    __m128i vJresetLo16 = _mm_set_epi16(-8,-9,-10,-11,-12,-13,-14,-15);
    __m128i vJresetHi16 = _mm_set_epi16(0,-1,-2,-3,-4,-5,-6,-7);
    __m128i vMaxScore = vNegInf;
    __m128i vMaxMatch = vNegInf;
    __m128i vMaxSimilar = vNegInf;
    __m128i vMaxLength = vNegInf;
    __m128i vEndILo = vNegInf;
    __m128i vEndIHi = vNegInf;
    __m128i vEndJLo = vNegInf;
    __m128i vEndJHi = vNegInf;
    __m128i vILimit16 = _mm_set1_epi16(s1Len);
    __m128i vJLimit16 = _mm_set1_epi16(s2Len);

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
        tbl_pr[j] = 0;
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
        __m128i vNscore = vNegInf0;
        __m128i vNmatch = vZero;
        __m128i vNsimilar = vZero;
        __m128i vNlength = vZero;
        __m128i vWscore = vNegInf0;
        __m128i vWmatch = vZero;
        __m128i vWsimilar = vZero;
        __m128i vWlength = vZero;
        __m128i vIns = vNegInf;
        __m128i vDel = vNegInf;
        __m128i vJLo16 = vJresetLo16;
        __m128i vJHi16 = vJresetHi16;
        __m128i vs1 = _mm_set_epi8(
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
        __m128i vs2 = vNegInf;
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
        __m128i vIltLimit = _mm_packs_epi16(
                    _mm_cmplt_epi16(vILo16, vILimit16),
                    _mm_cmplt_epi16(vIHi16, vILimit16));
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            __m128i vNWmatch = vNmatch;
            __m128i vNWsimilar = vNsimilar;
            __m128i vNWlength = vNlength;
            vNscore = _mm_srli_si128(vWscore, 1);
            vNscore = _mm_insert_epi8_rpl(vNscore, tbl_pr[j], 15);
            vNmatch = _mm_srli_si128(vWmatch, 1);
            vNmatch = _mm_insert_epi8_rpl(vNmatch, mch_pr[j], 15);
            vNsimilar = _mm_srli_si128(vWsimilar, 1);
            vNsimilar = _mm_insert_epi8_rpl(vNsimilar, sim_pr[j], 15);
            vNlength = _mm_srli_si128(vWlength, 1);
            vNlength = _mm_insert_epi8_rpl(vNlength, len_pr[j], 15);
            vDel = _mm_srli_si128(vDel, 1);
            vDel = _mm_insert_epi8_rpl(vDel, del_pr[j], 15);
            vDel = _mm_max_epi8_rpl(
                    _mm_subs_epi8(vNscore, vOpen),
                    _mm_subs_epi8(vDel, vGap));
            vIns = _mm_max_epi8_rpl(
                    _mm_subs_epi8(vWscore, vOpen),
                    _mm_subs_epi8(vIns, vGap));
            vs2 = _mm_srli_si128(vs2, 1);
            vs2 = _mm_insert_epi8_rpl(vs2, s2[j], 15);
            vMat = _mm_set_epi8(
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
            vNWscore = _mm_adds_epi8(vNWscore, vMat);
            vWscore = _mm_max_epi8_rpl(vNWscore, vIns);
            vWscore = _mm_max_epi8_rpl(vWscore, vDel);
            vWscore = _mm_max_epi8_rpl(vWscore, vZero);
            /* conditional block */
            {
                __m128i case1not;
                __m128i case2not;
                __m128i case2;
                __m128i case3;
                __m128i vCmatch;
                __m128i vCsimilar;
                __m128i vClength;
                case1not = _mm_or_si128(
                        _mm_cmplt_epi8(vNWscore,vDel),
                        _mm_cmplt_epi8(vNWscore,vIns));
                case2not = _mm_cmplt_epi8(vDel,vIns);
                case2 = _mm_andnot_si128(case2not,case1not);
                case3 = _mm_and_si128(case1not,case2not);
                vCmatch = _mm_andnot_si128(case1not,
                        _mm_adds_epi8(vNWmatch, _mm_and_si128(
                                _mm_cmpeq_epi8(vs1,vs2),vOne)));
                vCmatch = _mm_or_si128(vCmatch, _mm_and_si128(case2, vNmatch));
                vCmatch = _mm_or_si128(vCmatch, _mm_and_si128(case3, vWmatch));
                vCsimilar = _mm_andnot_si128(case1not,
                        _mm_adds_epi8(vNWsimilar, _mm_and_si128(
                                _mm_cmpgt_epi8(vMat,vZero),vOne)));
                vCsimilar = _mm_or_si128(vCsimilar, _mm_and_si128(case2, vNsimilar));
                vCsimilar = _mm_or_si128(vCsimilar, _mm_and_si128(case3, vWsimilar));
                vClength= _mm_andnot_si128(case1not,
                        _mm_adds_epi8(vNWlength, vOne));
                vClength= _mm_or_si128(vClength,_mm_and_si128(case2,
                            _mm_adds_epi8(vNlength, vOne)));
                vClength= _mm_or_si128(vClength,_mm_and_si128(case3,
                            _mm_adds_epi8(vWlength, vOne)));
                vWmatch = vCmatch;
                vWsimilar = vCsimilar;
                vWlength = vClength;
            }
            {
                __m128i cond_zero = _mm_cmpeq_epi8(vWscore, vZero);
                vWmatch = _mm_andnot_si128(cond_zero, vWmatch);
                vWsimilar = _mm_andnot_si128(cond_zero, vWsimilar);
                vWlength = _mm_andnot_si128(cond_zero, vWlength);
            }

            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m128i cond = _mm_packs_epi16(
                        _mm_cmpeq_epi16(vJLo16,vNegOne16),
                        _mm_cmpeq_epi16(vJHi16,vNegOne16));
                vWscore = _mm_andnot_si128(cond, vWscore);
                vWmatch = _mm_andnot_si128(cond, vWmatch);
                vWsimilar = _mm_andnot_si128(cond, vWsimilar);
                vWlength = _mm_andnot_si128(cond, vWlength);
                vDel = _mm_blendv_epi8_rpl(vDel, vNegInf, cond);
                vIns = _mm_blendv_epi8_rpl(vIns, vNegInf, cond);
            }

            /* check for saturation */
            {
                vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vWscore);
                vSaturationCheckMin = _mm_min_epi8_rpl(vSaturationCheckMin, vWscore);
                vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vWmatch);
                vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vWsimilar);
                vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vWlength);
            }
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vWscore, i, s1Len, j, s2Len);
            arr_store_si128(result->matches_table, vWmatch, i, s1Len, j, s2Len);
            arr_store_si128(result->similar_table, vWsimilar, i, s1Len, j, s2Len);
            arr_store_si128(result->length_table, vWlength, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->score_row, result->score_col, vWscore, i, s1Len, j, s2Len);
            arr_store_rowcol(result->matches_row, result->matches_col, vWmatch, i, s1Len, j, s2Len);
            arr_store_rowcol(result->similar_row, result->similar_col, vWsimilar, i, s1Len, j, s2Len);
            arr_store_rowcol(result->length_row, result->length_col, vWlength, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-15] = (int8_t)_mm_extract_epi8_rpl(vWscore,0);
            mch_pr[j-15] = (int8_t)_mm_extract_epi8_rpl(vWmatch,0);
            sim_pr[j-15] = (int8_t)_mm_extract_epi8_rpl(vWsimilar,0);
            len_pr[j-15] = (int8_t)_mm_extract_epi8_rpl(vWlength,0);
            del_pr[j-15] = (int8_t)_mm_extract_epi8_rpl(vDel,0);
            /* as minor diagonal vector passes across table, extract
             * max values within the i,j bounds */
            {
                __m128i cond_valid_J = _mm_and_si128(
                        _mm_packs_epi16(
                            _mm_cmpgt_epi16(vJLo16, vNegOne16),
                            _mm_cmpgt_epi16(vJHi16, vNegOne16)),
                        _mm_packs_epi16(
                            _mm_cmplt_epi16(vJLo16, vJLimit16),
                            _mm_cmplt_epi16(vJHi16, vJLimit16)));
                __m128i cond_valid_IJ = _mm_and_si128(cond_valid_J, vIltLimit);
                __m128i cond_eq = _mm_cmpeq_epi8(vWscore, vMaxScore);
                __m128i cond_max = _mm_cmpgt_epi8(vWscore, vMaxScore);
                __m128i cond_all = _mm_and_si128(cond_max, cond_valid_IJ);
                __m128i cond_Jlt = _mm_packs_epi16(
                        _mm_cmplt_epi16(vJLo16, vEndJLo),
                        _mm_cmplt_epi16(vJHi16, vEndJHi));
                __m128i cond_lo = _mm_unpacklo_epi8(cond_all, cond_all);
                __m128i cond_hi = _mm_unpackhi_epi8(cond_all, cond_all);
                vMaxScore = _mm_blendv_epi8_rpl(vMaxScore, vWscore, cond_all);
                vMaxMatch = _mm_blendv_epi8_rpl(vMaxMatch, vWmatch, cond_all);
                vMaxSimilar = _mm_blendv_epi8_rpl(vMaxSimilar, vWsimilar, cond_all);
                vMaxLength = _mm_blendv_epi8_rpl(vMaxLength, vWlength, cond_all);
                vEndILo = _mm_blendv_epi8_rpl(vEndILo, vILo16, cond_lo);
                vEndIHi = _mm_blendv_epi8_rpl(vEndIHi, vIHi16, cond_hi);
                vEndJLo = _mm_blendv_epi8_rpl(vEndJLo, vJLo16, cond_lo);
                vEndJHi = _mm_blendv_epi8_rpl(vEndJHi, vJHi16, cond_hi);
                cond_all = _mm_and_si128(cond_Jlt, cond_eq);
                cond_all = _mm_and_si128(cond_all, cond_valid_IJ);
                cond_lo = _mm_unpacklo_epi8(cond_all, cond_all);
                cond_hi = _mm_unpackhi_epi8(cond_all, cond_all);
                vMaxMatch = _mm_blendv_epi8_rpl(vMaxMatch, vWmatch, cond_all);
                vMaxSimilar = _mm_blendv_epi8_rpl(vMaxSimilar, vWsimilar, cond_all);
                vMaxLength = _mm_blendv_epi8_rpl(vMaxLength, vWlength, cond_all);
                vEndILo = _mm_blendv_epi8_rpl(vEndILo, vILo16, cond_lo);
                vEndIHi = _mm_blendv_epi8_rpl(vEndIHi, vIHi16, cond_hi);
                vEndJLo = _mm_blendv_epi8_rpl(vEndJLo, vJLo16, cond_lo);
                vEndJHi = _mm_blendv_epi8_rpl(vEndJHi, vJHi16, cond_hi);
            }
            vJLo16 = _mm_add_epi16(vJLo16, vOne16);
            vJHi16 = _mm_add_epi16(vJHi16, vOne16);
        }
        vILo16 = _mm_add_epi16(vILo16, vN16);
        vIHi16 = _mm_add_epi16(vIHi16, vN16);
    }

    /* alignment ending position */
    {
        int8_t *t = (int8_t*)&vMaxScore;
        int8_t *m = (int8_t*)&vMaxMatch;
        int8_t *s = (int8_t*)&vMaxSimilar;
        int8_t *l = (int8_t*)&vMaxLength;
        int16_t *ilo = (int16_t*)&vEndILo;
        int16_t *jlo = (int16_t*)&vEndJLo;
        int16_t *ihi = (int16_t*)&vEndIHi;
        int16_t *jhi = (int16_t*)&vEndJHi;
        int32_t k;
        for (k=0; k<N/2; ++k, ++t, ++m, ++s, ++l, ++ilo, ++jlo) {
            if (*t > score) {
                score = *t;
                matches = *m;
                similar = *s;
                length = *l;
                end_query = *ilo;
                end_ref = *jlo;
            }
            else if (*t == score) {
                if (*jlo < end_ref) {
                    matches = *m;
                    similar = *s;
                    length = *l;
                    end_query = *ilo;
                    end_ref = *jlo;
                }
                else if (*jlo == end_ref && *ilo < end_query) {
                    matches = *m;
                    similar = *s;
                    length = *l;
                    end_query = *ilo;
                    end_ref = *jlo;
                }
            }
        }
        for (k=N/2; k<N; ++k, ++t, ++m, ++s, ++l, ++ihi, ++jhi) {
            if (*t > score) {
                score = *t;
                matches = *m;
                similar = *s;
                length = *l;
                end_query = *ihi;
                end_ref = *jhi;
            }
            else if (*t == score) {
                if (*jhi < end_ref) {
                    matches = *m;
                    similar = *s;
                    length = *l;
                    end_query = *ihi;
                    end_ref = *jhi;
                }
                else if (*jhi == end_ref && *ihi < end_query) {
                    matches = *m;
                    similar = *s;
                    length = *l;
                    end_query = *ihi;
                    end_ref = *jhi;
                }
            }
        }
    }

    if (_mm_movemask_epi8(_mm_or_si128(
            _mm_cmpeq_epi8(vSaturationCheckMin, vNegLimit),
            _mm_cmpeq_epi8(vSaturationCheckMax, vPosLimit)))) {
        result->saturated = 1;
        score = INT8_MAX;
        matches = 0;
        similar = 0;
        length = 0;
        end_query = 0;
        end_ref = 0;
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


