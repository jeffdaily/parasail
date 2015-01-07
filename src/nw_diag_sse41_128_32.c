/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#include "config.h"

#include <assert.h>
#include <stdlib.h>

#include <emmintrin.h>
#include <smmintrin.h>

#include "parasail.h"
#include "parasail_internal.h"
#include "parasail_internal_sse.h"
#include "blosum/blosum_map.h"

#define NEG_INF_32 (INT32_MIN/(int32_t)(2))

/* shift given vector v, insert val, return shifted val */
static inline __m128i vshift32(const __m128i v, const int val)
{
    __m128i ret = _mm_srli_si128(v, 4);
    ret = _mm_insert_epi32(ret, val, 3);
    return ret;
}

#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        __m128i vWscore,
        int i,
        int s1Len,
        int j,
        int s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[(i+0)*s2Len + (j-0)] = (int32_t)_mm_extract_epi32(vWscore, 3);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[(i+1)*s2Len + (j-1)] = (int32_t)_mm_extract_epi32(vWscore, 2);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[(i+2)*s2Len + (j-2)] = (int32_t)_mm_extract_epi32(vWscore, 1);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[(i+3)*s2Len + (j-3)] = (int32_t)_mm_extract_epi32(vWscore, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME nw_table_diag_sse41_128_32
#else
#define FNAME nw_diag_sse41_128_32
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    const int N = 4; /* number of values in vector */
    const int PAD = N-1; /* N 32-byte values in vector, so N - 1 */
    const int PAD2 = PAD*2;
    int * const restrict s1 = parasail_memalign_int(16, s1Len+PAD);
    int * const restrict s2B= parasail_memalign_int(16, s2Len+PAD2);
    int * const restrict _tbl_pr = parasail_memalign_int(16, s2Len+PAD2);
    int * const restrict _del_pr = parasail_memalign_int(16, s2Len+PAD2);
    int * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int * const restrict tbl_pr = _tbl_pr+PAD;
    int * const restrict del_pr = _del_pr+PAD;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
    int i = 0;
    int j = 0;
    int score = NEG_INF_32;
    __m128i vNegInf = _mm_set1_epi32(NEG_INF_32);
    __m128i vOpen = _mm_set1_epi32(open);
    __m128i vGap  = _mm_set1_epi32(gap);
    __m128i vOne = _mm_set1_epi32(1);
    __m128i vN = _mm_set1_epi32(N);
    __m128i vGapN = _mm_set1_epi32(gap*N);
    __m128i vNegOne = _mm_set1_epi32(-1);
    __m128i vI = _mm_set_epi32(0,1,2,3);
    __m128i vJreset = _mm_set_epi32(0,-1,-2,-3);
    __m128i vMax = vNegInf;
    __m128i vILimit = _mm_set1_epi32(s1Len);
    __m128i vILimit1 = _mm_sub_epi32(vILimit, vOne);
    __m128i vJLimit = _mm_set1_epi32(s2Len);
    __m128i vJLimit1 = _mm_sub_epi32(vJLimit, vOne);
    __m128i vIBoundary = _mm_set_epi32(
            -open-0*gap,
            -open-1*gap,
            -open-2*gap,
            -open-3*gap);
    assert(s1Len > N);
    assert(s2Len > N);

    /* convert _s1 from char to int in range 0-23 */
    for (i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[(unsigned char)_s1[i]];
    }
    /* pad back of s1 with dummy values */
    for (i=s1Len; i<s1Len+PAD; ++i) {
        s1[i] = 0; /* point to first matrix row because we don't care */
    }

    /* convert _s2 from char to int in range 0-23 */
    for (j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[(unsigned char)_s2[j]];
    }
    /* pad front of s2 with dummy values */
    for (j=-PAD; j<0; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }
    /* pad back of s2 with dummy values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }

    /* set initial values for stored row */
    for (j=0; j<s2Len; ++j) {
        tbl_pr[j] = -open - j*gap;
        del_pr[j] = NEG_INF_32;
    }
    /* pad front of stored row values */
    for (j=-PAD; j<0; ++j) {
        tbl_pr[j] = NEG_INF_32;
        del_pr[j] = NEG_INF_32;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        tbl_pr[j] = NEG_INF_32;
        del_pr[j] = NEG_INF_32;
    }
    tbl_pr[-1] = 0; /* upper left corner */

    /* iterate over query sequence */
    for (i=0; i<s1Len-N; i+=N) {
        __m128i vNscore = vNegInf;
        __m128i vWscore = vNegInf;
        __m128i vIns = vNegInf;
        __m128i vDel = vNegInf;
        __m128i vJ = vJreset;
        const int * const restrict matrow0 = matrix[s1[i+0]];
        const int * const restrict matrow1 = matrix[s1[i+1]];
        const int * const restrict matrow2 = matrix[s1[i+2]];
        const int * const restrict matrow3 = matrix[s1[i+3]];
        vNscore = vshift32(vNscore, tbl_pr[-1]);
        vWscore = vshift32(vWscore, -open - i*gap);
        tbl_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<N; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            vNscore = vshift32(vWscore, tbl_pr[j]);
            vDel = vshift32(vDel, del_pr[j]);
            vDel = _mm_max_epi32(
                    _mm_sub_epi32(vNscore, vOpen),
                    _mm_sub_epi32(vDel, vGap));
            vIns = _mm_max_epi32(
                    _mm_sub_epi32(vWscore, vOpen),
                    _mm_sub_epi32(vIns, vGap));
            vMat = _mm_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]]
                    );
            vNWscore = _mm_add_epi32(vNWscore, vMat);
            vWscore = _mm_max_epi32(vNWscore, vIns);
            vWscore = _mm_max_epi32(vWscore, vDel);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m128i cond = _mm_cmpeq_epi32(vJ,vNegOne);
                vWscore = _mm_blendv_epi8(vWscore, vIBoundary, cond);
                vDel = _mm_blendv_epi8(vDel, vNegInf, cond);
                vIns = _mm_blendv_epi8(vIns, vNegInf, cond);
            }
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-3] = (int32_t)_mm_extract_epi32(vWscore,0);
            del_pr[j-3] = (int32_t)_mm_extract_epi32(vDel,0);
            vJ = _mm_add_epi32(vJ, vOne);
        }
        for (j=N; j<s2Len+PAD; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            vNscore = vshift32(vWscore, tbl_pr[j]);
            vDel = vshift32(vDel, del_pr[j]);
            vDel = _mm_max_epi32(
                    _mm_sub_epi32(vNscore, vOpen),
                    _mm_sub_epi32(vDel, vGap));
            vIns = _mm_max_epi32(
                    _mm_sub_epi32(vWscore, vOpen),
                    _mm_sub_epi32(vIns, vGap));
            vMat = _mm_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]]
                    );
            vNWscore = _mm_add_epi32(vNWscore, vMat);
            vWscore = _mm_max_epi32(vNWscore, vIns);
            vWscore = _mm_max_epi32(vWscore, vDel);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-3] = (int32_t)_mm_extract_epi32(vWscore,0);
            del_pr[j-3] = (int32_t)_mm_extract_epi32(vDel,0);
            vJ = _mm_add_epi32(vJ, vOne);
        }
        vI = _mm_add_epi32(vI, vN);
        vIBoundary = _mm_sub_epi32(vIBoundary, vGapN);
    }
    for (/*i=?*/; i<s1Len; i+=N) {
        __m128i vNscore = vNegInf;
        __m128i vWscore = vNegInf;
        __m128i vIns = vNegInf;
        __m128i vDel = vNegInf;
        __m128i vJ = vJreset;
        const int * const restrict matrow0 = matrix[s1[i+0]];
        const int * const restrict matrow1 = matrix[s1[i+1]];
        const int * const restrict matrow2 = matrix[s1[i+2]];
        const int * const restrict matrow3 = matrix[s1[i+3]];
        vNscore = vshift32(vNscore, tbl_pr[-1]);
        vWscore = vshift32(vWscore, -open - i*gap);
        tbl_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<N; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            vNscore = vshift32(vWscore, tbl_pr[j]);
            vDel = vshift32(vDel, del_pr[j]);
            vDel = _mm_max_epi32(
                    _mm_sub_epi32(vNscore, vOpen),
                    _mm_sub_epi32(vDel, vGap));
            vIns = _mm_max_epi32(
                    _mm_sub_epi32(vWscore, vOpen),
                    _mm_sub_epi32(vIns, vGap));
            vMat = _mm_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]]
                    );
            vNWscore = _mm_add_epi32(vNWscore, vMat);
            vWscore = _mm_max_epi32(vNWscore, vIns);
            vWscore = _mm_max_epi32(vWscore, vDel);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m128i cond = _mm_cmpeq_epi32(vJ,vNegOne);
                vWscore = _mm_blendv_epi8(vWscore, vIBoundary, cond);
                vDel = _mm_blendv_epi8(vDel, vNegInf, cond);
                vIns = _mm_blendv_epi8(vIns, vNegInf, cond);
            }
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-3] = (int32_t)_mm_extract_epi32(vWscore,0);
            del_pr[j-3] = (int32_t)_mm_extract_epi32(vDel,0);
            vJ = _mm_add_epi32(vJ, vOne);
        }
        for (j=N; j<s2Len-1; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            vNscore = vshift32(vWscore, tbl_pr[j]);
            vDel = vshift32(vDel, del_pr[j]);
            vDel = _mm_max_epi32(
                    _mm_sub_epi32(vNscore, vOpen),
                    _mm_sub_epi32(vDel, vGap));
            vIns = _mm_max_epi32(
                    _mm_sub_epi32(vWscore, vOpen),
                    _mm_sub_epi32(vIns, vGap));
            vMat = _mm_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]]
                    );
            vNWscore = _mm_add_epi32(vNWscore, vMat);
            vWscore = _mm_max_epi32(vNWscore, vIns);
            vWscore = _mm_max_epi32(vWscore, vDel);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-3] = (int32_t)_mm_extract_epi32(vWscore,0);
            del_pr[j-3] = (int32_t)_mm_extract_epi32(vDel,0);
            vJ = _mm_add_epi32(vJ, vOne);
        }
        for (j=s2Len-1; j<s2Len+PAD; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            vNscore = vshift32(vWscore, tbl_pr[j]);
            vDel = vshift32(vDel, del_pr[j]);
            vDel = _mm_max_epi32(
                    _mm_sub_epi32(vNscore, vOpen),
                    _mm_sub_epi32(vDel, vGap));
            vIns = _mm_max_epi32(
                    _mm_sub_epi32(vWscore, vOpen),
                    _mm_sub_epi32(vIns, vGap));
            vMat = _mm_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]]
                    );
            vNWscore = _mm_add_epi32(vNWscore, vMat);
            vWscore = _mm_max_epi32(vNWscore, vIns);
            vWscore = _mm_max_epi32(vWscore, vDel);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-3] = (int32_t)_mm_extract_epi32(vWscore,0);
            del_pr[j-3] = (int32_t)_mm_extract_epi32(vDel,0);
            /* as minor diagonal vector passes across table, extract
               last table value at the i,j bound */
            {
                __m128i cond_valid_I = _mm_cmpeq_epi32(vI, vILimit1);
                __m128i cond_valid_J = _mm_cmpeq_epi32(vJ, vJLimit1);
                __m128i cond_max = _mm_cmpgt_epi32(vWscore, vMax);
                __m128i cond_all = _mm_and_si128(cond_max,
                        _mm_and_si128(cond_valid_I, cond_valid_J));
                vMax = _mm_blendv_epi8(vMax, vWscore, cond_all);
            }
            vJ = _mm_add_epi32(vJ, vOne);
        }
        vI = _mm_add_epi32(vI, vN);
        vIBoundary = _mm_sub_epi32(vIBoundary, vGapN);
    }

    /* max in vMax */
    for (i=0; i<N; ++i) {
        int32_t value;
        value = (int32_t) _mm_extract_epi32(vMax, 3);
        if (value > score) {
            score = value;
        }
        vMax = _mm_slli_si128(vMax, 4);
    }

    result->score = score;

    parasail_free(_del_pr);
    parasail_free(_tbl_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}

