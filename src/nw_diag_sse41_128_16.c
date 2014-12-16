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
#include "blosum/blosum_map.h"

#define NEG_INF_16 (INT16_MIN/(int16_t)(2))
#define MAX(a,b) ((a)>(b)?(a):(b))

/* shift given vector v, insert val, return shifted val */
static inline __m128i vshift16(const __m128i v, const int val)
{
    __m128i ret = _mm_srli_si128(v, 2);
    ret = _mm_insert_epi16(ret, val, 7);
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
        array[(i+0)*s2Len + (j-0)] = (int16_t)_mm_extract_epi16(vWscore, 7);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[(i+1)*s2Len + (j-1)] = (int16_t)_mm_extract_epi16(vWscore, 6);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[(i+2)*s2Len + (j-2)] = (int16_t)_mm_extract_epi16(vWscore, 5);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[(i+3)*s2Len + (j-3)] = (int16_t)_mm_extract_epi16(vWscore, 4);
    }
    if (0 <= i+4 && i+4 < s1Len && 0 <= j-4 && j-4 < s2Len) {
        array[(i+4)*s2Len + (j-4)] = (int16_t)_mm_extract_epi16(vWscore, 3);
    }
    if (0 <= i+5 && i+5 < s1Len && 0 <= j-5 && j-5 < s2Len) {
        array[(i+5)*s2Len + (j-5)] = (int16_t)_mm_extract_epi16(vWscore, 2);
    }
    if (0 <= i+6 && i+6 < s1Len && 0 <= j-6 && j-6 < s2Len) {
        array[(i+6)*s2Len + (j-6)] = (int16_t)_mm_extract_epi16(vWscore, 1);
    }
    if (0 <= i+7 && i+7 < s1Len && 0 <= j-7 && j-7 < s2Len) {
        array[(i+7)*s2Len + (j-7)] = (int16_t)_mm_extract_epi16(vWscore, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME nw_table_diag_sse41_128_16
#else
#define FNAME nw_diag_sse41_128_16
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    const int N = 8; /* number of values in vector */
    const int PAD2 = N-1; /* N 16-byte values in vector, so N - 1 */
    const int PAD = PAD2*2;
    int * const restrict s1 = parasail_memalign_int(16, s1Len+PAD2);
    int * const restrict s2B= parasail_memalign_int(16, s2Len+PAD);
    int * const restrict _tbl_pr = parasail_memalign_int(16, s2Len+PAD);
    int * const restrict _del_pr = parasail_memalign_int(16, s2Len+PAD);
    int * const restrict s2 = s2B+PAD2; /* will allow later for negative indices */
    int * const restrict tbl_pr = _tbl_pr+PAD2;
    int * const restrict del_pr = _del_pr+PAD2;
#if PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
    int i = 0;
    int j = 0;
    int score = NEG_INF_16;
    __m128i vNegInf = _mm_set1_epi16(NEG_INF_16);
    __m128i vOpen = _mm_set1_epi16(open);
    __m128i vGap  = _mm_set1_epi16(gap);
    __m128i vOne = _mm_set1_epi16(1);
    __m128i vN = _mm_set1_epi16(N);
    __m128i vGapN = _mm_mullo_epi16(vN, vGap);
    __m128i vNegOne = _mm_set1_epi16(-1);
    __m128i vI = _mm_set_epi16(0,1,2,3,4,5,6,7);
    __m128i vJreset = _mm_set_epi16(0,-1,-2,-3,-4,-5,-6,-7);
    __m128i vMax = vNegInf;
    __m128i vILimit = _mm_set1_epi16(s1Len);
    __m128i vILimit1 = _mm_sub_epi16(vILimit, vOne);
    __m128i vJLimit = _mm_set1_epi16(s2Len);
    __m128i vJLimit1 = _mm_sub_epi16(vJLimit, vOne);
    __m128i vIBoundary = _mm_set_epi16(
            -open-0*gap,
            -open-1*gap,
            -open-2*gap,
            -open-3*gap,
            -open-4*gap,
            -open-5*gap,
            -open-6*gap,
            -open-7*gap);
    assert(s1Len > N);
    assert(s2Len > N);

    /* convert _s1 from char to int in range 0-23 */
    for (i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[(unsigned char)_s1[i]];
    }
    /* pad back of s1 with dummy values */
    for (i=s1Len; i<s1Len+PAD2; ++i) {
        s1[i] = 0; /* point to first matrix row because we don't care */
    }

    /* convert _s2 from char to int in range 0-23 */
    for (j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[(unsigned char)_s2[j]];
    }
    /* pad front of s2 with dummy values */
    for (j=-PAD2; j<0; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }
    /* pad back of s2 with dummy values */
    for (j=s2Len; j<s2Len+PAD2; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }

    /* set initial values for stored row */
    for (j=0; j<s2Len; ++j) {
        tbl_pr[j] = -open - j*gap;
        del_pr[j] = NEG_INF_16;
    }
    /* pad front of stored row values */
    for (j=-PAD2; j<0; ++j) {
        tbl_pr[j] = NEG_INF_16;
        del_pr[j] = NEG_INF_16;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD2; ++j) {
        tbl_pr[j] = NEG_INF_16;
        del_pr[j] = NEG_INF_16;
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
        const int * const restrict matrow4 = matrix[s1[i+4]];
        const int * const restrict matrow5 = matrix[s1[i+5]];
        const int * const restrict matrow6 = matrix[s1[i+6]];
        const int * const restrict matrow7 = matrix[s1[i+7]];
        vNscore = vshift16(vNscore, tbl_pr[-1]);
        vWscore = vshift16(vWscore, -open - i*gap);
        tbl_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<N; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            vNscore = vshift16(vWscore, tbl_pr[j]);
            vDel = vshift16(vDel, del_pr[j]);
            vDel = _mm_max_epi16(
                    _mm_sub_epi16(vNscore, vOpen),
                    _mm_sub_epi16(vDel, vGap));
            vIns = _mm_max_epi16(
                    _mm_sub_epi16(vWscore, vOpen),
                    _mm_sub_epi16(vIns, vGap));
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
            vNWscore = _mm_add_epi16(vNWscore, vMat);
            vWscore = _mm_max_epi16(vNWscore, vIns);
            vWscore = _mm_max_epi16(vWscore, vDel);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m128i cond = _mm_cmpeq_epi16(vJ,vNegOne);
                vWscore = _mm_andnot_si128(cond, vWscore); /* all but j=-1 */
                vWscore = _mm_or_si128(vWscore,
                        _mm_and_si128(cond, vIBoundary));
                vDel = _mm_andnot_si128(cond, vDel);
                vDel = _mm_or_si128(vDel, _mm_and_si128(cond, vNegInf));
                vIns = _mm_andnot_si128(cond, vIns);
                vIns = _mm_or_si128(vIns, _mm_and_si128(cond, vNegInf));
            }
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-7] = (int16_t)_mm_extract_epi16(vWscore,0);
            del_pr[j-7] = (int16_t)_mm_extract_epi16(vDel,0);
            vJ = _mm_add_epi16(vJ, vOne);
        }
        for (j=N; j<s2Len+PAD2; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            vNscore = vshift16(vWscore, tbl_pr[j]);
            vDel = vshift16(vDel, del_pr[j]);
            vDel = _mm_max_epi16(
                    _mm_sub_epi16(vNscore, vOpen),
                    _mm_sub_epi16(vDel, vGap));
            vIns = _mm_max_epi16(
                    _mm_sub_epi16(vWscore, vOpen),
                    _mm_sub_epi16(vIns, vGap));
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
            vNWscore = _mm_add_epi16(vNWscore, vMat);
            vWscore = _mm_max_epi16(vNWscore, vIns);
            vWscore = _mm_max_epi16(vWscore, vDel);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-7] = (int16_t)_mm_extract_epi16(vWscore,0);
            del_pr[j-7] = (int16_t)_mm_extract_epi16(vDel,0);
            vJ = _mm_add_epi16(vJ, vOne);
        }
        vI = _mm_add_epi16(vI, vN);
        vIBoundary = _mm_sub_epi16(vIBoundary, vGapN);
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
        const int * const restrict matrow4 = matrix[s1[i+4]];
        const int * const restrict matrow5 = matrix[s1[i+5]];
        const int * const restrict matrow6 = matrix[s1[i+6]];
        const int * const restrict matrow7 = matrix[s1[i+7]];
        vNscore = vshift16(vNscore, tbl_pr[-1]);
        vWscore = vshift16(vWscore, -open - i*gap);
        tbl_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<N; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            vNscore = vshift16(vWscore, tbl_pr[j]);
            vDel = vshift16(vDel, del_pr[j]);
            vDel = _mm_max_epi16(
                    _mm_sub_epi16(vNscore, vOpen),
                    _mm_sub_epi16(vDel, vGap));
            vIns = _mm_max_epi16(
                    _mm_sub_epi16(vWscore, vOpen),
                    _mm_sub_epi16(vIns, vGap));
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
            vNWscore = _mm_add_epi16(vNWscore, vMat);
            vWscore = _mm_max_epi16(vNWscore, vIns);
            vWscore = _mm_max_epi16(vWscore, vDel);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m128i cond = _mm_cmpeq_epi16(vJ,vNegOne);
                vWscore = _mm_andnot_si128(cond, vWscore); /* all but j=-1 */
                vWscore = _mm_or_si128(vWscore,
                        _mm_and_si128(cond, vIBoundary));
                vDel = _mm_andnot_si128(cond, vDel);
                vDel = _mm_or_si128(vDel, _mm_and_si128(cond, vNegInf));
                vIns = _mm_andnot_si128(cond, vIns);
                vIns = _mm_or_si128(vIns, _mm_and_si128(cond, vNegInf));
            }
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-7] = (int16_t)_mm_extract_epi16(vWscore,0);
            del_pr[j-7] = (int16_t)_mm_extract_epi16(vDel,0);
            vJ = _mm_add_epi16(vJ, vOne);
        }
        for (j=N; j<s2Len-1; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            vNscore = vshift16(vWscore, tbl_pr[j]);
            vDel = vshift16(vDel, del_pr[j]);
            vDel = _mm_max_epi16(
                    _mm_sub_epi16(vNscore, vOpen),
                    _mm_sub_epi16(vDel, vGap));
            vIns = _mm_max_epi16(
                    _mm_sub_epi16(vWscore, vOpen),
                    _mm_sub_epi16(vIns, vGap));
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
            vNWscore = _mm_add_epi16(vNWscore, vMat);
            vWscore = _mm_max_epi16(vNWscore, vIns);
            vWscore = _mm_max_epi16(vWscore, vDel);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-7] = (int16_t)_mm_extract_epi16(vWscore,0);
            del_pr[j-7] = (int16_t)_mm_extract_epi16(vDel,0);
            vJ = _mm_add_epi16(vJ, vOne);
        }
        for (j=s2Len-1; j<s2Len+PAD2; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            vNscore = vshift16(vWscore, tbl_pr[j]);
            vDel = vshift16(vDel, del_pr[j]);
            vDel = _mm_max_epi16(
                    _mm_sub_epi16(vNscore, vOpen),
                    _mm_sub_epi16(vDel, vGap));
            vIns = _mm_max_epi16(
                    _mm_sub_epi16(vWscore, vOpen),
                    _mm_sub_epi16(vIns, vGap));
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
            vNWscore = _mm_add_epi16(vNWscore, vMat);
            vWscore = _mm_max_epi16(vNWscore, vIns);
            vWscore = _mm_max_epi16(vWscore, vDel);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-7] = (int16_t)_mm_extract_epi16(vWscore,0);
            del_pr[j-7] = (int16_t)_mm_extract_epi16(vDel,0);
            /* as minor diagonal vector passes across table, extract
               last table value at the i,j bound */
            {
                __m128i cond_valid_I = _mm_cmpeq_epi16(vI, vILimit1);
                __m128i cond_valid_J = _mm_cmpeq_epi16(vJ, vJLimit1);
                __m128i cond_max = _mm_cmpgt_epi16(vWscore, vMax);
                __m128i cond_all = _mm_and_si128(cond_max,
                        _mm_and_si128(cond_valid_I, cond_valid_J));
                vMax = _mm_andnot_si128(cond_all, vMax); /* keep old */
                vMax = _mm_or_si128(vMax,
                        _mm_and_si128(cond_all, vWscore));
            }
            vJ = _mm_add_epi16(vJ, vOne);
        }
        vI = _mm_add_epi16(vI, vN);
        vIBoundary = _mm_sub_epi16(vIBoundary, vGapN);
    }

    /* max in vMax */
    for (i=0; i<N; ++i) {
        int16_t value;
        value = (int16_t) _mm_extract_epi16(vMax, 7);
        if (value > score) {
            score = value;
        }
        vMax = _mm_slli_si128(vMax, 2);
    }

    result->score = score;

    free(_del_pr);
    free(_tbl_pr);
    free(s2B);
    free(s1);

    return result;
}

