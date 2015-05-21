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
#include "parasail/internal_knc.h"
#include "parasail/matrices/blosum_map.h"

#define NEG_INF_32 (INT32_MIN/(int32_t)(2))

static inline __m512i insert(__m512i a, int32_t b, int imm) {
    __m512i_32_t tmp;
    tmp.m = a;
    tmp.v[imm] = b;
    return tmp.m;
}

static inline int32_t extract(__m512i a, int imm) {
    __m512i_32_t tmp;
    tmp.m = a;
    return tmp.v[imm];
}

/* shift given vector v, insert val, return shifted val */
static inline __m512i vshift8_(const __m512i permute2_idx, const __m512i v, const int val)
{
    __m512i ret = _mm512_permutevar_epi32(permute2_idx, v);
    ret = insert(ret, val, 15);
    return ret;
}

#define vshift8(v, val) vshift8_(permute2_idx, v, val)

#ifdef PARASAIL_TABLE
static inline void arr_store_si512(
        int *array,
        __m512i vWscore,
        int i,
        int s1Len,
        int j,
        int s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[(i+0)*s2Len + (j-0)] = (int32_t)extract(vWscore, 15);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[(i+1)*s2Len + (j-1)] = (int32_t)extract(vWscore, 14);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[(i+2)*s2Len + (j-2)] = (int32_t)extract(vWscore, 13);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[(i+3)*s2Len + (j-3)] = (int32_t)extract(vWscore, 12);
    }
    if (0 <= i+4 && i+4 < s1Len && 0 <= j-4 && j-4 < s2Len) {
        array[(i+4)*s2Len + (j-4)] = (int32_t)extract(vWscore, 11);
    }
    if (0 <= i+5 && i+5 < s1Len && 0 <= j-5 && j-5 < s2Len) {
        array[(i+5)*s2Len + (j-5)] = (int32_t)extract(vWscore, 10);
    }
    if (0 <= i+6 && i+6 < s1Len && 0 <= j-6 && j-6 < s2Len) {
        array[(i+6)*s2Len + (j-6)] = (int32_t)extract(vWscore, 9);
    }
    if (0 <= i+7 && i+7 < s1Len && 0 <= j-7 && j-7 < s2Len) {
        array[(i+7)*s2Len + (j-7)] = (int32_t)extract(vWscore, 8);
    }
    if (0 <= i+8 && i+8 < s1Len && 0 <= j-8 && j-8 < s2Len) {
        array[(i+8)*s2Len + (j-8)] = (int32_t)extract(vWscore, 7);
    }
    if (0 <= i+9 && i+9 < s1Len && 0 <= j-9 && j-9 < s2Len) {
        array[(i+9)*s2Len + (j-9)] = (int32_t)extract(vWscore, 6);
    }
    if (0 <= i+10 && i+10 < s1Len && 0 <= j-10 && j-10 < s2Len) {
        array[(i+10)*s2Len + (j-10)] = (int32_t)extract(vWscore, 5);
    }
    if (0 <= i+11 && i+11 < s1Len && 0 <= j-11 && j-11 < s2Len) {
        array[(i+11)*s2Len + (j-11)] = (int32_t)extract(vWscore, 4);
    }
    if (0 <= i+12 && i+12 < s1Len && 0 <= j-12 && j-12 < s2Len) {
        array[(i+12)*s2Len + (j-12)] = (int32_t)extract(vWscore, 3);
    }
    if (0 <= i+13 && i+13 < s1Len && 0 <= j-13 && j-13 < s2Len) {
        array[(i+13)*s2Len + (j-13)] = (int32_t)extract(vWscore, 2);
    }
    if (0 <= i+14 && i+14 < s1Len && 0 <= j-14 && j-14 < s2Len) {
        array[(i+14)*s2Len + (j-14)] = (int32_t)extract(vWscore, 1);
    }
    if (0 <= i+15 && i+15 < s1Len && 0 <= j-15 && j-15 < s2Len) {
        array[(i+15)*s2Len + (j-15)] = (int32_t)extract(vWscore, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME sw_table_diag_knc_512_32
#else
#define FNAME sw_diag_knc_512_32
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    const int N = 16; /* number of values in vector */
    const int PAD = N-1; /* N 8-bit values in vector, so N - 1 */
    const int PAD2 = PAD*2;
    int * const restrict s1 = parasail_memalign_int(64, s1Len+PAD);
    int * const restrict s2B= parasail_memalign_int(64, s2Len+PAD2);
    int * const restrict _tbl_pr = parasail_memalign_int(64, s2Len+PAD2);
    int * const restrict _del_pr = parasail_memalign_int(64, s2Len+PAD2);
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
    __m512i vNegInf = _mm512_set1_epi32(NEG_INF_32);
    __m512i vNegInf0 = insert(vNegInf, 0, 15);
    __m512i vOpen = _mm512_set1_epi32(open);
    __m512i vGap  = _mm512_set1_epi32(gap);
    __m512i vZero = _mm512_set1_epi32(0);
    __m512i vOne = _mm512_set1_epi32(1);
    __m512i vNegOne = _mm512_set1_epi32(-1);
    __m512i vN = _mm512_set1_epi32(N);
    __m512i vI = _mm512_set_epi32(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
    __m512i vJreset = _mm512_set_epi32(0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15);
    __m512i vMaxScore = vNegInf;
    __m512i vILimit = _mm512_set1_epi32(s1Len);
    __m512i vJLimit = _mm512_set1_epi32(s2Len);
    __m512i permute_idx = _mm512_set_16to16_pi(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15);
    __m512i permute2_idx = _mm512_set_16to16_pi(0,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1);

    /* convert _s1 from char to int in range 0-23 */
    for (i=0; i<s1Len; ++i) {
        s1[i] = parasail_blosum_map[(unsigned char)_s1[i]];
    }
    /* pad back of s1 with dummy values */
    for (i=s1Len; i<s1Len+PAD; ++i) {
        s1[i] = 0; /* point to first matrix row because we don't care */
    }

    /* convert _s2 from char to int in range 0-23 */
    for (j=0; j<s2Len; ++j) {
        s2[j] = parasail_blosum_map[(unsigned char)_s2[j]];
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
        tbl_pr[j] = 0;
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

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        __m512i vNscore = vNegInf0;
        __m512i vWscore = vZero;
        __m512i vIns = vNegInf;
        __m512i vDel = vNegInf;
        __m512i vJ = vJreset;
        const int * const restrict matrow0 = matrix[s1[i+0]];
        const int * const restrict matrow1 = matrix[s1[i+1]];
        const int * const restrict matrow2 = matrix[s1[i+2]];
        const int * const restrict matrow3 = matrix[s1[i+3]];
        const int * const restrict matrow4 = matrix[s1[i+4]];
        const int * const restrict matrow5 = matrix[s1[i+5]];
        const int * const restrict matrow6 = matrix[s1[i+6]];
        const int * const restrict matrow7 = matrix[s1[i+7]];
        const int * const restrict matrow8 = matrix[s1[i+8]];
        const int * const restrict matrow9 = matrix[s1[i+9]];
        const int * const restrict matrow10 = matrix[s1[i+10]];
        const int * const restrict matrow11 = matrix[s1[i+11]];
        const int * const restrict matrow12 = matrix[s1[i+12]];
        const int * const restrict matrow13 = matrix[s1[i+13]];
        const int * const restrict matrow14 = matrix[s1[i+14]];
        const int * const restrict matrow15 = matrix[s1[i+15]];
        __mmask16 vIltLimit = _mm512_cmplt_epi32_mask(vI, vILimit);
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m512i vMat;
            __m512i vNWscore = vNscore;
            vNscore = vshift8(vWscore, tbl_pr[j]);
            vDel = vshift8(vDel, del_pr[j]);
            vDel = _mm512_max_epi32(
                    _mm512_sub_epi32(vNscore, vOpen),
                    _mm512_sub_epi32(vDel, vGap));
            vIns = _mm512_max_epi32(
                    _mm512_sub_epi32(vWscore, vOpen),
                    _mm512_sub_epi32(vIns, vGap));
            vMat = _mm512_set_epi32(
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
            vNWscore = _mm512_add_epi32(vNWscore, vMat);
            vWscore = _mm512_max_epi32(vNWscore, vIns);
            vWscore = _mm512_max_epi32(vWscore, vDel);
            vWscore = _mm512_max_epi32(vWscore, vZero);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __mmask16 cond = _mm512_cmpeq_epi32_mask(vJ,vNegOne);
                vWscore = _mm512_mask_blend_epi32(cond, vWscore, vZero);
                vDel = _mm512_mask_blend_epi32(cond, vDel, vNegInf);
                vIns = _mm512_mask_blend_epi32(cond, vIns, vNegInf);
            }
#ifdef PARASAIL_TABLE
            arr_store_si512(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-15] = (int32_t)extract(vWscore,0);
            del_pr[j-15] = (int32_t)extract(vDel,0);
            /* as minor diagonal vector passes across table, extract
             * max values within the i,j bounds */
            {
                __mmask16 cond_valid_J = _mm512_kand(
                        _mm512_cmpgt_epi32_mask(vJ, vNegOne),
                        _mm512_cmplt_epi32_mask(vJ, vJLimit));
                __mmask16 cond_max = _mm512_cmpgt_epi32_mask(vWscore, vMaxScore);
                __mmask16 cond_all = _mm512_kand(cond_max,
                        _mm512_kand(vIltLimit, cond_valid_J));
                vMaxScore = _mm512_mask_blend_epi32(cond_all, vMaxScore, vWscore);
            }
            vJ = _mm512_add_epi32(vJ, vOne);
        }
        vI = _mm512_add_epi32(vI, vN);
    }

    /* max in vMaxScore */
    for (i=0; i<N; ++i) {
        int32_t value;
        value = (int32_t) extract(vMaxScore, 15);
        if (value > score) {
            score = value;
        }
        vMaxScore = _mm512_permutevar_epi32(permute_idx, vMaxScore);
    }

    result->score = score;

    parasail_free(_del_pr);
    parasail_free(_tbl_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}

