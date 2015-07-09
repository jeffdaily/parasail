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

#define NEG_INF (INT32_MIN/(int32_t)(2))

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

#ifdef PARASAIL_ROWCOL
static inline void arr_store_rowcol(
        int *row,
        int *col,
        __m512i vWscore,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (i+0 == s1Len-1 && 0 <= j-0 && j-0 < s2Len) {
        row[j-0] = (int32_t)extract(vWscore, 15);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int32_t)extract(vWscore, 15);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int32_t)extract(vWscore, 14);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int32_t)extract(vWscore, 14);
    }
    if (i+2 == s1Len-1 && 0 <= j-2 && j-2 < s2Len) {
        row[j-2] = (int32_t)extract(vWscore, 13);
    }
    if (j-2 == s2Len-1 && 0 <= i+2 && i+2 < s1Len) {
        col[(i+2)] = (int32_t)extract(vWscore, 13);
    }
    if (i+3 == s1Len-1 && 0 <= j-3 && j-3 < s2Len) {
        row[j-3] = (int32_t)extract(vWscore, 12);
    }
    if (j-3 == s2Len-1 && 0 <= i+3 && i+3 < s1Len) {
        col[(i+3)] = (int32_t)extract(vWscore, 12);
    }
    if (i+4 == s1Len-1 && 0 <= j-4 && j-4 < s2Len) {
        row[j-4] = (int32_t)extract(vWscore, 11);
    }
    if (j-4 == s2Len-1 && 0 <= i+4 && i+4 < s1Len) {
        col[(i+4)] = (int32_t)extract(vWscore, 11);
    }
    if (i+5 == s1Len-1 && 0 <= j-5 && j-5 < s2Len) {
        row[j-5] = (int32_t)extract(vWscore, 10);
    }
    if (j-5 == s2Len-1 && 0 <= i+5 && i+5 < s1Len) {
        col[(i+5)] = (int32_t)extract(vWscore, 10);
    }
    if (i+6 == s1Len-1 && 0 <= j-6 && j-6 < s2Len) {
        row[j-6] = (int32_t)extract(vWscore, 9);
    }
    if (j-6 == s2Len-1 && 0 <= i+6 && i+6 < s1Len) {
        col[(i+6)] = (int32_t)extract(vWscore, 9);
    }
    if (i+7 == s1Len-1 && 0 <= j-7 && j-7 < s2Len) {
        row[j-7] = (int32_t)extract(vWscore, 8);
    }
    if (j-7 == s2Len-1 && 0 <= i+7 && i+7 < s1Len) {
        col[(i+7)] = (int32_t)extract(vWscore, 8);
    }
    if (i+8 == s1Len-1 && 0 <= j-8 && j-8 < s2Len) {
        row[j-8] = (int32_t)extract(vWscore, 7);
    }
    if (j-8 == s2Len-1 && 0 <= i+8 && i+8 < s1Len) {
        col[(i+8)] = (int32_t)extract(vWscore, 7);
    }
    if (i+9 == s1Len-1 && 0 <= j-9 && j-9 < s2Len) {
        row[j-9] = (int32_t)extract(vWscore, 6);
    }
    if (j-9 == s2Len-1 && 0 <= i+9 && i+9 < s1Len) {
        col[(i+9)] = (int32_t)extract(vWscore, 6);
    }
    if (i+10 == s1Len-1 && 0 <= j-10 && j-10 < s2Len) {
        row[j-10] = (int32_t)extract(vWscore, 5);
    }
    if (j-10 == s2Len-1 && 0 <= i+10 && i+10 < s1Len) {
        col[(i+10)] = (int32_t)extract(vWscore, 5);
    }
    if (i+11 == s1Len-1 && 0 <= j-11 && j-11 < s2Len) {
        row[j-11] = (int32_t)extract(vWscore, 4);
    }
    if (j-11 == s2Len-1 && 0 <= i+11 && i+11 < s1Len) {
        col[(i+11)] = (int32_t)extract(vWscore, 4);
    }
    if (i+12 == s1Len-1 && 0 <= j-12 && j-12 < s2Len) {
        row[j-12] = (int32_t)extract(vWscore, 3);
    }
    if (j-12 == s2Len-1 && 0 <= i+12 && i+12 < s1Len) {
        col[(i+12)] = (int32_t)extract(vWscore, 3);
    }
    if (i+13 == s1Len-1 && 0 <= j-13 && j-13 < s2Len) {
        row[j-13] = (int32_t)extract(vWscore, 2);
    }
    if (j-13 == s2Len-1 && 0 <= i+13 && i+13 < s1Len) {
        col[(i+13)] = (int32_t)extract(vWscore, 2);
    }
    if (i+14 == s1Len-1 && 0 <= j-14 && j-14 < s2Len) {
        row[j-14] = (int32_t)extract(vWscore, 1);
    }
    if (j-14 == s2Len-1 && 0 <= i+14 && i+14 < s1Len) {
        col[(i+14)] = (int32_t)extract(vWscore, 1);
    }
    if (i+15 == s1Len-1 && 0 <= j-15 && j-15 < s2Len) {
        row[j-15] = (int32_t)extract(vWscore, 0);
    }
    if (j-15 == s2Len-1 && 0 <= i+15 && i+15 < s1Len) {
        col[(i+15)] = (int32_t)extract(vWscore, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sg_stats_table_diag_knc_512_32
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_stats_rowcol_diag_knc_512_32
#else
#define FNAME parasail_sg_stats_diag_knc_512_32
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    const int N = 16; /* number of values in vector */
    const int PAD = N-1;
    const int PAD2 = PAD*2;
    const int32_t s1Len_PAD = s1Len+PAD;
    const int32_t s2Len_PAD = s2Len+PAD;
    int * const restrict s1 = parasail_memalign_int(64, s1Len+PAD);
    int * const restrict s2B= parasail_memalign_int(64, s2Len+PAD2);
    int * const restrict _tbl_pr = parasail_memalign_int(64, s2Len+PAD2);
    int * const restrict _del_pr = parasail_memalign_int(64, s2Len+PAD2);
    int * const restrict _mch_pr = parasail_memalign_int(64, s2Len+PAD2);
    int * const restrict _sim_pr = parasail_memalign_int(64, s2Len+PAD2);
    int * const restrict _len_pr = parasail_memalign_int(64, s2Len+PAD2);
    int * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int * const restrict tbl_pr = _tbl_pr+PAD;
    int * const restrict del_pr = _del_pr+PAD;
    int * const restrict mch_pr = _mch_pr+PAD;
    int * const restrict sim_pr = _sim_pr+PAD;
    int * const restrict len_pr = _len_pr+PAD;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif
    int i = 0;
    int j = 0;
    int score = NEG_INF;
    int matches = NEG_INF;
    int similar = NEG_INF;
    int length = NEG_INF;
    __m512i vNegInf = _mm512_set1_epi32(NEG_INF);
    __m512i vNegInf0 = insert(vNegInf, 0, 15);
    __m512i vOpen = _mm512_set1_epi32(open);
    __m512i vGap  = _mm512_set1_epi32(gap);
    __m512i vZero = _mm512_set1_epi32(0);
    __m512i vOne = _mm512_set1_epi32(1);
    __m512i vN = _mm512_set1_epi32(N);
    __m512i vNegOne = _mm512_set1_epi32(-1);
    __m512i vI = _mm512_set_epi32(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
    __m512i vJreset = _mm512_set_epi32(0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15);
    __m512i vMaxScore = vNegInf;
    __m512i vMaxMatch = vNegInf;
    __m512i vMaxSimilar = vNegInf;
    __m512i vMaxLength = vNegInf;
    __m512i vILimit = _mm512_set1_epi32(s1Len);
    __m512i vJLimit = _mm512_set1_epi32(s2Len);
    __m512i vILimit1 = _mm512_set1_epi32(s1Len-1);
    __m512i vJLimit1 = _mm512_set1_epi32(s2Len-1);
    __m512i permute_idx = _mm512_set_16to16_pi(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15);
    __m512i permute2_idx = _mm512_set_16to16_pi(0,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1);

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
    for (j=s2Len; j<s2Len_PAD; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        sim_pr[j] = 0;
        len_pr[j] = 0;
    }

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        __m512i vNscore = vNegInf0;
        __m512i vNmatch = vZero;
        __m512i vNsimilar = vZero;
        __m512i vNlength = vZero;
        __m512i vWscore = vZero;
        __m512i vWmatch = vZero;
        __m512i vWsimilar = vZero;
        __m512i vWlength = vZero;
        __m512i vIns = vNegInf;
        __m512i vDel = vNegInf;
        __m512i vJ = vJreset;
        __m512i vs1 = _mm512_set_16to16_pi(
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
                s1[i+15]
                );
        __m512i vs2 = vNegInf;
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
        __mmask16 vIeqLimit1 = _mm512_cmpeq_epi32_mask(vI, vILimit1);
        /* iterate over database sequence */
        for (j=0; j<s2Len_PAD; ++j) {
            __m512i vMat;
            __m512i vNWscore = vNscore;
            __m512i vNWmatch = vNmatch;
            __m512i vNWsimilar = vNsimilar;
            __m512i vNWlength = vNlength;
            vNscore = vshift8(vWscore, tbl_pr[j]);
            vNmatch = vshift8(vWmatch, mch_pr[j]);
            vNsimilar = vshift8(vWsimilar, sim_pr[j]);
            vNlength = vshift8(vWlength, len_pr[j]);
            vDel = vshift8(vDel, del_pr[j]);
            vDel = _mm512_max_epi32(
                    _mm512_sub_epi32(vNscore, vOpen),
                    _mm512_sub_epi32(vDel, vGap));
            vIns = _mm512_max_epi32(
                    _mm512_sub_epi32(vWscore, vOpen),
                    _mm512_sub_epi32(vIns, vGap));
            vs2 = vshift8(vs2, s2[j]);
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
            /* conditional block */
            {
                __mmask16 case1not;
                __mmask16 case2not;
                __mmask16 case2;
                __mmask16 case3;
                __m512i vCmatch;
                __m512i vCsimilar;
                __m512i vClength;
                case1not = _mm512_kor(
                        _mm512_cmplt_epi32_mask(vNWscore,vDel),
                        _mm512_cmplt_epi32_mask(vNWscore,vIns));
                case2not = _mm512_cmplt_epi32_mask(vDel,vIns);
                case2 = _mm512_kandn(case2not,case1not);
                case3 = _mm512_kand(case1not,case2not);
                vCmatch = _mm512_mask_blend_epi32(case1not,
                        _mm512_mask_add_epi32(
                            vNWmatch,
                            _mm512_cmpeq_epi32_mask(vs1,vs2),
                            vNWmatch, vOne),
                        _mm512_mask_blend_epi32(case2, vWmatch, vNmatch));
                vCsimilar = _mm512_mask_blend_epi32(case1not,
                        _mm512_mask_add_epi32(
                            vNWsimilar,
                            _mm512_cmpgt_epi32_mask(vMat,vZero),
                            vNWsimilar, vOne),
                        _mm512_mask_blend_epi32(case2, vWsimilar, vNsimilar));
                vClength = _mm512_mask_blend_epi32(case1not,
                        _mm512_add_epi32(vNWlength, vOne),
                        _mm512_mask_blend_epi32(case2,
                            _mm512_add_epi32(vWlength, vOne),
                            _mm512_add_epi32(vNlength, vOne)));
                vWmatch = vCmatch;
                vWsimilar = vCsimilar;
                vWlength = vClength;
            }

            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __mmask16 cond = _mm512_cmpeq_epi32_mask(vJ,vNegOne);
                vWscore = _mm512_mask_blend_epi32(cond, vWscore, vZero);
                vWmatch = _mm512_mask_blend_epi32(cond, vWmatch, vZero);
                vWsimilar = _mm512_mask_blend_epi32(cond, vWsimilar, vZero);
                vWlength = _mm512_mask_blend_epi32(cond, vWlength, vZero);
                vDel = _mm512_mask_blend_epi32(cond, vDel, vNegInf);
                vIns = _mm512_mask_blend_epi32(cond, vIns, vNegInf);
            }
#ifdef PARASAIL_TABLE
            arr_store_si512(result->score_table, vWscore, i, s1Len, j, s2Len);
            arr_store_si512(result->matches_table, vWmatch, i, s1Len, j, s2Len);
            arr_store_si512(result->similar_table, vWsimilar, i, s1Len, j, s2Len);
            arr_store_si512(result->length_table, vWlength, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->score_row, result->score_col, vWscore, i, s1Len, j, s2Len);
            arr_store_rowcol(result->matches_row, result->matches_col, vWmatch, i, s1Len, j, s2Len);
            arr_store_rowcol(result->similar_row, result->similar_col, vWsimilar, i, s1Len, j, s2Len);
            arr_store_rowcol(result->length_row, result->length_col, vWlength, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-15] = (int32_t)extract(vWscore,0);
            mch_pr[j-15] = (int32_t)extract(vWmatch,0);
            sim_pr[j-15] = (int32_t)extract(vWsimilar,0);
            len_pr[j-15] = (int32_t)extract(vWlength,0);
            del_pr[j-15] = (int32_t)extract(vDel,0);
            /* as minor diagonal vector passes across table, extract
             * max values at the i,j bounds */
            {
                __mmask16 cond_valid_I =
                    _mm512_kand(vIeqLimit1,
                            _mm512_kand(
                                _mm512_cmpgt_epi32_mask(vJ, vNegOne),
                                _mm512_cmplt_epi32_mask(vJ, vJLimit)));
                __mmask16 cond_valid_J =
                    _mm512_kand(
                            _mm512_cmpeq_epi32_mask(vJ, vJLimit1),
                            _mm512_cmplt_epi32_mask(vI, vILimit));
                __mmask16 cond_max = _mm512_cmpgt_epi32_mask(vWscore, vMaxScore);
                __mmask16 cond_all = _mm512_kand(cond_max,
                        _mm512_kor(cond_valid_I, cond_valid_J));
                vMaxScore = _mm512_mask_blend_epi32(cond_all, vMaxScore, vWscore);
                vMaxMatch = _mm512_mask_blend_epi32(cond_all, vMaxMatch, vWmatch);
                vMaxSimilar = _mm512_mask_blend_epi32(cond_all, vMaxSimilar, vWsimilar);
                vMaxLength = _mm512_mask_blend_epi32(cond_all, vMaxLength, vWlength);
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
            matches = (int32_t) extract(vMaxMatch, 15);
            similar = (int32_t) extract(vMaxSimilar, 15);
            length = (int32_t) extract(vMaxLength, 15);
        }
        vMaxScore = _mm512_permutevar_epi32(permute_idx, vMaxScore);
        vMaxMatch = _mm512_permutevar_epi32(permute_idx, vMaxMatch);
        vMaxSimilar = _mm512_permutevar_epi32(permute_idx, vMaxSimilar);
        vMaxLength = _mm512_permutevar_epi32(permute_idx, vMaxLength);
    }

    result->score = score;
    result->matches = matches;
    result->similar = similar;
    result->length = length;

    parasail_free(_len_pr);
    parasail_free(_sim_pr);
    parasail_free(_mch_pr);
    parasail_free(_del_pr);
    parasail_free(_tbl_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}

