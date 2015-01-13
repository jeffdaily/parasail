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

#include <stdlib.h>

#include <immintrin.h>

#include "parasail.h"
#include "parasail_internal.h"
#include "parasail_internal_avx.h"
#include "blosum/blosum_map.h"

#define NEG_INF_32 (INT32_MIN/(int32_t)(2))

/* avx2 does not have _mm256_cmplt_epi32, emulate it */
static inline __m256i _mm256_cmplt_epi32(__m256i a, __m256i b) {
    return _mm256_cmpgt_epi32(b,a);
}

#if HAVE_AVX2_MM256_INSERT_EPI32
#else
static inline __m256i _mm256_insert_epi32(__m256i a, int32_t b, int imm) {
    __m256i_32_t tmp;
    tmp.m = a;
    tmp.v[imm] = b;
    return tmp.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI32
#else
static inline int32_t _mm256_extract_epi32(__m256i a, int imm) {
    __m256i_32_t tmp;
    tmp.m = a;
    return tmp.v[imm];
}
#endif

/* avx2 _mm256_srli_si256 does not shift across 128-bit lanes, emulate it */
static inline __m256i rshift32(__m256i a) {
    return _mm256_or_si256(
            _mm256_slli_si256(
                _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(3,0,0,1)),
                12),
            _mm256_srli_si256(a, 4));
}

static inline __m256i lshift32(__m256i a) {
    return _mm256_alignr_epi8(a,
            _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)),
            12);
}

/* shift given vector v, insert val, return shifted val */
static inline __m256i vshift32(const __m256i v, const int val)
{
    __m256i ret = rshift32(v);
    ret = _mm256_insert_epi32(ret, val, 7);
    return ret;
}

#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vWscore,
        int i,
        int s1Len,
        int j,
        int s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[(i+0)*s2Len + (j-0)] = (int32_t)_mm256_extract_epi32(vWscore, 7);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[(i+1)*s2Len + (j-1)] = (int32_t)_mm256_extract_epi32(vWscore, 6);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[(i+2)*s2Len + (j-2)] = (int32_t)_mm256_extract_epi32(vWscore, 5);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[(i+3)*s2Len + (j-3)] = (int32_t)_mm256_extract_epi32(vWscore, 4);
    }
    if (0 <= i+4 && i+4 < s1Len && 0 <= j-4 && j-4 < s2Len) {
        array[(i+4)*s2Len + (j-4)] = (int32_t)_mm256_extract_epi32(vWscore, 3);
    }
    if (0 <= i+5 && i+5 < s1Len && 0 <= j-5 && j-5 < s2Len) {
        array[(i+5)*s2Len + (j-5)] = (int32_t)_mm256_extract_epi32(vWscore, 2);
    }
    if (0 <= i+6 && i+6 < s1Len && 0 <= j-6 && j-6 < s2Len) {
        array[(i+6)*s2Len + (j-6)] = (int32_t)_mm256_extract_epi32(vWscore, 1);
    }
    if (0 <= i+7 && i+7 < s1Len && 0 <= j-7 && j-7 < s2Len) {
        array[(i+7)*s2Len + (j-7)] = (int32_t)_mm256_extract_epi32(vWscore, 0);
    }
}
#endif


#ifdef PARASAIL_TABLE
#define FNAME sg_table_diag_avx2_256_32
#else
#define FNAME sg_diag_avx2_256_32
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    const int N = 8; /* number of values in vector */
    const int PAD = N-1; /* N 32-byte values in vector, so N - 1 */
    const int PAD2 = PAD*2;
    int * const restrict s1 = parasail_memalign_int(32, s1Len+PAD);
    int * const restrict s2B= parasail_memalign_int(32, s2Len+PAD2);
    int * const restrict _tbl_pr = parasail_memalign_int(32, s2Len+PAD2);
    int * const restrict _del_pr = parasail_memalign_int(32, s2Len+PAD2);
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
    __m256i vNegInf = _mm256_set1_epi32(NEG_INF_32);
    __m256i vNegInf0 = rshift32(vNegInf); /* shift in a 0 */
    __m256i vOpen = _mm256_set1_epi32(open);
    __m256i vGap  = _mm256_set1_epi32(gap);
    __m256i vOne = _mm256_set1_epi32(1);
    __m256i vN = _mm256_set1_epi32(N);
    __m256i vNegOne = _mm256_set1_epi32(-1);
    __m256i vI = _mm256_set_epi32(0,1,2,3,4,5,6,7);
    __m256i vJreset = _mm256_set_epi32(0,-1,-2,-3,-4,-5,-6,-7);
    __m256i vMaxScore = vNegInf;
    __m256i vILimit = _mm256_set1_epi32(s1Len);
    __m256i vILimit1 = _mm256_sub_epi32(vILimit, vOne);
    __m256i vJLimit = _mm256_set1_epi32(s2Len);
    __m256i vJLimit1 = _mm256_sub_epi32(vJLimit, vOne);

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
    for (i=0; i<s1Len-N; i+=N) {
        __m256i vNscore = vNegInf0;
        __m256i vWscore = vNegInf0;
        __m256i vIns = vNegInf;
        __m256i vDel = vNegInf;
        __m256i vJ = vJreset;
        const int * const restrict matrow0 = matrix[s1[i+0]];
        const int * const restrict matrow1 = matrix[s1[i+1]];
        const int * const restrict matrow2 = matrix[s1[i+2]];
        const int * const restrict matrow3 = matrix[s1[i+3]];
        const int * const restrict matrow4 = matrix[s1[i+4]];
        const int * const restrict matrow5 = matrix[s1[i+5]];
        const int * const restrict matrow6 = matrix[s1[i+6]];
        const int * const restrict matrow7 = matrix[s1[i+7]];
        /* iterate over database sequence */
        for (j=0; j<N; ++j) {
            __m256i vMat;
            __m256i vNWscore = vNscore;
            vNscore = vshift32(vWscore, tbl_pr[j]);
            vDel = vshift32(vDel, del_pr[j]);
            vDel = _mm256_max_epi32(
                    _mm256_sub_epi32(vNscore, vOpen),
                    _mm256_sub_epi32(vDel, vGap));
            vIns = _mm256_max_epi32(
                    _mm256_sub_epi32(vWscore, vOpen),
                    _mm256_sub_epi32(vIns, vGap));
            vMat = _mm256_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]],
                    matrow4[s2[j-4]],
                    matrow5[s2[j-5]],
                    matrow6[s2[j-6]],
                    matrow7[s2[j-7]]
                    );
            vNWscore = _mm256_add_epi32(vNWscore, vMat);
            vWscore = _mm256_max_epi32(vNWscore, vIns);
            vWscore = _mm256_max_epi32(vWscore, vDel);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m256i cond = _mm256_cmpeq_epi32(vJ,vNegOne);
                vWscore = _mm256_andnot_si256(cond, vWscore);
                vDel = _mm256_blendv_epi8(vDel, vNegInf, cond);
                vIns = _mm256_blendv_epi8(vIns, vNegInf, cond);
            }
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-7] = (int32_t)_mm256_extract_epi32(vWscore,0);
            del_pr[j-7] = (int32_t)_mm256_extract_epi32(vDel,0);
            vJ = _mm256_add_epi32(vJ, vOne);
        }
        for (j=N; j<s2Len-1; ++j) {
            __m256i vMat;
            __m256i vNWscore = vNscore;
            vNscore = vshift32(vWscore, tbl_pr[j]);
            vDel = vshift32(vDel, del_pr[j]);
            vDel = _mm256_max_epi32(
                    _mm256_sub_epi32(vNscore, vOpen),
                    _mm256_sub_epi32(vDel, vGap));
            vIns = _mm256_max_epi32(
                    _mm256_sub_epi32(vWscore, vOpen),
                    _mm256_sub_epi32(vIns, vGap));
            vMat = _mm256_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]],
                    matrow4[s2[j-4]],
                    matrow5[s2[j-5]],
                    matrow6[s2[j-6]],
                    matrow7[s2[j-7]]
                    );
            vNWscore = _mm256_add_epi32(vNWscore, vMat);
            vWscore = _mm256_max_epi32(vNWscore, vIns);
            vWscore = _mm256_max_epi32(vWscore, vDel);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-7] = (int32_t)_mm256_extract_epi32(vWscore,0);
            del_pr[j-7] = (int32_t)_mm256_extract_epi32(vDel,0);
            vJ = _mm256_add_epi32(vJ, vOne);
        }
        for (j=s2Len-1; j<s2Len+PAD; ++j) {
            __m256i vMat;
            __m256i vNWscore = vNscore;
            vNscore = vshift32(vWscore, tbl_pr[j]);
            vDel = vshift32(vDel, del_pr[j]);
            vDel = _mm256_max_epi32(
                    _mm256_sub_epi32(vNscore, vOpen),
                    _mm256_sub_epi32(vDel, vGap));
            vIns = _mm256_max_epi32(
                    _mm256_sub_epi32(vWscore, vOpen),
                    _mm256_sub_epi32(vIns, vGap));
            vMat = _mm256_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]],
                    matrow4[s2[j-4]],
                    matrow5[s2[j-5]],
                    matrow6[s2[j-6]],
                    matrow7[s2[j-7]]
                    );
            vNWscore = _mm256_add_epi32(vNWscore, vMat);
            vWscore = _mm256_max_epi32(vNWscore, vIns);
            vWscore = _mm256_max_epi32(vWscore, vDel);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-7] = (int32_t)_mm256_extract_epi32(vWscore,0);
            del_pr[j-7] = (int32_t)_mm256_extract_epi32(vDel,0);
            /* as minor diagonal vector passes across the j limit
             * boundary, extract the last value of the row */
            {
                __m256i cond_j = _mm256_cmpeq_epi32(vJ, vJLimit1);
                __m256i cond_max = _mm256_cmpgt_epi32(vWscore, vMaxScore);
                __m256i cond_all = _mm256_and_si256(cond_max, cond_j);
                vMaxScore = _mm256_blendv_epi8(vMaxScore, vWscore, cond_all);
            }
            vJ = _mm256_add_epi32(vJ, vOne);
        }
        vI = _mm256_add_epi32(vI, vN);
    }
    for (/*i=?*/; i<s1Len; i+=N) {
        __m256i vNscore = vNegInf0;
        __m256i vWscore = vNegInf0;
        __m256i vIns = vNegInf;
        __m256i vDel = vNegInf;
        __m256i vJ = vJreset;
        const int * const restrict matrow0 = matrix[s1[i+0]];
        const int * const restrict matrow1 = matrix[s1[i+1]];
        const int * const restrict matrow2 = matrix[s1[i+2]];
        const int * const restrict matrow3 = matrix[s1[i+3]];
        const int * const restrict matrow4 = matrix[s1[i+4]];
        const int * const restrict matrow5 = matrix[s1[i+5]];
        const int * const restrict matrow6 = matrix[s1[i+6]];
        const int * const restrict matrow7 = matrix[s1[i+7]];
        __m256i vIltLimit = _mm256_cmplt_epi32(vI, vILimit);
        __m256i vIeqLimit1 = _mm256_cmpeq_epi32(vI, vILimit1);
        /* iterate over database sequence */
        for (j=0; j<N; ++j) {
            __m256i vMat;
            __m256i vNWscore = vNscore;
            vNscore = vshift32(vWscore, tbl_pr[j]);
            vDel = vshift32(vDel, del_pr[j]);
            vDel = _mm256_max_epi32(
                    _mm256_sub_epi32(vNscore, vOpen),
                    _mm256_sub_epi32(vDel, vGap));
            vIns = _mm256_max_epi32(
                    _mm256_sub_epi32(vWscore, vOpen),
                    _mm256_sub_epi32(vIns, vGap));
            vMat = _mm256_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]],
                    matrow4[s2[j-4]],
                    matrow5[s2[j-5]],
                    matrow6[s2[j-6]],
                    matrow7[s2[j-7]]
                    );
            vNWscore = _mm256_add_epi32(vNWscore, vMat);
            vWscore = _mm256_max_epi32(vNWscore, vIns);
            vWscore = _mm256_max_epi32(vWscore, vDel);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m256i cond = _mm256_cmpeq_epi32(vJ,vNegOne);
                vWscore = _mm256_andnot_si256(cond, vWscore);
                vDel = _mm256_blendv_epi8(vDel, vNegInf, cond);
                vIns = _mm256_blendv_epi8(vIns, vNegInf, cond);
            }
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-7] = (int32_t)_mm256_extract_epi32(vWscore,0);
            del_pr[j-7] = (int32_t)_mm256_extract_epi32(vDel,0);
            /* as minor diagonal vector passes across the i limit
             * boundary, extract the last value of the column */
            {
                __m256i cond_i = _mm256_and_si256(
                        vIeqLimit1,
                        _mm256_cmpgt_epi32(vJ, vNegOne));
                __m256i cond_max = _mm256_cmpgt_epi32(vWscore, vMaxScore);
                __m256i cond_all = _mm256_and_si256(cond_max, cond_i);
                vMaxScore = _mm256_blendv_epi8(vMaxScore, vWscore, cond_all);
            }
            vJ = _mm256_add_epi32(vJ, vOne);
        }
        for (j=N; j<s2Len-1; ++j) {
            __m256i vMat;
            __m256i vNWscore = vNscore;
            vNscore = vshift32(vWscore, tbl_pr[j]);
            vDel = vshift32(vDel, del_pr[j]);
            vDel = _mm256_max_epi32(
                    _mm256_sub_epi32(vNscore, vOpen),
                    _mm256_sub_epi32(vDel, vGap));
            vIns = _mm256_max_epi32(
                    _mm256_sub_epi32(vWscore, vOpen),
                    _mm256_sub_epi32(vIns, vGap));
            vMat = _mm256_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]],
                    matrow4[s2[j-4]],
                    matrow5[s2[j-5]],
                    matrow6[s2[j-6]],
                    matrow7[s2[j-7]]
                    );
            vNWscore = _mm256_add_epi32(vNWscore, vMat);
            vWscore = _mm256_max_epi32(vNWscore, vIns);
            vWscore = _mm256_max_epi32(vWscore, vDel);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-7] = (int32_t)_mm256_extract_epi32(vWscore,0);
            del_pr[j-7] = (int32_t)_mm256_extract_epi32(vDel,0);
            /* as minor diagonal vector passes across the i limit
             * boundary, extract the last value of the column */
            {
                __m256i cond_i = vIeqLimit1;
                __m256i cond_max = _mm256_cmpgt_epi32(vWscore, vMaxScore);
                __m256i cond_all = _mm256_and_si256(cond_max, cond_i);
                vMaxScore = _mm256_blendv_epi8(vMaxScore, vWscore, cond_all);
            }
            vJ = _mm256_add_epi32(vJ, vOne);
        }
        for (j=s2Len-1; j<s2Len+PAD; ++j) {
            __m256i vMat;
            __m256i vNWscore = vNscore;
            vNscore = vshift32(vWscore, tbl_pr[j]);
            vDel = vshift32(vDel, del_pr[j]);
            vDel = _mm256_max_epi32(
                    _mm256_sub_epi32(vNscore, vOpen),
                    _mm256_sub_epi32(vDel, vGap));
            vIns = _mm256_max_epi32(
                    _mm256_sub_epi32(vWscore, vOpen),
                    _mm256_sub_epi32(vIns, vGap));
            vMat = _mm256_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]],
                    matrow4[s2[j-4]],
                    matrow5[s2[j-5]],
                    matrow6[s2[j-6]],
                    matrow7[s2[j-7]]
                    );
            vNWscore = _mm256_add_epi32(vNWscore, vMat);
            vWscore = _mm256_max_epi32(vNWscore, vIns);
            vWscore = _mm256_max_epi32(vWscore, vDel);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-7] = (int32_t)_mm256_extract_epi32(vWscore,0);
            del_pr[j-7] = (int32_t)_mm256_extract_epi32(vDel,0);
            /* as minor diagonal vector passes across the i or j limit
             * boundary, extract the last value of the column or row */
            {
                __m256i cond_j = _mm256_and_si256(
                        vIltLimit,
                        _mm256_cmpeq_epi32(vJ, vJLimit1));
                __m256i cond_i = _mm256_and_si256(
                        vIeqLimit1,
                        _mm256_cmplt_epi32(vJ, vJLimit));
                __m256i cond_max = _mm256_cmpgt_epi32(vWscore, vMaxScore);
                __m256i cond_all = _mm256_and_si256(cond_max,
                        _mm256_or_si256(cond_i, cond_j));
                vMaxScore = _mm256_blendv_epi8(vMaxScore, vWscore, cond_all);
            }
            vJ = _mm256_add_epi32(vJ, vOne);
        }
        vI = _mm256_add_epi32(vI, vN);
    }

    /* max in vMaxScore */
    for (i=0; i<N; ++i) {
        int32_t value;
        value = (int32_t) _mm256_extract_epi32(vMaxScore, 7);
        if (value > score) {
            score = value;
        }
        vMaxScore = lshift32(vMaxScore);
    }

    result->score = score;

    parasail_free(_del_pr);
    parasail_free(_tbl_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}

