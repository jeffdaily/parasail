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

#include <stdint.h>
#include <stdlib.h>

#include <emmintrin.h>

#include "parasail.h"
#include "parasail_internal.h"
#include "parasail_internal_sse.h"
#include "blosum/blosum_map.h"

#define NEG_INF_32 (INT32_MIN/(int32_t)(2))
#define MAX(a,b) ((a)>(b)?(a):(b))

/* sse2 does not have _mm_insert_epi32, emulate it */
static inline __m128i _mm_insert_epi32(__m128i a, int32_t i, int imm) {
    __m128i_32_t tmp;
    tmp.m = a;
    tmp.v[imm] = i;
    return tmp.m;
}

/* sse2 does not have _mm_extract_epi32, emulate it */
static inline int32_t _mm_extract_epi32(__m128i a, int imm) {
    __m128i_32_t tmp;
    tmp.m = a;
    return tmp.v[imm];
}

/* sse2 does not have _mm_max_epi32, emulate it */
static inline __m128i _mm_max_epi32(__m128i a, __m128i b) {
    __m128i mask = _mm_cmpgt_epi32(a,b);
    a = _mm_and_si128(a,mask);
    b = _mm_andnot_si128(mask,b);
    return _mm_or_si128(a,b);
}

#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[(0*seglen+t)*dlen + d] = (int32_t)_mm_extract_epi32(vH, 0);
    array[(1*seglen+t)*dlen + d] = (int32_t)_mm_extract_epi32(vH, 1);
    array[(2*seglen+t)*dlen + d] = (int32_t)_mm_extract_epi32(vH, 2);
    array[(3*seglen+t)*dlen + d] = (int32_t)_mm_extract_epi32(vH, 3);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME nw_stats_table_scan_sse2_128_32
#else
#define FNAME nw_stats_scan_sse2_128_32
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    const int32_t n = 24; /* number of amino acids in table */
    const int32_t segWidth = 4; /* number of values in vector unit */
    int32_t segNum = 0;
    int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m128i* const restrict pvP  = parasail_memalign_m128i(16, n * segLen);
    __m128i* const restrict pvPm = parasail_memalign_m128i(16, n * segLen);
    __m128i* const restrict pvE  = parasail_memalign_m128i(16, segLen);
    __m128i* const restrict pvHt = parasail_memalign_m128i(16, segLen);
    __m128i* const restrict pvFt = parasail_memalign_m128i(16, segLen);
    __m128i* const restrict pvMt = parasail_memalign_m128i(16, segLen);
    __m128i* const restrict pvLt = parasail_memalign_m128i(16, segLen);
    __m128i* const restrict pvEx = parasail_memalign_m128i(16, segLen);
    __m128i* const restrict pvH  = parasail_memalign_m128i(16, segLen);
    __m128i* const restrict pvM  = parasail_memalign_m128i(16, segLen);
    __m128i* const restrict pvL  = parasail_memalign_m128i(16, segLen);
    int32_t* const restrict boundary = parasail_memalign_int32_t(16, s2Len+1);
    __m128i vGapO = _mm_set1_epi32(open);
    __m128i vGapE = _mm_set1_epi32(gap);
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi32(1);
    int32_t score = 0;
    int32_t matches = 0;
    int32_t length = 0;
#if PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif

    parasail_memset_m128i(pvM, vZero, segLen);
    parasail_memset_m128i(pvL, vZero, segLen);

    /* Generate query profile and match profile.
     * Rearrange query sequence & calculate the weight of match/mismatch.
     * Don't alias. */
    {
        int32_t index = 0;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                int32_t j = i;
                __m128i_32_t t;
                __m128i_32_t s;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = matrix[k][MAP_BLOSUM_[(unsigned char)s1[j]]];
                    s.v[segNum] = (k == MAP_BLOSUM_[(unsigned char)s1[j]]);
                    j += segLen;
                }
                _mm_store_si128(&pvP[index], t.m);
                _mm_store_si128(&pvPm[index], s.m);
                ++index;
            }
        }
    }

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m128i_32_t h;
            __m128i_32_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = -open-gap*(segNum*segLen+i);
                e.v[segNum] = NEG_INF_32;
            }
            _mm_store_si128(&pvH[index], h.m);
            _mm_store_si128(&pvE[index], e.m);
            ++index;
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            boundary[i] = -open-gap*(i-1);
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vE;
        __m128i vHt;
        __m128i vFt;
        __m128i vH;
        __m128i *pvW;
        __m128i vW;
        __m128i *pvC;
        __m128i vC;
        __m128i vM;
        __m128i vMp;
        __m128i vMt;
        __m128i vL;
        __m128i vLp;
        __m128i vLt;
        __m128i vEx;

        /* calculate E */
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vE = _mm_max_epi32(
                    _mm_sub_epi32(vE, vGapE),
                    _mm_sub_epi32(vH, vGapO));
            _mm_store_si128(pvE+i, vE);
        }

        /* calculate Ht */
        vH = _mm_slli_si128(_mm_load_si128(pvH+(segLen-1)), 4);
        vH = _mm_insert_epi32(vH, boundary[j], 0);
        vMp= _mm_slli_si128(_mm_load_si128(pvM+(segLen-1)), 4);
        vLp= _mm_slli_si128(_mm_load_si128(pvL+(segLen-1)), 4);
        vLp= _mm_add_epi32(vLp, vOne);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        pvC = pvPm+ MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            /* compute */
            vH = _mm_add_epi32(vH, vW);
            vHt = _mm_max_epi32(vH, vE);
            /* statistics */
            vC = _mm_load_si128(pvC+i);
            vMp = _mm_add_epi32(vMp, vC);
            vEx = _mm_cmpgt_epi32(vE, vH);
            vM = _mm_load_si128(pvM+i);
            vL = _mm_load_si128(pvL+i);
            vL = _mm_add_epi32(vL, vOne);
            vMt = _mm_and_si128(vEx, vM);
            vMt = _mm_or_si128(vMt, _mm_andnot_si128(vEx, vMp));
            vLt = _mm_and_si128(vEx, vL);
            vLt = _mm_or_si128(vLt, _mm_andnot_si128(vEx, vLp));
            /* store results */
            _mm_store_si128(pvHt+i, vHt);
            _mm_store_si128(pvEx+i, vEx);
            _mm_store_si128(pvMt+i, vMt);
            _mm_store_si128(pvLt+i, vLt);
            /* prep for next iteration */
            vH = _mm_load_si128(pvH+i);
            vMp = vM;
            vLp = vL;
        }

        /* calculate Ft */
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 4);
        vHt = _mm_insert_epi32(vHt, boundary[j+1], 0);
        vFt = _mm_set1_epi32(NEG_INF_32);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_sub_epi32(vFt, vGapE),
            vFt = _mm_max_epi32(vFt, vHt);
            vHt = _mm_load_si128(pvHt+i);
        }
        {
            __m128i_32_t tmp;
            tmp.m = vFt;
            tmp.v[1] = MAX(tmp.v[0]-segLen*gap, tmp.v[1]);
            tmp.v[2] = MAX(tmp.v[1]-segLen*gap, tmp.v[2]);
            tmp.v[3] = MAX(tmp.v[2]-segLen*gap, tmp.v[3]);
            vFt = tmp.m;
        }
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 4);
        vHt = _mm_insert_epi32(vHt, boundary[j+1], 0);
        vFt = _mm_slli_si128(vFt, 4);
        vFt = _mm_insert_epi32(vFt, NEG_INF_32, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_sub_epi32(vFt, vGapE),
            vFt = _mm_max_epi32(vFt, vHt);
            vHt = _mm_load_si128(pvHt+i);
            _mm_store_si128(pvFt+i, vFt);
        }

        /* calculate H,M,L */
        vMp = vZero;
        vLp = vOne;
        vC = _mm_cmpeq_epi32(vZero, vZero); /* check if prefix sum is needed */
        vC = _mm_srli_si128(vC, 4); /* zero out last value */
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vHt = _mm_load_si128(pvHt+i);
            vFt = _mm_load_si128(pvFt+i);
            /* compute */
            vFt = _mm_sub_epi32(vFt, vGapO);
            vH = _mm_max_epi32(vHt, vFt);
            /* statistics */
            vEx = _mm_load_si128(pvEx+i);
            vMt = _mm_load_si128(pvMt+i);
            vLt = _mm_load_si128(pvLt+i);
            vEx = _mm_or_si128(
                    _mm_and_si128(vEx, _mm_cmpeq_epi32(vHt, vFt)),
                    _mm_cmplt_epi32(vHt, vFt));
            vM = _mm_and_si128(vEx, vMp);
            vM = _mm_or_si128(vM, _mm_andnot_si128(vEx, vMt));
            vMp = vM;
            vL = _mm_and_si128(vEx, vLp);
            vL = _mm_or_si128(vL, _mm_andnot_si128(vEx, vLt));
            vLp = _mm_add_epi32(vL, vOne);
            vC = _mm_and_si128(vC, vEx);
            /* store results */
            _mm_store_si128(pvH+i, vH);
            _mm_store_si128(pvEx+i, vEx);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif
        }
#if PREFIX_SUM_CHECK
        /* check if local prefix sum for L is needed */
        if (_mm_movemask_epi8(vC))
#endif
        {
            vLp = _mm_sub_epi32(vLp, vOne);
            {
                __m128i_32_t uMp, uLp, uC;
                uC.m = vC;
                uMp.m = vMp;
                uMp.v[1] = uC.v[1] ? uMp.v[0] : uMp.v[1];
                uMp.v[2] = uC.v[2] ? uMp.v[1] : uMp.v[2];
                uMp.v[3] = uC.v[3] ? uMp.v[2] : uMp.v[3];
                vMp = uMp.m;
                uLp.m = vLp;
                uLp.v[1] = uC.v[1] ? uLp.v[1] + uLp.v[0] : uLp.v[1];
                uLp.v[2] = uC.v[2] ? uLp.v[2] + uLp.v[1] : uLp.v[2];
                uLp.v[3] = uC.v[3] ? uLp.v[3] + uLp.v[2] : uLp.v[3];
                vLp = uLp.m;
            }
            vLp = _mm_add_epi32(vLp, vOne);
        }
        /* final pass for M,L */
        vMp = _mm_slli_si128(vMp, 4);
        vLp = _mm_slli_si128(vLp, 4);
        for (i=0; i<segLen; ++i) {
            /* statistics */
            vEx = _mm_load_si128(pvEx+i);
            vMt = _mm_load_si128(pvMt+i);
            vLt = _mm_load_si128(pvLt+i);
            vM = _mm_and_si128(vEx, vMp);
            vM = _mm_or_si128(vM, _mm_andnot_si128(vEx, vMt));
            vMp = vM;
            vL = _mm_and_si128(vEx, vLp);
            vL = _mm_or_si128(vL, _mm_andnot_si128(vEx, vLt));
            vLp = _mm_add_epi32(vL, vOne);
            /* store results */
            _mm_store_si128(pvM+i, vM);
            _mm_store_si128(pvL+i, vL);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->matches_table, vM, i, segLen, j, s2Len);
            arr_store_si128(result->length_table, vL, i, segLen, j, s2Len);
#endif
        }
    }

    /* extract last value from the last column */
    {
        __m128i vH = _mm_load_si128(pvH + offset);
        __m128i vM = _mm_load_si128(pvM + offset);
        __m128i vL = _mm_load_si128(pvL + offset);
        for (k=0; k<position; ++k) {
            vH = _mm_slli_si128(vH, 4);
            vM = _mm_slli_si128(vM, 4);
            vL = _mm_slli_si128(vL, 4);
        }
        score = (int32_t) _mm_extract_epi32 (vH, 3);
        matches = (int32_t) _mm_extract_epi32 (vM, 3);
        length = (int32_t) _mm_extract_epi32 (vL, 3);
    }

    result->score = score;
    result->matches = matches;
    result->length = length;

    free(boundary);
    free(pvL);
    free(pvM);
    free(pvH);
    free(pvEx);
    free(pvLt);
    free(pvMt);
    free(pvFt);
    free(pvHt);
    free(pvE);
    free(pvPm);
    free(pvP);

    return result;
}

