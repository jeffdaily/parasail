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

#define NEG_INF_8 (INT8_MIN)
#define MAX(a,b) ((a)>(b)?(a):(b))

/* sse2 does not have _mm_insert_epi8, emulate it */
static inline __m128i _mm_insert_epi8(__m128i a, int8_t i, int imm) {
    __m128i_8_t tmp;
    tmp.m = a;
    tmp.v[imm] = i;
    return tmp.m;
}

/* sse2 does not have _mm_extract_epi8, emulate it */
static inline int8_t _mm_extract_epi8(__m128i a, int imm) {
    __m128i_8_t tmp;
    tmp.m = a;
    return tmp.v[imm];
}

/* sse2 does not have _mm_max_epi8, emulate it */
static inline __m128i _mm_max_epi8(__m128i a, __m128i b) {
    __m128i mask = _mm_cmpgt_epi8(a,b);
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
    array[(0*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 0);
    array[(1*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 1);
    array[(2*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 2);
    array[(3*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 3);
    array[(4*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 4);
    array[(5*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 5);
    array[(6*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 6);
    array[(7*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 7);
    array[(8*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 8);
    array[(9*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 9);
    array[(10*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 10);
    array[(11*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 11);
    array[(12*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 12);
    array[(13*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 13);
    array[(14*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 14);
    array[(15*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 15);

}

static inline void arr_store_si128_bias(
        int *array,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen,
        int32_t bias)
{
    array[(0*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 0);
    array[(1*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 1);
    array[(2*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 2);
    array[(3*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 3);
    array[(4*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 4);
    array[(5*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 5);
    array[(6*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 6);
    array[(7*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 7);
    array[(8*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 8);
    array[(9*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 9);
    array[(10*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 10);
    array[(11*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 11);
    array[(12*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 12);
    array[(13*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 13);
    array[(14*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 14);
    array[(15*seglen+t)*dlen + d] = bias + (int8_t)_mm_extract_epi8(vH, 15);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME sw_stats_table_scan_sse2_128_8
#else
#define FNAME sw_stats_scan_sse2_128_8
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
    const int32_t segWidth = 16; /* number of values in vector unit */
    int32_t segNum = 0;
    int32_t segLen = (s1Len + segWidth - 1) / segWidth;
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
    __m128i vGapO = _mm_set1_epi8(open);
    __m128i vGapE = _mm_set1_epi8(gap);
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi8(1);
    __m128i vOne16 = _mm_set1_epi16(1);
    int8_t bias = 127;
    int ibias = 127;
    __m128i vNegBias = _mm_set1_epi8(-bias);
    __m128i vSaturationCheck = _mm_setzero_si128();
    __m128i vNegLimit = _mm_set1_epi8(INT8_MIN);
    __m128i vPosLimit = _mm_set1_epi8(INT8_MAX);
    int8_t score = NEG_INF_8;
    int matches = 0;
    int length = 0;
    __m128i vMaxH = vZero;
    __m128i vMaxM = vZero;
    __m128i vMaxL = vZero;
    __m128i vQLimit16 = _mm_set1_epi16(s1Len);
    __m128i vQIndexHi16_reset = _mm_set_epi16(
            15*segLen,
            14*segLen,
            13*segLen,
            12*segLen,
            11*segLen,
            10*segLen,
            9*segLen,
            8*segLen);
    __m128i vQIndexLo16_reset = _mm_set_epi16(
            7*segLen,
            6*segLen,
            5*segLen,
            4*segLen,
            3*segLen,
            2*segLen,
            1*segLen,
            0*segLen);
#if PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif

    /* Generate query profile and match profile.
     * Rearrange query sequence & calculate the weight of match/mismatch.
     * Don't alias. */
    {
        int32_t index = 0;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                __m128i_8_t t;
                __m128i_8_t s;
                j = i;
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

    /* initialize H E M L */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m128i_8_t h;
            __m128i_8_t e;
            __m128i_8_t m;
            __m128i_8_t l;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = 0;
                e.v[segNum] = NEG_INF_8;
                m.v[segNum] = -bias;
                l.v[segNum] = -bias;
            }
            _mm_store_si128(&pvH[index], h.m);
            _mm_store_si128(&pvE[index], e.m);
            _mm_store_si128(&pvM[index], m.m);
            _mm_store_si128(&pvL[index], l.m);
            ++index;
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
        __m128i cond_max;
        __m128i vQIndexHi16;
        __m128i vQIndexLo16;

        /* calculate E */
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vE = _mm_max_epi8(
                    _mm_subs_epi8(vE, vGapE),
                    _mm_subs_epi8(vH, vGapO));
            _mm_store_si128(pvE+i, vE);
        }

        /* calculate Ht */
        vH = _mm_slli_si128(_mm_load_si128(pvH+(segLen-1)), 1);
        vMp= _mm_slli_si128(_mm_load_si128(pvM+(segLen-1)), 1);
        vMp= _mm_insert_epi8(vMp, -bias, 0);
        vLp= _mm_slli_si128(_mm_load_si128(pvL+(segLen-1)), 1);
        vLp= _mm_insert_epi8(vLp, -bias, 0);
        vLp= _mm_adds_epi8(vLp, vOne);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        pvC = pvPm+ MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            /* compute */
            vH = _mm_adds_epi8(vH, vW);
            vHt = _mm_max_epi8(vH, vE);
            vHt = _mm_max_epi8(vHt, vZero);
            /* statistics */
            vC = _mm_load_si128(pvC+i);
            vMp = _mm_adds_epi8(vMp, vC);
            vEx = _mm_cmpgt_epi8(vE, vH);
            vM = _mm_load_si128(pvM+i);
            vL = _mm_load_si128(pvL+i);
            vL = _mm_adds_epi8(vL, vOne);
            vMt = _mm_and_si128(vEx, vM);
            vMt = _mm_or_si128(vMt, _mm_andnot_si128(vEx, vMp));
            vLt = _mm_and_si128(vEx, vL);
            vLt = _mm_or_si128(vLt, _mm_andnot_si128(vEx, vLp));
            cond_max = _mm_cmpeq_epi8(vHt, vZero);
            vEx = _mm_andnot_si128(cond_max, vEx);
            vMt = _mm_andnot_si128(cond_max, vMt);
            vMt = _mm_or_si128(vMt, _mm_and_si128(cond_max, vNegBias));
            vLt = _mm_andnot_si128(cond_max, vLt);
            vLt = _mm_or_si128(vLt, _mm_and_si128(cond_max, vNegBias));
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
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 1);
        vFt = _mm_set1_epi8(NEG_INF_8);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_subs_epi8(vFt, vGapE);
            vFt = _mm_max_epi8(vFt, vHt);
            vHt = _mm_load_si128(pvHt+i);
        }
        {
            __m128i_8_t tmp;
            tmp.m = vFt;
            tmp.v[1] = MAX(tmp.v[0]-segLen*gap, tmp.v[1]);
            tmp.v[2] = MAX(tmp.v[1]-segLen*gap, tmp.v[2]);
            tmp.v[3] = MAX(tmp.v[2]-segLen*gap, tmp.v[3]);
            tmp.v[4] = MAX(tmp.v[3]-segLen*gap, tmp.v[4]);
            tmp.v[5] = MAX(tmp.v[4]-segLen*gap, tmp.v[5]);
            tmp.v[6] = MAX(tmp.v[5]-segLen*gap, tmp.v[6]);
            tmp.v[7] = MAX(tmp.v[6]-segLen*gap, tmp.v[7]);
            tmp.v[8]  = MAX(tmp.v[7] -segLen*gap, tmp.v[8]);
            tmp.v[9]  = MAX(tmp.v[8] -segLen*gap, tmp.v[9]);
            tmp.v[10] = MAX(tmp.v[9] -segLen*gap, tmp.v[10]);
            tmp.v[11] = MAX(tmp.v[10]-segLen*gap, tmp.v[11]);
            tmp.v[12] = MAX(tmp.v[11]-segLen*gap, tmp.v[12]);
            tmp.v[13] = MAX(tmp.v[12]-segLen*gap, tmp.v[13]);
            tmp.v[14] = MAX(tmp.v[13]-segLen*gap, tmp.v[14]);
            tmp.v[15] = MAX(tmp.v[14]-segLen*gap, tmp.v[15]);
            vFt = tmp.m;
        }
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 1);
        vFt = _mm_slli_si128(vFt, 1);
        vFt = _mm_insert_epi8(vFt, NEG_INF_8, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_subs_epi8(vFt, vGapE);
            vFt = _mm_max_epi8(vFt, vHt);
            vHt = _mm_load_si128(pvHt+i);
            _mm_store_si128(pvFt+i, vFt);
        }

        /* calculate H,M,L */
        vMp = vNegBias;
        vLp = _mm_adds_epi8(vNegBias, vOne);
        vC = _mm_cmpeq_epi8(vZero, vZero); /* check if prefix sum is needed */
        vC = _mm_srli_si128(vC, 1); /* zero out last value */
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vHt = _mm_load_si128(pvHt+i);
            vFt = _mm_load_si128(pvFt+i);
            /* compute */
            vFt = _mm_subs_epi8(vFt, vGapO);
            vH = _mm_max_epi8(vHt, vFt);
            /* statistics */
            vEx = _mm_load_si128(pvEx+i);
            vMt = _mm_load_si128(pvMt+i);
            vLt = _mm_load_si128(pvLt+i);
            vEx = _mm_or_si128(
                    _mm_and_si128(vEx, _mm_cmpeq_epi8(vHt, vFt)),
                    _mm_cmplt_epi8(vHt, vFt));
            vM = _mm_and_si128(vEx, vMp);
            vM = _mm_or_si128(vM, _mm_andnot_si128(vEx, vMt));
            vMp = vM;
            vL = _mm_and_si128(vEx, vLp);
            vL = _mm_or_si128(vL, _mm_andnot_si128(vEx, vLt));
            vLp = _mm_adds_epi8(vL, vOne);
            vC = _mm_and_si128(vC, vEx);
            /* store results */
            _mm_store_si128(pvH+i, vH);
            _mm_store_si128(pvEx+i, vEx);
            /* check for saturation */
            {
                vSaturationCheck = _mm_or_si128(vSaturationCheck,
                        _mm_or_si128(
                            _mm_cmpeq_epi8(vH, vNegLimit),
                            _mm_cmpeq_epi8(vH, vPosLimit)));
            }
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif
        }
#if PREFIX_SUM_CHECK
        /* check if local prefix sum for L is needed */
        if (_mm_movemask_epi8(vC))
#endif
        {
            vLp = _mm_subs_epi8(vLp, vOne);
            {
                __m128i_8_t uMp, uLp, uC;
                uC.m = vC;
                uMp.m = vMp;
                uMp.v[ 1] = uC.v[ 1] ? uMp.v[ 0] : uMp.v[ 1];
                uMp.v[ 2] = uC.v[ 2] ? uMp.v[ 1] : uMp.v[ 2];
                uMp.v[ 3] = uC.v[ 3] ? uMp.v[ 2] : uMp.v[ 3];
                uMp.v[ 4] = uC.v[ 4] ? uMp.v[ 3] : uMp.v[ 4];
                uMp.v[ 5] = uC.v[ 5] ? uMp.v[ 4] : uMp.v[ 5];
                uMp.v[ 6] = uC.v[ 6] ? uMp.v[ 5] : uMp.v[ 6];
                uMp.v[ 7] = uC.v[ 7] ? uMp.v[ 6] : uMp.v[ 7];
                uMp.v[ 8] = uC.v[ 8] ? uMp.v[ 7] : uMp.v[ 8];
                uMp.v[ 9] = uC.v[ 9] ? uMp.v[ 8] : uMp.v[ 9];
                uMp.v[10] = uC.v[10] ? uMp.v[ 9] : uMp.v[10];
                uMp.v[11] = uC.v[11] ? uMp.v[10] : uMp.v[11];
                uMp.v[12] = uC.v[12] ? uMp.v[11] : uMp.v[12];
                uMp.v[13] = uC.v[13] ? uMp.v[12] : uMp.v[13];
                uMp.v[14] = uC.v[14] ? uMp.v[13] : uMp.v[14];
                uMp.v[15] = uC.v[15] ? uMp.v[14] : uMp.v[15];
                vMp = uMp.m;
                uLp.m = vLp;
                uLp.v[ 1] = uC.v[ 1] ? uLp.v[ 1] + uLp.v[ 0] : uLp.v[ 1];
                uLp.v[ 2] = uC.v[ 2] ? uLp.v[ 2] + uLp.v[ 1] : uLp.v[ 2];
                uLp.v[ 3] = uC.v[ 3] ? uLp.v[ 3] + uLp.v[ 2] : uLp.v[ 3];
                uLp.v[ 4] = uC.v[ 4] ? uLp.v[ 4] + uLp.v[ 3] : uLp.v[ 4];
                uLp.v[ 5] = uC.v[ 5] ? uLp.v[ 5] + uLp.v[ 4] : uLp.v[ 5];
                uLp.v[ 6] = uC.v[ 6] ? uLp.v[ 6] + uLp.v[ 5] : uLp.v[ 6];
                uLp.v[ 7] = uC.v[ 7] ? uLp.v[ 7] + uLp.v[ 6] : uLp.v[ 7];
                uLp.v[ 8] = uC.v[ 8] ? uLp.v[ 8] + uLp.v[ 7] : uLp.v[ 8];
                uLp.v[ 9] = uC.v[ 9] ? uLp.v[ 9] + uLp.v[ 8] : uLp.v[ 9];
                uLp.v[10] = uC.v[10] ? uLp.v[10] + uLp.v[ 9] : uLp.v[10];
                uLp.v[11] = uC.v[11] ? uLp.v[11] + uLp.v[10] : uLp.v[11];
                uLp.v[12] = uC.v[12] ? uLp.v[12] + uLp.v[11] : uLp.v[12];
                uLp.v[13] = uC.v[13] ? uLp.v[13] + uLp.v[12] : uLp.v[13];
                uLp.v[14] = uC.v[14] ? uLp.v[14] + uLp.v[13] : uLp.v[14];
                uLp.v[15] = uC.v[15] ? uLp.v[15] + uLp.v[14] : uLp.v[15];
                vLp = uLp.m;
            }
            vLp = _mm_adds_epi8(vLp, vOne);
        }
        /* final pass for M,L */
        vQIndexLo16 = vQIndexLo16_reset;
        vQIndexHi16 = vQIndexHi16_reset;
        vMp = _mm_slli_si128(vMp, 1);
        vMp = _mm_insert_epi8(vMp, -bias, 0);
        vLp = _mm_slli_si128(vLp, 1);
        vLp = _mm_insert_epi8(vLp, -bias, 0);
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vH = _mm_load_si128(pvH+i);
            /* statistics */
            vEx = _mm_load_si128(pvEx+i);
            vMt = _mm_load_si128(pvMt+i);
            vLt = _mm_load_si128(pvLt+i);
            vM = _mm_and_si128(vEx, vMp);
            vM = _mm_or_si128(vM, _mm_andnot_si128(vEx, vMt));
            vMp = vM;
            vL = _mm_and_si128(vEx, vLp);
            vL = _mm_or_si128(vL, _mm_andnot_si128(vEx, vLt));
            vLp = _mm_adds_epi8(vL, vOne);
            /* store results */
            _mm_store_si128(pvM+i, vM);
            _mm_store_si128(pvL+i, vL);
            /* check for saturation */
            {
                vSaturationCheck = _mm_or_si128(vSaturationCheck,
                        _mm_or_si128(
                            _mm_cmpeq_epi8(vM, vNegLimit),
                            _mm_cmpeq_epi8(vM, vPosLimit)));
                vSaturationCheck = _mm_or_si128(vSaturationCheck,
                        _mm_or_si128(
                            _mm_cmpeq_epi8(vL, vNegLimit),
                            _mm_cmpeq_epi8(vL, vPosLimit)));
            }
#ifdef PARASAIL_TABLE
            arr_store_si128_bias(result->matches_table, vM, i, segLen, j, s2Len, bias);
            arr_store_si128_bias(result->length_table, vL, i, segLen, j, s2Len, bias);
#endif
            /* update max vector seen so far */
            {
                __m128i cond_lmt = _mm_packs_epi16(
                        _mm_cmplt_epi16(vQIndexLo16, vQLimit16),
                        _mm_cmplt_epi16(vQIndexHi16, vQLimit16));
                __m128i cond_max = _mm_cmpgt_epi8(vH, vMaxH);
                __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
                vMaxH = _mm_andnot_si128(cond_all, vMaxH);
                vMaxH = _mm_or_si128(vMaxH, _mm_and_si128(cond_all, vH));
                vMaxM = _mm_andnot_si128(cond_all, vMaxM);
                vMaxM = _mm_or_si128(vMaxM, _mm_and_si128(cond_all, vM));
                vMaxL = _mm_andnot_si128(cond_all, vMaxL);
                vMaxL = _mm_or_si128(vMaxL, _mm_and_si128(cond_all, vL));
                vQIndexLo16 = _mm_adds_epi16(vQIndexLo16, vOne16);
                vQIndexHi16 = _mm_adds_epi16(vQIndexHi16, vOne16);
            }
        }
    }

    /* max in vec */
    for (j=0; j<segWidth; ++j) {
        int8_t value = (int8_t) _mm_extract_epi8(vMaxH, 15);
        if (value > score) {
            score = value;
            matches = ibias + (int)((int8_t)_mm_extract_epi8(vMaxM, 15));
            length = ibias + (int)((int8_t)_mm_extract_epi8(vMaxL, 15));
        }
        vMaxH = _mm_slli_si128(vMaxH, 1);
        vMaxM = _mm_slli_si128(vMaxM, 1);
        vMaxL = _mm_slli_si128(vMaxL, 1);
    }

    result->score = score;
    result->matches = matches;
    result->length = length;

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
