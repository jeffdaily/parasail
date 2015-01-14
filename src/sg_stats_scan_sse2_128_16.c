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

#define NEG_INF_16 (INT16_MIN/(int16_t)(2))
#define MAX(a,b) ((a)>(b)?(a):(b))

#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[(0*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 0);
    array[(1*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 1);
    array[(2*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 2);
    array[(3*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 3);
    array[(4*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 4);
    array[(5*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 5);
    array[(6*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 6);
    array[(7*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 7);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME sg_stats_table_scan_sse2_128_16
#else
#define FNAME sg_stats_scan_sse2_128_16
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
    const int32_t segWidth = 8; /* number of values in vector unit */
    int32_t segNum = 0;
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
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
    __m128i vGapO = _mm_set1_epi16(open);
    __m128i vGapE = _mm_set1_epi16(gap);
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi16(1);
    __m128i vNegInf = _mm_set1_epi16(NEG_INF_16);
    int16_t score = NEG_INF_16;
    int16_t matches = 0;
    int16_t length = 0;
#ifdef PARASAIL_TABLE
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
                __m128i_16_t t;
                __m128i_16_t s;
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

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m128i_16_t h;
            __m128i_16_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = 0;
                e.v[segNum] = NEG_INF_16;
            }
            _mm_store_si128(&pvH[index], h.m);
            _mm_store_si128(&pvE[index], e.m);
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

        /* calculate E */
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vE = _mm_max_epi16(
                    _mm_sub_epi16(vE, vGapE),
                    _mm_sub_epi16(vH, vGapO));
            _mm_store_si128(pvE+i, vE);
        }

        /* calculate Ht */
        vH = _mm_slli_si128(_mm_load_si128(pvH+(segLen-1)), 2);
        vMp= _mm_slli_si128(_mm_load_si128(pvM+(segLen-1)), 2);
        vLp= _mm_slli_si128(_mm_load_si128(pvL+(segLen-1)), 2);
        vLp= _mm_add_epi16(vLp, vOne);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        pvC = pvPm+ MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            /* compute */
            vH = _mm_add_epi16(vH, vW);
            vHt = _mm_max_epi16(vH, vE);
            /* statistics */
            vC = _mm_load_si128(pvC+i);
            vMp = _mm_add_epi16(vMp, vC);
            vEx = _mm_cmpgt_epi16(vE, vH);
            vM = _mm_load_si128(pvM+i);
            vL = _mm_load_si128(pvL+i);
            vL = _mm_add_epi16(vL, vOne);
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
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 2);
        vFt = vNegInf;
        for (i=0; i<segLen; ++i) {
            vFt = _mm_sub_epi16(vFt, vGapE);
            vFt = _mm_max_epi16(vFt, vHt);
            vHt = _mm_load_si128(pvHt+i);
        }
        {
            __m128i_16_t tmp;
            tmp.m = vFt;
            tmp.v[1] = MAX(tmp.v[0]-segLen*gap, tmp.v[1]);
            tmp.v[2] = MAX(tmp.v[1]-segLen*gap, tmp.v[2]);
            tmp.v[3] = MAX(tmp.v[2]-segLen*gap, tmp.v[3]);
            tmp.v[4] = MAX(tmp.v[3]-segLen*gap, tmp.v[4]);
            tmp.v[5] = MAX(tmp.v[4]-segLen*gap, tmp.v[5]);
            tmp.v[6] = MAX(tmp.v[5]-segLen*gap, tmp.v[6]);
            tmp.v[7] = MAX(tmp.v[6]-segLen*gap, tmp.v[7]);
            vFt = tmp.m;
        }
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 2);
        vFt = _mm_slli_si128(vFt, 2);
        vFt = _mm_insert_epi16(vFt, NEG_INF_16, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_sub_epi16(vFt, vGapE);
            vFt = _mm_max_epi16(vFt, vHt);
            vHt = _mm_load_si128(pvHt+i);
            _mm_store_si128(pvFt+i, vFt);
        }

        /* calculate H,M,L */
        vMp = vZero;
        vLp = vOne;
        vC = _mm_cmpeq_epi16(vZero, vZero); /* check if prefix sum is needed */
        vC = _mm_srli_si128(vC, 2); /* zero out last value */
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vHt = _mm_load_si128(pvHt+i);
            vFt = _mm_load_si128(pvFt+i);
            /* compute */
            vFt = _mm_sub_epi16(vFt, vGapO);
            vH = _mm_max_epi16(vHt, vFt);
            /* statistics */
            vEx = _mm_load_si128(pvEx+i);
            vMt = _mm_load_si128(pvMt+i);
            vLt = _mm_load_si128(pvLt+i);
            vEx = _mm_or_si128(
                    _mm_and_si128(vEx, _mm_cmpeq_epi16(vHt, vFt)),
                    _mm_cmplt_epi16(vHt, vFt));
            vM = _mm_and_si128(vEx, vMp);
            vM = _mm_or_si128(vM, _mm_andnot_si128(vEx, vMt));
            vMp = vM;
            vL = _mm_and_si128(vEx, vLp);
            vL = _mm_or_si128(vL, _mm_andnot_si128(vEx, vLt));
            vLp = _mm_add_epi16(vL, vOne);
            vC = _mm_and_si128(vC, vEx);
            /* store results */
            _mm_store_si128(pvH+i, vH);
            _mm_store_si128(pvEx+i, vEx);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif
        }
        {
            vLp = _mm_sub_epi16(vLp, vOne);
            {
                __m128i_16_t uMp, uLp, uC;
                uC.m = vC;
                uMp.m = vMp;
                uMp.v[1] = uC.v[1] ? uMp.v[0] : uMp.v[1];
                uMp.v[2] = uC.v[2] ? uMp.v[1] : uMp.v[2];
                uMp.v[3] = uC.v[3] ? uMp.v[2] : uMp.v[3];
                uMp.v[4] = uC.v[4] ? uMp.v[3] : uMp.v[4];
                uMp.v[5] = uC.v[5] ? uMp.v[4] : uMp.v[5];
                uMp.v[6] = uC.v[6] ? uMp.v[5] : uMp.v[6];
                uMp.v[7] = uC.v[7] ? uMp.v[6] : uMp.v[7];
                vMp = uMp.m;
                uLp.m = vLp;
                uLp.v[1] = uC.v[1] ? uLp.v[1] + uLp.v[0] : uLp.v[1];
                uLp.v[2] = uC.v[2] ? uLp.v[2] + uLp.v[1] : uLp.v[2];
                uLp.v[3] = uC.v[3] ? uLp.v[3] + uLp.v[2] : uLp.v[3];
                uLp.v[4] = uC.v[4] ? uLp.v[4] + uLp.v[3] : uLp.v[4];
                uLp.v[5] = uC.v[5] ? uLp.v[5] + uLp.v[4] : uLp.v[5];
                uLp.v[6] = uC.v[6] ? uLp.v[6] + uLp.v[5] : uLp.v[6];
                uLp.v[7] = uC.v[7] ? uLp.v[7] + uLp.v[6] : uLp.v[7];
                vLp = uLp.m;
            }
            vLp = _mm_add_epi16(vLp, vOne);
        }
        /* final pass for M,L */
        vMp = _mm_slli_si128(vMp, 2);
        vLp = _mm_slli_si128(vLp, 2);
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
            vLp = _mm_add_epi16(vL, vOne);
            /* store results */
            _mm_store_si128(pvM+i, vM);
            _mm_store_si128(pvL+i, vL);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->matches_table, vM, i, segLen, j, s2Len);
            arr_store_si128(result->length_table, vL, i, segLen, j, s2Len);
#endif
        }

        /* extract last value from column */
        {
            int16_t value;
            vH = _mm_load_si128(pvH + offset);
            vM = _mm_load_si128(pvM + offset);
            vL = _mm_load_si128(pvL + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128(vH, 2);
                vM = _mm_slli_si128(vM, 2);
                vL = _mm_slli_si128(vL, 2);
            }
            value = (int16_t) _mm_extract_epi16(vH, 7);
            if (value > score) {
                score = value;
                matches = (int16_t) _mm_extract_epi16(vM, 7);
                length = (int16_t) _mm_extract_epi16(vL, 7);
            }
        }
    }

    /* max of last column */
    {
        __m128i vMaxH = vNegInf;
        __m128i vMaxM = vZero;
        __m128i vMaxL = vZero;
        __m128i vQLimit = _mm_set1_epi16(s1Len);
        __m128i vQIndex = _mm_set_epi16(
                7*segLen,
                6*segLen,
                5*segLen,
                4*segLen,
                3*segLen,
                2*segLen,
                1*segLen,
                0*segLen);

        for (i=0; i<segLen; ++i) {
            __m128i vH = _mm_load_si128(pvH + i);
            __m128i vM = _mm_load_si128(pvM + i);
            __m128i vL = _mm_load_si128(pvL + i);
            __m128i cond_lmt = _mm_cmplt_epi16(vQIndex, vQLimit);
            __m128i cond_max = _mm_cmpgt_epi16(vH, vMaxH);
            __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
            vMaxH = _mm_andnot_si128(cond_all, vMaxH);
            vMaxH = _mm_or_si128(vMaxH, _mm_and_si128(cond_all, vH));
            vMaxM = _mm_andnot_si128(cond_all, vMaxM);
            vMaxM = _mm_or_si128(vMaxM, _mm_and_si128(cond_all, vM));
            vMaxL = _mm_andnot_si128(cond_all, vMaxL);
            vMaxL = _mm_or_si128(vMaxL, _mm_and_si128(cond_all, vL));
            vQIndex = _mm_add_epi16(vQIndex, vOne);
        }

        /* max in vec */
        for (j=0; j<segWidth; ++j) {
            int16_t value = (int16_t) _mm_extract_epi16(vMaxH, 7);
            if (value > score) {
                score = value;
                matches = (int16_t) _mm_extract_epi16(vMaxM, 7);
                length = (int16_t) _mm_extract_epi16(vMaxL, 7);
            }
            vMaxH = _mm_slli_si128(vMaxH, 2);
            vMaxM = _mm_slli_si128(vMaxM, 2);
            vMaxL = _mm_slli_si128(vMaxL, 2);
        }
    }

    result->score = score;
    result->matches = matches;
    result->length = length;

    parasail_free(pvL);
    parasail_free(pvM);
    parasail_free(pvH);
    parasail_free(pvEx);
    parasail_free(pvLt);
    parasail_free(pvMt);
    parasail_free(pvFt);
    parasail_free(pvHt);
    parasail_free(pvE);
    parasail_free(pvPm);
    parasail_free(pvP);

    return result;
}

