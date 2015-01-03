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
#include <smmintrin.h>

#include "parasail.h"
#include "parasail_internal.h"
#include "parasail_internal_sse.h"
#include "blosum/blosum_map.h"

#define NEG_INF_8 (INT8_MIN)
#define MAX(a,b) ((a)>(b)?(a):(b))

#if PARASAIL_TABLE
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
#endif

#ifdef PARASAIL_TABLE
#define FNAME sg_table_scan_sse41_128_8
#else
#define FNAME sg_scan_sse41_128_8
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
    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m128i* const restrict pvP = parasail_memalign_m128i(16, n * segLen);
    __m128i* const restrict pvE = parasail_memalign_m128i(16, segLen);
    __m128i* const restrict pvHt= parasail_memalign_m128i(16, segLen);
    __m128i* const restrict pvFt= parasail_memalign_m128i(16, segLen);
    __m128i* const restrict pvH = parasail_memalign_m128i(16, segLen);
    __m128i vGapO = _mm_set1_epi8(open);
    __m128i vGapE = _mm_set1_epi8(gap);
    __m128i vSaturationCheck = _mm_setzero_si128();
    __m128i vNegLimit = _mm_set1_epi8(INT8_MIN);
    __m128i vPosLimit = _mm_set1_epi8(INT8_MAX);
    int8_t score = NEG_INF_8;
#if PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif

    /* Generate query profile.
     * Rearrange query sequence & calculate the weight of match/mismatch.
     * Don't alias. */
    {
        int32_t index = 0;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                __m128i_8_t t;
                j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = matrix[k][MAP_BLOSUM_[(unsigned char)s1[j]]];
                    j += segLen;
                }
                _mm_store_si128(&pvP[index], t.m);
                ++index;
            }
        }
    }

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m128i_8_t h;
            __m128i_8_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = 0;
                e.v[segNum] = NEG_INF_8;
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
        __m128i vHp;
        __m128i *pvW;
        __m128i vW;

#define LOOP_FUSION 1
#if LOOP_FUSION
        /* calculate E */
        /* calculate Ht */
        vHp = _mm_slli_si128(_mm_load_si128(pvH+(segLen-1)), 1);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            vE = _mm_max_epi8(
                    _mm_subs_epi8(vE, vGapE),
                    _mm_subs_epi8(vH, vGapO));
            vHt = _mm_max_epi8(
                    _mm_adds_epi8(vHp, vW),
                    vE);
            _mm_store_si128(pvE+i, vE);
            _mm_store_si128(pvHt+i, vHt);
            vHp = vH;
        }
#else
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
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            vHt = _mm_max_epi8(
                    _mm_adds_epi8(vH, vW),
                    vE);
            vH = _mm_load_si128(pvH+i);
            _mm_store_si128(pvHt+i, vHt);
        }
#endif

        /* calculate Ft */
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 1);
        vFt = _mm_set1_epi8(NEG_INF_8);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_max_epi8(
                    _mm_subs_epi8(vFt, vGapE),
                    vHt);
            vHt = _mm_load_si128(pvHt+i);
        }
        {
            __m128i_8_t tmp;
            tmp.m = vFt;
                tmp.v[1]  = MAX(tmp.v[0] -segLen*gap, tmp.v[1]);
                tmp.v[2]  = MAX(tmp.v[1] -segLen*gap, tmp.v[2]);
                tmp.v[3]  = MAX(tmp.v[2] -segLen*gap, tmp.v[3]);
                tmp.v[4]  = MAX(tmp.v[3] -segLen*gap, tmp.v[4]);
                tmp.v[5]  = MAX(tmp.v[4] -segLen*gap, tmp.v[5]);
                tmp.v[6]  = MAX(tmp.v[5] -segLen*gap, tmp.v[6]);
                tmp.v[7]  = MAX(tmp.v[6] -segLen*gap, tmp.v[7]);
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
            vFt = _mm_max_epi8(
                    _mm_subs_epi8(vFt, vGapE),
                    vHt);
            vHt = _mm_load_si128(pvHt+i);
            _mm_store_si128(pvFt+i, vFt);
        }

        /* calculate H */
        for (i=0; i<segLen; ++i) {
            vHt = _mm_load_si128(pvHt+i);
            vFt = _mm_load_si128(pvFt+i);
            vH = _mm_max_epi8(
                    vHt,
                    _mm_subs_epi8(vFt, vGapO));
            _mm_store_si128(pvH+i, vH);
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

        /* extract last value from column */
        {
            int8_t value;
            vH = _mm_load_si128(pvH + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128(vH, 1);
            }
            value = (int8_t) _mm_extract_epi8(vH, 15);
            if (value > score) {
                score = value;
            }
        }
    }

    /* max of last column */
    {
        __m128i vNegInf = _mm_set1_epi8(NEG_INF_8);
        __m128i vMaxH = vNegInf;
        __m128i vQIndexHi16 = _mm_set_epi16(
                15*segLen,
                14*segLen,
                13*segLen,
                12*segLen,
                11*segLen,
                10*segLen,
                9*segLen,
                8*segLen);
        __m128i vQIndexLo16 = _mm_set_epi16(
                7*segLen,
                6*segLen,
                5*segLen,
                4*segLen,
                3*segLen,
                2*segLen,
                1*segLen,
                0*segLen);
        __m128i vQLimit16 = _mm_set1_epi16(s1Len);
        __m128i vOne16 = _mm_set1_epi16(1);

        for (i=0; i<segLen; ++i) {
            /* load the last stored values */
            __m128i vH = _mm_load_si128(pvH + i);
            __m128i cond_lmt = _mm_packs_epi16(
                    _mm_cmplt_epi16(vQIndexLo16, vQLimit16),
                    _mm_cmplt_epi16(vQIndexHi16, vQLimit16));
            __m128i cond_max = _mm_cmpgt_epi8(vH, vMaxH);
            __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
            vMaxH = _mm_blendv_epi8(vMaxH, vH, cond_all);
            vQIndexLo16 = _mm_adds_epi16(vQIndexLo16, vOne16);
            vQIndexHi16 = _mm_adds_epi16(vQIndexHi16, vOne16);
        }

        /* max in vec */
        for (j=0; j<segWidth; ++j) {
            int8_t value = (int8_t) _mm_extract_epi8(vMaxH, 15);
            if (value > score) {
                score = value;
            }
            vMaxH = _mm_slli_si128(vMaxH, 1);
        }
    }

    if (_mm_movemask_epi8(vSaturationCheck)) {
        score = INT8_MAX;
    }

    result->score = score;

    free(pvH);
    free(pvFt);
    free(pvHt);
    free(pvE);
    free(pvP);

    return result;
}
