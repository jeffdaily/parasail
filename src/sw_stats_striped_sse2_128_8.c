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
    array[( 0*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  0);
    array[( 1*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  1);
    array[( 2*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  2);
    array[( 3*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  3);
    array[( 4*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  4);
    array[( 5*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  5);
    array[( 6*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  6);
    array[( 7*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  7);
    array[( 8*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  8);
    array[( 9*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  9);
    array[(10*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 10);
    array[(11*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 11);
    array[(12*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 12);
    array[(13*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 13);
    array[(14*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 14);
    array[(15*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 15);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME sw_stats_table_striped_sse2_128_8
#else
#define FNAME sw_stats_striped_sse2_128_8
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = 24; /* number of amino acids in table */
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m128i* const restrict vProfile  = parasail_memalign_m128i(16, n * segLen);
    __m128i* const restrict vProfileS = parasail_memalign_m128i(16, n * segLen);
    __m128i* restrict pvHStore        = parasail_memalign_m128i(16, segLen);
    __m128i* restrict pvHLoad         = parasail_memalign_m128i(16, segLen);
    __m128i* restrict pvHMStore       = parasail_memalign_m128i(16, segLen);
    __m128i* restrict pvHMLoad        = parasail_memalign_m128i(16, segLen);
    __m128i* restrict pvHLStore       = parasail_memalign_m128i(16, segLen);
    __m128i* restrict pvHLLoad        = parasail_memalign_m128i(16, segLen);
    __m128i* restrict pvEStore        = parasail_memalign_m128i(16, segLen);
    __m128i* restrict pvELoad         = parasail_memalign_m128i(16, segLen);
    __m128i* restrict pvEM            = parasail_memalign_m128i(16, segLen);
    __m128i* restrict pvEL            = parasail_memalign_m128i(16, segLen);
    __m128i vGapO = _mm_set1_epi8(open);
    __m128i vGapE = _mm_set1_epi8(gap);
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi8(1);
    int8_t score = NEG_INF_8;
    int8_t matches = NEG_INF_8;
    int8_t length = NEG_INF_8;
    /* for max calculation we don't want to include padded cells */
    __m128i vOne16 = _mm_set1_epi16(1);
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
    /* Trace the highest score of the whole SW matrix. */
    __m128i vMaxH = vZero;
    __m128i vMaxM = vZero;
    __m128i vMaxL = vZero;
    __m128i vSaturationCheck = _mm_setzero_si128();
    __m128i vNegLimit = _mm_set1_epi8(INT8_MIN);
    __m128i vPosLimit = _mm_set1_epi8(INT8_MAX);
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif

    parasail_memset_m128i(pvHMStore, vZero, segLen);
    parasail_memset_m128i(pvHLStore, vZero, segLen);

    /* Generate query profile.
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
                _mm_store_si128(&vProfile[index], t.m);
                _mm_store_si128(&vProfileS[index], s.m);
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
                e.v[segNum] = -open;
            }
            _mm_store_si128(&pvHStore[index], h.m);
            _mm_store_si128(&pvEStore[index], e.m);
            ++index;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vQIndexLo16 = vQIndexLo16_reset;
        __m128i vQIndexHi16 = vQIndexHi16_reset;
        __m128i vE;
        __m128i vEM;
        __m128i vEL;
        __m128i vF;
        __m128i vFM;
        __m128i vFL;
        __m128i vH;
        __m128i vHM;
        __m128i vHL;
        const __m128i* vP = NULL;
        const __m128i* vPS = NULL;
        __m128i* pv = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        vF = vZero;
        vFM = vZero;
        vFL = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = _mm_slli_si128(pvHStore[segLen - 1], 1);
        vHM = _mm_slli_si128(pvHMStore[segLen - 1], 1);
        vHL = _mm_slli_si128(pvHLStore[segLen - 1], 1);

        /* Correct part of the vProfile */
        vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;
        vPS = vProfileS + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        pv = pvHMLoad;
        pvHMLoad = pvHMStore;
        pvHMStore = pv;
        pv = pvHLLoad;
        pvHLLoad = pvHLStore;
        pvHLStore = pv;
        pv = pvELoad;
        pvELoad = pvEStore;
        pvEStore = pv;

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            __m128i case1not;
            __m128i case2not;
            __m128i case2;
            __m128i case3;
            __m128i cond_zero;
            vH = _mm_adds_epi8(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            case1not = _mm_or_si128(
                    _mm_cmplt_epi8(vH,vF),_mm_cmplt_epi8(vH,vE));
            case2not = _mm_cmplt_epi8(vF,vE);
            case2 = _mm_andnot_si128(case2not,case1not);
            case3 = _mm_and_si128(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi8(vH, vE);
            vH = _mm_max_epi8(vH, vF);
            vH = _mm_max_epi8(vH, vZero);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);
            cond_zero = _mm_cmpeq_epi8(vH, vZero);

            /* calculate vM */
            vEM = _mm_load_si128(pvEM + i);
            vHM = _mm_andnot_si128(case1not,
                    _mm_adds_epi8(vHM, _mm_load_si128(vPS + i)));
            vHM = _mm_or_si128(vHM, _mm_and_si128(case2, vFM));
            vHM = _mm_or_si128(vHM, _mm_and_si128(case3, vEM));
            vHM = _mm_andnot_si128(cond_zero, vHM);
            _mm_store_si128(pvHMStore + i, vHM);

            /* calculate vL */
            vEL = _mm_load_si128(pvEL + i);
            vHL = _mm_andnot_si128(case1not, _mm_adds_epi8(vHL, vOne));
            vHL = _mm_or_si128(vHL, _mm_and_si128(case2,
                        _mm_adds_epi8(vFL, vOne)));
            vHL = _mm_or_si128(vHL, _mm_and_si128(case3,
                        _mm_adds_epi8(vEL, vOne)));
            vHL = _mm_andnot_si128(cond_zero, vHL);
            _mm_store_si128(pvHLStore + i, vHL);

            /* check for saturation */
            {
                vSaturationCheck = _mm_or_si128(vSaturationCheck,
                        _mm_or_si128(
                            _mm_or_si128(
                                _mm_cmpeq_epi8(vH, vNegLimit),
                                _mm_cmpeq_epi8(vH, vPosLimit)),
                            _mm_or_si128(
                                _mm_cmpeq_epi8(vHM, vPosLimit),
                                _mm_cmpeq_epi8(vHL, vPosLimit))));
            }
#ifdef PARASAIL_TABLE
            arr_store_si128(result->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si128(result->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif
            /* update max vector seen so far */
            {
                __m128i cond_max = _mm_cmpgt_epi8(vH,vMaxH);
                __m128i cond_lmt = _mm_packs_epi16(
                        _mm_cmplt_epi16(vQIndexLo16, vQLimit16),
                        _mm_cmplt_epi16(vQIndexHi16, vQLimit16));
                __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
                vMaxH = _mm_andnot_si128(cond_all, vMaxH);
                vMaxH = _mm_or_si128(vMaxH, _mm_and_si128(cond_all, vH));
                vMaxM = _mm_andnot_si128(cond_all, vMaxM);
                vMaxM = _mm_or_si128(vMaxM, _mm_and_si128(cond_all, vHM));
                vMaxL = _mm_andnot_si128(cond_all, vMaxL);
                vMaxL = _mm_or_si128(vMaxL, _mm_and_si128(cond_all, vHL));
                vQIndexLo16 = _mm_adds_epi16(vQIndexLo16, vOne16);
                vQIndexHi16 = _mm_adds_epi16(vQIndexHi16, vOne16);
            }

            /* Update vE value. */
            vH = _mm_subs_epi8(vH, vGapO);
            vE = _mm_subs_epi8(vE, vGapE);
            vE = _mm_max_epi8(vE, vH);
            _mm_store_si128(pvEStore + i, vE);
            _mm_store_si128(pvEM + i, vHM);
            _mm_store_si128(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm_subs_epi8(vF, vGapE);
            vF = _mm_max_epi8(vF, vH);
            vFM = vHM;
            vFL = vHL;

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
            vHM = _mm_load_si128(pvHMLoad + i);
            vHL = _mm_load_si128(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            __m128i vHp = _mm_slli_si128(pvHLoad[segLen - 1], 1);
            vF = _mm_slli_si128(vF, 1);
            vFM = _mm_slli_si128(vFM, 1);
            vFL = _mm_slli_si128(vFL, 1);
            for (i=0; i<segLen; ++i) {
                __m128i case1not;
                __m128i case2not;
                __m128i case2;
                __m128i cond_zero;
                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm_adds_epi8(vHp, _mm_load_si128(vP + i));
                vE = _mm_load_si128(pvELoad + i);
                case1not = _mm_or_si128(
                        _mm_cmplt_epi8(vHp,vF),_mm_cmplt_epi8(vHp,vE));
                case2not = _mm_cmplt_epi8(vF,vE);
                case2 = _mm_andnot_si128(case2not,case1not);

                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi8(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
                cond_zero = _mm_cmpeq_epi8(vH, vZero);

                vHM = _mm_load_si128(pvHMStore + i);
                vHM = _mm_andnot_si128(case2, vHM);
                vHM = _mm_or_si128(vHM, _mm_and_si128(case2, vFM));
                vHM = _mm_andnot_si128(cond_zero, vHM);
                _mm_store_si128(pvHMStore + i, vHM);
                _mm_store_si128(pvEM + i, vHM);

                vHL = _mm_load_si128(pvHLStore + i);
                vHL = _mm_andnot_si128(case2, vHL);
                vHL = _mm_or_si128(vHL, _mm_and_si128(case2,
                            _mm_adds_epi8(vFL,vOne)));
                vHL = _mm_andnot_si128(cond_zero, vHL);
                _mm_store_si128(pvHLStore + i, vHL);
                _mm_store_si128(pvEL + i, vHL);

                /* check for saturation */
                {
                    vSaturationCheck = _mm_or_si128(vSaturationCheck,
                            _mm_or_si128(
                                _mm_or_si128(
                                    _mm_cmpeq_epi8(vH, vNegLimit),
                                    _mm_cmpeq_epi8(vH, vPosLimit)),
                                _mm_or_si128(
                                    _mm_cmpeq_epi8(vHM, vPosLimit),
                                    _mm_cmpeq_epi8(vHL, vPosLimit))));
                }
#ifdef PARASAIL_TABLE
                arr_store_si128(result->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si128(result->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm_subs_epi8(vH, vGapO);
                vF = _mm_subs_epi8(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi8(vF, vH))) goto end;
                vF = _mm_max_epi8(vF, vH);
                vFM = vHM;
                vFL = vHL;
                vHp = _mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
        }
    }

    /* max in vec */
    for (j=0; j<segWidth; ++j) {
        int8_t value = (int8_t) _mm_extract_epi8(vMaxH, 15);
        if (value > score) {
            score = value;
            matches = (int8_t) _mm_extract_epi8(vMaxM, 15);
            length = (int8_t) _mm_extract_epi8(vMaxL, 15);
        }
        vMaxH = _mm_slli_si128(vMaxH, 1);
        vMaxM = _mm_slli_si128(vMaxM, 1);
        vMaxL = _mm_slli_si128(vMaxL, 1);
    }

    if (_mm_movemask_epi8(vSaturationCheck)) {
        score = INT8_MAX;
        matches = 0;
        length = 0;
    }

    result->score = score;
    result->matches = matches;
    result->length = length;

    parasail_free(pvEL);
    parasail_free(pvEM);
    parasail_free(pvELoad);
    parasail_free(pvEStore);
    parasail_free(pvHLLoad);
    parasail_free(pvHLStore);
    parasail_free(pvHMLoad);
    parasail_free(pvHMStore);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);
    parasail_free(vProfileS);
    parasail_free(vProfile);

    return result;
}

