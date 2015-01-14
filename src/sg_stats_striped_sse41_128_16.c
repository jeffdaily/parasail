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

#define NEG_INF_16 (INT16_MIN/(int16_t)(2))

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
#define FNAME sg_stats_table_striped_sse41_128_16
#else
#define FNAME sg_stats_striped_sse41_128_16
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t nt = 0;
    int32_t segNum = 0;
    const int32_t n = 24; /* number of amino acids in table */
    const int32_t segWidth = 8; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
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
    __m128i* const restrict pvEM      = parasail_memalign_m128i(16, segLen);
    __m128i* const restrict pvEL      = parasail_memalign_m128i(16, segLen);
    __m128i vGapO = _mm_set1_epi16(open);
    __m128i vGapE = _mm_set1_epi16(gap);
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi16(1);
    int16_t score = NEG_INF_16;
    int16_t matches = NEG_INF_16;
    int16_t length = NEG_INF_16;
    __m128i initialF = _mm_set_epi16(
            -open-7*segLen*gap,
            -open-6*segLen*gap,
            -open-5*segLen*gap,
            -open-4*segLen*gap,
            -open-3*segLen*gap,
            -open-2*segLen*gap,
            -open-1*segLen*gap,
            -open-0*segLen*gap);
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
                __m128i_16_t t;
                __m128i_16_t s;
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
            __m128i_16_t h;
            __m128i_16_t e;
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
        vF = initialF;
        vFM = vZero;
        vFL = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = _mm_slli_si128(pvHStore[segLen - 1], 2);
        vHM = _mm_slli_si128(pvHMStore[segLen - 1], 2);
        vHL = _mm_slli_si128(pvHLStore[segLen - 1], 2);

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

            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            case1not = _mm_or_si128(
                    _mm_cmplt_epi16(vH,vF),_mm_cmplt_epi16(vH,vE));
            case2not = _mm_cmplt_epi16(vF,vE);
            case2 = _mm_andnot_si128(case2not,case1not);
            case3 = _mm_and_si128(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);

            /* calculate vM */
            vEM = _mm_load_si128(pvEM + i);
            vHM = _mm_andnot_si128(case1not,
                    _mm_add_epi16(vHM, _mm_load_si128(vPS + i)));
            vHM = _mm_or_si128(vHM, _mm_and_si128(case2, vFM));
            vHM = _mm_or_si128(vHM, _mm_and_si128(case3, vEM));
            _mm_store_si128(pvHMStore + i, vHM);

            /* calculate vL */
            vEL = _mm_load_si128(pvEL + i);
            vHL = _mm_andnot_si128(case1not, _mm_add_epi16(vHL, vOne));
            vHL = _mm_or_si128(vHL, _mm_and_si128(case2,
                        _mm_add_epi16(vFL, vOne)));
            vHL = _mm_or_si128(vHL, _mm_and_si128(case3,
                        _mm_add_epi16(vEL, vOne)));
            _mm_store_si128(pvHLStore + i, vHL);

#ifdef PARASAIL_TABLE
            arr_store_si128(result->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si128(result->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = _mm_subs_epi16(vH, vGapO);
            vE = _mm_subs_epi16(vE, vGapE);
            vE = _mm_max_epi16(vE, vH);
            _mm_store_si128(pvEStore + i, vE);
            _mm_store_si128(pvEM + i, vHM);
            _mm_store_si128(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm_subs_epi16(vF, vGapE);
            vF = _mm_max_epi16(vF, vH);
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
            __m128i vHp = _mm_slli_si128(pvHLoad[segLen - 1], 2);
            vF = _mm_slli_si128(vF, 2);
            vF = _mm_insert_epi16(vF, -open, 0);
            vFM = _mm_slli_si128(vFM, 2);
            vFL = _mm_slli_si128(vFL, 2);
            for (i=0; i<segLen; ++i) {
                __m128i case1not;
                __m128i case2not;
                __m128i case2;

                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm_adds_epi16(vHp, _mm_load_si128(vP + i));
                vE = _mm_load_si128(pvELoad + i);
                case1not = _mm_or_si128(
                        _mm_cmplt_epi16(vHp,vF),_mm_cmplt_epi16(vHp,vE));
                case2not = _mm_cmplt_epi16(vF,vE);
                case2 = _mm_andnot_si128(case2not,case1not);

                vHM = _mm_load_si128(pvHMStore + i);
                vHM = _mm_andnot_si128(case2, vHM);
                vHM = _mm_or_si128(vHM, _mm_and_si128(case2, vFM));
                _mm_store_si128(pvHMStore + i, vHM);
                _mm_store_si128(pvEM + i, vHM);

                vHL = _mm_load_si128(pvHLStore + i);
                vHL = _mm_andnot_si128(case2, vHL);
                vHL = _mm_or_si128(vHL, _mm_and_si128(case2,
                            _mm_add_epi16(vFL,vOne)));
                _mm_store_si128(pvHLStore + i, vHL);
                _mm_store_si128(pvEL + i, vHL);

                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi16(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
                arr_store_si128(result->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si128(result->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm_subs_epi16(vH, vGapO);
                vF = _mm_subs_epi16(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) goto end;
                vF = _mm_max_epi16(vF, vH);
                vFM = vHM;
                vFL = vHL;
                vHp = _mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
            int16_t tmp;
            /* extract last value from the column */
            vH = _mm_load_si128(pvHStore + offset);
            vHM = _mm_load_si128(pvHMStore + offset);
            vHL = _mm_load_si128(pvHLStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128 (vH, 2);
                vHM = _mm_slli_si128 (vHM, 2);
                vHL = _mm_slli_si128 (vHL, 2);
            }
            /* max of last value in each column */
            tmp = (int16_t) _mm_extract_epi16 (vH, 7);
            if (tmp > score) {
                score = tmp;
                matches = (int16_t)_mm_extract_epi16(vHM, 7);
                length = (int16_t)_mm_extract_epi16(vHL, 7);
            }
        }
    }

    /* max of last column */
    {
        __m128i vNegInf = _mm_set1_epi16(NEG_INF_16);
        __m128i vMaxLastColH = vNegInf;
        __m128i vMaxLastColHM = vNegInf;
        __m128i vMaxLastColHL = vNegInf;
        __m128i vQIndex = _mm_set_epi16(
                7*segLen,
                6*segLen,
                5*segLen,
                4*segLen,
                3*segLen,
                2*segLen,
                1*segLen,
                0*segLen);
        __m128i vQLimit = _mm_set1_epi16(s1Len);

        for (i=0; i<segLen; ++i) {
            /* load the last stored values */
            __m128i vH = _mm_load_si128(pvHStore + i);
            __m128i vHM = _mm_load_si128(pvHMStore + i);
            __m128i vHL = _mm_load_si128(pvHLStore + i);
            /* mask off the values that were padded */
            __m128i cond_lmt = _mm_cmplt_epi16(vQIndex, vQLimit);
            __m128i cond_max = _mm_cmpgt_epi16(vH, vMaxLastColH);
            __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
            vMaxLastColH = _mm_blendv_epi8(vMaxLastColH, vH, cond_all);
            vMaxLastColHM = _mm_blendv_epi8(vMaxLastColHM, vHM, cond_all);
            vMaxLastColHL = _mm_blendv_epi8(vMaxLastColHL, vHL, cond_all);
            vQIndex = _mm_adds_epi16(vQIndex, vOne);
        }

        /* max in vec */
        for (j=0; j<segWidth; ++j) {
            int16_t value = (int16_t) _mm_extract_epi16(vMaxLastColH, 7);
            if (value > score) {
                score = value;
                matches = (int16_t)_mm_extract_epi16(vMaxLastColHM, 7);
                length = (int16_t)_mm_extract_epi16(vMaxLastColHL, 7);
            }
            vMaxLastColH = _mm_slli_si128(vMaxLastColH, 2);
            vMaxLastColHM = _mm_slli_si128(vMaxLastColHM, 2);
            vMaxLastColHL = _mm_slli_si128(vMaxLastColHL, 2);
        }
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

