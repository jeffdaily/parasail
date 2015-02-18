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
#define FNAME sg_stats_table_striped_sse2_128_32
#else
#define FNAME sg_stats_striped_sse2_128_32
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
    const int32_t segWidth = 4; /* number of values in vector unit */
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
    __m128i vGapO = _mm_set1_epi32(open);
    __m128i vGapE = _mm_set1_epi32(gap);
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi32(1);
    int32_t score = NEG_INF_32;
    int32_t matches = NEG_INF_32;
    int32_t length = NEG_INF_32;
    __m128i vNegInf = _mm_set1_epi32(NEG_INF_32);
    __m128i vMaxH = vNegInf;
    __m128i vMaxHM = vNegInf;
    __m128i vMaxHL = vNegInf;
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
                __m128i_32_t t;
                __m128i_32_t s;
                j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = j >= s1Len ? 0 : matrix[k][MAP_BLOSUM_[(unsigned char)s1[j]]];
                    s.v[segNum] = j >= s1Len ? 0 : (k == MAP_BLOSUM_[(unsigned char)s1[j]]);
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
            __m128i_32_t h;
            __m128i_32_t e;
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

        /* Initialize F value to neg inf.  Any errors to vH values will
         * be corrected in the Lazy_F loop.  */
        vF = vNegInf;
        vFM = vZero;
        vFL = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = _mm_slli_si128(pvHStore[segLen - 1], 4);
        vHM = _mm_slli_si128(pvHMStore[segLen - 1], 4);
        vHL = _mm_slli_si128(pvHLStore[segLen - 1], 4);

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

            vH = _mm_add_epi32(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            case1not = _mm_or_si128(
                    _mm_cmplt_epi32(vH,vF),_mm_cmplt_epi32(vH,vE));
            case2not = _mm_cmplt_epi32(vF,vE);
            case2 = _mm_andnot_si128(case2not,case1not);
            case3 = _mm_and_si128(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi32(vH, vE);
            vH = _mm_max_epi32(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);

            /* calculate vM */
            vEM = _mm_load_si128(pvEM + i);
            vHM = _mm_andnot_si128(case1not,
                    _mm_add_epi32(vHM, _mm_load_si128(vPS + i)));
            vHM = _mm_or_si128(vHM, _mm_and_si128(case2, vFM));
            vHM = _mm_or_si128(vHM, _mm_and_si128(case3, vEM));
            _mm_store_si128(pvHMStore + i, vHM);

            /* calculate vL */
            vEL = _mm_load_si128(pvEL + i);
            vHL = _mm_andnot_si128(case1not, _mm_add_epi32(vHL, vOne));
            vHL = _mm_or_si128(vHL, _mm_and_si128(case2,
                        _mm_add_epi32(vFL, vOne)));
            vHL = _mm_or_si128(vHL, _mm_and_si128(case3,
                        _mm_add_epi32(vEL, vOne)));
            _mm_store_si128(pvHLStore + i, vHL);

#ifdef PARASAIL_TABLE
            arr_store_si128(result->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si128(result->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = _mm_sub_epi32(vH, vGapO);
            vE = _mm_sub_epi32(vE, vGapE);
            vE = _mm_max_epi32(vE, vH);
            _mm_store_si128(pvEStore + i, vE);
            _mm_store_si128(pvEM + i, vHM);
            _mm_store_si128(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm_sub_epi32(vF, vGapE);
            vF = _mm_max_epi32(vF, vH);
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
            __m128i vHp = _mm_slli_si128(pvHLoad[segLen - 1], 4);
            vF = _mm_slli_si128(vF, 4);
            vF = _mm_insert_epi32(vF, -open, 0);
            vFM = _mm_slli_si128(vFM, 4);
            vFL = _mm_slli_si128(vFL, 4);
            for (i=0; i<segLen; ++i) {
                __m128i case1not;
                __m128i case2not;
                __m128i case2;

                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm_add_epi32(vHp, _mm_load_si128(vP + i));
                vE = _mm_load_si128(pvELoad + i);
                case1not = _mm_or_si128(
                        _mm_cmplt_epi32(vHp,vF),_mm_cmplt_epi32(vHp,vE));
                case2not = _mm_cmplt_epi32(vF,vE);
                case2 = _mm_andnot_si128(case2not,case1not);

                vHM = _mm_load_si128(pvHMStore + i);
                vHM = _mm_andnot_si128(case2, vHM);
                vHM = _mm_or_si128(vHM, _mm_and_si128(case2, vFM));
                _mm_store_si128(pvHMStore + i, vHM);
                _mm_store_si128(pvEM + i, vHM);

                vHL = _mm_load_si128(pvHLStore + i);
                vHL = _mm_andnot_si128(case2, vHL);
                vHL = _mm_or_si128(vHL, _mm_and_si128(case2,
                            _mm_add_epi32(vFL,vOne)));
                _mm_store_si128(pvHLStore + i, vHL);
                _mm_store_si128(pvEL + i, vHL);

                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi32(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
                arr_store_si128(result->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si128(result->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm_sub_epi32(vH, vGapO);
                vF = _mm_sub_epi32(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi32(vF, vH))) goto end;
                /*vF = _mm_max_epi32(vF, vH);*/
                vFM = vHM;
                vFL = vHL;
                vHp = _mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
            /* extract vector containing last value from the column */
            __m128i cond_max;
            vH = _mm_load_si128(pvHStore + offset);
            vHM = _mm_load_si128(pvHMStore + offset);
            vHL = _mm_load_si128(pvHLStore + offset);
            cond_max = _mm_cmpgt_epi32(vH, vMaxH);
            vMaxH = _mm_andnot_si128(cond_max, vMaxH);
            vMaxH = _mm_or_si128(vMaxH,
                    _mm_and_si128(cond_max, vH));
            vMaxHM = _mm_andnot_si128(cond_max, vMaxHM);
            vMaxHM = _mm_or_si128(vMaxHM,
                    _mm_and_si128(cond_max, vHM));
            vMaxHL = _mm_andnot_si128(cond_max, vMaxHL);
            vMaxHL = _mm_or_si128(vMaxHL,
                    _mm_and_si128(cond_max, vHL));
        }
    }

    /* extract last value from the column */
    {
        int32_t tmp;
        for (k=0; k<position; ++k) {
            vMaxH  = _mm_slli_si128(vMaxH, 4);
            vMaxHM = _mm_slli_si128(vMaxHM, 4);
            vMaxHL = _mm_slli_si128(vMaxHL, 4);
        }
        tmp = (int32_t) _mm_extract_epi32 (vMaxH, 3);
        if (tmp > score) {
            score = tmp;
            matches = (int32_t)_mm_extract_epi32(vMaxHM, 3);
            length = (int32_t)_mm_extract_epi32(vMaxHL, 3);
        }
    }

    /* max of last column */
    {
        vMaxH = vNegInf;
        vMaxHM = vNegInf;
        vMaxHL = vNegInf;

        for (i=0; i<segLen; ++i) {
            /* load the last stored values */
            __m128i vH = _mm_load_si128(pvHStore + i);
            __m128i vHM = _mm_load_si128(pvHMStore + i);
            __m128i vHL = _mm_load_si128(pvHLStore + i);
            __m128i cond_max = _mm_cmpgt_epi32(vH, vMaxH);
            vMaxH = _mm_andnot_si128(cond_max, vMaxH);
            vMaxH = _mm_or_si128(vMaxH,
                    _mm_and_si128(cond_max, vH));
            vMaxHM = _mm_andnot_si128(cond_max, vMaxHM);
            vMaxHM = _mm_or_si128(vMaxHM,
                    _mm_and_si128(cond_max, vHM));
            vMaxHL = _mm_andnot_si128(cond_max, vMaxHL);
            vMaxHL = _mm_or_si128(vMaxHL,
                    _mm_and_si128(cond_max, vHL));
        }

        /* max in vec */
        for (j=0; j<segWidth; ++j) {
            int32_t value = (int32_t) _mm_extract_epi32(vMaxH, 3);
            if (value > score) {
                score = value;
                matches = (int32_t)_mm_extract_epi32(vMaxHM, 3);
                length = (int32_t)_mm_extract_epi32(vMaxHL, 3);
            }
            vMaxH = _mm_slli_si128(vMaxH, 4);
            vMaxHM = _mm_slli_si128(vMaxHM, 4);
            vMaxHL = _mm_slli_si128(vMaxHL, 4);
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

