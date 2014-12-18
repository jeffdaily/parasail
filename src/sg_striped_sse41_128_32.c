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
#include "blosum/blosum_map.h"

#define NEG_INF_32 (INT32_MIN/(int32_t)(2))

#if PARASAIL_TABLE
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
#define FNAME sg_table_striped_sse41_128_32
#else
#define FNAME sg_striped_sse41_128_32
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
    int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m128i* vProfile = parasail_memalign_m128i(16, n * segLen);
    __m128i* pvHStore = parasail_memalign_m128i(16, segLen);
    __m128i* pvHLoad =  parasail_memalign_m128i(16, segLen);
    __m128i* pvE = parasail_memalign_m128i(16, segLen);
    int score = NEG_INF_32;
    __m128i vGapO = _mm_set1_epi32(open);
    __m128i vGapE = _mm_set1_epi32(gap);
    __m128i vOne = _mm_set1_epi32(1);
    __m128i initialF = _mm_set_epi32(
            -open-open-3*segLen*gap,
            -open-open-2*segLen*gap,
            -open-open-1*segLen*gap,
            -open-open-0*segLen*gap);
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
                int32_t j = i;
                __m128i_32_t t;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = matrix[k][MAP_BLOSUM_[(unsigned char)s1[j]]];
                    j += segLen;
                }
                _mm_store_si128(&vProfile[index], t.m);
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
                e.v[segNum] = NEG_INF_32;
            }
            _mm_store_si128(&pvHStore[index], h.m);
            _mm_store_si128(&pvE[index], e.m);
            ++index;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vE;
        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        __m128i vF = initialF;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 4);

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vH = _mm_add_epi32(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi32(vH, vE);
            vH = _mm_max_epi32(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = _mm_sub_epi32(vH, vGapO);
            vE = _mm_sub_epi32(vE, vGapE);
            vE = _mm_max_epi32(vE, vH);
            _mm_store_si128(pvE + i, vE);

            /* Update vF value. */
            vF = _mm_sub_epi32(vF, vGapE);
            vF = _mm_max_epi32(vF, vH);

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<8; ++k) {
            vF = _mm_slli_si128(vF, 4);
            vF = _mm_insert_epi32(vF, -open, 0);
            for (i=0; i<segLen; ++i) {
                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi32(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
                arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm_sub_epi32(vH, vGapO);
                vF = _mm_sub_epi32(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi32(vF, vH))) goto end;
                vF = _mm_max_epi32(vF, vH);
            }
        }
end:
        {
            /* extract last value from the column */
            int32_t tmp;
            vH = _mm_load_si128(pvHStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128 (vH, 4);
            }
            /* max of last value in each column */
            tmp = (int32_t) _mm_extract_epi32 (vH, 3);
            if (tmp > score) {
                score = tmp;
            }
        }
    }

    /* max of last column */
    {
        __m128i vNegInf = _mm_set1_epi32(NEG_INF_32);
        __m128i vMaxLastColH = vNegInf;
        __m128i vQIndex = _mm_set_epi32(
                3*segLen,
                2*segLen,
                1*segLen,
                0*segLen);
        __m128i vQLimit = _mm_set1_epi32(s1Len);

        for (i=0; i<segLen; ++i) {
            __m128i vH = _mm_load_si128(pvHStore + i);
            __m128i cond_lmt = _mm_cmplt_epi32(vQIndex, vQLimit);
            __m128i cond_max = _mm_cmpgt_epi32(vH, vMaxLastColH);
            __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
            vMaxLastColH = _mm_blendv_epi8(vMaxLastColH, vH, cond_all);
            vQIndex = _mm_add_epi32(vQIndex, vOne);
        }

        /* max in vec */
        for (j=0; j<8; ++j) {
            int32_t value = (int32_t) _mm_extract_epi32(vMaxLastColH, 3);
            if (value > score) {
                score = value;
            }
            vMaxLastColH = _mm_slli_si128(vMaxLastColH, 4);
        }
    }

    result->score = score;

    free(pvE);
    free(pvHLoad);
    free(pvHStore);
    free(vProfile);

    return result;
}

