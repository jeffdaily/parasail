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

#define NEG_INF_32 (INT32_MIN/(int32_t)(2))

#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        __m128i vH,
        int32_t i,
        int32_t segWidth,
        int32_t j,
        int32_t s2Len)
{
    array[(i*segWidth+0)*s2Len + j] = (int32_t)_mm_extract_epi32(vH, 0);
    array[(i*segWidth+1)*s2Len + j] = (int32_t)_mm_extract_epi32(vH, 1);
    array[(i*segWidth+2)*s2Len + j] = (int32_t)_mm_extract_epi32(vH, 2);
    array[(i*segWidth+3)*s2Len + j] = (int32_t)_mm_extract_epi32(vH, 3);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME sw_table_blocked_sse41_128_32
#else
#define FNAME sw_blocked_sse41_128_32
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
    __m128i* const restrict vProfile = parasail_memalign_m128i(16, n * segLen);
    __m128i* restrict pvH = parasail_memalign_m128i(16, segLen);
    __m128i* const restrict pvE = parasail_memalign_m128i(16, segLen);
    int score = NEG_INF_32;
    __m128i vGapO = _mm_set1_epi32(open);
    __m128i vGapE = _mm_set1_epi32(gap);
    __m128i vZero = _mm_setzero_si128();
    __m128i vNegInf = _mm_set1_epi32(NEG_INF_32);
    __m128i vMaxH = vZero;
#ifdef PARASAIL_TABLE
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
                __m128i_32_t t;
                j = i*segWidth;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = j >= s1Len ? 0 : matrix[k][MAP_BLOSUM_[(unsigned char)s1[j]]];
                    j += 1;
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
                e.v[segNum] = -open;
            }
            _mm_store_si128(&pvH[index], h.m);
            _mm_store_si128(&pvE[index], e.m);
            ++index;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vX = vZero;
        __m128i vF = vNegInf;
        const __m128i* pvP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;
        for (i=0; i<segLen; ++i) {
            __m128i vP;
            __m128i vH;
            __m128i vE;
            __m128i vT1;

            vH = _mm_load_si128(pvH + i);
            vE = _mm_load_si128(pvE + i);

            vT1 = _mm_srli_si128(vH, 12); /* rshift 3 */
            vH = _mm_slli_si128(vH, 4); /* lshift 1 */
            vH = _mm_or_si128(vH, vX);
            vX = vT1;

            vP = _mm_load_si128(pvP + i);
            vH = _mm_add_epi32(vH, vP);
            vH = _mm_max_epi32(vH, vE);
            vH = _mm_max_epi32(vH, vZero);

            vF = _mm_srli_si128(vF, 12);
            vF = _mm_or_si128(vF, _mm_slli_si128(vH, 4));
            vF = _mm_sub_epi32(vF, vGapO);
            if (_mm_movemask_epi8(_mm_cmpgt_epi32(vF, vZero))) {
                __m128i vT2 = vF;
                while (_mm_movemask_epi8(_mm_cmpgt_epi32(vT2, vZero))) {
                    vT2 = _mm_slli_si128(vT2, 4);
                    vT2 = _mm_sub_epi32(vT2, vGapE);
                    vF = _mm_max_epi32(vF, vT2);
                }
                vH = _mm_max_epi32(vH, vF);
                vF = _mm_add_epi32(vF, vGapO);
                vF = _mm_sub_epi32(vF, vGapE);
                vF = _mm_max_epi32(vH, vF);
            }
            else {
                vF = vH;
            }

            _mm_store_si128(pvH + i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vH, i, segWidth, j, s2Len);
#endif
            vMaxH = _mm_max_epi32(vMaxH, vH);

            vH = _mm_sub_epi32(vH, vGapO);
            vE = _mm_sub_epi32(vE, vGapE);
            vE = _mm_max_epi32(vE, vH);
            _mm_store_si128(pvE + i, vE);

        }
    }

    /* max in vec */
    for (j=0; j<segWidth; ++j) {
        int32_t value = (int32_t) _mm_extract_epi32(vMaxH, 3);
        if (value > score) {
            score = value;
        }
        vMaxH = _mm_slli_si128(vMaxH, 4);
    }

    result->score = score;

    parasail_free(pvE);
    parasail_free(pvH);
    parasail_free(vProfile);

    return result;
}

