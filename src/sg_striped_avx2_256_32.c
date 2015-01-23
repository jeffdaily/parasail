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

#include <immintrin.h>

#include "parasail.h"
#include "parasail_internal.h"
#include "parasail_internal_avx.h"
#include "blosum/blosum_map.h"

#define NEG_INF_32 (INT32_MIN/(int32_t)(2))

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

/* avx2 _mm256_slli_si256 does not shift across 128-bit lanes, emulate it */
static inline __m256i shift(__m256i a) {
    return _mm256_alignr_epi8(a,
            _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)),
            12);
}

#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[(0*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 0);
    array[(1*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 1);
    array[(2*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 2);
    array[(3*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 3);
    array[(4*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 4);
    array[(5*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 5);
    array[(6*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 6);
    array[(7*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 7);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME sg_table_striped_avx2_256_32
#else
#define FNAME sg_striped_avx2_256_32
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
    const int32_t segWidth = 8; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m256i* const restrict vProfile = parasail_memalign_m256i(32, n * segLen);
    __m256i* restrict pvHStore = parasail_memalign_m256i(32, segLen);
    __m256i* restrict pvHLoad =  parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvE = parasail_memalign_m256i(32, segLen);
    int score = NEG_INF_32;
    __m256i vGapO = _mm256_set1_epi32(open);
    __m256i vGapE = _mm256_set1_epi32(gap);
    __m256i vOne = _mm256_set1_epi32(1);
    __m256i initialF = _mm256_set_epi32(
            -open-open-7*segLen*gap,
            -open-open-6*segLen*gap,
            -open-open-5*segLen*gap,
            -open-open-4*segLen*gap,
            -open-open-3*segLen*gap,
            -open-open-2*segLen*gap,
            -open-open-1*segLen*gap,
            -open-open-0*segLen*gap);
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
                __m256i_32_t t;
                j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = matrix[k][MAP_BLOSUM_[(unsigned char)s1[j]]];
                    j += segLen;
                    if (j >= s1Len) break;
                }
                _mm256_store_si256(&vProfile[index], t.m);
                ++index;
            }
        }
    }

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m256i_32_t h;
            __m256i_32_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = 0;
                e.v[segNum] = NEG_INF_32;
            }
            _mm256_store_si256(&pvHStore[index], h.m);
            _mm256_store_si256(&pvE[index], e.m);
            ++index;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m256i vE;
        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        __m256i vF = initialF;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m256i vH = shift(pvHStore[segLen - 1]);

        /* Correct part of the vProfile */
        const __m256i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m256i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vH = _mm256_add_epi32(vH, _mm256_load_si256(vP + i));
            vE = _mm256_load_si256(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = _mm256_max_epi32(vH, vE);
            vH = _mm256_max_epi32(vH, vF);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = _mm256_sub_epi32(vH, vGapO);
            vE = _mm256_sub_epi32(vE, vGapE);
            vE = _mm256_max_epi32(vE, vH);
            _mm256_store_si256(pvE + i, vE);

            /* Update vF value. */
            vF = _mm256_sub_epi32(vF, vGapE);
            vF = _mm256_max_epi32(vF, vH);

            /* Load the next vH. */
            vH = _mm256_load_si256(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            vF = shift(vF);
            vF = _mm256_insert_epi32(vF, -open, 0);
            for (i=0; i<segLen; ++i) {
                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi32(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
                arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm256_sub_epi32(vH, vGapO);
                vF = _mm256_sub_epi32(vF, vGapE);
                if (! _mm256_movemask_epi8(_mm256_cmpgt_epi32(vF, vH))) goto end;
                vF = _mm256_max_epi32(vF, vH);
            }
        }
end:
        {
            /* extract last value from the column */
            int32_t tmp;
            vH = _mm256_load_si256(pvHStore + offset);
            for (k=0; k<position; ++k) {
                vH = shift(vH);
            }
            /* max of last value in each column */
            tmp = (int32_t) _mm256_extract_epi32 (vH, 7);
            if (tmp > score) {
                score = tmp;
            }
        }
    }

    /* max of last column */
    {
        __m256i vNegInf = _mm256_set1_epi32(NEG_INF_32);
        __m256i vMaxLastColH = vNegInf;
        __m256i vQIndex = _mm256_set_epi32(
                7*segLen,
                6*segLen,
                5*segLen,
                4*segLen,
                3*segLen,
                2*segLen,
                1*segLen,
                0*segLen);
        __m256i vQLimit1 = _mm256_set1_epi32(s1Len-1);

        for (i=0; i<segLen; ++i) {
            __m256i vH = _mm256_load_si256(pvHStore + i);
            __m256i cond_lmt = _mm256_cmpgt_epi32(vQIndex, vQLimit1);
            __m256i cond_max = _mm256_cmpgt_epi32(vH, vMaxLastColH);
            __m256i cond_all = _mm256_andnot_si256(cond_lmt, cond_max);
            vMaxLastColH = _mm256_blendv_epi8(vMaxLastColH, vH, cond_all);
            vQIndex = _mm256_add_epi32(vQIndex, vOne);
        }

        /* max in vec */
        for (j=0; j<segWidth; ++j) {
            int32_t value = (int32_t) _mm256_extract_epi32(vMaxLastColH, 7);
            if (value > score) {
                score = value;
            }
            vMaxLastColH = shift(vMaxLastColH);
        }
    }

    result->score = score;

    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);
    parasail_free(vProfile);

    return result;
}

