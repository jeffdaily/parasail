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
#include "parasail/memory.h"
#include "parasail/internal_avx.h"
#include "parasail/matrices/blosum_map.h"

#define NEG_INF (INT16_MIN/(int16_t)(2))

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)


#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen,
        int32_t bias)
{
    array[( 0*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  0) - bias;
    array[( 1*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  1) - bias;
    array[( 2*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  2) - bias;
    array[( 3*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  3) - bias;
    array[( 4*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  4) - bias;
    array[( 5*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  5) - bias;
    array[( 6*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  6) - bias;
    array[( 7*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  7) - bias;
    array[( 8*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  8) - bias;
    array[( 9*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  9) - bias;
    array[(10*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 10) - bias;
    array[(11*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 11) - bias;
    array[(12*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 12) - bias;
    array[(13*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 13) - bias;
    array[(14*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 14) - bias;
    array[(15*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 15) - bias;
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME sw_table_striped_avx2_256_16
#else
#define FNAME sw_striped_avx2_256_16
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
    __m256i* const restrict vProfile = parasail_memalign___m256i(32, n * segLen);
    __m256i* restrict pvHStore = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLoad =  parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvE = parasail_memalign___m256i(32, segLen);
    __m256i vGapO = _mm256_set1_epi16(open);
    __m256i vGapE = _mm256_set1_epi16(gap);
    int16_t bias = INT16_MIN;
    int16_t score = bias;
    __m256i vBias = _mm256_set1_epi16(bias);
    __m256i vMaxH = vBias;
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
                __m256i_16_t t;
                j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = j >= s1Len ? 0 : matrix[k][parasail_blosum_map[(unsigned char)s1[j]]];
                    j += segLen;
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
            __m256i_16_t h;
            __m256i_16_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = bias;
                e.v[segNum] = bias;
            }
            _mm256_store_si256(&pvHStore[index], h.m);
            _mm256_store_si256(&pvE[index], e.m);
            ++index;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m256i vE;
        /* Initialize F value to 0.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        __m256i vF = vBias;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m256i vH = _mm256_slli_si256_rpl(pvHStore[segLen - 1], 2);
        vH = _mm256_insert_epi16(vH, bias, 0);

        /* Correct part of the vProfile */
        const __m256i* vP = vProfile + parasail_blosum_map[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m256i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vH = _mm256_adds_epi16(vH, _mm256_load_si256(vP + i));
            vE = _mm256_load_si256(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = _mm256_max_epi16(vH, vE);
            vH = _mm256_max_epi16(vH, vF);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len, bias);
#endif
            vMaxH = _mm256_max_epi16(vH, vMaxH);

            /* Update vE value. */
            vH = _mm256_subs_epi16(vH, vGapO);
            vE = _mm256_subs_epi16(vE, vGapE);
            vE = _mm256_max_epi16(vE, vH);
            _mm256_store_si256(pvE + i, vE);

            /* Update vF value. */
            vF = _mm256_subs_epi16(vF, vGapE);
            vF = _mm256_max_epi16(vF, vH);

            /* Load the next vH. */
            vH = _mm256_load_si256(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            vF = _mm256_slli_si256_rpl(vF, 2);
            vF = _mm256_insert_epi16(vF, bias, 0);
            for (i=0; i<segLen; ++i) {
                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi16(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
                arr_store_si256(result->score_table, vH, i, segLen, j, s2Len, bias);
#endif
                vMaxH = _mm256_max_epi16(vH, vMaxH);
                vH = _mm256_subs_epi16(vH, vGapO);
                vF = _mm256_subs_epi16(vF, vGapE);
                if (! _mm256_movemask_epi8(_mm256_cmpgt_epi16(vF, vH))) goto end;
                /*vF = _mm256_max_epi16(vF, vH);*/
            }
        }
end:
        {
        }
    }

    /* max in vec */
    for (j=0; j<segWidth; ++j) {
        int16_t value = (int16_t) _mm256_extract_epi16(vMaxH, 15);
        if (value > score) {
            score = value;
        }
        vMaxH = _mm256_slli_si256_rpl(vMaxH, 2);
    }

    if (score == INT16_MAX) {
        result->saturated = 1;
        score = INT16_MAX;
    }

    result->score = score - bias;

    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);
    parasail_free(vProfile);

    return result;
}

