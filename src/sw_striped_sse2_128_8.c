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
#include "blosum/blosum_map.h"

#define NEG_INF_8 (INT8_MIN)

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
#define FNAME sw_table_striped_sse2_128_8
#else
#define FNAME sw_striped_sse2_128_8
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
    int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m128i* const restrict vProfile = parasail_memalign_m128i(16, n * segLen);
    __m128i* restrict pvHStore = parasail_memalign_m128i(16, segLen);
    __m128i* restrict pvHLoad =  parasail_memalign_m128i(16, segLen);
    __m128i* const restrict pvE = parasail_memalign_m128i(16, segLen);
    int score = NEG_INF_8;
    __m128i vGapO = _mm_set1_epi8(open);
    __m128i vGapE = _mm_set1_epi8(gap);
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne16 = _mm_set1_epi16(1);
    /* for max calculation we don't want to include padded cells */
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
    __m128i vQLimit16 = _mm_set1_epi16(s1Len);
    __m128i vMaxH = vZero;
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
                __m128i_8_t t;
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
            __m128i_8_t h;
            __m128i_8_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = 0;
                e.v[segNum] = NEG_INF_8;
            }
            _mm_store_si128(&pvHStore[index], h.m);
            _mm_store_si128(&pvE[index], e.m);
            ++index;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vQIndexLo16 = vQIndexLo16_reset;
        __m128i vQIndexHi16 = vQIndexHi16_reset;
        __m128i vE;
        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        __m128i vF = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 1);

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vH = _mm_adds_epi8(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi8(vH, vE);
            vH = _mm_max_epi8(vH, vF);
            vH = _mm_max_epi8(vH, vZero);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
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
                vQIndexLo16 = _mm_add_epi16(vQIndexLo16, vOne16);
                vQIndexHi16 = _mm_add_epi16(vQIndexHi16, vOne16);
            }

            /* Update vE value. */
            vH = _mm_subs_epi8(vH, vGapO);
            vE = _mm_subs_epi8(vE, vGapE);
            vE = _mm_max_epi8(vE, vH);
            _mm_store_si128(pvE + i, vE);

            /* Update vF value. */
            vF = _mm_subs_epi8(vF, vGapE);
            vF = _mm_max_epi8(vF, vH);

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<8; ++k) {
            vF = _mm_slli_si128(vF, 1);
            for (i=0; i<segLen; ++i) {
                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi8(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
                arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm_subs_epi8(vH, vGapO);
                vF = _mm_subs_epi8(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi8(vF, vH))) goto end;
                vF = _mm_max_epi8(vF, vH);
            }
        }
end:
        {
        }
    }

    /* max in vec */
    for (j=0; j<8; ++j) {
        int8_t value = (int8_t) _mm_extract_epi8(vMaxH, 15);
        if (value > score) {
            score = value;
        }
        vMaxH = _mm_slli_si128(vMaxH, 1);
    }

    result->score = score;

    free(pvE);
    free(pvHLoad);
    free(pvHStore);
    free(vProfile);

    return result;
}

