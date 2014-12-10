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

#ifdef PARASAIL_TABLE
#include "align_striped_128_8_table.h"
#else
#include "align_striped_128_8.h"
#endif
#include "blosum/blosum_map.h"


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
#define FNAME sg_striped_128_8_table
#else
#define FNAME sg_striped_128_8
#endif

int FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix
#ifdef PARASAIL_TABLE
        , int * const restrict score_table
#endif
        )
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t nt = 0;
    int32_t segNum = 0;
    const int32_t n = 24; /* number of amino acids in table */
    int32_t segLen = (s1Len + 15) / 16;
    __m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));

    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = 15 - (s1Len - 1) / segLen;

    /* the max alignment score */
    int max = NEG_INF_8;

    /* Define 16 byte 0 vector. */
    __m128i vNegInf = _mm_set1_epi8(NEG_INF_8);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));

    __m128i vSaturationCheck = _mm_setzero_si128();
    __m128i vNegLimit = _mm_set1_epi8(INT8_MIN);
    __m128i vPosLimit = _mm_set1_epi8(INT8_MAX);

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi8(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi8(gap);

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    {
        int8_t *t = (int8_t*)vProfile;
        for (nt=0; nt<n; ++nt) {
            for (i=0; i<segLen; ++i) {
                int32_t j = i;
                for (segNum=0; segNum<16; ++segNum) {
                    *t++ = matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    j += segLen;
                }
            }
        }
    }

    /* initialize E */
    {
        int8_t *e = (int8_t*)pvE;
        for (i=0; i<segLen; ++i) {
            for (segNum=0; segNum<16; ++segNum) {
                *e++ = -open;
            }
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vE;
        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        __m128i vF = vNegInf;

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
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);
            /* check for saturation */
            {
                vSaturationCheck = _mm_or_si128(vSaturationCheck,
                        _mm_or_si128(
                            _mm_cmpeq_epi8(vH, vNegLimit),
                            _mm_cmpeq_epi8(vH, vPosLimit)));
            }
#ifdef PARASAIL_TABLE
            arr_store_si128(score_table, vH, i, segLen, j, s2Len);
#endif

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
        for (k=0; k<16; ++k) {
            vF = _mm_slli_si128(vF, 1);
            vF = _mm_insert_epi8(vF, -open, 0);
            for (i=0; i<segLen; ++i) {
                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi8(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
                /* check for saturation */
                {
                    vSaturationCheck = _mm_or_si128(vSaturationCheck,
                            _mm_or_si128(
                                _mm_cmpeq_epi8(vH, vNegLimit),
                                _mm_cmpeq_epi8(vH, vPosLimit)));
                }
#ifdef PARASAIL_TABLE
                arr_store_si128(score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm_subs_epi8(vH, vGapO);
                vF = _mm_subs_epi8(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi8(vF, vH))) goto end;
                vF = _mm_max_epi8(vF, vH);
            }
        }
end:
        {
            /* extract last value from the column */
            int8_t tmp;
            vH = _mm_load_si128(pvHStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128 (vH, 1);
            }
            /* max of last value in each column */
            tmp = (int8_t) _mm_extract_epi8 (vH, 15);
            if (tmp > max) {
                max = tmp;
            }
        }
    }

    /* max of last column */
    {
        __m128i vMaxLastColH = vNegInf;
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
            __m128i vH = _mm_load_si128(pvHStore + i);
            __m128i cond_lmt = _mm_packs_epi16(
                    _mm_cmplt_epi16(vQIndexLo16, vQLimit16),
                    _mm_cmplt_epi16(vQIndexHi16, vQLimit16));
            __m128i cond_max = _mm_cmpgt_epi8(vH, vMaxLastColH);
            __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
            vMaxLastColH = _mm_andnot_si128(cond_all, vMaxLastColH);
            vMaxLastColH = _mm_or_si128(vMaxLastColH, _mm_and_si128(cond_all, vH));
            vQIndexLo16 = _mm_adds_epi8(vQIndexLo16, vOne16);
            vQIndexHi16 = _mm_adds_epi8(vQIndexHi16, vOne16);
        }

        /* max in vec */
        for (j=0; j<16; ++j) {
            int8_t value = (int8_t) _mm_extract_epi8(vMaxLastColH, 15);
            if (value > max) {
                max = value;
            }
            vMaxLastColH = _mm_slli_si128(vMaxLastColH, 1);
        }
    }

    if (_mm_movemask_epi8(vSaturationCheck)) {
        max = INT8_MAX;
    }

    free(vProfile);
    free(pvHStore);
    free(pvHLoad);
    free(pvE);

    return max;
}

