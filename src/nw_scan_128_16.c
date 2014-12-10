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

#ifdef PARASAIL_TABLE
#include "align_scan_128_16_table.h"
#else
#include "align_scan_128_16.h"
#endif
#include "blosum/blosum_map.h"


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

#define PARALLEL_PREFIX_OP(vFt, gap, segLen)            \
{                                                       \
    union {                                             \
        __m128i m;                                      \
        int16_t v[8];                                   \
    } tmp;                                              \
    tmp.m = vFt;                                        \
    tmp.v[1] = MAX(tmp.v[0]-segLen*gap, tmp.v[1]);      \
    tmp.v[2] = MAX(tmp.v[1]-segLen*gap, tmp.v[2]);      \
    tmp.v[3] = MAX(tmp.v[2]-segLen*gap, tmp.v[3]);      \
    tmp.v[4] = MAX(tmp.v[3]-segLen*gap, tmp.v[4]);      \
    tmp.v[5] = MAX(tmp.v[4]-segLen*gap, tmp.v[5]);      \
    tmp.v[6] = MAX(tmp.v[5]-segLen*gap, tmp.v[6]);      \
    tmp.v[7] = MAX(tmp.v[6]-segLen*gap, tmp.v[7]);      \
    vFt = tmp.m;                                        \
}

#ifdef PARASAIL_TABLE
#define FNAME nw_scan_128_16_table
#else
#define FNAME nw_scan_128_16
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
    const int32_t n = 24; /* number of amino acids in table */
    const int32_t segWidth = 8; /* number of values in vector unit */
    int32_t segNum = 0;
    int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m128i* pvP = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    __m128i* pvE = (__m128i*)malloc(segLen * sizeof(__m128i));
    __m128i* pvHt = (__m128i*)malloc(segLen * sizeof(__m128i));
    __m128i* pvFt = (__m128i*)malloc(segLen * sizeof(__m128i));
    __m128i* pvH = (__m128i*)malloc(segLen * sizeof(__m128i));
    int16_t* boundary = (int16_t*)malloc((s2Len+1) * sizeof(int16_t));
    __m128i vGapO = _mm_set1_epi16(open);
    __m128i vGapE = _mm_set1_epi16(gap);
    int16_t score = 0;

    /* Generate query profile.
     * Rearrange query sequence & calculate the weight of match/mismatch */
    {
        int16_t *t = (int16_t*)pvP;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                int32_t j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    *t++ = matrix[k*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    j += segLen;
                }
            }
        }
    }

    /* initialize H and E */
    {
        int16_t *h = (int16_t*)pvH;
        int16_t *e = (int16_t*)pvE;
        for (i=0; i<segLen; ++i) {
            for (segNum=0; segNum<segWidth; ++segNum) {
                *h = -open-gap*(segNum*segLen+i);
                *e = NEG_INF_16;
                ++h;
                ++e;
            }
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            boundary[i] = -open-gap*(i-1);
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
        vHp = _mm_slli_si128(_mm_load_si128(pvH+(segLen-1)), 2);
        vHp = _mm_insert_epi16(vHp, boundary[j], 0);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            vE = _mm_max_epi16(
                    _mm_sub_epi16(vE, vGapE),
                    _mm_sub_epi16(vH, vGapO));
            vHt = _mm_max_epi16(
                    _mm_add_epi16(vHp, vW),
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
            vE = _mm_max_epi16(
                    _mm_sub_epi16(vE, vGapE),
                    _mm_sub_epi16(vH, vGapO));
            _mm_store_si128(pvE+i, vE);
        }

        /* calculate Ht */
        vH = _mm_slli_si128(_mm_load_si128(pvH+(segLen-1)), 2);
        vH = _mm_insert_epi16(vH, boundary[j], 0);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            vHt = _mm_max_epi16(
                    _mm_add_epi16(vH, vW),
                    vE);
            vH = _mm_load_si128(pvH+i);
            _mm_store_si128(pvHt+i, vHt);
        }
#endif

        /* calculate Ft */
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 2);
        vHt = _mm_insert_epi16(vHt, boundary[j+1], 0);
        vFt = _mm_set1_epi16(NEG_INF_16);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_max_epi16(
                    _mm_sub_epi16(vFt, vGapE),
                    vHt);
            vHt = _mm_load_si128(pvHt+i);
        }
        PARALLEL_PREFIX_OP(vFt, gap, segLen)
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 2);
        vHt = _mm_insert_epi16(vHt, boundary[j+1], 0);
        vFt = _mm_slli_si128(vFt, 2);
        vFt = _mm_insert_epi16(vFt, NEG_INF_16, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_max_epi16(
                    _mm_sub_epi16(vFt, vGapE),
                    vHt);
            vHt = _mm_load_si128(pvHt+i);
            _mm_store_si128(pvFt+i, vFt);
        }

        /* calculate H */
        for (i=0; i<segLen; ++i) {
            vHt = _mm_load_si128(pvHt+i);
            vFt = _mm_load_si128(pvFt+i);
            vH = _mm_max_epi16(
                    vHt,
                    _mm_sub_epi16(vFt, vGapO));
            _mm_store_si128(pvH+i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si128(score_table, vH, i, segLen, j, s2Len);
#endif
        }
    }

    /* extract last value from the last column */
    {
        __m128i vH = _mm_load_si128(pvH + offset);
        for (k=0; k<position; ++k) {
            vH = _mm_slli_si128(vH, 2);
        }
        score = (int16_t) _mm_extract_epi16 (vH, 7);
    }

    free(pvP);
    free(pvE);
    free(pvHt);
    free(pvFt);
    free(pvH);
    free(boundary);

    return score;
}
