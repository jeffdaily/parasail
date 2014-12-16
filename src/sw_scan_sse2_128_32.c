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

#define NEG_INF_32 (INT32_MIN/(int32_t)(2))
#define MAX(a,b) ((a)>(b)?(a):(b))

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
#define FNAME sw_table_scan_sse2_128_32
#else
#define FNAME sw_scan_sse2_128_32
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    const int32_t n = 24; /* number of amino acids in table */
    const int32_t segWidth = 4; /* number of values in vector unit */
    int32_t segNum = 0;
    int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m128i* restrict pvP = parasail_memalign_m128i(16, n * segLen);
    __m128i* restrict pvE = parasail_memalign_m128i(16, segLen);
    __m128i* restrict pvHt= parasail_memalign_m128i(16, segLen);
    __m128i* restrict pvFt= parasail_memalign_m128i(16, segLen);
    __m128i* restrict pvH = parasail_memalign_m128i(16, segLen);
    __m128i vGapO = _mm_set1_epi32(open);
    __m128i vGapE = _mm_set1_epi32(gap);
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi32(1);
    int32_t score = NEG_INF_32;
#if PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
    __m128i vMaxH = _mm_set1_epi32(NEG_INF_32);
    __m128i vQLimit = _mm_set1_epi32(s1Len);
    __m128i vQIndex_reset = _mm_set_epi32(
            3*segLen,
            2*segLen,
            1*segLen,
            0*segLen);

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
                _mm_store_si128(&pvP[index], t.m);
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
            _mm_store_si128(&pvH[index], h.m);
            _mm_store_si128(&pvE[index], e.m);
            ++index;
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
        __m128i vQIndex = vQIndex_reset;

        /* calculate E */
        /* calculate Ht */
        vHp = _mm_slli_si128(_mm_load_si128(pvH+(segLen-1)), 4);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            vE = _mm_max_epi32(
                    _mm_sub_epi32(vE, vGapE),
                    _mm_sub_epi32(vH, vGapO));
            vHt = _mm_max_epi32(
                    _mm_add_epi32(vHp, vW),
                    vE);
            _mm_store_si128(pvE+i, vE);
            _mm_store_si128(pvHt+i, vHt);
            vHp = vH;
        }

        /* calculate Ft */
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 4);
        vFt = _mm_set1_epi32(NEG_INF_32);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_max_epi32(
                    _mm_sub_epi32(vFt, vGapE),
                    vHt);
            vHt = _mm_load_si128(pvHt+i);
        }
        {
            __m128i_32_t tmp;
            tmp.m = vFt;
            tmp.v[1] = MAX(tmp.v[0]-segLen*gap, tmp.v[1]);
            tmp.v[2] = MAX(tmp.v[1]-segLen*gap, tmp.v[2]);
            tmp.v[3] = MAX(tmp.v[2]-segLen*gap, tmp.v[3]);
            vFt = tmp.m;
        }
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 4);
        vFt = _mm_slli_si128(vFt, 4);
        vFt = _mm_insert_epi32(vFt, NEG_INF_32, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_max_epi32(
                    _mm_sub_epi32(vFt, vGapE),
                    vHt);
            vHt = _mm_load_si128(pvHt+i);
            _mm_store_si128(pvFt+i, vFt);
        }

        /* calculate H */
        for (i=0; i<segLen; ++i) {
            vHt = _mm_load_si128(pvHt+i);
            vFt = _mm_load_si128(pvFt+i);
            vH = _mm_max_epi32(
                    vHt,
                    _mm_sub_epi32(vFt, vGapO));
            vH = _mm_max_epi32(vH, vZero);
            _mm_store_si128(pvH+i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif
            /* update max vector seen so far */
            {
                __m128i cond_max = _mm_cmpgt_epi32(vH,vMaxH);
                __m128i cond_lmt = _mm_cmplt_epi32(vQIndex,vQLimit);
                __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
                vMaxH = _mm_andnot_si128(cond_all, vMaxH);
                vMaxH = _mm_or_si128(vMaxH, _mm_and_si128(cond_all, vH));
                vQIndex = _mm_add_epi32(vQIndex, vOne);
            }
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

    free(pvH);
    free(pvFt);
    free(pvHt);
    free(pvE);
    free(pvP);

    return result;
}
