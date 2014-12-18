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
#define MAX(a,b) ((a)>(b)?(a):(b))

/* avx2 does not have _mm256_insert_epi32, emulate it */
static inline __m256i _mm256_insert_epi32(__m256i a, int32_t i, int imm) {
    __m256i_32_t tmp;
    tmp.m = a;
    tmp.v[imm] = i;
    return tmp.m;
}

/* avx2 does not have _mm256_extract_epi32, emulate it */
static inline int32_t _mm256_extract_epi32(__m256i a, int imm) {
    __m256i_32_t tmp;
    tmp.m = a;
    return tmp.v[imm];
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
#define FNAME nw_table_scan_avx2_256_32
#else
#define FNAME nw_scan_avx2_256_32
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
    const int32_t segWidth = 8; /* number of values in vector unit */
    int32_t segNum = 0;
    int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m256i* const restrict pvP = parasail_memalign_m256i(32, n * segLen);
    __m256i* const restrict pvE = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvHt= parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvFt= parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvH = parasail_memalign_m256i(32, segLen);
    int32_t* const restrict boundary = parasail_memalign_int32_t(32, s2Len+1);
    __m256i vGapO = _mm256_set1_epi32(open);
    __m256i vGapE = _mm256_set1_epi32(gap);
    int32_t score = 0;
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
                __m256i_32_t t;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = matrix[k][MAP_BLOSUM_[(unsigned char)s1[j]]];
                    j += segLen;
                }
                _mm256_store_si256(&pvP[index], t.m);
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
                h.v[segNum] = -open-gap*(segNum*segLen+i);
                e.v[segNum] = NEG_INF_32;
            }
            _mm256_store_si256(&pvH[index], h.m);
            _mm256_store_si256(&pvE[index], e.m);
            ++index;
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
        __m256i vE;
        __m256i vHt;
        __m256i vFt;
        __m256i vH;
        __m256i vHp;
        __m256i *pvW;
        __m256i vW;

        /* calculate E */
        /* calculate Ht */
        vHp = _mm256_slli_si256(_mm256_load_si256(pvH+(segLen-1)), 2);
        vHp = _mm256_insert_epi32(vHp, boundary[j], 0);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            vH = _mm256_load_si256(pvH+i);
            vE = _mm256_load_si256(pvE+i);
            vW = _mm256_load_si256(pvW+i);
            vE = _mm256_max_epi32(
                    _mm256_sub_epi32(vE, vGapE),
                    _mm256_sub_epi32(vH, vGapO));
            vHt = _mm256_max_epi32(
                    _mm256_add_epi32(vHp, vW),
                    vE);
            _mm256_store_si256(pvE+i, vE);
            _mm256_store_si256(pvHt+i, vHt);
            vHp = vH;
        }

        /* calculate Ft */
        vHt = _mm256_slli_si256(_mm256_load_si256(pvHt+(segLen-1)), 2);
        vHt = _mm256_insert_epi32(vHt, boundary[j+1], 0);
        vFt = _mm256_set1_epi32(NEG_INF_32);
        for (i=0; i<segLen; ++i) {
            vFt = _mm256_max_epi32(
                    _mm256_sub_epi32(vFt, vGapE),
                    vHt);
            vHt = _mm256_load_si256(pvHt+i);
        }
        {
            __m256i_32_t tmp;
            tmp.m = vFt;
            tmp.v[1] = MAX(tmp.v[0]-segLen*gap, tmp.v[1]);
            tmp.v[2] = MAX(tmp.v[1]-segLen*gap, tmp.v[2]);
            tmp.v[3] = MAX(tmp.v[2]-segLen*gap, tmp.v[3]);
            tmp.v[4] = MAX(tmp.v[3]-segLen*gap, tmp.v[4]);
            tmp.v[5] = MAX(tmp.v[4]-segLen*gap, tmp.v[5]);
            tmp.v[6] = MAX(tmp.v[5]-segLen*gap, tmp.v[6]);
            tmp.v[7] = MAX(tmp.v[6]-segLen*gap, tmp.v[7]);
            vFt = tmp.m;
        }
        vHt = _mm256_slli_si256(_mm256_load_si256(pvHt+(segLen-1)), 2);
        vHt = _mm256_insert_epi32(vHt, boundary[j+1], 0);
        vFt = _mm256_slli_si256(vFt, 2);
        vFt = _mm256_insert_epi32(vFt, NEG_INF_32, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm256_max_epi32(
                    _mm256_sub_epi32(vFt, vGapE),
                    vHt);
            vHt = _mm256_load_si256(pvHt+i);
            _mm256_store_si256(pvFt+i, vFt);
        }

        /* calculate H */
        for (i=0; i<segLen; ++i) {
            vHt = _mm256_load_si256(pvHt+i);
            vFt = _mm256_load_si256(pvFt+i);
            vH = _mm256_max_epi32(
                    vHt,
                    _mm256_sub_epi32(vFt, vGapO));
            _mm256_store_si256(pvH+i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif
        }
    }

    /* extract last value from the last column */
    {
        __m256i vH = _mm256_load_si256(pvH + offset);
        for (k=0; k<position; ++k) {
            vH = _mm256_slli_si256(vH, 2);
        }
        score = (int32_t) _mm256_extract_epi32 (vH, 7);
    }

    result->score = score;

    free(boundary);
    free(pvH);
    free(pvFt);
    free(pvHt);
    free(pvE);
    free(pvP);

    return result;
}

