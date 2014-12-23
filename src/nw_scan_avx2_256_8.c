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

#define NEG_INF_8 (INT8_MIN)
#define MAX(a,b) ((a)>(b)?(a):(b))

#if 0
/* avx2 does not have _mm256_insert_epi8, emulate it */
static inline __m256i _mm256_insert_epi8(__m256i a, int8_t i, int imm) {
    __m256i_8_t tmp;
    tmp.m = a;
    tmp.v[imm] = i;
    return tmp.m;
}

/* avx2 does not have _mm256_extract_epi8, emulate it */
static inline int8_t _mm256_extract_epi8(__m256i a, int imm) {
    __m256i_8_t tmp;
    tmp.m = a;
    return tmp.v[imm];
}
#endif

/* avx2 _mm256_slli_si256 does not shift across 128-bit lanes, emulate it */
static inline __m256i shift(__m256i a, __m256i idx) {
#if 0
    return  _mm256_permutevar8x8_epi8(a, idx);
#else
    __m256i_8_t tmp;
    tmp.m = a;
    tmp.v[31] = tmp.v[30];
    tmp.v[30] = tmp.v[29];
    tmp.v[29] = tmp.v[28];
    tmp.v[28] = tmp.v[27];
    tmp.v[27] = tmp.v[26];
    tmp.v[26] = tmp.v[25];
    tmp.v[25] = tmp.v[24];
    tmp.v[24] = tmp.v[23];
    tmp.v[23] = tmp.v[22];
    tmp.v[22] = tmp.v[21];
    tmp.v[21] = tmp.v[20];
    tmp.v[20] = tmp.v[19];
    tmp.v[19] = tmp.v[18];
    tmp.v[18] = tmp.v[17];
    tmp.v[17] = tmp.v[16];
    tmp.v[16] = tmp.v[15];
    tmp.v[15] = tmp.v[14];
    tmp.v[14] = tmp.v[13];
    tmp.v[13] = tmp.v[12];
    tmp.v[12] = tmp.v[11];
    tmp.v[11] = tmp.v[10];
    tmp.v[10] = tmp.v[ 9];
    tmp.v[ 9] = tmp.v[ 8];
    tmp.v[ 8] = tmp.v[ 7];
    tmp.v[ 7] = tmp.v[ 6];
    tmp.v[ 6] = tmp.v[ 5];
    tmp.v[ 5] = tmp.v[ 4];
    tmp.v[ 4] = tmp.v[ 3];
    tmp.v[ 3] = tmp.v[ 2];
    tmp.v[ 2] = tmp.v[ 1];
    tmp.v[ 1] = tmp.v[ 0];
    tmp.v[ 0] = 0;
    return tmp.m;
#endif
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
    array[( 0*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  0);
    array[( 1*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  1);
    array[( 2*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  2);
    array[( 3*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  3);
    array[( 4*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  4);
    array[( 5*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  5);
    array[( 6*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  6);
    array[( 7*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  7);
    array[( 8*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  8);
    array[( 9*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  9);
    array[(10*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 10);
    array[(11*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 11);
    array[(12*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 12);
    array[(13*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 13);
    array[(14*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 14);
    array[(15*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 15);
    array[(16*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 16);
    array[(17*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 17);
    array[(18*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 18);
    array[(19*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 19);
    array[(20*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 20);
    array[(21*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 21);
    array[(22*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 22);
    array[(23*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 23);
    array[(24*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 24);
    array[(25*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 25);
    array[(26*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 26);
    array[(27*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 27);
    array[(28*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 28);
    array[(29*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 29);
    array[(30*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 30);
    array[(31*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 31);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME nw_table_scan_avx2_256_8
#else
#define FNAME nw_scan_avx2_256_8
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
    const int32_t segWidth = 32; /* number of values in vector unit */
    int32_t segNum = 0;
    int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m256i* const restrict pvP = parasail_memalign_m256i(32, n * segLen);
    __m256i* const restrict pvE = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvHt= parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvFt= parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvH = parasail_memalign_m256i(32, segLen);
    int8_t* const restrict boundary = parasail_memalign_int8_t(32, s2Len+1);
    __m256i vGapO = _mm256_set1_epi8(open);
    __m256i vGapE = _mm256_set1_epi8(gap);
    __m256i vSaturationCheck = _mm256_setzero_si256();
    __m256i vNegLimit = _mm256_set1_epi8(INT8_MIN);
    __m256i vPosLimit = _mm256_set1_epi8(INT8_MAX);
    __m256i idx = _mm256_set_epi8(30,29,28,27,26,25,24,23,
                                  22,21,20,19,18,17,16,15,
                                  14,13,12,11,10, 9, 8, 7,
                                   6, 5, 4, 3, 2, 1, 0,31);
    int8_t score = 0;
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
                __m256i_8_t t;
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
            __m256i_8_t h;
            __m256i_8_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = -open-gap*(segNum*segLen+i);
                e.v[segNum] = NEG_INF_8;
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
        vHp = _mm256_load_si256(pvH+(segLen-1));
        vHp = shift(vHp, idx);
        vHp = _mm256_insert_epi8(vHp, boundary[j], 0);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            vH = _mm256_load_si256(pvH+i);
            vE = _mm256_load_si256(pvE+i);
            vW = _mm256_load_si256(pvW+i);
            vE = _mm256_max_epi8(
                    _mm256_subs_epi8(vE, vGapE),
                    _mm256_subs_epi8(vH, vGapO));
            vHt = _mm256_max_epi8(
                    _mm256_adds_epi8(vHp, vW),
                    vE);
            _mm256_store_si256(pvE+i, vE);
            _mm256_store_si256(pvHt+i, vHt);
            vHp = vH;
        }

        /* calculate Ft */
        vHt = _mm256_load_si256(pvHt+(segLen-1));
        vHt = shift(vHt, idx);
        vHt = _mm256_insert_epi8(vHt, boundary[j+1], 0);
        vFt = _mm256_set1_epi8(NEG_INF_8);
        for (i=0; i<segLen; ++i) {
            vFt = _mm256_max_epi8(
                    _mm256_subs_epi8(vFt, vGapE),
                    vHt);
            vHt = _mm256_load_si256(pvHt+i);
        }
        {
            __m256i_8_t tmp;
            tmp.m = vFt;
            tmp.v[ 1] = MAX(tmp.v[ 0]-segLen*gap, tmp.v[ 1]);
            tmp.v[ 2] = MAX(tmp.v[ 1]-segLen*gap, tmp.v[ 2]);
            tmp.v[ 3] = MAX(tmp.v[ 2]-segLen*gap, tmp.v[ 3]);
            tmp.v[ 4] = MAX(tmp.v[ 3]-segLen*gap, tmp.v[ 4]);
            tmp.v[ 5] = MAX(tmp.v[ 4]-segLen*gap, tmp.v[ 5]);
            tmp.v[ 6] = MAX(tmp.v[ 5]-segLen*gap, tmp.v[ 6]);
            tmp.v[ 7] = MAX(tmp.v[ 6]-segLen*gap, tmp.v[ 7]);
            tmp.v[ 8] = MAX(tmp.v[ 7]-segLen*gap, tmp.v[ 8]);
            tmp.v[ 9] = MAX(tmp.v[ 8]-segLen*gap, tmp.v[ 9]);
            tmp.v[10] = MAX(tmp.v[ 9]-segLen*gap, tmp.v[10]);
            tmp.v[11] = MAX(tmp.v[10]-segLen*gap, tmp.v[11]);
            tmp.v[12] = MAX(tmp.v[11]-segLen*gap, tmp.v[12]);
            tmp.v[13] = MAX(tmp.v[12]-segLen*gap, tmp.v[13]);
            tmp.v[14] = MAX(tmp.v[13]-segLen*gap, tmp.v[14]);
            tmp.v[15] = MAX(tmp.v[14]-segLen*gap, tmp.v[15]);
            tmp.v[16] = MAX(tmp.v[15]-segLen*gap, tmp.v[16]);
            tmp.v[17] = MAX(tmp.v[16]-segLen*gap, tmp.v[17]);
            tmp.v[18] = MAX(tmp.v[17]-segLen*gap, tmp.v[18]);
            tmp.v[19] = MAX(tmp.v[18]-segLen*gap, tmp.v[19]);
            tmp.v[20] = MAX(tmp.v[19]-segLen*gap, tmp.v[20]);
            tmp.v[21] = MAX(tmp.v[20]-segLen*gap, tmp.v[21]);
            tmp.v[22] = MAX(tmp.v[21]-segLen*gap, tmp.v[22]);
            tmp.v[23] = MAX(tmp.v[22]-segLen*gap, tmp.v[23]);
            tmp.v[24] = MAX(tmp.v[23]-segLen*gap, tmp.v[24]);
            tmp.v[25] = MAX(tmp.v[24]-segLen*gap, tmp.v[25]);
            tmp.v[26] = MAX(tmp.v[25]-segLen*gap, tmp.v[26]);
            tmp.v[27] = MAX(tmp.v[26]-segLen*gap, tmp.v[27]);
            tmp.v[28] = MAX(tmp.v[27]-segLen*gap, tmp.v[28]);
            tmp.v[29] = MAX(tmp.v[28]-segLen*gap, tmp.v[29]);
            tmp.v[30] = MAX(tmp.v[29]-segLen*gap, tmp.v[30]);
            tmp.v[31] = MAX(tmp.v[30]-segLen*gap, tmp.v[31]);
            vFt = tmp.m;
        }
        vHt = _mm256_load_si256(pvHt+(segLen-1));
        vHt = shift(vHt, idx);
        vHt = _mm256_insert_epi8(vHt, boundary[j+1], 0);
        vFt = shift(vFt, idx);
        vFt = _mm256_insert_epi8(vFt, NEG_INF_8, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm256_max_epi8(
                    _mm256_subs_epi8(vFt, vGapE),
                    vHt);
            vHt = _mm256_load_si256(pvHt+i);
            _mm256_store_si256(pvFt+i, vFt);
        }

        /* calculate H */
        for (i=0; i<segLen; ++i) {
            vHt = _mm256_load_si256(pvHt+i);
            vFt = _mm256_load_si256(pvFt+i);
            vH = _mm256_max_epi8(
                    vHt,
                    _mm256_subs_epi8(vFt, vGapO));
            _mm256_store_si256(pvH+i, vH);
            /* check for saturation */
            {
                vSaturationCheck = _mm256_or_si256(vSaturationCheck,
                        _mm256_or_si256(
                            _mm256_cmpeq_epi8(vH, vNegLimit),
                            _mm256_cmpeq_epi8(vH, vPosLimit)));
            }
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif
        }
    }

    /* extract last value from the last column */
    {
        __m256i vH = _mm256_load_si256(pvH + offset);
        for (k=0; k<position; ++k) {
            vH = shift(vH, idx);
        }
        score = (int8_t) _mm256_extract_epi8 (vH, 31);
    }

    if (_mm256_movemask_epi8(vSaturationCheck)) {
        score = INT8_MAX;
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

