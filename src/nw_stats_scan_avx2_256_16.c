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

#define NEG_INF_16 (INT16_MIN/(int16_t)(2))
#define MAX(a,b) ((a)>(b)?(a):(b))

/* avx2 does not have _mm256_cmplt_epi16, emulate it */
static inline __m256i _mm256_cmplt_epi16(__m256i a, __m256i b) {
    return _mm256_cmpgt_epi16(b,a);
}

#if HAVE_AVX2_MM256_INSERT_EPI16
#else
static inline __m256i _mm256_insert_epi16(__m256i a, int16_t b, int imm) {
    __m256i_16_t tmp;
    tmp.m = a;
    tmp.v[imm] = b;
    return tmp.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI16
#else
static inline int16_t _mm256_extract_epi16(__m256i a, int imm) {
    __m256i_16_t tmp;
    tmp.m = a;
    return tmp.v[imm];
}
#endif

/* avx2 _mm256_slli_si256 does not shift across 128-bit lanes, emulate it */
static inline __m256i shift(__m256i a) {
    return _mm256_alignr_epi8(a,
            _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)),
            14);
}

/* avx2 _mm256_srli_si256 does not shift across 128-bit lanes, emulate it */
static inline __m256i rshift(__m256i a) {
    return _mm256_or_si256(
            _mm256_slli_si256(
                _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(3,0,0,1)),
                14),
            _mm256_srli_si256(a, 2));
}

static inline __m256i lrotate16(__m256i a) {
    return _mm256_alignr_epi8(a,
            _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,1)),
            14);
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
    array[( 0*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 0);
    array[( 1*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 1);
    array[( 2*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 2);
    array[( 3*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 3);
    array[( 4*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 4);
    array[( 5*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 5);
    array[( 6*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 6);
    array[( 7*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 7);
    array[( 8*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 8);
    array[( 9*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 9);
    array[(10*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 10);
    array[(11*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 11);
    array[(12*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 12);
    array[(13*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 13);
    array[(14*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 14);
    array[(15*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 15);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME nw_stats_table_scan_avx2_256_16
#else
#define FNAME nw_stats_scan_avx2_256_16
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
    const int32_t segWidth = 16; /* number of values in vector unit */
    int32_t segNum = 0;
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m256i* const restrict pvP  = parasail_memalign_m256i(32, n * segLen);
    __m256i* const restrict pvPm = parasail_memalign_m256i(32, n * segLen);
    __m256i* const restrict pvE  = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvHt = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvFt = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvMt = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvLt = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvEx = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvH  = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvM  = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvL  = parasail_memalign_m256i(32, segLen);
    int16_t* const restrict boundary = parasail_memalign_int16_t(32, s2Len+1);
    __m256i vGapO = _mm256_set1_epi16(open);
    __m256i vGapE = _mm256_set1_epi16(gap);
    __m256i vZero = _mm256_setzero_si256();
    __m256i vOne = _mm256_set1_epi16(1);
    int16_t score = 0;
    int16_t matches = 0;
    int16_t length = 0;
    __m256i segLenXgap_reset = _mm256_set_epi16(
            NEG_INF_16, NEG_INF_16, NEG_INF_16, NEG_INF_16,
            NEG_INF_16, NEG_INF_16, NEG_INF_16, NEG_INF_16,
            NEG_INF_16, NEG_INF_16, NEG_INF_16, NEG_INF_16,
            NEG_INF_16, NEG_INF_16, NEG_INF_16, -segLen*gap);
    __m256i insert_mask = _mm256_cmpeq_epi16(_mm256_setzero_si256(),
            _mm256_set_epi16(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1));
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif

    parasail_memset_m256i(pvM, vZero, segLen);
    parasail_memset_m256i(pvL, vZero, segLen);

    /* Generate query profile and match profile.
     * Rearrange query sequence & calculate the weight of match/mismatch.
     * Don't alias. */
    {
        int32_t index = 0;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                __m256i_16_t t;
                __m256i_16_t s;
                j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = matrix[k][MAP_BLOSUM_[(unsigned char)s1[j]]];
                    s.v[segNum] = (k == MAP_BLOSUM_[(unsigned char)s1[j]]);
                    j += segLen;
                    if (j >= s1Len) break;
                }
                _mm256_store_si256(&pvP[index], t.m);
                _mm256_store_si256(&pvPm[index], s.m);
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
                h.v[segNum] = -open-gap*(segNum*segLen+i);
                e.v[segNum] = NEG_INF_16;
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
        __m256i *pvW;
        __m256i vW;
        __m256i *pvC;
        __m256i vC;
        __m256i vM;
        __m256i vMp;
        __m256i vMt;
        __m256i vL;
        __m256i vLp;
        __m256i vLt;
        __m256i vEx;

        /* calculate E */
        for (i=0; i<segLen; ++i) {
            vH = _mm256_load_si256(pvH+i);
            vE = _mm256_load_si256(pvE+i);
            vE = _mm256_max_epi16(
                    _mm256_sub_epi16(vE, vGapE),
                    _mm256_sub_epi16(vH, vGapO));
            _mm256_store_si256(pvE+i, vE);
        }

        /* calculate Ht */
        vH = shift(_mm256_load_si256(pvH+(segLen-1)));
        vH = _mm256_insert_epi16(vH, boundary[j], 0);
        vMp= shift(_mm256_load_si256(pvM+(segLen-1)));
        vLp= shift(_mm256_load_si256(pvL+(segLen-1)));
        vLp= _mm256_add_epi16(vLp, vOne);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        pvC = pvPm+ MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vE = _mm256_load_si256(pvE+i);
            vW = _mm256_load_si256(pvW+i);
            /* compute */
            vH = _mm256_add_epi16(vH, vW);
            vHt = _mm256_max_epi16(vH, vE);
            /* statistics */
            vC = _mm256_load_si256(pvC+i);
            vMp = _mm256_add_epi16(vMp, vC);
            vEx = _mm256_cmpgt_epi16(vE, vH);
            vM = _mm256_load_si256(pvM+i);
            vL = _mm256_load_si256(pvL+i);
            vL = _mm256_add_epi16(vL, vOne);
            vMt = _mm256_blendv_epi8(vMp, vM, vEx);
            vLt = _mm256_blendv_epi8(vLp, vL, vEx);
            /* store results */
            _mm256_store_si256(pvHt+i, vHt);
            _mm256_store_si256(pvEx+i, vEx);
            _mm256_store_si256(pvMt+i, vMt);
            _mm256_store_si256(pvLt+i, vLt);
            /* prep for next iteration */
            vH = _mm256_load_si256(pvH+i);
            vMp = vM;
            vLp = vL;
        }

        /* calculate Ft */
        vHt = shift(_mm256_load_si256(pvHt+(segLen-1)));
        vHt = _mm256_insert_epi16(vHt, boundary[j+1], 0);
        vFt = _mm256_set1_epi16(NEG_INF_16);
        for (i=0; i<segLen; ++i) {
            vFt = _mm256_sub_epi16(vFt, vGapE),
            vFt = _mm256_max_epi16(vFt, vHt);
            vHt = _mm256_load_si256(pvHt+i);
        }
#if 0
        {
            __m256i_16_t tmp;
            tmp.m = vFt;
            tmp.v[1] = MAX(tmp.v[0]-segLen*gap, tmp.v[1]);
            tmp.v[2] = MAX(tmp.v[1]-segLen*gap, tmp.v[2]);
            tmp.v[3] = MAX(tmp.v[2]-segLen*gap, tmp.v[3]);
            tmp.v[4] = MAX(tmp.v[3]-segLen*gap, tmp.v[4]);
            tmp.v[5] = MAX(tmp.v[4]-segLen*gap, tmp.v[5]);
            tmp.v[6] = MAX(tmp.v[5]-segLen*gap, tmp.v[6]);
            tmp.v[7] = MAX(tmp.v[6]-segLen*gap, tmp.v[7]);
            tmp.v[8]  = MAX(tmp.v[7] -segLen*gap, tmp.v[8]);
            tmp.v[9]  = MAX(tmp.v[8] -segLen*gap, tmp.v[9]);
            tmp.v[10] = MAX(tmp.v[9] -segLen*gap, tmp.v[10]);
            tmp.v[11] = MAX(tmp.v[10]-segLen*gap, tmp.v[11]);
            tmp.v[12] = MAX(tmp.v[11]-segLen*gap, tmp.v[12]);
            tmp.v[13] = MAX(tmp.v[12]-segLen*gap, tmp.v[13]);
            tmp.v[14] = MAX(tmp.v[13]-segLen*gap, tmp.v[14]);
            tmp.v[15] = MAX(tmp.v[14]-segLen*gap, tmp.v[15]);
            vFt = tmp.m;
        }
#else
        {
            __m256i vFt_save = vFt;
            __m256i segLenXgap = segLenXgap_reset;
            for (i=0; i<segWidth-1; ++i) {
                __m256i vFtt = shift(vFt);
                segLenXgap = lrotate16(segLenXgap);
                vFtt = _mm256_add_epi16(vFtt, segLenXgap);
                vFt = _mm256_max_epi16(vFt, vFtt);
            }
            vFt = _mm256_blendv_epi8(vFt_save, vFt, insert_mask);
        }
#endif
        vHt = shift(_mm256_load_si256(pvHt+(segLen-1)));
        vHt = _mm256_insert_epi16(vHt, boundary[j+1], 0);
        vFt = shift(vFt);
        vFt = _mm256_insert_epi16(vFt, NEG_INF_16, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm256_sub_epi16(vFt, vGapE),
            vFt = _mm256_max_epi16(vFt, vHt);
            vHt = _mm256_load_si256(pvHt+i);
            _mm256_store_si256(pvFt+i, vFt);
        }

        /* calculate H,M,L */
        vMp = vZero;
        vLp = vOne;
        vC = _mm256_cmpeq_epi16(vZero, vZero); /* check if prefix sum is needed */
        vC = rshift(vC); /* zero out last value */
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vHt = _mm256_load_si256(pvHt+i);
            vFt = _mm256_load_si256(pvFt+i);
            /* compute */
            vFt = _mm256_sub_epi16(vFt, vGapO);
            vH = _mm256_max_epi16(vHt, vFt);
            /* statistics */
            vEx = _mm256_load_si256(pvEx+i);
            vMt = _mm256_load_si256(pvMt+i);
            vLt = _mm256_load_si256(pvLt+i);
            vEx = _mm256_or_si256(
                    _mm256_and_si256(vEx, _mm256_cmpeq_epi16(vHt, vFt)),
                    _mm256_cmplt_epi16(vHt, vFt));
            vM = _mm256_blendv_epi8(vMt, vMp, vEx);
            vL = _mm256_blendv_epi8(vLt, vLp, vEx);
            vMp = vM;
            vLp = _mm256_add_epi16(vL, vOne);
            vC = _mm256_and_si256(vC, vEx);
            /* store results */
            _mm256_store_si256(pvH+i, vH);
            _mm256_store_si256(pvEx+i, vEx);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif
        }
        {
            vLp = _mm256_sub_epi16(vLp, vOne);
            {
                __m256i_16_t uMp, uLp, uC;
                uC.m = vC;
                uMp.m = vMp;
                uMp.v[ 1] = uC.v[ 1] ? uMp.v[ 0] : uMp.v[ 1];
                uMp.v[ 2] = uC.v[ 2] ? uMp.v[ 1] : uMp.v[ 2];
                uMp.v[ 3] = uC.v[ 3] ? uMp.v[ 2] : uMp.v[ 3];
                uMp.v[ 4] = uC.v[ 4] ? uMp.v[ 3] : uMp.v[ 4];
                uMp.v[ 5] = uC.v[ 5] ? uMp.v[ 4] : uMp.v[ 5];
                uMp.v[ 6] = uC.v[ 6] ? uMp.v[ 5] : uMp.v[ 6];
                uMp.v[ 7] = uC.v[ 7] ? uMp.v[ 6] : uMp.v[ 7];
                uMp.v[ 8] = uC.v[ 8] ? uMp.v[ 7] : uMp.v[ 8];
                uMp.v[ 9] = uC.v[ 9] ? uMp.v[ 8] : uMp.v[ 9];
                uMp.v[10] = uC.v[10] ? uMp.v[ 9] : uMp.v[10];
                uMp.v[11] = uC.v[11] ? uMp.v[10] : uMp.v[11];
                uMp.v[12] = uC.v[12] ? uMp.v[11] : uMp.v[12];
                uMp.v[13] = uC.v[13] ? uMp.v[12] : uMp.v[13];
                uMp.v[14] = uC.v[14] ? uMp.v[13] : uMp.v[14];
                uMp.v[15] = uC.v[15] ? uMp.v[14] : uMp.v[15];
                vMp = uMp.m;
                uLp.m = vLp;
                uLp.v[ 1] = uC.v[ 1] ? uLp.v[ 1] + uLp.v[ 0] : uLp.v[ 1];
                uLp.v[ 2] = uC.v[ 2] ? uLp.v[ 2] + uLp.v[ 1] : uLp.v[ 2];
                uLp.v[ 3] = uC.v[ 3] ? uLp.v[ 3] + uLp.v[ 2] : uLp.v[ 3];
                uLp.v[ 4] = uC.v[ 4] ? uLp.v[ 4] + uLp.v[ 3] : uLp.v[ 4];
                uLp.v[ 5] = uC.v[ 5] ? uLp.v[ 5] + uLp.v[ 4] : uLp.v[ 5];
                uLp.v[ 6] = uC.v[ 6] ? uLp.v[ 6] + uLp.v[ 5] : uLp.v[ 6];
                uLp.v[ 7] = uC.v[ 7] ? uLp.v[ 7] + uLp.v[ 6] : uLp.v[ 7];
                uLp.v[ 8] = uC.v[ 8] ? uLp.v[ 8] + uLp.v[ 7] : uLp.v[ 8];
                uLp.v[ 9] = uC.v[ 9] ? uLp.v[ 9] + uLp.v[ 8] : uLp.v[ 9];
                uLp.v[10] = uC.v[10] ? uLp.v[10] + uLp.v[ 9] : uLp.v[10];
                uLp.v[11] = uC.v[11] ? uLp.v[11] + uLp.v[10] : uLp.v[11];
                uLp.v[12] = uC.v[12] ? uLp.v[12] + uLp.v[11] : uLp.v[12];
                uLp.v[13] = uC.v[13] ? uLp.v[13] + uLp.v[12] : uLp.v[13];
                uLp.v[14] = uC.v[14] ? uLp.v[14] + uLp.v[13] : uLp.v[14];
                uLp.v[15] = uC.v[15] ? uLp.v[15] + uLp.v[14] : uLp.v[15];
                vLp = uLp.m;
            }
            vLp = _mm256_add_epi16(vLp, vOne);
        }
        /* final pass for M,L */
        vMp = shift(vMp);
        vLp = shift(vLp);
        for (i=0; i<segLen; ++i) {
            /* statistics */
            vEx = _mm256_load_si256(pvEx+i);
            vMt = _mm256_load_si256(pvMt+i);
            vLt = _mm256_load_si256(pvLt+i);
            vM = _mm256_blendv_epi8(vMt, vMp, vEx);
            vL = _mm256_blendv_epi8(vLt, vLp, vEx);
            vMp = vM;
            vLp = _mm256_add_epi16(vL, vOne);
            /* store results */
            _mm256_store_si256(pvM+i, vM);
            _mm256_store_si256(pvL+i, vL);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->matches_table, vM, i, segLen, j, s2Len);
            arr_store_si256(result->length_table, vL, i, segLen, j, s2Len);
#endif
        }
    }

    /* extract last value from the last column */
    {
        __m256i vH = _mm256_load_si256(pvH + offset);
        __m256i vM = _mm256_load_si256(pvM + offset);
        __m256i vL = _mm256_load_si256(pvL + offset);
        for (k=0; k<position; ++k) {
            vH = shift(vH);
            vM = shift(vM);
            vL = shift(vL);
        }
        score = (int16_t) _mm256_extract_epi16 (vH, 15);
        matches = (int16_t) _mm256_extract_epi16 (vM, 15);
        length = (int16_t) _mm256_extract_epi16 (vL, 15);
    }

    result->score = score;
    result->matches = matches;
    result->length = length;

    parasail_free(boundary);
    parasail_free(pvL);
    parasail_free(pvM);
    parasail_free(pvH);
    parasail_free(pvEx);
    parasail_free(pvLt);
    parasail_free(pvMt);
    parasail_free(pvFt);
    parasail_free(pvHt);
    parasail_free(pvE);
    parasail_free(pvPm);
    parasail_free(pvP);

    return result;
}

