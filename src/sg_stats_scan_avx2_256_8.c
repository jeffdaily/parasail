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

/* avx2 _mm256_packs_epi16 does not pack across 128-bit lanes, emulate it */
static inline __m256i pack16(__m256i a, __m256i b) {
    return _mm256_permute4x64_epi64(
            _mm256_packs_epi16(a, b),
            _MM_SHUFFLE(3,1,2,0));
}

/* avx2 does not have _mm256_cmplt_epi8, emulate it */
static inline __m256i _mm256_cmplt_epi8(__m256i a, __m256i b) {
    return _mm256_cmpgt_epi8(b,a);
}

/* avx2 does not have _mm256_cmplt_epi16, emulate it */
static inline __m256i _mm256_cmplt_epi16(__m256i a, __m256i b) {
    return _mm256_cmpgt_epi16(b,a);
}

#if HAVE_AVX2_MM256_INSERT_EPI8
#else
static inline __m256i _mm256_insert_epi8(__m256i a, int8_t b, int imm) {
    __m256i_8_t tmp;
    tmp.m = a;
    tmp.v[imm] = b;
    return tmp.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI8
#else
static inline int8_t _mm256_extract_epi8(__m256i a, int imm) {
    __m256i_8_t tmp;
    tmp.m = a;
    return tmp.v[imm];
}
#endif

/* avx2 _mm256_slli_si256 does not shift across 128-bit lanes, emulate it */
static inline __m256i shift(__m256i a) {
    return _mm256_alignr_epi8(a,
            _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)),
            15);
}

/* avx2 _mm256_srli_si256 does not shift across 128-bit lanes, emulate it */
static inline __m256i rshift(__m256i a) {
    return _mm256_or_si256(
            _mm256_slli_si256(
                _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(3,0,0,1)),
                15),
            _mm256_srli_si256(a, 1));
}

static inline __m256i lrotate8(__m256i a) {
    return _mm256_alignr_epi8(a,
            _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,1)),
            15);
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
    array[( 0*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 0);
    array[( 1*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 1);
    array[( 2*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 2);
    array[( 3*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 3);
    array[( 4*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 4);
    array[( 5*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 5);
    array[( 6*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 6);
    array[( 7*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 7);
    array[( 8*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 8);
    array[( 9*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 9);
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

static inline void arr_store_si256_bias(
        int *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen,
        int32_t bias)
{
    array[( 0*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 0);
    array[( 1*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 1);
    array[( 2*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 2);
    array[( 3*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 3);
    array[( 4*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 4);
    array[( 5*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 5);
    array[( 6*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 6);
    array[( 7*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 7);
    array[( 8*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 8);
    array[( 9*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 9);
    array[(10*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 10);
    array[(11*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 11);
    array[(12*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 12);
    array[(13*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 13);
    array[(14*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 14);
    array[(15*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 15);
    array[(16*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 16);
    array[(17*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 17);
    array[(18*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 18);
    array[(19*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 19);
    array[(20*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 20);
    array[(21*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 21);
    array[(22*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 22);
    array[(23*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 23);
    array[(24*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 24);
    array[(25*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 25);
    array[(26*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 26);
    array[(27*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 27);
    array[(28*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 28);
    array[(29*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 29);
    array[(30*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 30);
    array[(31*seglen+t)*dlen + d] = bias + (int8_t)_mm256_extract_epi8(vH, 31);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME sg_stats_table_scan_avx2_256_8
#else
#define FNAME sg_stats_scan_avx2_256_8
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
    __m256i vGapO = _mm256_set1_epi8(open);
    __m256i vGapE = _mm256_set1_epi8(gap);
    __m256i vZero = _mm256_setzero_si256();
    __m256i vOne = _mm256_set1_epi8(1);
    int8_t bias = 127;
    int ibias = 127;
    __m256i vNegBias = _mm256_set1_epi8(-bias);
    __m256i vNegInf = _mm256_set1_epi8(NEG_INF_8);
    __m256i vSaturationCheck = _mm256_setzero_si256();
    __m256i vNegLimit = _mm256_set1_epi8(INT8_MIN);
    __m256i vPosLimit = _mm256_set1_epi8(INT8_MAX);
    int8_t score = NEG_INF_8;
    int8_t matches = 0;
    int8_t length = 0;
    __m256i vMaxH = vNegInf;
    __m256i vMaxM = vZero;
    __m256i vMaxL = vZero;
    __m256i segLenXgap_reset = _mm256_set_epi8(
            NEG_INF_8, NEG_INF_8, NEG_INF_8, NEG_INF_8,
            NEG_INF_8, NEG_INF_8, NEG_INF_8, NEG_INF_8,
            NEG_INF_8, NEG_INF_8, NEG_INF_8, NEG_INF_8,
            NEG_INF_8, NEG_INF_8, NEG_INF_8, NEG_INF_8,
            NEG_INF_8, NEG_INF_8, NEG_INF_8, NEG_INF_8,
            NEG_INF_8, NEG_INF_8, NEG_INF_8, NEG_INF_8,
            NEG_INF_8, NEG_INF_8, NEG_INF_8, NEG_INF_8,
            NEG_INF_8, NEG_INF_8, NEG_INF_8, -segLen*gap);
    __m256i insert_mask = _mm256_cmpeq_epi8(_mm256_setzero_si256(),
            _mm256_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1));
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif

    parasail_memset_m256i(pvM, vNegBias, segLen);
    parasail_memset_m256i(pvL, vNegBias, segLen);

    /* Generate query profile and match profile.
     * Rearrange query sequence & calculate the weight of match/mismatch.
     * Don't alias. */
    {
        int32_t index = 0;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                __m256i_8_t t;
                __m256i_8_t s;
                j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = j >= s1Len ? 0 : matrix[k][MAP_BLOSUM_[(unsigned char)s1[j]]];
                    s.v[segNum] = (k == MAP_BLOSUM_[(unsigned char)s1[j]]);
                    j += segLen;
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
            __m256i_8_t h;
            __m256i_8_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = 0;
                e.v[segNum] = NEG_INF_8;
            }
            _mm256_store_si256(&pvH[index], h.m);
            _mm256_store_si256(&pvE[index], e.m);
            ++index;
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
            vE = _mm256_max_epi8(
                    _mm256_subs_epi8(vE, vGapE),
                    _mm256_subs_epi8(vH, vGapO));
            _mm256_store_si256(pvE+i, vE);
        }

        /* calculate Ht */
        vH = shift(_mm256_load_si256(pvH+(segLen-1)));
        vMp= shift(_mm256_load_si256(pvM+(segLen-1)));
        vMp= _mm256_insert_epi8(vMp, -bias, 0);
        vLp= shift(_mm256_load_si256(pvL+(segLen-1)));
        vLp= _mm256_insert_epi8(vLp, -bias, 0);
        vLp= _mm256_adds_epi8(vLp, vOne);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        pvC = pvPm+ MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vE = _mm256_load_si256(pvE+i);
            vW = _mm256_load_si256(pvW+i);
            /* compute */
            vH = _mm256_adds_epi8(vH, vW);
            vHt = _mm256_max_epi8(vH, vE);
            /* statistics */
            vC = _mm256_load_si256(pvC+i);
            vMp = _mm256_adds_epi8(vMp, vC);
            vEx = _mm256_cmpgt_epi8(vE, vH);
            vM = _mm256_load_si256(pvM+i);
            vL = _mm256_load_si256(pvL+i);
            vL = _mm256_adds_epi8(vL, vOne);
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
        vFt = vNegInf;
        for (i=0; i<segLen; ++i) {
            vFt = _mm256_subs_epi8(vFt, vGapE);
            vFt = _mm256_max_epi8(vFt, vHt);
            vHt = _mm256_load_si256(pvHt+i);
        }
#if 0
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
#else
        {
            __m256i vFt_save = vFt;
            __m256i segLenXgap = segLenXgap_reset;
            for (i=0; i<segWidth-1; ++i) {
                __m256i vFtt = shift(vFt);
                segLenXgap = lrotate8(segLenXgap);
                vFtt = _mm256_adds_epi8(vFtt, segLenXgap);
                vFt = _mm256_max_epi8(vFtt, vFt);
            }
            vFt = _mm256_blendv_epi8(vFt_save, vFt, insert_mask);
        }
#endif
        vHt = shift(_mm256_load_si256(pvHt+(segLen-1)));
        vFt = shift(vFt);
        vFt = _mm256_insert_epi8(vFt, NEG_INF_8, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm256_subs_epi8(vFt, vGapE);
            vFt = _mm256_max_epi8(vFt, vHt);
            vHt = _mm256_load_si256(pvHt+i);
            _mm256_store_si256(pvFt+i, vFt);
        }

        /* calculate H,M,L */
        vMp = vNegBias;
        vLp = _mm256_adds_epi8(vNegBias, vOne);
        vC = _mm256_cmpeq_epi8(vZero, vZero); /* check if prefix sum is needed */
        vC = rshift(vC); /* zero out last value */
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vHt = _mm256_load_si256(pvHt+i);
            vFt = _mm256_load_si256(pvFt+i);
            /* compute */
            vFt = _mm256_subs_epi8(vFt, vGapO);
            vH = _mm256_max_epi8(vHt, vFt);
            /* statistics */
            vEx = _mm256_load_si256(pvEx+i);
            vMt = _mm256_load_si256(pvMt+i);
            vLt = _mm256_load_si256(pvLt+i);
            vEx = _mm256_or_si256(
                    _mm256_and_si256(vEx, _mm256_cmpeq_epi8(vHt, vFt)),
                    _mm256_cmplt_epi8(vHt, vFt));
            vM = _mm256_blendv_epi8(vMt, vMp, vEx);
            vL = _mm256_blendv_epi8(vLt, vLp, vEx);
            vMp = vM;
            vLp = _mm256_adds_epi8(vL, vOne);
            vC = _mm256_and_si256(vC, vEx);
            /* store results */
            _mm256_store_si256(pvH+i, vH);
            _mm256_store_si256(pvEx+i, vEx);
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
        {
            vLp = _mm256_subs_epi8(vLp, vOne);
            {
                __m256i_8_t uMp, uLp, uC;
                uC.m = vC;
                uMp.m = vMp;
                uMp.v[1] = uC.v[1] ? uMp.v[0] : uMp.v[1];
                uMp.v[2] = uC.v[2] ? uMp.v[1] : uMp.v[2];
                uMp.v[3] = uC.v[3] ? uMp.v[2] : uMp.v[3];
                uMp.v[4] = uC.v[4] ? uMp.v[3] : uMp.v[4];
                uMp.v[5] = uC.v[5] ? uMp.v[4] : uMp.v[5];
                uMp.v[6] = uC.v[6] ? uMp.v[5] : uMp.v[6];
                uMp.v[7] = uC.v[7] ? uMp.v[6] : uMp.v[7];
                vMp = uMp.m;
                uLp.m = vLp;
                uLp.v[1] = uC.v[1] ? uLp.v[1] + uLp.v[0] : uLp.v[1];
                uLp.v[2] = uC.v[2] ? uLp.v[2] + uLp.v[1] : uLp.v[2];
                uLp.v[3] = uC.v[3] ? uLp.v[3] + uLp.v[2] : uLp.v[3];
                uLp.v[4] = uC.v[4] ? uLp.v[4] + uLp.v[3] : uLp.v[4];
                uLp.v[5] = uC.v[5] ? uLp.v[5] + uLp.v[4] : uLp.v[5];
                uLp.v[6] = uC.v[6] ? uLp.v[6] + uLp.v[5] : uLp.v[6];
                uLp.v[7] = uC.v[7] ? uLp.v[7] + uLp.v[6] : uLp.v[7];
                vLp = uLp.m;
            }
            vLp = _mm256_adds_epi8(vLp, vOne);
        }
        /* final pass for M,L */
        vMp = shift(vMp);
        vMp = _mm256_insert_epi8(vMp, -bias, 0);
        vLp = shift(vLp);
        vLp = _mm256_insert_epi8(vLp, -bias, 0);
        for (i=0; i<segLen; ++i) {
            /* statistics */
            vEx = _mm256_load_si256(pvEx+i);
            vMt = _mm256_load_si256(pvMt+i);
            vLt = _mm256_load_si256(pvLt+i);
            vM = _mm256_blendv_epi8(vMt, vMp, vEx);
            vL = _mm256_blendv_epi8(vLt, vLp, vEx);
            vMp = vM;
            vLp = _mm256_adds_epi8(vL, vOne);
            /* store results */
            _mm256_store_si256(pvM+i, vM);
            _mm256_store_si256(pvL+i, vL);
            /* check for saturation */
            {
                vSaturationCheck = _mm256_or_si256(vSaturationCheck,
                        _mm256_or_si256(
                            _mm256_cmpeq_epi8(vM, vPosLimit),
                            _mm256_cmpeq_epi8(vL, vPosLimit)));
            }
#ifdef PARASAIL_TABLE
            arr_store_si256_bias(result->matches_table, vM, i, segLen, j, s2Len, bias);
            arr_store_si256_bias(result->length_table, vL, i, segLen, j, s2Len, bias);
#endif
        }

        /* extract vector containing last value from column */
        {
            __m256i cond_max;
            vH = _mm256_load_si256(pvH + offset);
            vM = _mm256_load_si256(pvM + offset);
            vL = _mm256_load_si256(pvL + offset);
            cond_max = _mm256_cmpgt_epi8(vH, vMaxH);
            vMaxH = _mm256_blendv_epi8(vMaxH, vH, cond_max);
            vMaxM = _mm256_blendv_epi8(vMaxM, vM, cond_max);
            vMaxL = _mm256_blendv_epi8(vMaxL, vL, cond_max);
        }
    }

    /* max last value from all columns */
    {
        int8_t value;
        for (k=0; k<position; ++k) {
            vMaxH = shift(vMaxH);
            vMaxM = shift(vMaxM);
            vMaxL = shift(vMaxL);
        }
        value = (int8_t) _mm256_extract_epi8(vMaxH, 31);
        if (value > score) {
            score = value;
            matches = ibias + (int)((int8_t)_mm256_extract_epi8(vMaxM, 31));
            length = ibias + (int)((int8_t)_mm256_extract_epi8(vMaxL, 31));
        }
    }

    /* max of last column */
    {
        __m256i vQLimit16 = _mm256_set1_epi16(s1Len);
        __m256i vQIndexHi16 = _mm256_set_epi16(
                31*segLen,
                30*segLen,
                29*segLen,
                28*segLen,
                27*segLen,
                26*segLen,
                25*segLen,
                24*segLen,
                23*segLen,
                22*segLen,
                21*segLen,
                20*segLen,
                19*segLen,
                18*segLen,
                17*segLen,
                16*segLen);
        __m256i vQIndexLo16 = _mm256_set_epi16(
                15*segLen,
                14*segLen,
                13*segLen,
                12*segLen,
                11*segLen,
                10*segLen,
                9*segLen,
                8*segLen,
                7*segLen,
                6*segLen,
                5*segLen,
                4*segLen,
                3*segLen,
                2*segLen,
                1*segLen,
                0*segLen);
        __m256i vOne16 = _mm256_set1_epi16(1);
        vMaxH = vNegInf;
        vMaxM = vZero;
        vMaxL = vZero;

        for (i=0; i<segLen; ++i) {
            __m256i vH = _mm256_load_si256(pvH + i);
            __m256i vM = _mm256_load_si256(pvM + i);
            __m256i vL = _mm256_load_si256(pvL + i);
            __m256i cond_lmt = pack16(
                _mm256_cmplt_epi16(vQIndexLo16, vQLimit16),
                _mm256_cmplt_epi16(vQIndexHi16, vQLimit16));
            __m256i cond_max = _mm256_cmpgt_epi8(vH, vMaxH);
            __m256i cond_all = _mm256_and_si256(cond_max, cond_lmt);
            vMaxH = _mm256_blendv_epi8(vMaxH, vH, cond_all);
            vMaxM = _mm256_blendv_epi8(vMaxM, vM, cond_all);
            vMaxL = _mm256_blendv_epi8(vMaxL, vL, cond_all);
            vQIndexLo16 = _mm256_adds_epi8(vQIndexLo16, vOne16);
            vQIndexHi16 = _mm256_adds_epi8(vQIndexHi16, vOne16);
        }

        /* max in vec */
        for (j=0; j<segWidth; ++j) {
            int8_t value = (int8_t) _mm256_extract_epi8(vMaxH, 31);
            if (value > score) {
                score = value;
                matches = ibias + (int)((int8_t)_mm256_extract_epi8(vMaxM, 31));
                length = ibias + (int)((int8_t)_mm256_extract_epi8(vMaxL, 31));
            }
            vMaxH = shift(vMaxH);
            vMaxM = shift(vMaxM);
            vMaxL = shift(vMaxL);
        }
    }

    if (_mm256_movemask_epi8(vSaturationCheck)) {
        result->saturated = 1;
        score = INT8_MAX;
        matches = 0;
        length = 0;
    }

    result->score = score;
    result->matches = matches;
    result->length = length;

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

