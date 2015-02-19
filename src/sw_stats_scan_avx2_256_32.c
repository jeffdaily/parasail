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

/* avx2 does not have _mm256_cmplt_epi32, emulate it */
static inline __m256i _mm256_cmplt_epi32(__m256i a, __m256i b) {
    return _mm256_cmpgt_epi32(b,a);
}

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

/* avx2 _mm256_srli_si256 does not shift across 128-bit lanes, emulate it */
static inline __m256i rshift(__m256i a) {
    return _mm256_or_si256(
            _mm256_slli_si256(
                _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(3,0,0,1)),
                12),
            _mm256_srli_si256(a, 4));
}

static inline __m256i lrotate32(__m256i a) {
    return _mm256_alignr_epi8(a,
            _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,1)),
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
#define FNAME sw_stats_table_scan_avx2_256_32
#else
#define FNAME sw_stats_scan_avx2_256_32
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
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m256i* const restrict pvP  = parasail_memalign_m256i(32, n * segLen);
    __m256i* const restrict pvPm = parasail_memalign_m256i(32, n * segLen);
    __m256i* const restrict pvPs = parasail_memalign_m256i(32, n * segLen);
    __m256i* const restrict pvE  = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvHt = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvFt = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvMt = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvSt = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvLt = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvEx = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvH  = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvM  = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvS  = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvL  = parasail_memalign_m256i(32, segLen);
    __m256i vGapO = _mm256_set1_epi32(open);
    __m256i vGapE = _mm256_set1_epi32(gap);
    __m256i vZero = _mm256_setzero_si256();
    __m256i vOne = _mm256_set1_epi32(1);
    int32_t score = 0;
    int32_t matches = 0;
    int32_t similar = 0;
    int32_t length = 0;
    __m256i vMaxH = vZero;
    __m256i vMaxM = vZero;
    __m256i vMaxS = vZero;
    __m256i vMaxL = vZero;
    __m256i segLenXgap_reset = _mm256_set_epi32(
            NEG_INF_32, NEG_INF_32, NEG_INF_32, NEG_INF_32,
            NEG_INF_32, NEG_INF_32, NEG_INF_32, -segLen*gap);
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif

    parasail_memset_m256i(pvM, vZero, segLen);
    parasail_memset_m256i(pvS, vZero, segLen);
    parasail_memset_m256i(pvL, vZero, segLen);

    /* Generate query profile and match profile.
     * Rearrange query sequence & calculate the weight of match/mismatch.
     * Don't alias. */
    {
        int32_t index = 0;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                __m256i_32_t p;
                __m256i_32_t m;
                __m256i_32_t s;
                j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    p.v[segNum] = j >= s1Len ? 0 : matrix[k][MAP_BLOSUM_[(unsigned char)s1[j]]];
                    m.v[segNum] = j >= s1Len ? 0 : (k == MAP_BLOSUM_[(unsigned char)s1[j]]);
                    s.v[segNum] = p.v[segNum] > 0;
                    j += segLen;
                }
                _mm256_store_si256(&pvP[index], p.m);
                _mm256_store_si256(&pvPm[index], m.m);
                _mm256_store_si256(&pvPs[index], s.m);
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
        __m256i *pvD;
        __m256i vC;
        __m256i vD;
        __m256i vM;
        __m256i vMp;
        __m256i vMt;
        __m256i vS;
        __m256i vSp;
        __m256i vSt;
        __m256i vL;
        __m256i vLp;
        __m256i vLt;
        __m256i vEx;

        /* calculate E */
        for (i=0; i<segLen; ++i) {
            vH = _mm256_load_si256(pvH+i);
            vE = _mm256_load_si256(pvE+i);
            vE = _mm256_max_epi32(
                    _mm256_sub_epi32(vE, vGapE),
                    _mm256_sub_epi32(vH, vGapO));
            _mm256_store_si256(pvE+i, vE);
        }

        /* calculate Ht */
        vH = shift(_mm256_load_si256(pvH+(segLen-1)));
        vMp= shift(_mm256_load_si256(pvM+(segLen-1)));
        vSp= shift(_mm256_load_si256(pvS+(segLen-1)));
        vLp= shift(_mm256_load_si256(pvL+(segLen-1)));
        vLp= _mm256_add_epi32(vLp, vOne);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        pvC = pvPm+ MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        pvD = pvPs+ MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            __m256i cond_max;
            /* load values we need */
            vE = _mm256_load_si256(pvE+i);
            vW = _mm256_load_si256(pvW+i);
            /* compute */
            vH = _mm256_add_epi32(vH, vW);
            vHt = _mm256_max_epi32(vH, vE);
            vHt = _mm256_max_epi32(vHt, vZero);
            /* statistics */
            vC = _mm256_load_si256(pvC+i);
            vD = _mm256_load_si256(pvD+i);
            vMp = _mm256_add_epi32(vMp, vC);
            vSp = _mm256_add_epi32(vSp, vD);
            vEx = _mm256_cmpgt_epi32(vE, vH);
            vM = _mm256_load_si256(pvM+i);
            vS = _mm256_load_si256(pvS+i);
            vL = _mm256_load_si256(pvL+i);
            vL = _mm256_add_epi32(vL, vOne);
            vMt = _mm256_blendv_epi8(vMp, vM, vEx);
            vSt = _mm256_blendv_epi8(vSp, vS, vEx);
            vLt = _mm256_blendv_epi8(vLp, vL, vEx);
            cond_max = _mm256_cmpeq_epi32(vHt, vZero);
            vEx = _mm256_andnot_si256(cond_max, vEx);
            vMt = _mm256_andnot_si256(cond_max, vMt);
            vSt = _mm256_andnot_si256(cond_max, vSt);
            vLt = _mm256_andnot_si256(cond_max, vLt);
            /* store results */
            _mm256_store_si256(pvHt+i, vHt);
            _mm256_store_si256(pvEx+i, vEx);
            _mm256_store_si256(pvMt+i, vMt);
            _mm256_store_si256(pvSt+i, vSt);
            _mm256_store_si256(pvLt+i, vLt);
            /* prep for next iteration */
            vH = _mm256_load_si256(pvH+i);
            vMp = vM;
            vSp = vS;
            vLp = vL;
        }

        /* calculate Ft */
        vHt = shift(_mm256_load_si256(pvHt+(segLen-1)));
        vFt = _mm256_set1_epi32(NEG_INF_32);
        for (i=0; i<segLen; ++i) {
            vFt = _mm256_sub_epi32(vFt, vGapE);
            vFt = _mm256_max_epi32(vFt, vHt);
            vHt = _mm256_load_si256(pvHt+i);
        }
#if 0
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
#else
        {
            __m256i vFt_save = vFt;
            __m256i segLenXgap = segLenXgap_reset;
            for (i=0; i<segWidth-1; ++i) {
                __m256i vFtt = shift(vFt);
                segLenXgap = lrotate32(segLenXgap);
                vFtt = _mm256_add_epi32(vFtt, segLenXgap);
                vFt = _mm256_max_epi32(vFt, vFtt);
            }
            vFt = _mm256_blend_epi32(vFt_save, vFt, 0xFE);
        }
#endif
            vHt = shift(_mm256_load_si256(pvHt+(segLen-1)));
        vFt = shift(vFt);
        vFt = _mm256_insert_epi32(vFt, NEG_INF_32, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm256_sub_epi32(vFt, vGapE);
            vFt = _mm256_max_epi32(vFt, vHt);
            vHt = _mm256_load_si256(pvHt+i);
            _mm256_store_si256(pvFt+i, vFt);
        }

        /* calculate H,M,L */
        vMp = vZero;
        vSp = vZero;
        vLp = vOne;
        vC = _mm256_cmpeq_epi32(vZero, vZero); /* check if prefix sum is needed */
        vC = rshift(vC); /* zero out last value */
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vHt = _mm256_load_si256(pvHt+i);
            vFt = _mm256_load_si256(pvFt+i);
            /* compute */
            vFt = _mm256_sub_epi32(vFt, vGapO);
            vH = _mm256_max_epi32(vHt, vFt);
            /* statistics */
            vEx = _mm256_load_si256(pvEx+i);
            vMt = _mm256_load_si256(pvMt+i);
            vSt = _mm256_load_si256(pvSt+i);
            vLt = _mm256_load_si256(pvLt+i);
            vEx = _mm256_or_si256(
                    _mm256_and_si256(vEx, _mm256_cmpeq_epi32(vHt, vFt)),
                    _mm256_cmplt_epi32(vHt, vFt));
            vM = _mm256_blendv_epi8(vMt, vMp, vEx);
            vS = _mm256_blendv_epi8(vSt, vSp, vEx);
            vL = _mm256_blendv_epi8(vLt, vLp, vEx);
            vMp = vM;
            vSp = vS;
            vLp = _mm256_add_epi32(vL, vOne);
            vC = _mm256_and_si256(vC, vEx);
            /* store results */
            _mm256_store_si256(pvH+i, vH);
            _mm256_store_si256(pvEx+i, vEx);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif
        }
        {
            vLp = _mm256_sub_epi32(vLp, vOne);
            {
                __m256i_32_t uMp, uSp, uLp, uC;
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
                uSp.m = vSp;
                uSp.v[1] = uC.v[1] ? uSp.v[0] : uSp.v[1];
                uSp.v[2] = uC.v[2] ? uSp.v[1] : uSp.v[2];
                uSp.v[3] = uC.v[3] ? uSp.v[2] : uSp.v[3];
                uSp.v[4] = uC.v[4] ? uSp.v[3] : uSp.v[4];
                uSp.v[5] = uC.v[5] ? uSp.v[4] : uSp.v[5];
                uSp.v[6] = uC.v[6] ? uSp.v[5] : uSp.v[6];
                uSp.v[7] = uC.v[7] ? uSp.v[6] : uSp.v[7];
                vSp = uSp.m;
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
            vLp = _mm256_add_epi32(vLp, vOne);
        }
        /* final pass for M,L */
        vMp = shift(vMp);
        vSp = shift(vSp);
        vLp = shift(vLp);
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vH = _mm256_load_si256(pvH+i);
            /* statistics */
            vEx = _mm256_load_si256(pvEx+i);
            vMt = _mm256_load_si256(pvMt+i);
            vSt = _mm256_load_si256(pvSt+i);
            vLt = _mm256_load_si256(pvLt+i);
            vM = _mm256_blendv_epi8(vMt, vMp, vEx);
            vS = _mm256_blendv_epi8(vSt, vSp, vEx);
            vL = _mm256_blendv_epi8(vLt, vLp, vEx);
            vMp = vM;
            vSp = vS;
            vLp = _mm256_add_epi32(vL, vOne);
            /* store results */
            _mm256_store_si256(pvM+i, vM);
            _mm256_store_si256(pvS+i, vS);
            _mm256_store_si256(pvL+i, vL);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->matches_table, vM, i, segLen, j, s2Len);
            arr_store_si256(result->similar_table, vS, i, segLen, j, s2Len);
            arr_store_si256(result->length_table, vL, i, segLen, j, s2Len);
#endif
            /* update max vector seen so far */
            {
                __m256i cond_max = _mm256_cmpgt_epi32(vH, vMaxH);
                vMaxH = _mm256_blendv_epi8(vMaxH, vH, cond_max);
                vMaxM = _mm256_blendv_epi8(vMaxM, vM, cond_max);
                vMaxS = _mm256_blendv_epi8(vMaxS, vS, cond_max);
                vMaxL = _mm256_blendv_epi8(vMaxL, vL, cond_max);
            }
        }
    }

    /* max in vec */
    for (j=0; j<segWidth; ++j) {
        int32_t value = (int32_t) _mm256_extract_epi32(vMaxH, 7);
        if (value > score) {
            matches = (int32_t) _mm256_extract_epi32(vMaxM, 7);
            similar = (int32_t) _mm256_extract_epi32(vMaxS, 7);
            length = (int32_t) _mm256_extract_epi32(vMaxL, 7);
            score = value;
        }
        vMaxH = shift(vMaxH);
        vMaxM = shift(vMaxM);
        vMaxS = shift(vMaxS);
        vMaxL = shift(vMaxL);
    }

    result->score = score;
    result->matches = matches;
    result->similar = similar;
    result->length = length;

    parasail_free(pvL);
    parasail_free(pvS);
    parasail_free(pvM);
    parasail_free(pvH);
    parasail_free(pvEx);
    parasail_free(pvLt);
    parasail_free(pvSt);
    parasail_free(pvMt);
    parasail_free(pvFt);
    parasail_free(pvHt);
    parasail_free(pvE);
    parasail_free(pvPs);
    parasail_free(pvPm);
    parasail_free(pvP);

    return result;
}
