/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include <emmintrin.h>

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_sse.h"

#define NEG_INF INT8_MIN
#define MAX(a,b) ((a)>(b)?(a):(b))

static inline __m128i _mm_blendv_epi8_rpl(__m128i a, __m128i b, __m128i mask) {
    a = _mm_andnot_si128(mask, a);
    a = _mm_or_si128(a, _mm_and_si128(mask, b));
    return a;
}

static inline __m128i _mm_insert_epi8_rpl(__m128i a, int8_t i, const int imm) {
    __m128i_8_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}

static inline __m128i _mm_max_epi8_rpl(__m128i a, __m128i b) {
    __m128i mask = _mm_cmpgt_epi8(a, b);
    a = _mm_and_si128(a, mask);
    b = _mm_andnot_si128(mask, b);
    return _mm_or_si128(a, b);
}

static inline int8_t _mm_extract_epi8_rpl(__m128i a, const int imm) {
    __m128i_8_t A;
    A.m = a;
    return A.v[imm];
}

#define _mm_rlli_si128_rpl(a,imm) _mm_or_si128(_mm_slli_si128(a,imm),_mm_srli_si128(a,16-imm))

static inline int8_t _mm_hmax_epi8_rpl(__m128i a) {
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 8));
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 4));
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 2));
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 1));
    return _mm_extract_epi8_rpl(a, 0);
}

static inline __m128i _mm_min_epi8_rpl(__m128i a, __m128i b) {
    __m128i mask = _mm_cmpgt_epi8(b, a);
    a = _mm_and_si128(a, mask);
    b = _mm_andnot_si128(mask, b);
    return _mm_or_si128(a, b);
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
    array[( 0*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  0);
    array[( 1*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  1);
    array[( 2*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  2);
    array[( 3*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  3);
    array[( 4*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  4);
    array[( 5*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  5);
    array[( 6*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  6);
    array[( 7*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  7);
    array[( 8*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  8);
    array[( 9*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  9);
    array[(10*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 10);
    array[(11*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 11);
    array[(12*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 12);
    array[(13*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 13);
    array[(14*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 14);
    array[(15*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 15);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m128i vH,
        int32_t t,
        int32_t seglen)
{
    col[ 0*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  0);
    col[ 1*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  1);
    col[ 2*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  2);
    col[ 3*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  3);
    col[ 4*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  4);
    col[ 5*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  5);
    col[ 6*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  6);
    col[ 7*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  7);
    col[ 8*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  8);
    col[ 9*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  9);
    col[10*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 10);
    col[11*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 11);
    col[12*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 12);
    col[13*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 13);
    col[14*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 14);
    col[15*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 15);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sg_stats_table_scan_sse2_128_8
#define PNAME parasail_sg_stats_table_scan_profile_sse2_128_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_stats_rowcol_scan_sse2_128_8
#define PNAME parasail_sg_stats_rowcol_scan_profile_sse2_128_8
#else
#define FNAME parasail_sg_stats_scan_sse2_128_8
#define PNAME parasail_sg_stats_scan_profile_sse2_128_8
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_sse_128_8(s1, s1Len, matrix);
    parasail_result_t *result = PNAME(profile, s2, s2Len, open, gap);
    parasail_profile_free(profile);
    return result;
}

parasail_result_t* PNAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int32_t segNum = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const int32_t segWidth = 16;
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m128i* const restrict pvP  = (__m128i*)profile->profile8.score;
    __m128i* const restrict pvPm = (__m128i*)profile->profile8.matches;
    __m128i* const restrict pvPs = (__m128i*)profile->profile8.similar;
    __m128i* const restrict pvE  = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHt = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvFt = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvMt = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvSt = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvLt = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvEx = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvH  = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvM  = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvS  = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvL  = parasail_memalign___m128i(16, segLen);
    __m128i vGapO = _mm_set1_epi8(open);
    __m128i vGapE = _mm_set1_epi8(gap);
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi8(1);
    __m128i vNegInf = _mm_set1_epi8(NEG_INF);
    int8_t score = NEG_INF;
    int8_t matches = 0;
    int8_t similar = 0;
    int8_t length = 0;
    __m128i vMaxH = vNegInf;
    __m128i vMaxM = vZero;
    __m128i vMaxS = vZero;
    __m128i vMaxL = vZero;
    __m128i vPosMask = _mm_cmpeq_epi8(_mm_set1_epi8(position),
            _mm_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15));
    const int8_t segLenXgap = -segLen*gap;
    __m128i insert_mask = _mm_cmpeq_epi8(_mm_setzero_si128(),
            _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1));
    __m128i vSegLenXgap_reset = _mm_blendv_epi8_rpl(vNegInf,
            _mm_set1_epi8(segLenXgap),
            insert_mask);
    __m128i vNegLimit = _mm_set1_epi8(INT8_MIN);
    __m128i vPosLimit = _mm_set1_epi8(INT8_MAX);
    __m128i vSaturationCheckMin = vPosLimit;
    __m128i vSaturationCheckMax = vNegLimit;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif

    parasail_memset___m128i(pvM, vZero, segLen);
    parasail_memset___m128i(pvS, vZero, segLen);
    parasail_memset___m128i(pvL, vZero, segLen);

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m128i_8_t h;
            __m128i_8_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = 0;
                e.v[segNum] = NEG_INF;
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
        __m128i *pvW;
        __m128i vW;
        __m128i *pvC;
        __m128i *pvD;
        __m128i vC;
        __m128i vD;
        __m128i vM;
        __m128i vMp;
        __m128i vMt;
        __m128i vS;
        __m128i vSp;
        __m128i vSt;
        __m128i vL;
        __m128i vLp;
        __m128i vLt;
        __m128i vEx;

        /* calculate E */
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vE = _mm_max_epi8_rpl(
                    _mm_subs_epi8(vE, vGapE),
                    _mm_subs_epi8(vH, vGapO));
            _mm_store_si128(pvE+i, vE);
        }

        /* calculate Ht */
        vH = _mm_slli_si128(_mm_load_si128(pvH+(segLen-1)), 1);
        vMp= _mm_slli_si128(_mm_load_si128(pvM+(segLen-1)), 1);
        vSp= _mm_slli_si128(_mm_load_si128(pvS+(segLen-1)), 1);
        vLp= _mm_slli_si128(_mm_load_si128(pvL+(segLen-1)), 1);
        vLp= _mm_adds_epi8(vLp, vOne);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvC = pvPm+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvD = pvPs+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            /* compute */
            vH = _mm_adds_epi8(vH, vW);
            vHt = _mm_max_epi8_rpl(vH, vE);
            /* statistics */
            vC = _mm_load_si128(pvC+i);
            vD = _mm_load_si128(pvD+i);
            vMp = _mm_adds_epi8(vMp, vC);
            vSp = _mm_adds_epi8(vSp, vD);
            vEx = _mm_cmpgt_epi8(vE, vH);
            vM = _mm_load_si128(pvM+i);
            vS = _mm_load_si128(pvS+i);
            vL = _mm_load_si128(pvL+i);
            vL = _mm_adds_epi8(vL, vOne);
            vMt = _mm_blendv_epi8_rpl(vMp, vM, vEx);
            vSt = _mm_blendv_epi8_rpl(vSp, vS, vEx);
            vLt = _mm_blendv_epi8_rpl(vLp, vL, vEx);
            /* store results */
            _mm_store_si128(pvHt+i, vHt);
            _mm_store_si128(pvEx+i, vEx);
            _mm_store_si128(pvMt+i, vMt);
            _mm_store_si128(pvSt+i, vSt);
            _mm_store_si128(pvLt+i, vLt);
            /* prep for next iteration */
            vH = _mm_load_si128(pvH+i);
            vMp = vM;
            vSp = vS;
            vLp = vL;
        }

        /* calculate Ft */
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 1);
        vFt = vNegInf;
        for (i=0; i<segLen; ++i) {
            vFt = _mm_subs_epi8(vFt, vGapE);
            vFt = _mm_max_epi8_rpl(vFt, vHt);
            vHt = _mm_load_si128(pvHt+i);
        }
        {
            __m128i vFt_save = vFt;
            __m128i segLenXgap = vSegLenXgap_reset;
            for (i=0; i<segWidth-1; ++i) {
                __m128i vFtt = _mm_slli_si128(vFt, 1);
                segLenXgap = _mm_rlli_si128_rpl(segLenXgap, 1);
                vFtt = _mm_adds_epi8(vFtt, segLenXgap);
                vFt = _mm_max_epi8_rpl(vFt, vFtt);
            }
            vFt = _mm_blendv_epi8_rpl(vFt_save, vFt, insert_mask);
        }
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 1);
        vFt = _mm_slli_si128(vFt, 1);
        vFt = _mm_insert_epi8_rpl(vFt, NEG_INF, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_subs_epi8(vFt, vGapE);
            vFt = _mm_max_epi8_rpl(vFt, vHt);
            vHt = _mm_load_si128(pvHt+i);
            _mm_store_si128(pvFt+i, vFt);
        }

        /* calculate H,M,L */
        vMp = vZero;
        vSp = vZero;
        vLp = vOne;
        vC = _mm_cmpeq_epi8(vZero, vZero); /* check if prefix sum is needed */
        vC = _mm_srli_si128(vC, 1); /* zero out last value */
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vHt = _mm_load_si128(pvHt+i);
            vFt = _mm_load_si128(pvFt+i);
            /* compute */
            vFt = _mm_subs_epi8(vFt, vGapO);
            vH = _mm_max_epi8_rpl(vHt, vFt);
            /* statistics */
            vEx = _mm_load_si128(pvEx+i);
            vMt = _mm_load_si128(pvMt+i);
            vSt = _mm_load_si128(pvSt+i);
            vLt = _mm_load_si128(pvLt+i);
            vEx = _mm_or_si128(
                    _mm_and_si128(vEx, _mm_cmpeq_epi8(vHt, vFt)),
                    _mm_cmplt_epi8(vHt, vFt));
            vM = _mm_blendv_epi8_rpl(vMt, vMp, vEx);
            vS = _mm_blendv_epi8_rpl(vSt, vSp, vEx);
            vL = _mm_blendv_epi8_rpl(vLt, vLp, vEx);
            vMp = vM;
            vSp = vS;
            vLp = _mm_adds_epi8(vL, vOne);
            vC = _mm_and_si128(vC, vEx);
            /* store results */
            _mm_store_si128(pvH+i, vH);
            _mm_store_si128(pvEx+i, vEx);
            /* check for saturation */
            {
                vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vH);
                vSaturationCheckMin = _mm_min_epi8_rpl(vSaturationCheckMin, vH);
            }
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif
        }
        {
            vLp = _mm_subs_epi8(vLp, vOne);
            {
                __m128i_8_t uMp, uSp, uLp, uC;
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
                uSp.m = vSp;
                uSp.v[ 1] = uC.v[ 1] ? uSp.v[ 0] : uSp.v[ 1];
                uSp.v[ 2] = uC.v[ 2] ? uSp.v[ 1] : uSp.v[ 2];
                uSp.v[ 3] = uC.v[ 3] ? uSp.v[ 2] : uSp.v[ 3];
                uSp.v[ 4] = uC.v[ 4] ? uSp.v[ 3] : uSp.v[ 4];
                uSp.v[ 5] = uC.v[ 5] ? uSp.v[ 4] : uSp.v[ 5];
                uSp.v[ 6] = uC.v[ 6] ? uSp.v[ 5] : uSp.v[ 6];
                uSp.v[ 7] = uC.v[ 7] ? uSp.v[ 6] : uSp.v[ 7];
                uSp.v[ 8] = uC.v[ 8] ? uSp.v[ 7] : uSp.v[ 8];
                uSp.v[ 9] = uC.v[ 9] ? uSp.v[ 8] : uSp.v[ 9];
                uSp.v[10] = uC.v[10] ? uSp.v[ 9] : uSp.v[10];
                uSp.v[11] = uC.v[11] ? uSp.v[10] : uSp.v[11];
                uSp.v[12] = uC.v[12] ? uSp.v[11] : uSp.v[12];
                uSp.v[13] = uC.v[13] ? uSp.v[12] : uSp.v[13];
                uSp.v[14] = uC.v[14] ? uSp.v[13] : uSp.v[14];
                uSp.v[15] = uC.v[15] ? uSp.v[14] : uSp.v[15];
                vSp = uSp.m;
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
        }
        /* final pass for M,L */
        vMp = _mm_slli_si128(vMp, 1);
        vSp = _mm_slli_si128(vSp, 1);
        vLp = _mm_slli_si128(vLp, 1);
        vLp = _mm_adds_epi8(vLp, vOne);
        for (i=0; i<segLen; ++i) {
            /* statistics */
            vEx = _mm_load_si128(pvEx+i);
            vMt = _mm_load_si128(pvMt+i);
            vSt = _mm_load_si128(pvSt+i);
            vLt = _mm_load_si128(pvLt+i);
            vM = _mm_blendv_epi8_rpl(vMt, vMp, vEx);
            vS = _mm_blendv_epi8_rpl(vSt, vSp, vEx);
            vL = _mm_blendv_epi8_rpl(vLt, vLp, vEx);
            vMp = vM;
            vSp = vS;
            vLp = _mm_adds_epi8(vL, vOne);
            /* store results */
            _mm_store_si128(pvM+i, vM);
            _mm_store_si128(pvS+i, vS);
            _mm_store_si128(pvL+i, vL);
            /* check for saturation */
            {
                vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vM);
                vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vS);
                vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vL);
            }
#ifdef PARASAIL_TABLE
            arr_store_si128(result->matches_table, vM, i, segLen, j, s2Len);
            arr_store_si128(result->similar_table, vS, i, segLen, j, s2Len);
            arr_store_si128(result->length_table, vL, i, segLen, j, s2Len);
#endif
        }

        /* extract vector containing last value from column */
        {
            __m128i cond_max;
            vH = _mm_load_si128(pvH + offset);
            vM = _mm_load_si128(pvM + offset);
            vS = _mm_load_si128(pvS + offset);
            vL = _mm_load_si128(pvL + offset);
            cond_max = _mm_cmpgt_epi8(vH, vMaxH);
            vMaxH = _mm_blendv_epi8_rpl(vMaxH, vH, cond_max);
            vMaxM = _mm_blendv_epi8_rpl(vMaxM, vM, cond_max);
            vMaxS = _mm_blendv_epi8_rpl(vMaxS, vS, cond_max);
            vMaxL = _mm_blendv_epi8_rpl(vMaxL, vL, cond_max);
            if (_mm_movemask_epi8(_mm_and_si128(vPosMask, cond_max))) {
                end_ref = j;
                end_query = s1Len - 1;
            }
#ifdef PARASAIL_ROWCOL
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128(vH, 1);
                vM = _mm_slli_si128(vM, 1);
                vS = _mm_slli_si128(vS, 1);
                vL = _mm_slli_si128(vL, 1);
            }
            result->score_row[j] = (int8_t) _mm_extract_epi8_rpl (vH, 15);
            result->matches_row[j] = (int8_t) _mm_extract_epi8_rpl (vM, 15);
            result->similar_row[j] = (int8_t) _mm_extract_epi8_rpl (vS, 15);
            result->length_row[j] = (int8_t) _mm_extract_epi8_rpl (vL, 15);
#endif
        }
    }

    /* max last value from all columns */
    {
        for (k=0; k<position; ++k) {
            vMaxH = _mm_slli_si128(vMaxH, 1);
            vMaxM = _mm_slli_si128(vMaxM, 1);
            vMaxS = _mm_slli_si128(vMaxS, 1);
            vMaxL = _mm_slli_si128(vMaxL, 1);
        }
        score = (int8_t) _mm_extract_epi8_rpl(vMaxH, 15);
        matches = (int8_t) _mm_extract_epi8_rpl(vMaxM, 15);
        similar = (int8_t) _mm_extract_epi8_rpl(vMaxS, 15);
        length = (int8_t) _mm_extract_epi8_rpl(vMaxL, 15);
    }

    /* max of last column */
    {
        int8_t score_last;
        vMaxH = vNegInf;

        for (i=0; i<segLen; ++i) {
            /* load the last stored values */
            __m128i vH = _mm_load_si128(pvH + i);
#ifdef PARASAIL_ROWCOL
            __m128i vM = _mm_load_si128(pvM + i);
            __m128i vS = _mm_load_si128(pvS + i);
            __m128i vL = _mm_load_si128(pvL + i);
            arr_store_col(result->score_col, vH, i, segLen);
            arr_store_col(result->matches_col, vM, i, segLen);
            arr_store_col(result->similar_col, vS, i, segLen);
            arr_store_col(result->length_col, vL, i, segLen);
#endif
            vMaxH = _mm_max_epi8_rpl(vH, vMaxH);
        }

        /* max in vec */
        score_last = _mm_hmax_epi8_rpl(vMaxH);
        if (score_last > score) {
            score = score_last;
            end_ref = s2Len - 1;
            end_query = s1Len;
            /* Trace the alignment ending position on read. */
            {
                int8_t *t = (int8_t*)pvH;
                int8_t *m = (int8_t*)pvM;
                int8_t *s = (int8_t*)pvS;
                int8_t *l = (int8_t*)pvL;
                int32_t column_len = segLen * segWidth;
                for (i = 0; i<column_len; ++i, ++t, ++m, ++s, ++l) {
                    if (*t == score) {
                        int32_t temp = i / segWidth + i % segWidth * segLen;
                        if (temp < end_query) {
                            end_query = temp;
                            matches = *m;
                            similar = *s;
                            length = *l;
                        }
                    }
                }
            }
        }
    }

    if (_mm_movemask_epi8(_mm_or_si128(
            _mm_cmpeq_epi8(vSaturationCheckMin, vNegLimit),
            _mm_cmpeq_epi8(vSaturationCheckMax, vPosLimit)))) {
        result->saturated = 1;
        score = INT8_MAX;
        matches = 0;
        similar = 0;
        length = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->matches = matches;
    result->similar = similar;
    result->length = length;
    result->end_query = end_query;
    result->end_ref = end_ref;

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

    return result;
}


