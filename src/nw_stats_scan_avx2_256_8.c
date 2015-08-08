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

#include <immintrin.h>

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_avx.h"

#define NEG_INF INT8_MIN
#define MAX(a,b) ((a)>(b)?(a):(b))

#if HAVE_AVX2_MM256_INSERT_EPI8
#define _mm256_insert_epi8_rpl _mm256_insert_epi8
#else
static inline __m256i _mm256_insert_epi8_rpl(__m256i a, int8_t i, int imm) {
    __m256i_8_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI8
#define _mm256_extract_epi8_rpl _mm256_extract_epi8
#else
static inline int8_t _mm256_extract_epi8_rpl(__m256i a, int imm) {
    __m256i_8_t A;
    A.m = a;
    return A.v[imm];
}
#endif

#define _mm256_cmplt_epi8_rpl(a,b) _mm256_cmpgt_epi8(b,a)

#define _mm256_srli_si256_rpl(a,imm) _mm256_or_si256(_mm256_slli_si256(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(3,0,0,1)), 16-imm), _mm256_srli_si256(a, imm))

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)


#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[( 0*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  0);
    array[( 1*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  1);
    array[( 2*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  2);
    array[( 3*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  3);
    array[( 4*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  4);
    array[( 5*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  5);
    array[( 6*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  6);
    array[( 7*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  7);
    array[( 8*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  8);
    array[( 9*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  9);
    array[(10*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 10);
    array[(11*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 11);
    array[(12*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 12);
    array[(13*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 13);
    array[(14*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 14);
    array[(15*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 15);
    array[(16*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 16);
    array[(17*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 17);
    array[(18*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 18);
    array[(19*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 19);
    array[(20*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 20);
    array[(21*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 21);
    array[(22*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 22);
    array[(23*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 23);
    array[(24*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 24);
    array[(25*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 25);
    array[(26*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 26);
    array[(27*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 27);
    array[(28*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 28);
    array[(29*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 29);
    array[(30*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 30);
    array[(31*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 31);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m256i vH,
        int32_t t,
        int32_t seglen)
{
    col[ 0*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  0);
    col[ 1*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  1);
    col[ 2*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  2);
    col[ 3*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  3);
    col[ 4*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  4);
    col[ 5*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  5);
    col[ 6*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  6);
    col[ 7*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  7);
    col[ 8*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  8);
    col[ 9*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  9);
    col[10*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 10);
    col[11*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 11);
    col[12*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 12);
    col[13*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 13);
    col[14*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 14);
    col[15*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 15);
    col[16*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 16);
    col[17*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 17);
    col[18*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 18);
    col[19*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 19);
    col[20*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 20);
    col[21*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 21);
    col[22*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 22);
    col[23*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 23);
    col[24*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 24);
    col[25*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 25);
    col[26*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 26);
    col[27*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 27);
    col[28*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 28);
    col[29*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 29);
    col[30*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 30);
    col[31*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 31);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_nw_stats_table_scan_avx2_256_8
#define PNAME parasail_nw_stats_table_scan_profile_avx2_256_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_stats_rowcol_scan_avx2_256_8
#define PNAME parasail_nw_stats_rowcol_scan_profile_avx2_256_8
#else
#define FNAME parasail_nw_stats_scan_avx2_256_8
#define PNAME parasail_nw_stats_scan_profile_avx2_256_8
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_avx_256_8(s1, s1Len, matrix);
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
    const int32_t segWidth = 32;
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m256i* const restrict pvP  = (__m256i*)profile->profile8.score;
    __m256i* const restrict pvPm = (__m256i*)profile->profile8.matches;
    __m256i* const restrict pvPs = (__m256i*)profile->profile8.similar;
    __m256i* const restrict pvE  = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvHt = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvFt = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvMt = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvSt = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvLt = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvEx = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvH  = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvM  = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvS  = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvL  = parasail_memalign___m256i(32, segLen);
    int8_t* const restrict boundary = parasail_memalign_int8_t(32, s2Len+1);
    __m256i vGapO = _mm256_set1_epi8(open);
    __m256i vGapE = _mm256_set1_epi8(gap);
    __m256i vZero = _mm256_setzero_si256();
    __m256i vOne = _mm256_set1_epi8(1);
    __m256i vNegInf = _mm256_set1_epi8(NEG_INF);
    int8_t score = NEG_INF;
    int8_t matches = 0;
    int8_t similar = 0;
    int8_t length = 0;
    __m256i vNegLimit = _mm256_set1_epi8(INT8_MIN);
    __m256i vPosLimit = _mm256_set1_epi8(INT8_MAX);
    __m256i vSaturationCheckMin = vPosLimit;
    __m256i vSaturationCheckMax = vNegLimit;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif

    parasail_memset___m256i(pvM, vZero, segLen);
    parasail_memset___m256i(pvS, vZero, segLen);
    parasail_memset___m256i(pvL, vZero, segLen);

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m256i_8_t h;
            __m256i_8_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.v[segNum] = tmp < INT8_MIN ? INT8_MIN : tmp;
                e.v[segNum] = NEG_INF;
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
            int64_t tmp = -open-gap*(i-1);
            boundary[i] = tmp < INT8_MIN ? INT8_MIN : tmp;
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
            vE = _mm256_max_epi8(
                    _mm256_subs_epi8(vE, vGapE),
                    _mm256_subs_epi8(vH, vGapO));
            _mm256_store_si256(pvE+i, vE);
        }

        /* calculate Ht */
        vH = _mm256_slli_si256_rpl(_mm256_load_si256(pvH+(segLen-1)), 1);
        vH = _mm256_insert_epi8_rpl(vH, boundary[j], 0);
        vMp= _mm256_slli_si256_rpl(_mm256_load_si256(pvM+(segLen-1)), 1);
        vSp= _mm256_slli_si256_rpl(_mm256_load_si256(pvS+(segLen-1)), 1);
        vLp= _mm256_slli_si256_rpl(_mm256_load_si256(pvL+(segLen-1)), 1);
        vLp= _mm256_adds_epi8(vLp, vOne);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvC = pvPm+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvD = pvPs+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vE = _mm256_load_si256(pvE+i);
            vW = _mm256_load_si256(pvW+i);
            /* compute */
            vH = _mm256_adds_epi8(vH, vW);
            vHt = _mm256_max_epi8(vH, vE);
            /* statistics */
            vC = _mm256_load_si256(pvC+i);
            vD = _mm256_load_si256(pvD+i);
            vMp = _mm256_adds_epi8(vMp, vC);
            vSp = _mm256_adds_epi8(vSp, vD);
            vEx = _mm256_cmpgt_epi8(vE, vH);
            vM = _mm256_load_si256(pvM+i);
            vS = _mm256_load_si256(pvS+i);
            vL = _mm256_load_si256(pvL+i);
            vL = _mm256_adds_epi8(vL, vOne);
            vMt = _mm256_blendv_epi8(vMp, vM, vEx);
            vSt = _mm256_blendv_epi8(vSp, vS, vEx);
            vLt = _mm256_blendv_epi8(vLp, vL, vEx);
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
        vHt = _mm256_slli_si256_rpl(_mm256_load_si256(pvHt+(segLen-1)), 1);
        vHt = _mm256_insert_epi8_rpl(vHt, boundary[j+1], 0);
        vFt = vNegInf;
        for (i=0; i<segLen; ++i) {
            vFt = _mm256_subs_epi8(vFt, vGapE);
            vFt = _mm256_max_epi8(vFt, vHt);
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
        vHt = _mm256_slli_si256_rpl(_mm256_load_si256(pvHt+(segLen-1)), 1);
        vHt = _mm256_insert_epi8_rpl(vHt, boundary[j+1], 0);
        vFt = _mm256_slli_si256_rpl(vFt, 1);
        vFt = _mm256_insert_epi8_rpl(vFt, NEG_INF, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm256_subs_epi8(vFt, vGapE);
            vFt = _mm256_max_epi8(vFt, vHt);
            vHt = _mm256_load_si256(pvHt+i);
            _mm256_store_si256(pvFt+i, vFt);
        }

        /* calculate H,M,L */
        vMp = vZero;
        vSp = vZero;
        vLp = vOne;
        vC = _mm256_cmpeq_epi8(vZero, vZero); /* check if prefix sum is needed */
        vC = _mm256_srli_si256_rpl(vC, 1); /* zero out last value */
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
            vSt = _mm256_load_si256(pvSt+i);
            vLt = _mm256_load_si256(pvLt+i);
            vEx = _mm256_or_si256(
                    _mm256_and_si256(vEx, _mm256_cmpeq_epi8(vHt, vFt)),
                    _mm256_cmplt_epi8_rpl(vHt, vFt));
            vM = _mm256_blendv_epi8(vMt, vMp, vEx);
            vS = _mm256_blendv_epi8(vSt, vSp, vEx);
            vL = _mm256_blendv_epi8(vLt, vLp, vEx);
            vMp = vM;
            vSp = vS;
            vLp = _mm256_adds_epi8(vL, vOne);
            vC = _mm256_and_si256(vC, vEx);
            /* store results */
            _mm256_store_si256(pvH+i, vH);
            _mm256_store_si256(pvEx+i, vEx);
            /* check for saturation */
            {
                vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vH);
                vSaturationCheckMin = _mm256_min_epi8(vSaturationCheckMin, vH);
            }
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif
        }
        {
            vLp = _mm256_subs_epi8(vLp, vOne);
            {
                __m256i_8_t uMp, uSp, uLp, uC;
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
                uMp.v[16] = uC.v[16] ? uMp.v[15] : uMp.v[16];
                uMp.v[17] = uC.v[17] ? uMp.v[16] : uMp.v[17];
                uMp.v[18] = uC.v[18] ? uMp.v[17] : uMp.v[18];
                uMp.v[19] = uC.v[19] ? uMp.v[18] : uMp.v[19];
                uMp.v[20] = uC.v[20] ? uMp.v[19] : uMp.v[20];
                uMp.v[21] = uC.v[21] ? uMp.v[20] : uMp.v[21];
                uMp.v[22] = uC.v[22] ? uMp.v[21] : uMp.v[22];
                uMp.v[23] = uC.v[23] ? uMp.v[22] : uMp.v[23];
                uMp.v[24] = uC.v[24] ? uMp.v[23] : uMp.v[24];
                uMp.v[25] = uC.v[25] ? uMp.v[24] : uMp.v[25];
                uMp.v[26] = uC.v[26] ? uMp.v[25] : uMp.v[26];
                uMp.v[27] = uC.v[27] ? uMp.v[26] : uMp.v[27];
                uMp.v[28] = uC.v[28] ? uMp.v[27] : uMp.v[28];
                uMp.v[29] = uC.v[29] ? uMp.v[28] : uMp.v[29];
                uMp.v[30] = uC.v[30] ? uMp.v[29] : uMp.v[30];
                uMp.v[31] = uC.v[31] ? uMp.v[30] : uMp.v[31];
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
                uSp.v[16] = uC.v[16] ? uSp.v[15] : uSp.v[16];
                uSp.v[17] = uC.v[17] ? uSp.v[16] : uSp.v[17];
                uSp.v[18] = uC.v[18] ? uSp.v[17] : uSp.v[18];
                uSp.v[19] = uC.v[19] ? uSp.v[18] : uSp.v[19];
                uSp.v[20] = uC.v[20] ? uSp.v[19] : uSp.v[20];
                uSp.v[21] = uC.v[21] ? uSp.v[20] : uSp.v[21];
                uSp.v[22] = uC.v[22] ? uSp.v[21] : uSp.v[22];
                uSp.v[23] = uC.v[23] ? uSp.v[22] : uSp.v[23];
                uSp.v[24] = uC.v[24] ? uSp.v[23] : uSp.v[24];
                uSp.v[25] = uC.v[25] ? uSp.v[24] : uSp.v[25];
                uSp.v[26] = uC.v[26] ? uSp.v[25] : uSp.v[26];
                uSp.v[27] = uC.v[27] ? uSp.v[26] : uSp.v[27];
                uSp.v[28] = uC.v[28] ? uSp.v[27] : uSp.v[28];
                uSp.v[29] = uC.v[29] ? uSp.v[28] : uSp.v[29];
                uSp.v[30] = uC.v[30] ? uSp.v[29] : uSp.v[30];
                uSp.v[31] = uC.v[31] ? uSp.v[30] : uSp.v[31];
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
                uLp.v[16] = uC.v[16] ? uLp.v[16] + uLp.v[15] : uLp.v[16];
                uLp.v[17] = uC.v[17] ? uLp.v[17] + uLp.v[16] : uLp.v[17];
                uLp.v[18] = uC.v[18] ? uLp.v[18] + uLp.v[17] : uLp.v[18];
                uLp.v[19] = uC.v[19] ? uLp.v[19] + uLp.v[18] : uLp.v[19];
                uLp.v[20] = uC.v[20] ? uLp.v[20] + uLp.v[19] : uLp.v[20];
                uLp.v[21] = uC.v[21] ? uLp.v[21] + uLp.v[20] : uLp.v[21];
                uLp.v[22] = uC.v[22] ? uLp.v[22] + uLp.v[21] : uLp.v[22];
                uLp.v[23] = uC.v[23] ? uLp.v[23] + uLp.v[22] : uLp.v[23];
                uLp.v[24] = uC.v[24] ? uLp.v[24] + uLp.v[23] : uLp.v[24];
                uLp.v[25] = uC.v[25] ? uLp.v[25] + uLp.v[24] : uLp.v[25];
                uLp.v[26] = uC.v[26] ? uLp.v[26] + uLp.v[25] : uLp.v[26];
                uLp.v[27] = uC.v[27] ? uLp.v[27] + uLp.v[26] : uLp.v[27];
                uLp.v[28] = uC.v[28] ? uLp.v[28] + uLp.v[27] : uLp.v[28];
                uLp.v[29] = uC.v[29] ? uLp.v[29] + uLp.v[28] : uLp.v[29];
                uLp.v[30] = uC.v[30] ? uLp.v[30] + uLp.v[29] : uLp.v[30];
                uLp.v[31] = uC.v[31] ? uLp.v[31] + uLp.v[30] : uLp.v[31];
                vLp = uLp.m;
            }
            vLp = _mm256_adds_epi8(vLp, vOne);
        }
        /* final pass for M,L */
        vMp = _mm256_slli_si256_rpl(vMp, 1);
        vSp = _mm256_slli_si256_rpl(vSp, 1);
        vLp = _mm256_slli_si256_rpl(vLp, 1);
        for (i=0; i<segLen; ++i) {
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
            vLp = _mm256_adds_epi8(vL, vOne);
            /* store results */
            _mm256_store_si256(pvM+i, vM);
            _mm256_store_si256(pvS+i, vS);
            _mm256_store_si256(pvL+i, vL);
            /* check for saturation */
            {
                vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vM);
                vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vS);
                vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vL);
            }
#ifdef PARASAIL_TABLE
            arr_store_si256(result->matches_table, vM, i, segLen, j, s2Len);
            arr_store_si256(result->similar_table, vS, i, segLen, j, s2Len);
            arr_store_si256(result->length_table, vL, i, segLen, j, s2Len);
#endif
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = _mm256_load_si256(pvH + offset);
            vM = _mm256_load_si256(pvM + offset);
            vS = _mm256_load_si256(pvS + offset);
            vL = _mm256_load_si256(pvL + offset);
            for (k=0; k<position; ++k) {
                vH = _mm256_slli_si256_rpl(vH, 1);
                vM = _mm256_slli_si256_rpl(vM, 1);
                vS = _mm256_slli_si256_rpl(vS, 1);
                vL = _mm256_slli_si256_rpl(vL, 1);
            }
            result->score_row[j] = (int8_t) _mm256_extract_epi8_rpl (vH, 31);
            result->matches_row[j] = (int8_t) _mm256_extract_epi8_rpl (vM, 31);
            result->similar_row[j] = (int8_t) _mm256_extract_epi8_rpl (vS, 31);
            result->length_row[j] = (int8_t) _mm256_extract_epi8_rpl (vL, 31);
        }
#endif
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        __m256i vH = _mm256_load_si256(pvH+i);
        __m256i vM = _mm256_load_si256(pvM+i);
        __m256i vS = _mm256_load_si256(pvS+i);
        __m256i vL = _mm256_load_si256(pvL+i);
        arr_store_col(result->score_col, vH, i, segLen);
        arr_store_col(result->matches_col, vM, i, segLen);
        arr_store_col(result->similar_col, vS, i, segLen);
        arr_store_col(result->length_col, vL, i, segLen);
    }
#endif

    /* extract last value from the last column */
    {
        __m256i vH = _mm256_load_si256(pvH + offset);
        __m256i vM = _mm256_load_si256(pvM + offset);
        __m256i vS = _mm256_load_si256(pvS + offset);
        __m256i vL = _mm256_load_si256(pvL + offset);
        for (k=0; k<position; ++k) {
            vH = _mm256_slli_si256_rpl(vH, 1);
            vM = _mm256_slli_si256_rpl(vM, 1);
            vS = _mm256_slli_si256_rpl(vS, 1);
            vL = _mm256_slli_si256_rpl(vL, 1);
        }
        score = (int8_t) _mm256_extract_epi8_rpl (vH, 31);
        matches = (int8_t) _mm256_extract_epi8_rpl (vM, 31);
        similar = (int8_t) _mm256_extract_epi8_rpl (vS, 31);
        length = (int8_t) _mm256_extract_epi8_rpl (vL, 31);
    }

    if (_mm256_movemask_epi8(_mm256_or_si256(
            _mm256_cmpeq_epi8(vSaturationCheckMin, vNegLimit),
            _mm256_cmpeq_epi8(vSaturationCheckMax, vPosLimit)))) {
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

    parasail_free(boundary);
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


