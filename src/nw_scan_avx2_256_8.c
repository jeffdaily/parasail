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

#define _mm256_rlli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,1)), 16-imm)

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
#define FNAME parasail_nw_table_scan_avx2_256_8
#define PNAME parasail_nw_table_scan_profile_avx2_256_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_rowcol_scan_avx2_256_8
#define PNAME parasail_nw_rowcol_scan_profile_avx2_256_8
#else
#define FNAME parasail_nw_scan_avx2_256_8
#define PNAME parasail_nw_scan_profile_avx2_256_8
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_avx_256_8(s1, s1Len, matrix);
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
    const int32_t segWidth = 32; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m256i* const restrict pvP = (__m256i*)profile->profile8.score;
    __m256i* const restrict pvE = parasail_memalign___m256i(32, segLen);
    int8_t* const restrict boundary = parasail_memalign_int8_t(32, s2Len+1);
    __m256i* const restrict pvHt= parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvH = parasail_memalign___m256i(32, segLen);
    __m256i vGapO = _mm256_set1_epi8(open);
    __m256i vGapE = _mm256_set1_epi8(gap);
    __m256i vNegInf = _mm256_set1_epi8(NEG_INF);
    int8_t score = NEG_INF;
    const int8_t segLenXgap = -segLen*gap;
    __m256i insert_mask = _mm256_cmpeq_epi8(_mm256_setzero_si256(),
            _mm256_set_epi8(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0));
    __m256i vSegLenXgap1 = _mm256_set1_epi8((segLen-1)*gap);
    __m256i vSegLenXgap = _mm256_blendv_epi8(vNegInf,
            _mm256_set1_epi8(segLenXgap),
            insert_mask);
    __m256i vNegLimit = _mm256_set1_epi8(INT8_MIN);
    __m256i vPosLimit = _mm256_set1_epi8(INT8_MAX);
    __m256i vSaturationCheckMin = vPosLimit;
    __m256i vSaturationCheckMax = vNegLimit;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m256i_8_t h;
            __m256i_8_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.v[segNum] = tmp < INT8_MIN ? INT8_MIN : tmp;
                tmp = tmp - open;
                e.v[segNum] = tmp < INT8_MIN ? INT8_MIN : tmp;
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
        __m256i vHp;
        __m256i *pvW;
        __m256i vW;

        /* calculate E */
        /* calculate Ht */
        /* calculate Ft first pass */
        vHp = _mm256_load_si256(pvH+(segLen-1));
        vHp = _mm256_slli_si256_rpl(vHp, 1);
        vHp = _mm256_insert_epi8_rpl(vHp, boundary[j], 0);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = vNegInf;
        vFt = vNegInf;
        for (i=0; i<segLen; ++i) {
            vH = _mm256_load_si256(pvH+i);
            vE = _mm256_load_si256(pvE+i);
            vW = _mm256_load_si256(pvW+i);
            vE = _mm256_max_epi8(
                    _mm256_subs_epi8(vE, vGapE),
                    _mm256_subs_epi8(vH, vGapO));
            vFt = _mm256_subs_epi8(vFt, vGapE);
            vFt = _mm256_max_epi8(vFt, vHt);
            vHt = _mm256_max_epi8(
                    _mm256_adds_epi8(vHp, vW),
                    vE);
            _mm256_store_si256(pvE+i, vE);
            _mm256_store_si256(pvHt+i, vHt);
            vHp = vH;
        }

        /* adjust Ft before local prefix scan */
        vHt = _mm256_slli_si256_rpl(vHt, 1);
        vHt = _mm256_insert_epi8_rpl(vHt, boundary[j+1], 0);
        vFt = _mm256_max_epi8(vFt,
                _mm256_subs_epi8(vHt, vSegLenXgap1));
        /* local prefix scan */
        vFt = _mm256_blendv_epi8(vNegInf, vFt, insert_mask);
            for (i=0; i<segWidth-1; ++i) {
                __m256i vFtt = _mm256_rlli_si256_rpl(vFt, 1);
                vFtt = _mm256_adds_epi8(vFtt, vSegLenXgap);
                vFt = _mm256_max_epi8(vFt, vFtt);
            }
        vFt = _mm256_rlli_si256_rpl(vFt, 1);

        /* second Ft pass */
        /* calculate vH */
        for (i=0; i<segLen; ++i) {
            vFt = _mm256_subs_epi8(vFt, vGapE);
            vFt = _mm256_max_epi8(vFt, vHt);
            vHt = _mm256_load_si256(pvHt+i);
            vH = _mm256_max_epi8(vHt, _mm256_subs_epi8(vFt, vGapO));
            _mm256_store_si256(pvH+i, vH);
            /* check for saturation */
            {
                vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vH);
                vSaturationCheckMin = _mm256_min_epi8(vSaturationCheckMin, vH);
            }
#ifdef PARASAIL_TABLE
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = _mm256_load_si256(pvH + offset);
            for (k=0; k<position; ++k) {
                vH = _mm256_slli_si256_rpl(vH, 1);
            }
            result->score_row[j] = (int8_t) _mm256_extract_epi8_rpl (vH, 31);
        }
#endif
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        __m256i vH = _mm256_load_si256(pvH+i);
        arr_store_col(result->score_col, vH, i, segLen);
    }
#endif

    /* extract last value from the last column */
    {
        __m256i vH = _mm256_load_si256(pvH + offset);
        for (k=0; k<position; ++k) {
            vH = _mm256_slli_si256_rpl(vH, 1);
        }
        score = (int8_t) _mm256_extract_epi8_rpl (vH, 31);
    }

    if (_mm256_movemask_epi8(_mm256_or_si256(
            _mm256_cmpeq_epi8(vSaturationCheckMin, vNegLimit),
            _mm256_cmpeq_epi8(vSaturationCheckMax, vPosLimit)))) {
        result->saturated = 1;
        score = INT8_MAX;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(boundary);
    parasail_free(pvH);
    parasail_free(pvHt);
    parasail_free(pvE);

    return result;
}


