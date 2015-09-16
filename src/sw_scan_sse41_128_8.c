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
#include <string.h>

#include <emmintrin.h>
#include <smmintrin.h>

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_sse.h"

#define NEG_INF INT8_MIN

#define _mm_rlli_si128_rpl(a,imm) _mm_alignr_epi8(a, a, 16-imm)

static inline int8_t _mm_hmax_epi8_rpl(__m128i a) {
    a = _mm_max_epi8(a, _mm_srli_si128(a, 8));
    a = _mm_max_epi8(a, _mm_srli_si128(a, 4));
    a = _mm_max_epi8(a, _mm_srli_si128(a, 2));
    a = _mm_max_epi8(a, _mm_srli_si128(a, 1));
    return _mm_extract_epi8(a, 0);
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
    array[( 0*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  0);
    array[( 1*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  1);
    array[( 2*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  2);
    array[( 3*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  3);
    array[( 4*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  4);
    array[( 5*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  5);
    array[( 6*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  6);
    array[( 7*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  7);
    array[( 8*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  8);
    array[( 9*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  9);
    array[(10*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 10);
    array[(11*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 11);
    array[(12*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 12);
    array[(13*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 13);
    array[(14*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 14);
    array[(15*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 15);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m128i vH,
        int32_t t,
        int32_t seglen)
{
    col[ 0*seglen+t] = (int8_t)_mm_extract_epi8(vH,  0);
    col[ 1*seglen+t] = (int8_t)_mm_extract_epi8(vH,  1);
    col[ 2*seglen+t] = (int8_t)_mm_extract_epi8(vH,  2);
    col[ 3*seglen+t] = (int8_t)_mm_extract_epi8(vH,  3);
    col[ 4*seglen+t] = (int8_t)_mm_extract_epi8(vH,  4);
    col[ 5*seglen+t] = (int8_t)_mm_extract_epi8(vH,  5);
    col[ 6*seglen+t] = (int8_t)_mm_extract_epi8(vH,  6);
    col[ 7*seglen+t] = (int8_t)_mm_extract_epi8(vH,  7);
    col[ 8*seglen+t] = (int8_t)_mm_extract_epi8(vH,  8);
    col[ 9*seglen+t] = (int8_t)_mm_extract_epi8(vH,  9);
    col[10*seglen+t] = (int8_t)_mm_extract_epi8(vH, 10);
    col[11*seglen+t] = (int8_t)_mm_extract_epi8(vH, 11);
    col[12*seglen+t] = (int8_t)_mm_extract_epi8(vH, 12);
    col[13*seglen+t] = (int8_t)_mm_extract_epi8(vH, 13);
    col[14*seglen+t] = (int8_t)_mm_extract_epi8(vH, 14);
    col[15*seglen+t] = (int8_t)_mm_extract_epi8(vH, 15);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_table_scan_sse41_128_8
#define PNAME parasail_sw_table_scan_profile_sse41_128_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_rowcol_scan_sse41_128_8
#define PNAME parasail_sw_rowcol_scan_profile_sse41_128_8
#else
#define FNAME parasail_sw_scan_sse41_128_8
#define PNAME parasail_sw_scan_profile_sse41_128_8
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_sse_128_8(s1, s1Len, matrix);
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
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int32_t segNum = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m128i* const restrict pvP = (__m128i*)profile->profile8.score;
    __m128i* const restrict pvE = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHt= parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvH = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHMax = parasail_memalign___m128i(16, segLen);
    __m128i vGapO = _mm_set1_epi8(open);
    __m128i vGapE = _mm_set1_epi8(gap);
    __m128i vNegInf = _mm_set1_epi8(NEG_INF);
    __m128i vZero = _mm_setzero_si128();
    int8_t score = NEG_INF;
    __m128i vMaxH = vNegInf;
    __m128i vMaxHUnit = vNegInf;
    const int8_t segLenXgap = -segLen*gap;
    __m128i insert_mask = _mm_cmpeq_epi8(vZero,
            _mm_set_epi8(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0));
    __m128i vSegLenXgap1 = _mm_set1_epi8((segLen-1)*gap);
    __m128i vSegLenXgap = _mm_blendv_epi8(vNegInf,
            _mm_set1_epi8(segLenXgap),
            insert_mask);
    __m128i vNegLimit = _mm_set1_epi8(INT8_MIN);
    __m128i vPosLimit = _mm_set1_epi8(INT8_MAX);
    __m128i vSaturationCheckMin = vPosLimit;
    __m128i vSaturationCheckMax = vNegLimit;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(segLen*segWidth, s2Len);
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif

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
        __m128i vHp;
        __m128i *pvW;
        __m128i vW;

        /* calculate E */
        /* calculate Ht */
        /* calculate Ft first pass */
        vHp = _mm_load_si128(pvH+(segLen-1));
        vHp = _mm_slli_si128(vHp, 1);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = vNegInf;
        vFt = vNegInf;
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            vE = _mm_max_epi8(
                    _mm_subs_epi8(vE, vGapE),
                    _mm_subs_epi8(vH, vGapO));
            vFt = _mm_subs_epi8(vFt, vGapE);
            vFt = _mm_max_epi8(vFt, vHt);
            vHt = _mm_max_epi8(
                    _mm_adds_epi8(vHp, vW),
                    vE);
            _mm_store_si128(pvE+i, vE);
            _mm_store_si128(pvHt+i, vHt);
            vHp = vH;
        }

        /* adjust Ft before local prefix scan */
        vHt = _mm_slli_si128(vHt, 1);
        vFt = _mm_max_epi8(vFt,
                _mm_subs_epi8(vHt, vSegLenXgap1));
        /* local prefix scan */
        vFt = _mm_blendv_epi8(vNegInf, vFt, insert_mask);
        for (i=0; i<segWidth-1; ++i) {
                __m128i vFtt = _mm_rlli_si128_rpl(vFt, 1);
                vFtt = _mm_adds_epi8(vFtt, vSegLenXgap);
                vFt = _mm_max_epi8(vFt, vFtt);
        }
        vFt = _mm_rlli_si128_rpl(vFt, 1);

        /* second Ft pass */
        /* calculate vH */
        for (i=0; i<segLen; ++i) {
            vFt = _mm_subs_epi8(vFt, vGapE);
            vFt = _mm_max_epi8(vFt, vHt);
            vHt = _mm_load_si128(pvHt+i);
            vH = _mm_max_epi8(vHt, _mm_subs_epi8(vFt, vGapO));
            vH = _mm_max_epi8(vH, vZero);
            _mm_store_si128(pvH+i, vH);
            /* check for saturation */
            {
                vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vH);
                vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vH);
            }
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif
            vMaxH = _mm_max_epi8(vH, vMaxH);
        }

        {
            __m128i vCompare = _mm_cmpgt_epi8(vMaxH, vMaxHUnit);
            if (_mm_movemask_epi8(vCompare)) {
                score = _mm_hmax_epi8_rpl(vMaxH);
                vMaxHUnit = _mm_set1_epi8(score);
                end_ref = j;
                (void)memcpy(pvHMax, pvH, sizeof(__m128i)*segLen);
            }
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            int32_t k = 0;
            vH = _mm_load_si128(pvH + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128(vH, 1);
            }
            result->score_row[j] = (int8_t) _mm_extract_epi8 (vH, 15);
        }
#endif
    }

    /* Trace the alignment ending position on read. */
    {
        int8_t *t = (int8_t*)pvHMax;
        int32_t column_len = segLen * segWidth;
        end_query = s1Len;
        for (i = 0; i<column_len; ++i, ++t) {
            if (*t == score) {
                int32_t temp = i / segWidth + i % segWidth * segLen;
                if (temp < end_query) {
                    end_query = temp;
                }
            }
        }
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        __m128i vH = _mm_load_si128(pvH+i);
        arr_store_col(result->score_col, vH, i, segLen);
    }
#endif

    if (_mm_movemask_epi8(_mm_or_si128(
            _mm_cmpeq_epi8(vSaturationCheckMin, vNegLimit),
            _mm_cmpeq_epi8(vSaturationCheckMax, vPosLimit)))) {
        result->saturated = 1;
        score = INT8_MAX;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(pvHMax);
    parasail_free(pvH);
    parasail_free(pvHt);
    parasail_free(pvE);

    return result;
}

