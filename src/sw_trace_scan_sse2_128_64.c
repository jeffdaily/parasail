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

#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <emmintrin.h>
#endif

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_sse.h"

#define NEG_INF (INT64_MIN/(int64_t)(2))

static inline int64_t _mm_extract_epi64_rpl(__m128i a, const int imm) {
    __m128i_64_t A;
    A.m = a;
    return A.v[imm];
}

static inline __m128i _mm_blendv_epi8_rpl(__m128i a, __m128i b, __m128i mask) {
    a = _mm_andnot_si128(mask, a);
    a = _mm_or_si128(a, _mm_and_si128(mask, b));
    return a;
}

static inline __m128i _mm_cmpgt_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? 0xFFFFFFFFFFFFFFFF : 0;
    A.v[1] = (A.v[1]>B.v[1]) ? 0xFFFFFFFFFFFFFFFF : 0;
    return A.m;
}

static inline __m128i _mm_cmpeq_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]==B.v[0]) ? 0xFFFFFFFFFFFFFFFF : 0;
    A.v[1] = (A.v[1]==B.v[1]) ? 0xFFFFFFFFFFFFFFFF : 0;
    return A.m;
}

static inline __m128i _mm_max_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]>B.v[1]) ? A.v[1] : B.v[1];
    return A.m;
}

#if HAVE_SSE2_MM_SET_EPI64X
#define _mm_set_epi64x_rpl _mm_set_epi64x
#else
static inline __m128i _mm_set_epi64x_rpl(int64_t e1, int64_t e0) {
    __m128i_64_t A;
    A.v[0] = e0;
    A.v[1] = e1;
    return A.m;
}
#endif

#define _mm_rlli_si128_rpl(a,imm) _mm_or_si128(_mm_slli_si128(a,imm),_mm_srli_si128(a,16-imm))

#if HAVE_SSE2_MM_SET1_EPI64X
#define _mm_set1_epi64x_rpl _mm_set1_epi64x
#else
static inline __m128i _mm_set1_epi64x_rpl(int64_t i) {
    __m128i_64_t A;
    A.v[0] = i;
    A.v[1] = i;
    return A.m;
}
#endif

static inline int64_t _mm_hmax_epi64_rpl(__m128i a) {
    a = _mm_max_epi64_rpl(a, _mm_srli_si128(a, 8));
    return _mm_extract_epi64_rpl(a, 0);
}


static inline void arr_store_si128(
        int *array,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[(0*seglen+t)*dlen + d] = (int64_t)_mm_extract_epi64_rpl(vH, 0);
    array[(1*seglen+t)*dlen + d] = (int64_t)_mm_extract_epi64_rpl(vH, 1);
}

#define FNAME parasail_sw_trace_scan_sse2_128_64
#define PNAME parasail_sw_trace_scan_profile_sse2_128_64

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_sse_128_64(s1, s1Len, matrix);
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
    const int32_t segWidth = 2; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m128i* const restrict pvP = (__m128i*)profile->profile64.score;
    __m128i* const restrict pvE = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHt= parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvH = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHMax = parasail_memalign___m128i(16, segLen);
    __m128i vGapO = _mm_set1_epi64x_rpl(open);
    __m128i vGapE = _mm_set1_epi64x_rpl(gap);
    __m128i vNegInf = _mm_set1_epi64x_rpl(NEG_INF);
    __m128i vZero = _mm_setzero_si128();
    int64_t score = NEG_INF;
    __m128i vMaxH = vNegInf;
    __m128i vMaxHUnit = vNegInf;
    const int64_t segLenXgap = -segLen*gap;
    __m128i insert_mask = _mm_cmpeq_epi64_rpl(vZero,
            _mm_set_epi64x_rpl(1,0));
    __m128i vSegLenXgap1 = _mm_set1_epi64x_rpl((segLen-1)*gap);
    __m128i vSegLenXgap = _mm_blendv_epi8_rpl(vNegInf,
            _mm_set1_epi64x_rpl(segLenXgap),
            insert_mask);
    
    parasail_result_t *result = parasail_result_new_trace(segLen*segWidth, s2Len, 8);

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m128i_64_t h;
            __m128i_64_t e;
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
        vHp = _mm_slli_si128(vHp, 8);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = vNegInf;
        vFt = vNegInf;
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            vE = _mm_max_epi64_rpl(
                    _mm_sub_epi64(vE, vGapE),
                    _mm_sub_epi64(vH, vGapO));
            vFt = _mm_sub_epi64(vFt, vGapE);
            vFt = _mm_max_epi64_rpl(vFt, vHt);
            vHt = _mm_max_epi64_rpl(
                    _mm_add_epi64(vHp, vW),
                    vE);
            _mm_store_si128(pvE+i, vE);
            _mm_store_si128(pvHt+i, vHt);
            vHp = vH;
        }

        /* adjust Ft before local prefix scan */
        vHt = _mm_slli_si128(vHt, 8);
        vFt = _mm_max_epi64_rpl(vFt,
                _mm_sub_epi64(vHt, vSegLenXgap1));
        /* local prefix scan */
        vFt = _mm_blendv_epi8_rpl(vNegInf, vFt, insert_mask);
        for (i=0; i<segWidth-1; ++i) {
                __m128i vFtt = _mm_rlli_si128_rpl(vFt, 8);
                vFtt = _mm_add_epi64(vFtt, vSegLenXgap);
                vFt = _mm_max_epi64_rpl(vFt, vFtt);
        }
        vFt = _mm_rlli_si128_rpl(vFt, 8);

        /* second Ft pass */
        /* calculate vH */
        for (i=0; i<segLen; ++i) {
            vFt = _mm_sub_epi64(vFt, vGapE);
            vFt = _mm_max_epi64_rpl(vFt, vHt);
            vHt = _mm_load_si128(pvHt+i);
            vH = _mm_max_epi64_rpl(vHt, _mm_sub_epi64(vFt, vGapO));
            vH = _mm_max_epi64_rpl(vH, vZero);
            _mm_store_si128(pvH+i, vH);
            
            vMaxH = _mm_max_epi64_rpl(vH, vMaxH);
        }

        {
            __m128i vCompare = _mm_cmpgt_epi64_rpl(vMaxH, vMaxHUnit);
            if (_mm_movemask_epi8(vCompare)) {
                score = _mm_hmax_epi64_rpl(vMaxH);
                vMaxHUnit = _mm_set1_epi64x_rpl(score);
                end_ref = j;
                (void)memcpy(pvHMax, pvH, sizeof(__m128i)*segLen);
            }
        }
    }

    /* Trace the alignment ending position on read. */
    {
        int64_t *t = (int64_t*)pvHMax;
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

    

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag = PARASAIL_FLAG_SW | PARASAIL_FLAG_SCAN
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_2;

    parasail_free(pvHMax);
    parasail_free(pvH);
    parasail_free(pvHt);
    parasail_free(pvE);

    return result;
}


