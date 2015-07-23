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

#define NEG_INF (INT16_MIN/(int16_t)(2))
#define MAX(a,b) ((a)>(b)?(a):(b))

#define _mm_rlli_si128_rpl(a,imm) _mm_alignr_epi8(a, a, 16-imm)

static inline int16_t _mm_hmax_epi16_rpl(__m128i a) {
    a = _mm_max_epi16(a, _mm_srli_si128(a, 8));
    a = _mm_max_epi16(a, _mm_srli_si128(a, 4));
    a = _mm_max_epi16(a, _mm_srli_si128(a, 2));
    return _mm_extract_epi16(a, 0);
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
    array[(0*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 0);
    array[(1*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 1);
    array[(2*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 2);
    array[(3*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 3);
    array[(4*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 4);
    array[(5*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 5);
    array[(6*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 6);
    array[(7*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 7);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m128i vH,
        int32_t t,
        int32_t seglen)
{
    col[0*seglen+t] = (int16_t)_mm_extract_epi16(vH, 0);
    col[1*seglen+t] = (int16_t)_mm_extract_epi16(vH, 1);
    col[2*seglen+t] = (int16_t)_mm_extract_epi16(vH, 2);
    col[3*seglen+t] = (int16_t)_mm_extract_epi16(vH, 3);
    col[4*seglen+t] = (int16_t)_mm_extract_epi16(vH, 4);
    col[5*seglen+t] = (int16_t)_mm_extract_epi16(vH, 5);
    col[6*seglen+t] = (int16_t)_mm_extract_epi16(vH, 6);
    col[7*seglen+t] = (int16_t)_mm_extract_epi16(vH, 7);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_stats_table_scan_sse41_128_16
#define PNAME parasail_sw_stats_table_scan_profile_sse41_128_16
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_stats_rowcol_scan_sse41_128_16
#define PNAME parasail_sw_stats_rowcol_scan_profile_sse41_128_16
#else
#define FNAME parasail_sw_stats_scan_sse41_128_16
#define PNAME parasail_sw_stats_scan_profile_sse41_128_16
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_sse_128_16(s1, s1Len, matrix);
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
    const int32_t segWidth = 8;
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m128i* const restrict pvP  = (__m128i*)profile->profile16.score;
    __m128i* const restrict pvPm = (__m128i*)profile->profile16.matches;
    __m128i* const restrict pvPs = (__m128i*)profile->profile16.similar;
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
    __m128i* const restrict pvHMax = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHMMax = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHSMax = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvHLMax = parasail_memalign___m128i(16, segLen);
    __m128i vGapO = _mm_set1_epi16(open);
    __m128i vGapE = _mm_set1_epi16(gap);
    __m128i vNegInf = _mm_set1_epi16(NEG_INF);
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi16(1);
    int16_t score = NEG_INF;
    int16_t matches = 0;
    int16_t similar = 0;
    int16_t length = 0;
    __m128i vMaxH = vNegInf;
    __m128i vMaxM = vNegInf;
    __m128i vMaxS = vNegInf;
    __m128i vMaxL = vNegInf;
    __m128i vMaxHUnit = vNegInf;
    const int16_t segLenXgap = -segLen*gap;
    __m128i insert_mask = _mm_cmpeq_epi16(vZero,
            _mm_set_epi16(0,0,0,0,0,0,0,1));
    __m128i vSegLenXgap_reset = _mm_blendv_epi8(vNegInf,
            _mm_set1_epi16(segLenXgap),
            insert_mask);
    
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
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
            __m128i_16_t h;
            __m128i_16_t e;
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
            vE = _mm_max_epi16(
                    _mm_sub_epi16(vE, vGapE),
                    _mm_sub_epi16(vH, vGapO));
            _mm_store_si128(pvE+i, vE);
        }

        /* calculate Ht */
        vH = _mm_slli_si128(_mm_load_si128(pvH+(segLen-1)), 2);
        vMp= _mm_slli_si128(_mm_load_si128(pvM+(segLen-1)), 2);
        vSp= _mm_slli_si128(_mm_load_si128(pvS+(segLen-1)), 2);
        vLp= _mm_slli_si128(_mm_load_si128(pvL+(segLen-1)), 2);
        vLp= _mm_add_epi16(vLp, vOne);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvC = pvPm+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvD = pvPs+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            __m128i cond_max;
            /* load values we need */
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            /* compute */
            vH = _mm_add_epi16(vH, vW);
            vHt = _mm_max_epi16(vH, vE);
            vHt = _mm_max_epi16(vHt, vZero);
            /* statistics */
            vC = _mm_load_si128(pvC+i);
            vD = _mm_load_si128(pvD+i);
            vMp = _mm_add_epi16(vMp, vC);
            vSp = _mm_add_epi16(vSp, vD);
            vEx = _mm_cmpgt_epi16(vE, vH);
            vM = _mm_load_si128(pvM+i);
            vS = _mm_load_si128(pvS+i);
            vL = _mm_load_si128(pvL+i);
            vL = _mm_add_epi16(vL, vOne);
            vMt = _mm_blendv_epi8(vMp, vM, vEx);
            vSt = _mm_blendv_epi8(vSp, vS, vEx);
            vLt = _mm_blendv_epi8(vLp, vL, vEx);
            cond_max = _mm_cmpeq_epi16(vHt, vZero);
            vEx = _mm_andnot_si128(cond_max, vEx);
            vMt = _mm_andnot_si128(cond_max, vMt);
            vSt = _mm_andnot_si128(cond_max, vSt);
            vLt = _mm_andnot_si128(cond_max, vLt);
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
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 2);
        vFt = vNegInf;
        for (i=0; i<segLen; ++i) {
            vFt = _mm_sub_epi16(vFt, vGapE);
            vFt = _mm_max_epi16(vFt, vHt);
            vHt = _mm_load_si128(pvHt+i);
        }
        {
            __m128i vFt_save = vFt;
            __m128i segLenXgap = vSegLenXgap_reset;
            for (i=0; i<segWidth-1; ++i) {
                __m128i vFtt = _mm_slli_si128(vFt, 2);
                segLenXgap = _mm_rlli_si128_rpl(segLenXgap, 2);
                vFtt = _mm_add_epi16(vFtt, segLenXgap);
                vFt = _mm_max_epi16(vFt, vFtt);
            }
            vFt = _mm_blendv_epi8(vFt_save, vFt, insert_mask);
        }
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 2);
        vFt = _mm_slli_si128(vFt, 2);
        vFt = _mm_insert_epi16(vFt, NEG_INF, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_sub_epi16(vFt, vGapE);
            vFt = _mm_max_epi16(vFt, vHt);
            vHt = _mm_load_si128(pvHt+i);
            _mm_store_si128(pvFt+i, vFt);
        }

        /* calculate H,M,L */
        vMp = vZero;
        vSp = vZero;
        vLp = vOne;
        vC = _mm_cmpeq_epi16(vZero, vZero); /* check if prefix sum is needed */
        vC = _mm_srli_si128(vC, 2); /* zero out last value */
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vHt = _mm_load_si128(pvHt+i);
            vFt = _mm_load_si128(pvFt+i);
            /* compute */
            vFt = _mm_sub_epi16(vFt, vGapO);
            vH = _mm_max_epi16(vHt, vFt);
            /* statistics */
            vEx = _mm_load_si128(pvEx+i);
            vMt = _mm_load_si128(pvMt+i);
            vSt = _mm_load_si128(pvSt+i);
            vLt = _mm_load_si128(pvLt+i);
            vEx = _mm_or_si128(
                    _mm_and_si128(vEx, _mm_cmpeq_epi16(vHt, vFt)),
                    _mm_cmplt_epi16(vHt, vFt));
            vM = _mm_blendv_epi8(vMt, vMp, vEx);
            vS = _mm_blendv_epi8(vSt, vSp, vEx);
            vL = _mm_blendv_epi8(vLt, vLp, vEx);
            vMp = vM;
            vSp = vS;
            vLp = _mm_add_epi16(vL, vOne);
            vC = _mm_and_si128(vC, vEx);
            /* store results */
            _mm_store_si128(pvH+i, vH);
            _mm_store_si128(pvEx+i, vEx);
            
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif
        }
        {
            vLp = _mm_sub_epi16(vLp, vOne);
            {
                __m128i_16_t uMp, uSp, uLp, uC;
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
            vLp = _mm_add_epi16(vLp, vOne);
        }
        /* final pass for M,L */
        vMp = _mm_slli_si128(vMp, 2);
        vSp = _mm_slli_si128(vSp, 2);
        vLp = _mm_slli_si128(vLp, 2);
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vH = _mm_load_si128(pvH+i);
            /* statistics */
            vEx = _mm_load_si128(pvEx+i);
            vMt = _mm_load_si128(pvMt+i);
            vSt = _mm_load_si128(pvSt+i);
            vLt = _mm_load_si128(pvLt+i);
            vM = _mm_blendv_epi8(vMt, vMp, vEx);
            vS = _mm_blendv_epi8(vSt, vSp, vEx);
            vL = _mm_blendv_epi8(vLt, vLp, vEx);
            vMp = vM;
            vSp = vS;
            vLp = _mm_add_epi16(vL, vOne);
            /* store results */
            _mm_store_si128(pvM+i, vM);
            _mm_store_si128(pvS+i, vS);
            _mm_store_si128(pvL+i, vL);
            
#ifdef PARASAIL_TABLE
            arr_store_si128(result->matches_table, vM, i, segLen, j, s2Len);
            arr_store_si128(result->similar_table, vS, i, segLen, j, s2Len);
            arr_store_si128(result->length_table, vL, i, segLen, j, s2Len);
#endif
            /* update max vector seen so far */
            {
                __m128i cond_max = _mm_cmpgt_epi16(vH, vMaxH);
                vMaxH = _mm_blendv_epi8(vMaxH, vH, cond_max);
                vMaxM = _mm_blendv_epi8(vMaxM, vM, cond_max);
                vMaxS = _mm_blendv_epi8(vMaxS, vS, cond_max);
                vMaxL = _mm_blendv_epi8(vMaxL, vL, cond_max);
            }
        }

        {
            __m128i vCompare = _mm_cmpgt_epi16(vMaxH, vMaxHUnit);
            if (_mm_movemask_epi8(vCompare)) {
                score = _mm_hmax_epi16_rpl(vMaxH);
                vMaxHUnit = _mm_set1_epi16(score);
                end_ref = j;
                (void)memcpy(pvHMax, pvH, sizeof(__m128i)*segLen);
                (void)memcpy(pvHMMax, pvM, sizeof(__m128i)*segLen);
                (void)memcpy(pvHSMax, pvS, sizeof(__m128i)*segLen);
                (void)memcpy(pvHLMax, pvL, sizeof(__m128i)*segLen);
            }
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            int32_t k = 0;
            vH = _mm_load_si128(pvH + offset);
            vM = _mm_load_si128(pvM + offset);
            vS = _mm_load_si128(pvS + offset);
            vL = _mm_load_si128(pvL + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128(vH, 2);
                vM = _mm_slli_si128(vM, 2);
                vS = _mm_slli_si128(vS, 2);
                vL = _mm_slli_si128(vL, 2);
            }
            result->score_row[j] = (int16_t) _mm_extract_epi16 (vH, 7);
            result->matches_row[j] = (int16_t) _mm_extract_epi16 (vM, 7);
            result->similar_row[j] = (int16_t) _mm_extract_epi16 (vS, 7);
            result->length_row[j] = (int16_t) _mm_extract_epi16 (vL, 7);
        }
#endif
    }

    /* Trace the alignment ending position on read. */
    {
        int16_t *t = (int16_t*)pvHMax;
        int16_t *m = (int16_t*)pvHMMax;
        int16_t *s = (int16_t*)pvHSMax;
        int16_t *l = (int16_t*)pvHLMax;
        int32_t column_len = segLen * segWidth;
        end_query = s1Len;
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

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        __m128i vH = _mm_load_si128(pvH+i);
        __m128i vM = _mm_load_si128(pvM+i);
        __m128i vS = _mm_load_si128(pvS+i);
        __m128i vL = _mm_load_si128(pvL+i);
        arr_store_col(result->score_col, vH, i, segLen);
        arr_store_col(result->matches_col, vM, i, segLen);
        arr_store_col(result->similar_col, vS, i, segLen);
        arr_store_col(result->length_col, vL, i, segLen);
    }
#endif

    

    result->score = score;
    result->matches = matches;
    result->similar = similar;
    result->length = length;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(pvHLMax);
    parasail_free(pvHSMax);
    parasail_free(pvHMMax);
    parasail_free(pvHMax);
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


