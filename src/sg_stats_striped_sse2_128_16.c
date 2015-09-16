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

#define FASTSTATS

#define NEG_INF (INT16_MIN/(int16_t)(2))

static inline __m128i _mm_blendv_epi8_rpl(__m128i a, __m128i b, __m128i mask) {
    a = _mm_andnot_si128(mask, a);
    a = _mm_or_si128(a, _mm_and_si128(mask, b));
    return a;
}

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
#define FNAME parasail_sg_stats_table_striped_sse2_128_16
#define PNAME parasail_sg_stats_table_striped_profile_sse2_128_16
#define INAME PNAME
#define STATIC
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_stats_rowcol_striped_sse2_128_16
#define PNAME parasail_sg_stats_rowcol_striped_profile_sse2_128_16
#define INAME PNAME
#define STATIC
#else
#define FNAME parasail_sg_stats_striped_sse2_128_16
#ifdef FASTSTATS
#define PNAME parasail_sg_stats_striped_profile_sse2_128_16_internal
#define INAME parasail_sg_stats_striped_profile_sse2_128_16
#define STATIC static
#else
#define PNAME parasail_sg_stats_striped_profile_sse2_128_16
#define INAME PNAME
#define STATIC
#endif
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_sse_128_16(s1, s1Len, matrix);
    parasail_result_t *result = INAME(profile, s2, s2Len, open, gap);
    parasail_profile_free(profile);
    return result;
}

STATIC parasail_result_t* PNAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const int32_t segWidth = 8; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m128i* const restrict vProfile  = (__m128i*)profile->profile16.score;
    __m128i* const restrict vProfileM = (__m128i*)profile->profile16.matches;
    __m128i* const restrict vProfileS = (__m128i*)profile->profile16.similar;
    __m128i* restrict pvHStore        = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHLoad         = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHMStore       = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHMLoad        = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHSStore       = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHSLoad        = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHLStore       = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHLLoad        = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvEStore        = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvELoad         = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvEM      = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvES      = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvEL      = parasail_memalign___m128i(16, segLen);
    __m128i vGapO = _mm_set1_epi16(open);
    __m128i vGapE = _mm_set1_epi16(gap);
    __m128i vNegInf = _mm_set1_epi16(NEG_INF);
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi16(1);
    int16_t score = NEG_INF;
    int16_t matches = NEG_INF;
    int16_t similar = NEG_INF;
    int16_t length = NEG_INF;
    
    __m128i vMaxH = vNegInf;
    __m128i vMaxHM = vNegInf;
    __m128i vMaxHS = vNegInf;
    __m128i vMaxHL = vNegInf;
    __m128i vPosMask = _mm_cmpeq_epi16(_mm_set1_epi16(position),
            _mm_set_epi16(0,1,2,3,4,5,6,7));
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif

    parasail_memset___m128i(pvHMStore, vZero, segLen);
    parasail_memset___m128i(pvHSStore, vZero, segLen);
    parasail_memset___m128i(pvHLStore, vZero, segLen);
    parasail_memset___m128i(pvEM, vZero, segLen);
    parasail_memset___m128i(pvES, vZero, segLen);
    parasail_memset___m128i(pvEL, vZero, segLen);
    parasail_memset___m128i(pvHStore, vZero, segLen);
    parasail_memset___m128i(pvEStore, _mm_set1_epi16(-open), segLen);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vE;
        __m128i vEM;
        __m128i vES;
        __m128i vEL;
        __m128i vF;
        __m128i vFM;
        __m128i vFS;
        __m128i vFL;
        __m128i vH;
        __m128i vHM;
        __m128i vHS;
        __m128i vHL;
        const __m128i* vP = NULL;
        const __m128i* vPM = NULL;
        const __m128i* vPS = NULL;
        __m128i* pv = NULL;

        /* Initialize F value to neg inf.  Any errors to vH values will
         * be corrected in the Lazy_F loop.  */
        vF = vNegInf;
        vFM = vZero;
        vFS = vZero;
        vFL = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = _mm_slli_si128(pvHStore[segLen - 1], 2);
        vHM = _mm_slli_si128(pvHMStore[segLen - 1], 2);
        vHS = _mm_slli_si128(pvHSStore[segLen - 1], 2);
        vHL = _mm_slli_si128(pvHLStore[segLen - 1], 2);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;
        vPM = vProfileM + matrix->mapper[(unsigned char)s2[j]] * segLen;
        vPS = vProfileS + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        pv = pvHMLoad;
        pvHMLoad = pvHMStore;
        pvHMStore = pv;
        pv = pvHSLoad;
        pvHSLoad = pvHSStore;
        pvHSStore = pv;
        pv = pvHLLoad;
        pvHLLoad = pvHLStore;
        pvHLStore = pv;
        pv = pvELoad;
        pvELoad = pvEStore;
        pvEStore = pv;

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            __m128i case1not;
            __m128i case2not;
            __m128i case2;
            __m128i case3;

            vH = _mm_add_epi16(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            case1not = _mm_or_si128(
                    _mm_cmplt_epi16(vH,vF),_mm_cmplt_epi16(vH,vE));
            case2not = _mm_cmplt_epi16(vF,vE);
            case2 = _mm_andnot_si128(case2not,case1not);
            case3 = _mm_and_si128(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);

            /* calculate vM */
            vEM = _mm_load_si128(pvEM + i);
            vHM = _mm_blendv_epi8_rpl(
                    _mm_add_epi16(vHM, _mm_load_si128(vPM + i)),
                    _mm_or_si128(
                        _mm_and_si128(case2, vFM),
                        _mm_and_si128(case3, vEM)),
                    case1not);
            _mm_store_si128(pvHMStore + i, vHM);

            /* calculate vS */
            vES = _mm_load_si128(pvES + i);
            vHS = _mm_blendv_epi8_rpl(
                    _mm_add_epi16(vHS, _mm_load_si128(vPS + i)),
                    _mm_or_si128(
                        _mm_and_si128(case2, vFS),
                        _mm_and_si128(case3, vES)),
                    case1not);
            _mm_store_si128(pvHSStore + i, vHS);

            /* calculate vL */
            vEL = _mm_load_si128(pvEL + i);
            vHL = _mm_blendv_epi8_rpl(
                    _mm_add_epi16(vHL, vOne),
                    _mm_or_si128(
                        _mm_and_si128(case2, _mm_add_epi16(vFL, vOne)),
                        _mm_and_si128(case3, _mm_add_epi16(vEL, vOne))),
                    case1not);
            _mm_store_si128(pvHLStore + i, vHL);
            
#ifdef PARASAIL_TABLE
            arr_store_si128(result->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si128(result->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si128(result->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = _mm_sub_epi16(vH, vGapO);
            vE = _mm_sub_epi16(vE, vGapE);
            vE = _mm_max_epi16(vE, vH);
            _mm_store_si128(pvEStore + i, vE);
            _mm_store_si128(pvEM + i, vHM);
            _mm_store_si128(pvES + i, vHS);
            _mm_store_si128(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm_sub_epi16(vF, vGapE);
            vF = _mm_max_epi16(vF, vH);
            vFM = vHM;
            vFS = vHS;
            vFL = vHL;

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
            vHM = _mm_load_si128(pvHMLoad + i);
            vHS = _mm_load_si128(pvHSLoad + i);
            vHL = _mm_load_si128(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            __m128i vHp = _mm_slli_si128(pvHLoad[segLen - 1], 2);
            vF = _mm_slli_si128(vF, 2);
            vF = _mm_insert_epi16(vF, -open, 0);
            vFM = _mm_slli_si128(vFM, 2);
            vFS = _mm_slli_si128(vFS, 2);
            vFL = _mm_slli_si128(vFL, 2);
            for (i=0; i<segLen; ++i) {
                __m128i case1not;
                __m128i case2not;
                __m128i case2;

                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm_add_epi16(vHp, _mm_load_si128(vP + i));
                vE = _mm_load_si128(pvELoad + i);
                case1not = _mm_or_si128(
                        _mm_cmplt_epi16(vHp,vF),_mm_cmplt_epi16(vHp,vE));
                case2not = _mm_cmplt_epi16(vF,vE);
                case2 = _mm_andnot_si128(case2not,case1not);

                vHM = _mm_load_si128(pvHMStore + i);
                vHM = _mm_blendv_epi8_rpl(vHM, vFM, case2);
                _mm_store_si128(pvHMStore + i, vHM);
                _mm_store_si128(pvEM + i, vHM);

                vHS = _mm_load_si128(pvHSStore + i);
                vHS = _mm_blendv_epi8_rpl(vHS, vFS, case2);
                _mm_store_si128(pvHSStore + i, vHS);
                _mm_store_si128(pvES + i, vHS);

                vHL = _mm_load_si128(pvHLStore + i);
                vHL = _mm_blendv_epi8_rpl(vHL, _mm_add_epi16(vFL,vOne), case2);
                _mm_store_si128(pvHLStore + i, vHL);
                _mm_store_si128(pvEL + i, vHL);

                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi16(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
                arr_store_si128(result->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si128(result->similar_table, vHS, i, segLen, j, s2Len);
                arr_store_si128(result->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si128(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm_sub_epi16(vH, vGapO);
                vF = _mm_sub_epi16(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) goto end;
                /*vF = _mm_max_epi16(vF, vH);*/
                vFM = vHM;
                vFS = vHS;
                vFL = vHL;
                vHp = _mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
            /* extract vector containing last value from the column */
            __m128i cond_max;
            vH = _mm_load_si128(pvHStore + offset);
            vHM = _mm_load_si128(pvHMStore + offset);
            vHS = _mm_load_si128(pvHSStore + offset);
            vHL = _mm_load_si128(pvHLStore + offset);
            cond_max = _mm_cmpgt_epi16(vH, vMaxH);
            vMaxH = _mm_blendv_epi8_rpl(vMaxH, vH, cond_max);
            vMaxHM = _mm_blendv_epi8_rpl(vMaxHM, vHM, cond_max);
            vMaxHS = _mm_blendv_epi8_rpl(vMaxHS, vHS, cond_max);
            vMaxHL = _mm_blendv_epi8_rpl(vMaxHL, vHL, cond_max);
            if (_mm_movemask_epi8(_mm_and_si128(vPosMask, cond_max))) {
                end_ref = j;
                end_query = s1Len - 1;
            }
#ifdef PARASAIL_ROWCOL
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128(vH, 2);
                vHM = _mm_slli_si128(vHM, 2);
                vHS = _mm_slli_si128(vHS, 2);
                vHL = _mm_slli_si128(vHL, 2);
            }
            result->score_row[j] = (int16_t) _mm_extract_epi16 (vH, 7);
            result->matches_row[j] = (int16_t) _mm_extract_epi16 (vHM, 7);
            result->similar_row[j] = (int16_t) _mm_extract_epi16 (vHS, 7);
            result->length_row[j] = (int16_t) _mm_extract_epi16 (vHL, 7);
#endif
        }
    }

    {
        /* extract last value from the column */
        for (k=0; k<position; ++k) {
            vMaxH  = _mm_slli_si128 (vMaxH, 2);
            vMaxHM = _mm_slli_si128 (vMaxHM, 2);
            vMaxHS = _mm_slli_si128 (vMaxHS, 2);
            vMaxHL = _mm_slli_si128 (vMaxHL, 2);
        }
        score = (int16_t) _mm_extract_epi16 (vMaxH, 7);
        matches = (int16_t)_mm_extract_epi16(vMaxHM, 7);
        similar = (int16_t)_mm_extract_epi16(vMaxHS, 7);
        length = (int16_t)_mm_extract_epi16(vMaxHL, 7);
    }

    /* max of last column */
    if (INT32_MAX == profile->stop || 0 == profile->stop)
    {
        int16_t score_last;
        vMaxH = vNegInf;

        if (0 == profile->stop) {
            /* ignore last row contributions */
            score = NEG_INF;
            matches = NEG_INF;
            similar = NEG_INF;
            length = NEG_INF;
            end_query = s1Len;
            end_ref = s2Len - 1;
        }

        for (i=0; i<segLen; ++i) {
            /* load the last stored values */
            __m128i vH = _mm_load_si128(pvHStore + i);
#ifdef PARASAIL_ROWCOL
            __m128i vHM = _mm_load_si128(pvHMStore + i);
            __m128i vHS = _mm_load_si128(pvHSStore + i);
            __m128i vHL = _mm_load_si128(pvHLStore + i);
            arr_store_col(result->score_col, vH, i, segLen);
            arr_store_col(result->matches_col, vHM, i, segLen);
            arr_store_col(result->similar_col, vHS, i, segLen);
            arr_store_col(result->length_col, vHL, i, segLen);
#endif
            vMaxH = _mm_max_epi16(vH, vMaxH);
        }

        /* max in vec */
        score_last = _mm_hmax_epi16_rpl(vMaxH);
        if (score_last > score) {
            end_query = s1Len;
            end_ref = s2Len - 1;
            /* Trace the alignment ending position on read. */
            {
                int16_t *t = (int16_t*)pvHStore;
                int16_t *m = (int16_t*)pvHMStore;
                int16_t *s = (int16_t*)pvHSStore;
                int16_t *l = (int16_t*)pvHLStore;
                int32_t column_len = segLen * segWidth;
                for (i = 0; i<column_len; ++i, ++t, ++m, ++s, ++l) {
                    int32_t temp = i / segWidth + i % segWidth * segLen;
                    if (temp < s1Len) {
                        if (*t > score || (*t == score && temp < end_query)) {
                            score = *t;
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

    

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->matches = matches;
    result->similar = similar;
    result->length = length;

    parasail_free(pvEL);
    parasail_free(pvES);
    parasail_free(pvEM);
    parasail_free(pvELoad);
    parasail_free(pvEStore);
    parasail_free(pvHLLoad);
    parasail_free(pvHLStore);
    parasail_free(pvHSLoad);
    parasail_free(pvHSStore);
    parasail_free(pvHMLoad);
    parasail_free(pvHMStore);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

#ifdef FASTSTATS
#ifdef PARASAIL_TABLE
#else
#ifdef PARASAIL_ROWCOL
#else
#include <assert.h>
parasail_result_t* INAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    const char *s1 = profile->s1;
    const parasail_matrix_t *matrix = profile->matrix;

    /* find the end loc first with the faster implementation */
    parasail_result_t *result = parasail_sg_striped_profile_sse2_128_16(profile, s2, s2Len, open, gap);
    if (!result->saturated) {
        int s1Len_new = 0;
        int s2Len_new = 0;
        parasail_result_t *result_final = NULL;

        /* using the end loc, call the original stats function */
        s1Len_new = result->end_query+1;
        s2Len_new = result->end_ref+1;

        if (s1Len_new == profile->s1Len) {
            /* special 'stop' value tells stats function not to
             * consider last column results */
            int stop_save = profile->stop;
            ((parasail_profile_t*)profile)->stop = 1;
            result_final = PNAME(
                    profile, s2, s2Len_new, open, gap);
            ((parasail_profile_t*)profile)->stop = stop_save;
        }
        else {
            parasail_profile_t *profile_final = NULL;
            profile_final = parasail_profile_create_stats_sse_128_16(
                    s1, s1Len_new, matrix);
            /* special 'stop' value tells stats function not to
             * consider last row results */
            profile_final->stop = 0;
            result_final = PNAME(
                    profile_final, s2, s2Len_new, open, gap);

            parasail_profile_free(profile_final);
        }

        parasail_result_free(result);

        /* correct the end locations before returning */
        result_final->end_query = s1Len_new-1;
        result_final->end_ref = s2Len_new-1;
        return result_final;
    }
    else {
        return result;
    }
}
#endif
#endif
#endif


