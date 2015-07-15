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

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_sse.h"

#define NEG_INF INT8_MIN

static inline __m128i _mm_blendv_epi8_rpl(__m128i a, __m128i b, __m128i mask) {
    a = _mm_andnot_si128(mask, a);
    a = _mm_or_si128(a, _mm_and_si128(mask, b));
    return a;
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

static inline int8_t _mm_hmax_epi8_rpl(__m128i a) {
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 8));
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 4));
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 2));
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 1));
    return _mm_extract_epi8_rpl(a, 0);
}


#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen,
        int32_t bias)
{
    array[( 0*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  0) - bias;
    array[( 1*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  1) - bias;
    array[( 2*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  2) - bias;
    array[( 3*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  3) - bias;
    array[( 4*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  4) - bias;
    array[( 5*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  5) - bias;
    array[( 6*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  6) - bias;
    array[( 7*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  7) - bias;
    array[( 8*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  8) - bias;
    array[( 9*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH,  9) - bias;
    array[(10*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 10) - bias;
    array[(11*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 11) - bias;
    array[(12*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 12) - bias;
    array[(13*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 13) - bias;
    array[(14*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 14) - bias;
    array[(15*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8_rpl(vH, 15) - bias;
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t bias)
{
    col[ 0*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  0) - bias;
    col[ 1*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  1) - bias;
    col[ 2*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  2) - bias;
    col[ 3*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  3) - bias;
    col[ 4*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  4) - bias;
    col[ 5*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  5) - bias;
    col[ 6*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  6) - bias;
    col[ 7*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  7) - bias;
    col[ 8*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  8) - bias;
    col[ 9*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH,  9) - bias;
    col[10*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 10) - bias;
    col[11*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 11) - bias;
    col[12*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 12) - bias;
    col[13*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 13) - bias;
    col[14*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 14) - bias;
    col[15*seglen+t] = (int8_t)_mm_extract_epi8_rpl(vH, 15) - bias;
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_stats_table_striped_sse2_128_8
#define PNAME parasail_sw_stats_table_striped_profile_sse2_128_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_stats_rowcol_striped_sse2_128_8
#define PNAME parasail_sw_stats_rowcol_striped_profile_sse2_128_8
#else
#define FNAME parasail_sw_stats_striped_sse2_128_8
#define PNAME parasail_sw_stats_striped_profile_sse2_128_8
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
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m128i* const restrict vProfile  = (__m128i*)profile->profile8.score;
    __m128i* const restrict vProfileM = (__m128i*)profile->profile8.matches;
    __m128i* const restrict vProfileS = (__m128i*)profile->profile8.similar;
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
    __m128i* restrict pvHMax          = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHMMax          = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHSMax          = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHLMax          = parasail_memalign___m128i(16, segLen);
    __m128i vGapO = _mm_set1_epi8(open);
    __m128i vGapE = _mm_set1_epi8(gap);
    __m128i vZero = _mm_set1_epi8(0);
    __m128i vOne = _mm_set1_epi8(1);
    int8_t bias = INT8_MIN;
    int8_t score = bias;
    int8_t matches = bias;
    int8_t similar = bias;
    int8_t length = bias;
    __m128i vBias = _mm_set1_epi8(bias);
    __m128i vMaxH = vBias;
    __m128i vMaxHUnit = vBias;
    __m128i insert_mask = _mm_cmpgt_epi8(
            _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
            vZero);
    __m128i vSaturationCheckMax = vBias;
    __m128i vPosLimit = _mm_set1_epi8(INT8_MAX);
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

    parasail_memset___m128i(pvHMStore, vBias, segLen);
    parasail_memset___m128i(pvHSStore, vBias, segLen);
    parasail_memset___m128i(pvHLStore, vBias, segLen);

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m128i_8_t h;
            __m128i_8_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = bias;
                e.v[segNum] = bias;
            }
            _mm_store_si128(&pvHStore[index], h.m);
            _mm_store_si128(&pvEStore[index], e.m);
            ++index;
        }
    }

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

        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        vF = vBias;
        vFM = vBias;
        vFS = vBias;
        vFL = vBias;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = _mm_slli_si128(pvHStore[segLen - 1], 1);
        vHM = _mm_slli_si128(pvHMStore[segLen - 1], 1);
        vHS = _mm_slli_si128(pvHSStore[segLen - 1], 1);
        vHL = _mm_slli_si128(pvHLStore[segLen - 1], 1);
        vH = _mm_blendv_epi8_rpl(vH, vBias, insert_mask);
        vHM = _mm_blendv_epi8_rpl(vHM, vBias, insert_mask);
        vHS = _mm_blendv_epi8_rpl(vHS, vBias, insert_mask);
        vHL = _mm_blendv_epi8_rpl(vHL, vBias, insert_mask);

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
            __m128i cond_zero;

            vH = _mm_adds_epi8(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            case1not = _mm_or_si128(
                    _mm_cmplt_epi8(vH,vF),_mm_cmplt_epi8(vH,vE));
            case2not = _mm_cmplt_epi8(vF,vE);
            case2 = _mm_andnot_si128(case2not,case1not);
            case3 = _mm_and_si128(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi8_rpl(vH, vE);
            vH = _mm_max_epi8_rpl(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);
            cond_zero = _mm_cmpeq_epi8(vH, vBias);

            /* calculate vM */
            vEM = _mm_load_si128(pvEM + i);
            vHM = _mm_blendv_epi8_rpl(
                    _mm_adds_epi8(vHM, _mm_load_si128(vPM + i)),
                    _mm_or_si128(
                        _mm_and_si128(case2, vFM),
                        _mm_and_si128(case3, vEM)),
                    case1not);
            vHM = _mm_blendv_epi8_rpl(vHM, vBias, cond_zero);
            _mm_store_si128(pvHMStore + i, vHM);

            /* calculate vS */
            vES = _mm_load_si128(pvES + i);
            vHS = _mm_blendv_epi8_rpl(
                    _mm_adds_epi8(vHS, _mm_load_si128(vPS + i)),
                    _mm_or_si128(
                        _mm_and_si128(case2, vFS),
                        _mm_and_si128(case3, vES)),
                    case1not);
            vHS = _mm_blendv_epi8_rpl(vHS, vBias, cond_zero);
            _mm_store_si128(pvHSStore + i, vHS);

            /* calculate vL */
            vEL = _mm_load_si128(pvEL + i);
            vHL = _mm_blendv_epi8_rpl(
                    _mm_adds_epi8(vHL, vOne),
                    _mm_or_si128(
                        _mm_and_si128(case2, _mm_adds_epi8(vFL, vOne)),
                        _mm_and_si128(case3, _mm_adds_epi8(vEL, vOne))),
                    case1not);
            vHL = _mm_blendv_epi8_rpl(vHL, vBias, cond_zero);
            _mm_store_si128(pvHLStore + i, vHL);
            vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHM);
            vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHS);
            vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->matches_table, vHM, i, segLen, j, s2Len, bias);
            arr_store_si128(result->similar_table, vHS, i, segLen, j, s2Len, bias);
            arr_store_si128(result->length_table, vHL, i, segLen, j, s2Len, bias);
            arr_store_si128(result->score_table, vH, i, segLen, j, s2Len, bias);
#endif
            vMaxH = _mm_max_epi8_rpl(vH, vMaxH);

            /* Update vE value. */
            vH = _mm_subs_epi8(vH, vGapO);
            vE = _mm_subs_epi8(vE, vGapE);
            vE = _mm_max_epi8_rpl(vE, vH);
            _mm_store_si128(pvEStore + i, vE);
            _mm_store_si128(pvEM + i, vHM);
            _mm_store_si128(pvES + i, vHS);
            _mm_store_si128(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm_subs_epi8(vF, vGapE);
            vF = _mm_max_epi8_rpl(vF, vH);
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
            __m128i vHp = _mm_slli_si128(pvHLoad[segLen - 1], 1);
            vF = _mm_slli_si128(vF, 1);
            vFM = _mm_slli_si128(vFM, 1);
            vFS = _mm_slli_si128(vFS, 1);
            vFL = _mm_slli_si128(vFL, 1);
            vHp = _mm_blendv_epi8_rpl(vHp, vBias, insert_mask);
            vF = _mm_blendv_epi8_rpl(vF, vBias, insert_mask);
            vFM = _mm_blendv_epi8_rpl(vFM, vBias, insert_mask);
            vFS = _mm_blendv_epi8_rpl(vFS, vBias, insert_mask);
            vFL = _mm_blendv_epi8_rpl(vFL, vBias, insert_mask);
            for (i=0; i<segLen; ++i) {
                __m128i case1not;
                __m128i case2not;
                __m128i case2;
                __m128i cond_zero;

                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm_adds_epi8(vHp, _mm_load_si128(vP + i));
                vE = _mm_load_si128(pvELoad + i);
                case1not = _mm_or_si128(
                        _mm_cmplt_epi8(vHp,vF),_mm_cmplt_epi8(vHp,vE));
                case2not = _mm_cmplt_epi8(vF,vE);
                case2 = _mm_andnot_si128(case2not,case1not);

                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi8_rpl(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
                cond_zero = _mm_cmpeq_epi8(vH, vBias);

                vHM = _mm_load_si128(pvHMStore + i);
                vHM = _mm_blendv_epi8_rpl(vHM, vFM, case2);
                vHM = _mm_blendv_epi8_rpl(vHM, vBias, cond_zero);
                _mm_store_si128(pvHMStore + i, vHM);
                _mm_store_si128(pvEM + i, vHM);

                vHS = _mm_load_si128(pvHSStore + i);
                vHS = _mm_blendv_epi8_rpl(vHS, vFS, case2);
                vHS = _mm_blendv_epi8_rpl(vHS, vBias, cond_zero);
                _mm_store_si128(pvHSStore + i, vHS);
                _mm_store_si128(pvES + i, vHS);

                vHL = _mm_load_si128(pvHLStore + i);
                vHL = _mm_blendv_epi8_rpl(vHL, _mm_adds_epi8(vFL,vOne), case2);
                vHL = _mm_blendv_epi8_rpl(vHL, vBias, cond_zero);
                _mm_store_si128(pvHLStore + i, vHL);
                _mm_store_si128(pvEL + i, vHL);

                vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHM);
                vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHS);
                vSaturationCheckMax = _mm_max_epi8_rpl(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
                arr_store_si128(result->matches_table, vHM, i, segLen, j, s2Len, bias);
                arr_store_si128(result->similar_table, vHS, i, segLen, j, s2Len, bias);
                arr_store_si128(result->length_table, vHL, i, segLen, j, s2Len, bias);
                arr_store_si128(result->score_table, vH, i, segLen, j, s2Len, bias);
#endif
                vMaxH = _mm_max_epi8_rpl(vH, vMaxH);
                vH = _mm_subs_epi8(vH, vGapO);
                vF = _mm_subs_epi8(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi8(vF, vH))) goto end;
                /*vF = _mm_max_epi8_rpl(vF, vH);*/
                vFM = vHM;
                vFS = vHS;
                vFL = vHL;
                vHp = _mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
        }

        {
            __m128i vCompare = _mm_cmpgt_epi8(vMaxH, vMaxHUnit);
            if (_mm_movemask_epi8(vCompare)) {
                score = _mm_hmax_epi8_rpl(vMaxH);
                vMaxHUnit = _mm_set1_epi8(score);
                end_ref = j;
                (void)memcpy(pvHMax, pvHStore, sizeof(__m128i)*segLen);
                (void)memcpy(pvHMMax, pvHMStore, sizeof(__m128i)*segLen);
                (void)memcpy(pvHSMax, pvHSStore, sizeof(__m128i)*segLen);
                (void)memcpy(pvHLMax, pvHLStore, sizeof(__m128i)*segLen);
            }
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = _mm_load_si128(pvHStore + offset);
            vHM = _mm_load_si128(pvHMStore + offset);
            vHS = _mm_load_si128(pvHSStore + offset);
            vHL = _mm_load_si128(pvHLStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128(vH, 1);
                vHM = _mm_slli_si128(vHM, 1);
                vHS = _mm_slli_si128(vHS, 1);
                vHL = _mm_slli_si128(vHL, 1);
            }
            result->score_row[j] = (int8_t) _mm_extract_epi8_rpl (vH, 15) - bias;
            result->matches_row[j] = (int8_t) _mm_extract_epi8_rpl (vHM, 15) - bias;
            result->similar_row[j] = (int8_t) _mm_extract_epi8_rpl (vHS, 15) - bias;
            result->length_row[j] = (int8_t) _mm_extract_epi8_rpl (vHL, 15) - bias;
        }
#endif
    }

    /* Trace the alignment ending position on read. */
    {
        int8_t *t = (int8_t*)pvHMax;
        int8_t *m = (int8_t*)pvHMMax;
        int8_t *s = (int8_t*)pvHSMax;
        int8_t *l = (int8_t*)pvHLMax;
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
        __m128i vH = _mm_load_si128(pvHStore+i);
        __m128i vHM = _mm_load_si128(pvHMStore+i);
        __m128i vHS = _mm_load_si128(pvHSStore+i);
        __m128i vHL = _mm_load_si128(pvHLStore+i);
        arr_store_col(result->score_col, vH, i, segLen, bias);
        arr_store_col(result->matches_col, vHM, i, segLen, bias);
        arr_store_col(result->similar_col, vHS, i, segLen, bias);
        arr_store_col(result->length_col, vHL, i, segLen, bias);
    }
#endif

    if (score == INT8_MAX
            || _mm_movemask_epi8(_mm_cmpeq_epi8(vSaturationCheckMax,vPosLimit))) {
        result->saturated = 1;
        score = INT8_MAX;
        matches = INT8_MIN;
        similar = INT8_MIN;
        length = INT8_MIN;
    }

    result->score = score - bias;
    result->matches = matches - bias;
    result->similar = similar - bias;
    result->length = length - bias;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(pvHLMax);
    parasail_free(pvHSMax);
    parasail_free(pvHMMax);
    parasail_free(pvHMax);
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

