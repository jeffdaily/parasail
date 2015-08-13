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

#include <immintrin.h>

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_avx.h"

#define FASTSTATS

#define NEG_INF (INT16_MIN/(int16_t)(2))

#define _mm256_cmplt_epi16_rpl(a,b) _mm256_cmpgt_epi16(b,a)

#if HAVE_AVX2_MM256_EXTRACT_EPI16
#define _mm256_extract_epi16_rpl _mm256_extract_epi16
#else
static inline int16_t _mm256_extract_epi16_rpl(__m256i a, int imm) {
    __m256i_16_t A;
    A.m = a;
    return A.v[imm];
}
#endif

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)

static inline int16_t _mm256_hmax_epi16_rpl(__m256i a) {
    a = _mm256_max_epi16(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
    a = _mm256_max_epi16(a, _mm256_slli_si256(a, 8));
    a = _mm256_max_epi16(a, _mm256_slli_si256(a, 4));
    a = _mm256_max_epi16(a, _mm256_slli_si256(a, 2));
    return _mm256_extract_epi16_rpl(a, 15);
}


#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen,
        int32_t bias)
{
    array[( 0*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  0) - bias;
    array[( 1*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  1) - bias;
    array[( 2*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  2) - bias;
    array[( 3*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  3) - bias;
    array[( 4*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  4) - bias;
    array[( 5*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  5) - bias;
    array[( 6*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  6) - bias;
    array[( 7*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  7) - bias;
    array[( 8*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  8) - bias;
    array[( 9*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  9) - bias;
    array[(10*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 10) - bias;
    array[(11*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 11) - bias;
    array[(12*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 12) - bias;
    array[(13*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 13) - bias;
    array[(14*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 14) - bias;
    array[(15*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 15) - bias;
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t bias)
{
    col[ 0*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  0) - bias;
    col[ 1*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  1) - bias;
    col[ 2*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  2) - bias;
    col[ 3*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  3) - bias;
    col[ 4*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  4) - bias;
    col[ 5*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  5) - bias;
    col[ 6*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  6) - bias;
    col[ 7*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  7) - bias;
    col[ 8*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  8) - bias;
    col[ 9*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  9) - bias;
    col[10*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 10) - bias;
    col[11*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 11) - bias;
    col[12*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 12) - bias;
    col[13*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 13) - bias;
    col[14*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 14) - bias;
    col[15*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 15) - bias;
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_stats_table_striped_avx2_256_16
#define PNAME parasail_sw_stats_table_striped_profile_avx2_256_16
#define INAME PNAME
#define STATIC
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_stats_rowcol_striped_avx2_256_16
#define PNAME parasail_sw_stats_rowcol_striped_profile_avx2_256_16
#define INAME PNAME
#define STATIC
#else
#define FNAME parasail_sw_stats_striped_avx2_256_16
#ifdef FASTSTATS
#define PNAME parasail_sw_stats_striped_profile_avx2_256_16_internal
#define INAME parasail_sw_stats_striped_profile_avx2_256_16
#define STATIC static
#else
#define PNAME parasail_sw_stats_striped_profile_avx2_256_16
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
    parasail_profile_t *profile = parasail_profile_create_stats_avx_256_16(s1, s1Len, matrix);
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
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m256i* const restrict vProfile  = (__m256i*)profile->profile16.score;
    __m256i* const restrict vProfileM = (__m256i*)profile->profile16.matches;
    __m256i* const restrict vProfileS = (__m256i*)profile->profile16.similar;
    __m256i* restrict pvHStore        = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLoad         = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHMStore       = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHMLoad        = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHSStore       = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHSLoad        = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLStore       = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLLoad        = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvEStore        = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvELoad         = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvEM      = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvES      = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvEL      = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHMax          = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHMMax          = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHSMax          = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLMax          = parasail_memalign___m256i(32, segLen);
    __m256i vGapO = _mm256_set1_epi16(open);
    __m256i vGapE = _mm256_set1_epi16(gap);
    __m256i vZero = _mm256_set1_epi16(0);
    __m256i vOne = _mm256_set1_epi16(1);
    int16_t bias = INT16_MIN;
    int16_t score = bias;
    int16_t matches = bias;
    int16_t similar = bias;
    int16_t length = bias;
    __m256i vBias = _mm256_set1_epi16(bias);
    __m256i vMaxH = vBias;
    __m256i vMaxHUnit = vBias;
    __m256i insert_mask = _mm256_cmpgt_epi16(
            _mm256_set_epi16(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
            vZero);
    __m256i vSaturationCheckMax = vBias;
    __m256i vPosLimit = _mm256_set1_epi16(INT16_MAX);
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

    parasail_memset___m256i(pvHMStore, vBias, segLen);
    parasail_memset___m256i(pvHSStore, vBias, segLen);
    parasail_memset___m256i(pvHLStore, vBias, segLen);
    parasail_memset___m256i(pvEM, vBias, segLen);
    parasail_memset___m256i(pvES, vBias, segLen);
    parasail_memset___m256i(pvEL, vBias, segLen);
    parasail_memset___m256i(pvHStore, vBias, segLen);
    parasail_memset___m256i(pvEStore, vBias, segLen);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m256i vE;
        __m256i vEM;
        __m256i vES;
        __m256i vEL;
        __m256i vF;
        __m256i vFM;
        __m256i vFS;
        __m256i vFL;
        __m256i vH;
        __m256i vHM;
        __m256i vHS;
        __m256i vHL;
        const __m256i* vP = NULL;
        const __m256i* vPM = NULL;
        const __m256i* vPS = NULL;
        __m256i* pv = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        vF = vBias;
        vFM = vBias;
        vFS = vBias;
        vFL = vBias;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = _mm256_slli_si256_rpl(pvHStore[segLen - 1], 2);
        vHM = _mm256_slli_si256_rpl(pvHMStore[segLen - 1], 2);
        vHS = _mm256_slli_si256_rpl(pvHSStore[segLen - 1], 2);
        vHL = _mm256_slli_si256_rpl(pvHLStore[segLen - 1], 2);
        vH = _mm256_blendv_epi8(vH, vBias, insert_mask);
        vHM = _mm256_blendv_epi8(vHM, vBias, insert_mask);
        vHS = _mm256_blendv_epi8(vHS, vBias, insert_mask);
        vHL = _mm256_blendv_epi8(vHL, vBias, insert_mask);

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
            __m256i case1not;
            __m256i case2not;
            __m256i case2;
            __m256i case3;
            __m256i cond_zero;

            vH = _mm256_adds_epi16(vH, _mm256_load_si256(vP + i));
            vE = _mm256_load_si256(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            case1not = _mm256_or_si256(
                    _mm256_cmplt_epi16_rpl(vH,vF),_mm256_cmplt_epi16_rpl(vH,vE));
            case2not = _mm256_cmplt_epi16_rpl(vF,vE);
            case2 = _mm256_andnot_si256(case2not,case1not);
            case3 = _mm256_and_si256(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm256_max_epi16(vH, vE);
            vH = _mm256_max_epi16(vH, vF);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);
            cond_zero = _mm256_cmpeq_epi16(vH, vBias);

            /* calculate vM */
            vEM = _mm256_load_si256(pvEM + i);
            vHM = _mm256_blendv_epi8(
                    _mm256_adds_epi16(vHM, _mm256_load_si256(vPM + i)),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, vFM),
                        _mm256_and_si256(case3, vEM)),
                    case1not);
            vHM = _mm256_blendv_epi8(vHM, vBias, cond_zero);
            _mm256_store_si256(pvHMStore + i, vHM);

            /* calculate vS */
            vES = _mm256_load_si256(pvES + i);
            vHS = _mm256_blendv_epi8(
                    _mm256_adds_epi16(vHS, _mm256_load_si256(vPS + i)),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, vFS),
                        _mm256_and_si256(case3, vES)),
                    case1not);
            vHS = _mm256_blendv_epi8(vHS, vBias, cond_zero);
            _mm256_store_si256(pvHSStore + i, vHS);

            /* calculate vL */
            vEL = _mm256_load_si256(pvEL + i);
            vHL = _mm256_blendv_epi8(
                    _mm256_adds_epi16(vHL, vOne),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, _mm256_adds_epi16(vFL, vOne)),
                        _mm256_and_si256(case3, _mm256_adds_epi16(vEL, vOne))),
                    case1not);
            vHL = _mm256_blendv_epi8(vHL, vBias, cond_zero);
            _mm256_store_si256(pvHLStore + i, vHL);
            vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vHM);
            vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vHS);
            vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->matches_table, vHM, i, segLen, j, s2Len, bias);
            arr_store_si256(result->similar_table, vHS, i, segLen, j, s2Len, bias);
            arr_store_si256(result->length_table, vHL, i, segLen, j, s2Len, bias);
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len, bias);
#endif
            vMaxH = _mm256_max_epi16(vH, vMaxH);

            /* Update vE value. */
            vH = _mm256_subs_epi16(vH, vGapO);
            vE = _mm256_subs_epi16(vE, vGapE);
            vE = _mm256_max_epi16(vE, vH);
            _mm256_store_si256(pvEStore + i, vE);
            _mm256_store_si256(pvEM + i, vHM);
            _mm256_store_si256(pvES + i, vHS);
            _mm256_store_si256(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm256_subs_epi16(vF, vGapE);
            vF = _mm256_max_epi16(vF, vH);
            vFM = vHM;
            vFS = vHS;
            vFL = vHL;

            /* Load the next vH. */
            vH = _mm256_load_si256(pvHLoad + i);
            vHM = _mm256_load_si256(pvHMLoad + i);
            vHS = _mm256_load_si256(pvHSLoad + i);
            vHL = _mm256_load_si256(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            __m256i vHp = _mm256_slli_si256_rpl(pvHLoad[segLen - 1], 2);
            vF = _mm256_slli_si256_rpl(vF, 2);
            vFM = _mm256_slli_si256_rpl(vFM, 2);
            vFS = _mm256_slli_si256_rpl(vFS, 2);
            vFL = _mm256_slli_si256_rpl(vFL, 2);
            vHp = _mm256_blendv_epi8(vHp, vBias, insert_mask);
            vF = _mm256_blendv_epi8(vF, vBias, insert_mask);
            vFM = _mm256_blendv_epi8(vFM, vBias, insert_mask);
            vFS = _mm256_blendv_epi8(vFS, vBias, insert_mask);
            vFL = _mm256_blendv_epi8(vFL, vBias, insert_mask);
            for (i=0; i<segLen; ++i) {
                __m256i case1not;
                __m256i case2not;
                __m256i case2;
                __m256i cond_zero;

                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm256_adds_epi16(vHp, _mm256_load_si256(vP + i));
                vE = _mm256_load_si256(pvELoad + i);
                case1not = _mm256_or_si256(
                        _mm256_cmplt_epi16_rpl(vHp,vF),_mm256_cmplt_epi16_rpl(vHp,vE));
                case2not = _mm256_cmplt_epi16_rpl(vF,vE);
                case2 = _mm256_andnot_si256(case2not,case1not);

                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi16(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
                cond_zero = _mm256_cmpeq_epi16(vH, vBias);

                vHM = _mm256_load_si256(pvHMStore + i);
                vHM = _mm256_blendv_epi8(vHM, vFM, case2);
                vHM = _mm256_blendv_epi8(vHM, vBias, cond_zero);
                _mm256_store_si256(pvHMStore + i, vHM);
                _mm256_store_si256(pvEM + i, vHM);

                vHS = _mm256_load_si256(pvHSStore + i);
                vHS = _mm256_blendv_epi8(vHS, vFS, case2);
                vHS = _mm256_blendv_epi8(vHS, vBias, cond_zero);
                _mm256_store_si256(pvHSStore + i, vHS);
                _mm256_store_si256(pvES + i, vHS);

                vHL = _mm256_load_si256(pvHLStore + i);
                vHL = _mm256_blendv_epi8(vHL, _mm256_adds_epi16(vFL,vOne), case2);
                vHL = _mm256_blendv_epi8(vHL, vBias, cond_zero);
                _mm256_store_si256(pvHLStore + i, vHL);
                _mm256_store_si256(pvEL + i, vHL);

                vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vHM);
                vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vHS);
                vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
                arr_store_si256(result->matches_table, vHM, i, segLen, j, s2Len, bias);
                arr_store_si256(result->similar_table, vHS, i, segLen, j, s2Len, bias);
                arr_store_si256(result->length_table, vHL, i, segLen, j, s2Len, bias);
                arr_store_si256(result->score_table, vH, i, segLen, j, s2Len, bias);
#endif
                vMaxH = _mm256_max_epi16(vH, vMaxH);
                vH = _mm256_subs_epi16(vH, vGapO);
                vF = _mm256_subs_epi16(vF, vGapE);
                if (! _mm256_movemask_epi8(_mm256_cmpgt_epi16(vF, vH))) goto end;
                /*vF = _mm256_max_epi16(vF, vH);*/
                vFM = vHM;
                vFS = vHS;
                vFL = vHL;
                vHp = _mm256_load_si256(pvHLoad + i);
            }
        }
end:
        {
        }

        {
            __m256i vCompare = _mm256_cmpgt_epi16(vMaxH, vMaxHUnit);
            if (_mm256_movemask_epi8(vCompare)) {
                score = _mm256_hmax_epi16_rpl(vMaxH);
                vMaxHUnit = _mm256_set1_epi16(score);
                end_ref = j;
                (void)memcpy(pvHMax, pvHStore, sizeof(__m256i)*segLen);
                (void)memcpy(pvHMMax, pvHMStore, sizeof(__m256i)*segLen);
                (void)memcpy(pvHSMax, pvHSStore, sizeof(__m256i)*segLen);
                (void)memcpy(pvHLMax, pvHLStore, sizeof(__m256i)*segLen);
            }
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = _mm256_load_si256(pvHStore + offset);
            vHM = _mm256_load_si256(pvHMStore + offset);
            vHS = _mm256_load_si256(pvHSStore + offset);
            vHL = _mm256_load_si256(pvHLStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm256_slli_si256_rpl(vH, 2);
                vHM = _mm256_slli_si256_rpl(vHM, 2);
                vHS = _mm256_slli_si256_rpl(vHS, 2);
                vHL = _mm256_slli_si256_rpl(vHL, 2);
            }
            result->score_row[j] = (int16_t) _mm256_extract_epi16_rpl (vH, 15) - bias;
            result->matches_row[j] = (int16_t) _mm256_extract_epi16_rpl (vHM, 15) - bias;
            result->similar_row[j] = (int16_t) _mm256_extract_epi16_rpl (vHS, 15) - bias;
            result->length_row[j] = (int16_t) _mm256_extract_epi16_rpl (vHL, 15) - bias;
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
        __m256i vH = _mm256_load_si256(pvHStore+i);
        __m256i vHM = _mm256_load_si256(pvHMStore+i);
        __m256i vHS = _mm256_load_si256(pvHSStore+i);
        __m256i vHL = _mm256_load_si256(pvHLStore+i);
        arr_store_col(result->score_col, vH, i, segLen, bias);
        arr_store_col(result->matches_col, vHM, i, segLen, bias);
        arr_store_col(result->similar_col, vHS, i, segLen, bias);
        arr_store_col(result->length_col, vHL, i, segLen, bias);
    }
#endif

    if (score == INT16_MAX
            || _mm256_movemask_epi8(_mm256_cmpeq_epi16(vSaturationCheckMax,vPosLimit))) {
        result->saturated = 1;
        score = INT16_MAX;
        matches = INT16_MIN;
        similar = INT16_MIN;
        length = INT16_MIN;
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
    parasail_result_t *result = parasail_sw_striped_profile_avx2_256_16(profile, s2, s2Len, open, gap);
    if (!result->saturated) {
#if 0
        int s1Len_new = 0;
        int s2Len_new = 0;
        char *s1_new = NULL;
        char *s2_new = NULL;
        parasail_profile_t *profile_new = NULL;
        parasail_result_t *result_new = NULL;
        int s1_begin = 0;
        int s2_begin = 0;
        int s1Len_final = 0;
        int s2Len_final = 0;
        parasail_profile_t *profile_final = NULL;
        parasail_result_t *result_final = NULL;

        /* using the end loc and the non-stats version of the function,
         * reverse the inputs and find the beg loc */
        s1Len_new = result->end_query+1;
        s2Len_new = result->end_ref+1;
        s1_new = parasail_reverse(s1, s1Len_new);
        s2_new = parasail_reverse(s2, s2Len_new);
        profile_new = parasail_profile_create_avx_256_16(
                s1_new, s1Len_new, matrix);
        profile_new->stop = result->score;
        result_new = parasail_sw_striped_profile_avx2_256_16(
                profile_new, s2_new, s2Len_new, open, gap);

        /* using both the beg and end loc, call the original stats func */
        s1_begin = s1Len_new - result_new->end_query - 1;
        s2_begin = s2Len_new - result_new->end_ref - 1;
        s1Len_final = s1Len_new - s1_begin;
        s2Len_final = s2Len_new - s2_begin;
        assert(s1_begin >= 0);
        assert(s2_begin >= 0);
        assert(s1Len_new > s1_begin);
        assert(s2Len_new > s2_begin);
        profile_final = parasail_profile_create_stats_avx_256_16(
                &s1[s1_begin], s1Len_final, matrix);
        result_final = PNAME(
                profile_final, &s2[s2_begin], s2Len_final, open, gap);

        /* clean up all the temporary profiles, sequences, and results */
        free(s1_new);
        free(s2_new);
        parasail_profile_free(profile_new);
        parasail_profile_free(profile_final);
        parasail_result_free(result);
        parasail_result_free(result_new);

        /* correct the end locations before returning */
        result_final->end_query = s1Len_new-1;
        result_final->end_ref = s2Len_new-1;
        return result_final;
#else
        int s1Len_new = 0;
        int s2Len_new = 0;
        parasail_profile_t *profile_final = NULL;
        parasail_result_t *result_final = NULL;

        /* using the end loc, call the original stats function */
        s1Len_new = result->end_query+1;
        s2Len_new = result->end_ref+1;
        profile_final = parasail_profile_create_stats_avx_256_16(
                s1, s1Len_new, matrix);
        result_final = PNAME(
                profile_final, s2, s2Len_new, open, gap);

        /* clean up all the temporary profiles, sequences, and results */
        parasail_profile_free(profile_final);
        parasail_result_free(result);

        /* correct the end locations before returning */
        result_final->end_query = s1Len_new-1;
        result_final->end_ref = s2Len_new-1;
        return result_final;
#endif
    }
    else {
        return result;
    }
}
#endif
#endif
#endif


