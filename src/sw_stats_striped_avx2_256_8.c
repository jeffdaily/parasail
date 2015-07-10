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

#define _mm256_cmplt_epi8_rpl(a,b) _mm256_cmpgt_epi8(b,a)

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)


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
    array[( 0*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  0) - bias;
    array[( 1*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  1) - bias;
    array[( 2*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  2) - bias;
    array[( 3*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  3) - bias;
    array[( 4*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  4) - bias;
    array[( 5*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  5) - bias;
    array[( 6*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  6) - bias;
    array[( 7*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  7) - bias;
    array[( 8*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  8) - bias;
    array[( 9*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  9) - bias;
    array[(10*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 10) - bias;
    array[(11*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 11) - bias;
    array[(12*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 12) - bias;
    array[(13*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 13) - bias;
    array[(14*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 14) - bias;
    array[(15*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 15) - bias;
    array[(16*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 16) - bias;
    array[(17*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 17) - bias;
    array[(18*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 18) - bias;
    array[(19*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 19) - bias;
    array[(20*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 20) - bias;
    array[(21*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 21) - bias;
    array[(22*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 22) - bias;
    array[(23*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 23) - bias;
    array[(24*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 24) - bias;
    array[(25*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 25) - bias;
    array[(26*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 26) - bias;
    array[(27*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 27) - bias;
    array[(28*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 28) - bias;
    array[(29*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 29) - bias;
    array[(30*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 30) - bias;
    array[(31*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 31) - bias;
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
    col[ 0*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  0) - bias;
    col[ 1*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  1) - bias;
    col[ 2*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  2) - bias;
    col[ 3*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  3) - bias;
    col[ 4*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  4) - bias;
    col[ 5*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  5) - bias;
    col[ 6*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  6) - bias;
    col[ 7*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  7) - bias;
    col[ 8*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  8) - bias;
    col[ 9*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  9) - bias;
    col[10*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 10) - bias;
    col[11*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 11) - bias;
    col[12*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 12) - bias;
    col[13*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 13) - bias;
    col[14*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 14) - bias;
    col[15*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 15) - bias;
    col[16*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 16) - bias;
    col[17*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 17) - bias;
    col[18*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 18) - bias;
    col[19*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 19) - bias;
    col[20*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 20) - bias;
    col[21*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 21) - bias;
    col[22*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 22) - bias;
    col[23*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 23) - bias;
    col[24*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 24) - bias;
    col[25*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 25) - bias;
    col[26*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 26) - bias;
    col[27*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 27) - bias;
    col[28*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 28) - bias;
    col[29*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 29) - bias;
    col[30*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 30) - bias;
    col[31*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 31) - bias;
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_stats_table_striped_avx2_256_8
#define PNAME parasail_sw_stats_table_striped_profile_avx2_256_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_stats_rowcol_striped_avx2_256_8
#define PNAME parasail_sw_stats_rowcol_striped_profile_avx2_256_8
#else
#define FNAME parasail_sw_stats_striped_avx2_256_8
#define PNAME parasail_sw_stats_striped_profile_avx2_256_8
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
    int32_t segNum = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const int32_t segWidth = 32; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m256i* const restrict vProfile  = (__m256i*)profile->profile8.score;
    __m256i* const restrict vProfileM = (__m256i*)profile->profile8.similar;
    __m256i* const restrict vProfileS = (__m256i*)profile->profile8.matches;
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
    __m256i vGapO = _mm256_set1_epi8(open);
    __m256i vGapE = _mm256_set1_epi8(gap);
    __m256i vOne = _mm256_set1_epi8(1);
    int8_t bias = INT8_MIN;
    int8_t score = bias;
    int8_t matches = bias;
    int8_t similar = bias;
    int8_t length = bias;
    __m256i vBias = _mm256_set1_epi8(bias);
    __m256i vMaxH = vBias;
    __m256i vMaxHM = vBias;
    __m256i vMaxHS = vBias;
    __m256i vMaxHL = vBias;
    __m256i vSaturationCheckMax = vBias;
    __m256i vPosLimit = _mm256_set1_epi8(INT8_MAX);
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

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m256i_8_t h;
            __m256i_8_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = bias;
                e.v[segNum] = bias;
            }
            _mm256_store_si256(&pvHStore[index], h.m);
            _mm256_store_si256(&pvEStore[index], e.m);
            ++index;
        }
    }

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
        vH = _mm256_slli_si256_rpl(pvHStore[segLen - 1], 1);
        vHM = _mm256_slli_si256_rpl(pvHMStore[segLen - 1], 1);
        vHS = _mm256_slli_si256_rpl(pvHSStore[segLen - 1], 1);
        vHL = _mm256_slli_si256_rpl(pvHLStore[segLen - 1], 1);
        vH = _mm256_insert_epi8_rpl(vH, bias, 0);
        vHM = _mm256_insert_epi8_rpl(vHM, bias, 0);
        vHS = _mm256_insert_epi8_rpl(vHS, bias, 0);
        vHL = _mm256_insert_epi8_rpl(vHL, bias, 0);

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

            vH = _mm256_adds_epi8(vH, _mm256_load_si256(vP + i));
            vE = _mm256_load_si256(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            case1not = _mm256_or_si256(
                    _mm256_cmplt_epi8_rpl(vH,vF),_mm256_cmplt_epi8_rpl(vH,vE));
            case2not = _mm256_cmplt_epi8_rpl(vF,vE);
            case2 = _mm256_andnot_si256(case2not,case1not);
            case3 = _mm256_and_si256(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm256_max_epi8(vH, vE);
            vH = _mm256_max_epi8(vH, vF);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);
            cond_zero = _mm256_cmpeq_epi8(vH, vBias);

            /* calculate vM */
            vEM = _mm256_load_si256(pvEM + i);
            vHM = _mm256_blendv_epi8(
                    _mm256_adds_epi8(vHM, _mm256_load_si256(vPM + i)),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, vFM),
                        _mm256_and_si256(case3, vEM)),
                    case1not);
            vHM = _mm256_blendv_epi8(vHM, vBias, cond_zero);
            _mm256_store_si256(pvHMStore + i, vHM);

            /* calculate vS */
            vES = _mm256_load_si256(pvES + i);
            vHS = _mm256_blendv_epi8(
                    _mm256_adds_epi8(vHS, _mm256_load_si256(vPS + i)),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, vFS),
                        _mm256_and_si256(case3, vES)),
                    case1not);
            vHS = _mm256_blendv_epi8(vHS, vBias, cond_zero);
            _mm256_store_si256(pvHSStore + i, vHS);

            /* calculate vL */
            vEL = _mm256_load_si256(pvEL + i);
            vHL = _mm256_blendv_epi8(
                    _mm256_adds_epi8(vHL, vOne),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, _mm256_adds_epi8(vFL, vOne)),
                        _mm256_and_si256(case3, _mm256_adds_epi8(vEL, vOne))),
                    case1not);
            vHL = _mm256_blendv_epi8(vHL, vBias, cond_zero);
            _mm256_store_si256(pvHLStore + i, vHL);
            vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vHM);
            vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vHS);
            vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->matches_table, vHM, i, segLen, j, s2Len, bias);
            arr_store_si256(result->similar_table, vHS, i, segLen, j, s2Len, bias);
            arr_store_si256(result->length_table, vHL, i, segLen, j, s2Len, bias);
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len, bias);
#endif
            /* update max vector seen so far */
            {
                __m256i cond_max = _mm256_cmpgt_epi8(vH, vMaxH);
                vMaxH = _mm256_blendv_epi8(vMaxH, vH,  cond_max);
                vMaxHM = _mm256_blendv_epi8(vMaxHM, vHM, cond_max);
                vMaxHS = _mm256_blendv_epi8(vMaxHS, vHS, cond_max);
                vMaxHL = _mm256_blendv_epi8(vMaxHL, vHL, cond_max);
            }

            /* Update vE value. */
            vH = _mm256_subs_epi8(vH, vGapO);
            vE = _mm256_subs_epi8(vE, vGapE);
            vE = _mm256_max_epi8(vE, vH);
            _mm256_store_si256(pvEStore + i, vE);
            _mm256_store_si256(pvEM + i, vHM);
            _mm256_store_si256(pvES + i, vHS);
            _mm256_store_si256(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm256_subs_epi8(vF, vGapE);
            vF = _mm256_max_epi8(vF, vH);
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
            __m256i vHp = _mm256_slli_si256_rpl(pvHLoad[segLen - 1], 1);
            vF = _mm256_slli_si256_rpl(vF, 1);
            vFM = _mm256_slli_si256_rpl(vFM, 1);
            vFS = _mm256_slli_si256_rpl(vFS, 1);
            vFL = _mm256_slli_si256_rpl(vFL, 1);
            vHp = _mm256_insert_epi8_rpl(vHp, bias, 0);
            vF  = _mm256_insert_epi8_rpl(vF,  bias, 0);
            vFM = _mm256_insert_epi8_rpl(vFM, bias, 0);
            vFS = _mm256_insert_epi8_rpl(vFS, bias, 0);
            vFL = _mm256_insert_epi8_rpl(vFL, bias, 0);
            for (i=0; i<segLen; ++i) {
                __m256i case1not;
                __m256i case2not;
                __m256i case2;
                __m256i cond_zero;

                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm256_adds_epi8(vHp, _mm256_load_si256(vP + i));
                vE = _mm256_load_si256(pvELoad + i);
                case1not = _mm256_or_si256(
                        _mm256_cmplt_epi8_rpl(vHp,vF),_mm256_cmplt_epi8_rpl(vHp,vE));
                case2not = _mm256_cmplt_epi8_rpl(vF,vE);
                case2 = _mm256_andnot_si256(case2not,case1not);

                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi8(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
                cond_zero = _mm256_cmpeq_epi8(vH, vBias);

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
                vHL = _mm256_blendv_epi8(vHL, _mm256_adds_epi8(vFL,vOne), case2);
                vHL = _mm256_blendv_epi8(vHL, vBias, cond_zero);
                _mm256_store_si256(pvHLStore + i, vHL);
                _mm256_store_si256(pvEL + i, vHL);

                vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vHM);
                vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vHS);
                vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
                arr_store_si256(result->matches_table, vHM, i, segLen, j, s2Len, bias);
                arr_store_si256(result->similar_table, vHS, i, segLen, j, s2Len, bias);
                arr_store_si256(result->length_table, vHL, i, segLen, j, s2Len, bias);
                arr_store_si256(result->score_table, vH, i, segLen, j, s2Len, bias);
#endif
                vH = _mm256_subs_epi8(vH, vGapO);
                vF = _mm256_subs_epi8(vF, vGapE);
                if (! _mm256_movemask_epi8(_mm256_cmpgt_epi8(vF, vH))) goto end;
                /*vF = _mm256_max_epi8(vF, vH);*/
                vFM = vHM;
                vFS = vHS;
                vFL = vHL;
                vHp = _mm256_load_si256(pvHLoad + i);
            }
        }
end:
        {
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = _mm256_load_si256(pvHStore + offset);
            vHM = _mm256_load_si256(pvHMStore + offset);
            vHS = _mm256_load_si256(pvHSStore + offset);
            vHL = _mm256_load_si256(pvHLStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm256_slli_si256_rpl(vH, 1);
                vHM = _mm256_slli_si256_rpl(vHM, 1);
                vHS = _mm256_slli_si256_rpl(vHS, 1);
                vHL = _mm256_slli_si256_rpl(vHL, 1);
            }
            result->score_row[j] = (int8_t) _mm256_extract_epi8_rpl (vH, 31) - bias;
            result->matches_row[j] = (int8_t) _mm256_extract_epi8_rpl (vHM, 31) - bias;
            result->similar_row[j] = (int8_t) _mm256_extract_epi8_rpl (vHS, 31) - bias;
            result->length_row[j] = (int8_t) _mm256_extract_epi8_rpl (vHL, 31) - bias;
        }
#endif
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

    /* max in vec */
    for (j=0; j<segWidth; ++j) {
        int8_t value = (int8_t) _mm256_extract_epi8_rpl(vMaxH, 31);
        if (value > score) {
            score = value;
            matches = (int8_t)_mm256_extract_epi8_rpl(vMaxHM, 31);
            similar = (int8_t)_mm256_extract_epi8_rpl(vMaxHS, 31);
            length = (int8_t)_mm256_extract_epi8_rpl(vMaxHL, 31);
        }
        vMaxH = _mm256_slli_si256_rpl(vMaxH, 1);
        vMaxHM = _mm256_slli_si256_rpl(vMaxHM, 1);
        vMaxHS = _mm256_slli_si256_rpl(vMaxHS, 1);
        vMaxHL = _mm256_slli_si256_rpl(vMaxHL, 1);
    }

    if (score == INT8_MAX
            || _mm256_movemask_epi8(_mm256_cmpeq_epi8(vSaturationCheckMax,vPosLimit))) {
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


