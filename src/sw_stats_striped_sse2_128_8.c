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
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_stats_rowcol_striped_sse2_128_8
#else
#define FNAME parasail_sw_stats_striped_sse2_128_8
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m128i* const restrict vProfile  = parasail_memalign___m128i(16, n * segLen);
    __m128i* const restrict vProfileM = parasail_memalign___m128i(16, n * segLen);
    __m128i* const restrict vProfileS = parasail_memalign___m128i(16, n * segLen);
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
    __m128i vGapO = _mm_set1_epi8(open);
    __m128i vGapE = _mm_set1_epi8(gap);
    __m128i vOne = _mm_set1_epi8(1);
    int8_t bias = INT8_MIN;
    int8_t score = bias;
    int8_t matches = bias;
    int8_t similar = bias;
    int8_t length = bias;
    __m128i vBias = _mm_set1_epi8(bias);
    __m128i vMaxH = vBias;
    __m128i vMaxHM = vBias;
    __m128i vMaxHS = vBias;
    __m128i vMaxHL = vBias;
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

    /* Generate query profile.
     * Rearrange query sequence & calculate the weight of match/mismatch.
     * Don't alias. */
    {
        int32_t index = 0;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                __m128i_8_t p;
                __m128i_8_t m;
                __m128i_8_t s;
                j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    p.v[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                    m.v[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                    s.v[segNum] = p.v[segNum] > 0;
                    j += segLen;
                }
                _mm_store_si128(&vProfile[index], p.m);
                _mm_store_si128(&vProfileM[index], m.m);
                _mm_store_si128(&vProfileS[index], s.m);
                ++index;
            }
        }
    }

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
        vH = _mm_insert_epi8_rpl(vH, bias, 0);
        vHM = _mm_insert_epi8_rpl(vHM, bias, 0);
        vHS = _mm_insert_epi8_rpl(vHS, bias, 0);
        vHL = _mm_insert_epi8_rpl(vHL, bias, 0);

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
            /* update max vector seen so far */
            {
                __m128i cond_max = _mm_cmpgt_epi8(vH, vMaxH);
                vMaxH = _mm_blendv_epi8_rpl(vMaxH, vH,  cond_max);
                vMaxHM = _mm_blendv_epi8_rpl(vMaxHM, vHM, cond_max);
                vMaxHS = _mm_blendv_epi8_rpl(vMaxHS, vHS, cond_max);
                vMaxHL = _mm_blendv_epi8_rpl(vMaxHL, vHL, cond_max);
            }

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
            vHp = _mm_insert_epi8_rpl(vHp, bias, 0);
            vF  = _mm_insert_epi8_rpl(vF,  bias, 0);
            vFM = _mm_insert_epi8_rpl(vFM, bias, 0);
            vFS = _mm_insert_epi8_rpl(vFS, bias, 0);
            vFL = _mm_insert_epi8_rpl(vFL, bias, 0);
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

    /* max in vec */
    for (j=0; j<segWidth; ++j) {
        int8_t value = (int8_t) _mm_extract_epi8_rpl(vMaxH, 15);
        if (value > score) {
            score = value;
            matches = (int8_t)_mm_extract_epi8_rpl(vMaxHM, 15);
            similar = (int8_t)_mm_extract_epi8_rpl(vMaxHS, 15);
            length = (int8_t)_mm_extract_epi8_rpl(vMaxHL, 15);
        }
        vMaxH = _mm_slli_si128(vMaxH, 1);
        vMaxHM = _mm_slli_si128(vMaxHM, 1);
        vMaxHS = _mm_slli_si128(vMaxHS, 1);
        vMaxHL = _mm_slli_si128(vMaxHL, 1);
    }

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
    parasail_free(vProfileS);
    parasail_free(vProfileM);
    parasail_free(vProfile);

    return result;
}


