/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include <immintrin.h>

#include "parasail.h"
#include "parasail_internal.h"
#include "parasail_internal_avx.h"
#include "blosum/blosum_map.h"

#define NEG_INF_8 (INT8_MIN)

/* avx2 does not have _mm256_cmplt_epi8, emulate it */
static inline __m256i _mm256_cmplt_epi8(__m256i a, __m256i b) {
    return _mm256_cmpgt_epi8(b,a);
}

#if HAVE_AVX2_MM256_INSERT_EPI8
#else
static inline __m256i _mm256_insert_epi8(__m256i a, int i, int imm) {
    __m256i_8_t tmp;
    tmp.m = a;
    tmp.v[imm] = i;
    return tmp.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI8
#else
static inline int8_t _mm256_extract_epi8(__m256i a, int imm) {
    __m256i_8_t tmp;
    tmp.m = a;
    return tmp.v[imm];
}
#endif

/* avx2 _mm256_slli_si256 does not shift across 128-bit lanes, emulate it */
static inline __m256i shift(__m256i a) {
    return _mm256_alignr_epi8(a,
            _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)),
            15);
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
    array[( 0*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  0) - bias;
    array[( 1*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  1) - bias;
    array[( 2*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  2) - bias;
    array[( 3*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  3) - bias;
    array[( 4*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  4) - bias;
    array[( 5*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  5) - bias;
    array[( 6*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  6) - bias;
    array[( 7*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  7) - bias;
    array[( 8*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  8) - bias;
    array[( 9*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  9) - bias;
    array[(10*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 10) - bias;
    array[(11*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 11) - bias;
    array[(12*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 12) - bias;
    array[(13*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 13) - bias;
    array[(14*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 14) - bias;
    array[(15*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 15) - bias;
    array[(16*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 16) - bias;
    array[(17*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 17) - bias;
    array[(18*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 18) - bias;
    array[(19*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 19) - bias;
    array[(20*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 20) - bias;
    array[(21*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 21) - bias;
    array[(22*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 22) - bias;
    array[(23*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 23) - bias;
    array[(24*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 24) - bias;
    array[(25*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 25) - bias;
    array[(26*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 26) - bias;
    array[(27*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 27) - bias;
    array[(28*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 28) - bias;
    array[(29*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 29) - bias;
    array[(30*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 30) - bias;
    array[(31*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 31) - bias;
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME sw_stats_table_striped_avx2_256_8
#else
#define FNAME sw_stats_striped_avx2_256_8
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = 24; /* number of amino acids in table */
    const int32_t segWidth = 32; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m256i* const restrict vProfile  = parasail_memalign___m256i(32, n * segLen);
    __m256i* const restrict vProfileM = parasail_memalign___m256i(32, n * segLen);
    __m256i* const restrict vProfileS = parasail_memalign___m256i(32, n * segLen);
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
    __m256i* restrict pvEM            = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvES            = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvEL            = parasail_memalign___m256i(32, segLen);
    int score = NEG_INF_8;
    int matches = NEG_INF_8;
    int similar = NEG_INF_8;
    int length = NEG_INF_8;
    __m256i vGapO = _mm256_set1_epi8(open);
    __m256i vGapE = _mm256_set1_epi8(gap);
    __m256i vNegInf = _mm256_set1_epi8(NEG_INF_8);
    __m256i vOne = _mm256_set1_epi8(1);
    int8_t bias = INT8_MIN;
    /* Trace the highest score of the whole SW matrix. */
    __m256i vMaxH = vNegInf;
    __m256i vMaxM = vNegInf;
    __m256i vMaxS = vNegInf;
    __m256i vMaxL = vNegInf;
    __m256i vSaturationCheck = _mm256_setzero_si256();
    __m256i vPosLimit = _mm256_set1_epi8(INT8_MAX);
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif

    parasail_memset___m256i(pvHMStore, vNegInf, segLen);
    parasail_memset___m256i(pvHSStore, vNegInf, segLen);
    parasail_memset___m256i(pvHLStore, vNegInf, segLen);

    /* Generate query profile.
     * Rearrange query sequence & calculate the weight of match/mismatch.
     * Don't alias. */
    {
        int32_t index = 0;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                __m256i_8_t p;
                __m256i_8_t m;
                __m256i_8_t s;
                j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    p.v[segNum] = j >= s1Len ? 0 : matrix[k][MAP_BLOSUM_[(unsigned char)s1[j]]];
                    m.v[segNum] = j >= s1Len ? 0 : (k == MAP_BLOSUM_[(unsigned char)s1[j]]);
                    s.v[segNum] = p.v[segNum] > 0;
                    j += segLen;
                }
                _mm256_store_si256(&vProfile[index], p.m);
                _mm256_store_si256(&vProfileM[index], m.m);
                _mm256_store_si256(&vProfileS[index], s.m);
                ++index;
            }
        }
    }

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m256i_8_t h;
            __m256i_8_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = bias;
                e.v[segNum] = NEG_INF_8;
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
        vF  = vNegInf;
        vFM = vNegInf;
        vFS = vNegInf;
        vFL = vNegInf;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH  = shift(pvHStore[segLen - 1]);
        vHM = shift(pvHMStore[segLen - 1]);
        vHS = shift(pvHSStore[segLen - 1]);
        vHL = shift(pvHLStore[segLen - 1]);
        vH  = _mm256_insert_epi8(vH,  bias, 0);
        vHM = _mm256_insert_epi8(vHM, bias, 0);
        vHS = _mm256_insert_epi8(vHS, bias, 0);
        vHL = _mm256_insert_epi8(vHL, bias, 0);

        /* Correct part of the vProfile */
        vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;
        vPM = vProfileM + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;
        vPS = vProfileS + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

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
                    _mm256_cmplt_epi8(vH,vF),_mm256_cmplt_epi8(vH,vE));
            case2not = _mm256_cmplt_epi8(vF,vE);
            case2 = _mm256_andnot_si256(case2not,case1not);
            case3 = _mm256_and_si256(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm256_max_epi8(vH, vE);
            vH = _mm256_max_epi8(vH, vF);
            /*vH = _mm256_max_epi8(vH, vZero);*/
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);
            cond_zero = _mm256_cmpeq_epi8(vH, vNegInf);

            /* calculate vM */
            vEM = _mm256_load_si256(pvEM + i);
            vHM = _mm256_andnot_si256(case1not,
                    _mm256_adds_epi8(vHM, _mm256_load_si256(vPM + i)));
            vHM = _mm256_or_si256(vHM, _mm256_and_si256(case2, vFM));
            vHM = _mm256_or_si256(vHM, _mm256_and_si256(case3, vEM));
            vHM = _mm256_andnot_si256(cond_zero, vHM);
            vHM = _mm256_or_si256(vHM, _mm256_and_si256(cond_zero, vNegInf));
            _mm256_store_si256(pvHMStore + i, vHM);

            /* calculate vS */
            vES = _mm256_load_si256(pvES + i);
            vHS = _mm256_andnot_si256(case1not,
                    _mm256_adds_epi8(vHS, _mm256_load_si256(vPS + i)));
            vHS = _mm256_or_si256(vHS, _mm256_and_si256(case2, vFS));
            vHS = _mm256_or_si256(vHS, _mm256_and_si256(case3, vES));
            vHS = _mm256_andnot_si256(cond_zero, vHS);
            vHS = _mm256_or_si256(vHS, _mm256_and_si256(cond_zero, vNegInf));
            _mm256_store_si256(pvHSStore + i, vHS);

            /* calculate vL */
            vEL = _mm256_load_si256(pvEL + i);
            vHL = _mm256_andnot_si256(case1not, _mm256_adds_epi8(vHL, vOne));
            vHL = _mm256_or_si256(vHL, _mm256_and_si256(case2,
                        _mm256_adds_epi8(vFL, vOne)));
            vHL = _mm256_or_si256(vHL, _mm256_and_si256(case3,
                        _mm256_adds_epi8(vEL, vOne)));
            vHL = _mm256_andnot_si256(cond_zero, vHL);
            vHL = _mm256_or_si256(vHL, _mm256_and_si256(cond_zero, vNegInf));
            _mm256_store_si256(pvHLStore + i, vHL);

            /* check for saturation */
            {
                vSaturationCheck = _mm256_or_si256(vSaturationCheck,
                        _mm256_or_si256(
                            _mm256_cmpeq_epi8(vH, vPosLimit),
                            _mm256_or_si256(
                                _mm256_cmpeq_epi8(vHM, vPosLimit),
                                _mm256_or_si256(
                                    _mm256_cmpeq_epi8(vHS, vPosLimit),
                                    _mm256_cmpeq_epi8(vHL, vPosLimit)))));
            }
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
                vMaxM = _mm256_blendv_epi8(vMaxM, vHM, cond_max);
                vMaxS = _mm256_blendv_epi8(vMaxS, vHS, cond_max);
                vMaxL = _mm256_blendv_epi8(vMaxL, vHL, cond_max);
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
            __m256i vHp = shift(pvHLoad[segLen - 1]);
            vHp = _mm256_insert_epi8(vHp, bias, 0);
            vF  = shift(vF);
            vFM = shift(vFM);
            vFS = shift(vFS);
            vFL = shift(vFL);
            vF  = _mm256_insert_epi8(vF,  bias, 0);
            vFM = _mm256_insert_epi8(vFM, bias, 0);
            vFS = _mm256_insert_epi8(vFS, bias, 0);
            vFL = _mm256_insert_epi8(vFL, bias, 0);
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
                        _mm256_cmplt_epi8(vHp,vF),_mm256_cmplt_epi8(vHp,vE));
                case2not = _mm256_cmplt_epi8(vF,vE);
                case2 = _mm256_andnot_si256(case2not,case1not);

                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi8(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
                cond_zero = _mm256_cmpeq_epi8(vH, vNegInf);

                vHM = _mm256_load_si256(pvHMStore + i);
                vHM = _mm256_andnot_si256(case2, vHM);
                vHM = _mm256_or_si256(vHM, _mm256_and_si256(case2, vFM));
                vHM = _mm256_andnot_si256(cond_zero, vHM);
                vHM = _mm256_or_si256(vHM, _mm256_and_si256(cond_zero, vNegInf));
                _mm256_store_si256(pvHMStore + i, vHM);
                _mm256_store_si256(pvEM + i, vHM);

                vHS = _mm256_load_si256(pvHSStore + i);
                vHS = _mm256_andnot_si256(case2, vHS);
                vHS = _mm256_or_si256(vHS, _mm256_and_si256(case2, vFS));
                vHS = _mm256_andnot_si256(cond_zero, vHS);
                vHS = _mm256_or_si256(vHS, _mm256_and_si256(cond_zero, vNegInf));
                _mm256_store_si256(pvHSStore + i, vHS);
                _mm256_store_si256(pvES + i, vHS);

                vHL = _mm256_load_si256(pvHLStore + i);
                vHL = _mm256_andnot_si256(case2, vHL);
                vHL = _mm256_or_si256(vHL, _mm256_and_si256(case2,
                            _mm256_adds_epi8(vFL,vOne)));
                vHL = _mm256_andnot_si256(cond_zero, vHL);
                vHL = _mm256_or_si256(vHL, _mm256_and_si256(cond_zero, vNegInf));
                _mm256_store_si256(pvHLStore + i, vHL);
                _mm256_store_si256(pvEL + i, vHL);

                /* check for saturation */
                {
                    vSaturationCheck = _mm256_or_si256(vSaturationCheck,
                            _mm256_or_si256(
                                _mm256_cmpeq_epi8(vH, vPosLimit),
                                _mm256_or_si256(
                                    _mm256_cmpeq_epi8(vHM, vPosLimit),
                                    _mm256_or_si256(
                                        _mm256_cmpeq_epi8(vHS, vPosLimit),
                                        _mm256_cmpeq_epi8(vHL, vPosLimit)))));
                }
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
    }

    /* max in vec */
    for (j=0; j<segWidth; ++j) {
        int value = (int8_t) _mm256_extract_epi8(vMaxH, 31) - (int)bias;
        if (value > score) {
            score = value;
            matches = (int8_t) _mm256_extract_epi8(vMaxM, 31) - (int)bias;
            similar = (int8_t) _mm256_extract_epi8(vMaxS, 31) - (int)bias;
            length = (int8_t) _mm256_extract_epi8(vMaxL, 31) - (int)bias;
        }
        vMaxH = shift(vMaxH);
        vMaxM = shift(vMaxM);
        vMaxS = shift(vMaxS);
        vMaxL = shift(vMaxL);
    }

    if (_mm256_movemask_epi8(vSaturationCheck)) {
        result->saturated = 1;
        score = INT8_MAX;
        matches = 0;
        similar = 0;
        length = 0;
    }

    result->score = score;
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
    parasail_free(vProfileS);
    parasail_free(vProfileM);
    parasail_free(vProfile);

    return result;
}

