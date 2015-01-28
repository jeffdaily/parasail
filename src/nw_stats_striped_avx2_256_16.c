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

#define NEG_INF_16 (INT16_MIN/(int16_t)(2))

/* avx2 does not have _mm256_cmplt_epi16, emulate it */
static inline __m256i _mm256_cmplt_epi16(__m256i a, __m256i b) {
    return _mm256_cmpgt_epi16(b,a);
}

#if HAVE_AVX2_MM256_INSERT_EPI16
#else
static inline __m256i _mm256_insert_epi16(__m256i a, int16_t b, int imm) {
    __m256i_16_t tmp;
    tmp.m = a;
    tmp.v[imm] = b;
    return tmp.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI16
#else
static inline int16_t _mm256_extract_epi16(__m256i a, int imm) {
    __m256i_16_t tmp;
    tmp.m = a;
    return tmp.v[imm];
}
#endif

/* avx2 _mm256_slli_si256 does not shift across 128-bit lanes, emulate it */
static inline __m256i shift(__m256i a) {
    return _mm256_alignr_epi8(a,
            _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)),
            14);
}

#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[( 0*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  0);
    array[( 1*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  1);
    array[( 2*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  2);
    array[( 3*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  3);
    array[( 4*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  4);
    array[( 5*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  5);
    array[( 6*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  6);
    array[( 7*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  7);
    array[( 8*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  8);
    array[( 9*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH,  9);
    array[(10*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 10);
    array[(11*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 11);
    array[(12*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 12);
    array[(13*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 13);
    array[(14*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 14);
    array[(15*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16(vH, 15);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME nw_stats_table_striped_avx2_256_16
#else
#define FNAME nw_stats_striped_avx2_256_16
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
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m256i* const restrict vProfile  = parasail_memalign_m256i(32, n * segLen);
    __m256i* const restrict vProfileS = parasail_memalign_m256i(32, n * segLen);
    __m256i* restrict pvHStore        = parasail_memalign_m256i(32, segLen);
    __m256i* restrict pvHLoad         = parasail_memalign_m256i(32, segLen);
    __m256i* restrict pvHMStore       = parasail_memalign_m256i(32, segLen);
    __m256i* restrict pvHMLoad        = parasail_memalign_m256i(32, segLen);
    __m256i* restrict pvHLStore       = parasail_memalign_m256i(32, segLen);
    __m256i* restrict pvHLLoad        = parasail_memalign_m256i(32, segLen);
    __m256i* restrict pvEStore        = parasail_memalign_m256i(32, segLen);
    __m256i* restrict pvELoad         = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvEM      = parasail_memalign_m256i(32, segLen);
    __m256i* const restrict pvEL      = parasail_memalign_m256i(32, segLen);
    int16_t* const restrict boundary  = parasail_memalign_int16_t(32, s2Len+1);
    __m256i vGapO = _mm256_set1_epi16(open);
    __m256i vGapE = _mm256_set1_epi16(gap);
    __m256i vZero = _mm256_setzero_si256();
    __m256i vOne = _mm256_set1_epi16(1);
    int16_t score;
    int16_t matches;
    int16_t length;
    __m256i initialF = _mm256_set_epi16(
            -open-open-15*segLen*gap,
            -open-open-14*segLen*gap,
            -open-open-13*segLen*gap,
            -open-open-12*segLen*gap,
            -open-open-11*segLen*gap,
            -open-open-10*segLen*gap,
            -open-open-9*segLen*gap,
            -open-open-8*segLen*gap,
            -open-open-7*segLen*gap,
            -open-open-6*segLen*gap,
            -open-open-5*segLen*gap,
            -open-open-4*segLen*gap,
            -open-open-3*segLen*gap,
            -open-open-2*segLen*gap,
            -open-open-1*segLen*gap,
            -open-open-0*segLen*gap);
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif

    parasail_memset_m256i(pvHMStore, vZero, segLen);
    parasail_memset_m256i(pvHLStore, vZero, segLen);

    /* Generate query profile.
     * Rearrange query sequence & calculate the weight of match/mismatch.
     * Don't alias. */
    {
        int32_t index = 0;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                __m256i_16_t t;
                __m256i_16_t s;
                j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = j >= s1Len ? 0 : matrix[k][MAP_BLOSUM_[(unsigned char)s1[j]]];
                    s.v[segNum] = (k == MAP_BLOSUM_[(unsigned char)s1[j]]);
                    j += segLen;
                }
                _mm256_store_si256(&vProfile[index], t.m);
                _mm256_store_si256(&vProfileS[index], s.m);
                ++index;
            }
        }
    }

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m256i_16_t h;
            __m256i_16_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = -open-gap*(segNum*segLen+i);
                e.v[segNum] = NEG_INF_16;
            }
            _mm256_store_si256(&pvHStore[index], h.m);
            _mm256_store_si256(&pvEStore[index], e.m);
            ++index;
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            boundary[i] = -open-gap*(i-1);
        }
    }

    initialF = _mm256_add_epi16(initialF, vGapE);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m256i vE;
        __m256i vEM;
        __m256i vEL;
        __m256i vF;
        __m256i vFM;
        __m256i vFL;
        __m256i vH;
        __m256i vHM;
        __m256i vHL;
        const __m256i* vP = NULL;
        const __m256i* vPS = NULL;
        __m256i* pv = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        initialF = _mm256_sub_epi16(initialF, vGapE);
        vF = initialF;
        vFM = vZero;
        vFL = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = shift(pvHStore[segLen - 1]);
        vHM = shift(pvHMStore[segLen - 1]);
        vHL = shift(pvHLStore[segLen - 1]);

        /* insert upper boundary condition */
        vH = _mm256_insert_epi16(vH, boundary[j], 0);

        /* Correct part of the vProfile */
        vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;
        vPS = vProfileS + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        pv = pvHMLoad;
        pvHMLoad = pvHMStore;
        pvHMStore = pv;
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

            vH = _mm256_add_epi16(vH, _mm256_load_si256(vP + i));
            vE = _mm256_load_si256(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            case1not = _mm256_or_si256(
                    _mm256_cmplt_epi16(vH,vF),_mm256_cmplt_epi16(vH,vE));
            case2not = _mm256_cmplt_epi16(vF,vE);
            case2 = _mm256_andnot_si256(case2not,case1not);
            case3 = _mm256_and_si256(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm256_max_epi16(vH, vE);
            vH = _mm256_max_epi16(vH, vF);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);

            /* calculate vM */
            vEM = _mm256_load_si256(pvEM + i);
            vHM = _mm256_andnot_si256(case1not,
                    _mm256_add_epi16(vHM, _mm256_load_si256(vPS + i)));
            vHM = _mm256_or_si256(vHM, _mm256_and_si256(case2, vFM));
            vHM = _mm256_or_si256(vHM, _mm256_and_si256(case3, vEM));
            _mm256_store_si256(pvHMStore + i, vHM);

            /* calculate vL */
            vEL = _mm256_load_si256(pvEL + i);
            vHL = _mm256_andnot_si256(case1not, _mm256_add_epi16(vHL, vOne));
            vHL = _mm256_or_si256(vHL, _mm256_and_si256(case2,
                        _mm256_add_epi16(vFL, vOne)));
            vHL = _mm256_or_si256(vHL, _mm256_and_si256(case3,
                        _mm256_add_epi16(vEL, vOne)));
            _mm256_store_si256(pvHLStore + i, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si256(result->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = _mm256_sub_epi16(vH, vGapO);
            vE = _mm256_sub_epi16(vE, vGapE);
            vE = _mm256_max_epi16(vE, vH);
            _mm256_store_si256(pvEStore + i, vE);
            _mm256_store_si256(pvEM + i, vHM);
            _mm256_store_si256(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm256_sub_epi16(vF, vGapE);
            vF = _mm256_max_epi16(vF, vH);
            vFM = vHM;
            vFL = vHL;

            /* Load the next vH. */
            vH = _mm256_load_si256(pvHLoad + i);
            vHM = _mm256_load_si256(pvHMLoad + i);
            vHL = _mm256_load_si256(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            __m256i vHp = shift(pvHLoad[segLen - 1]);
            vHp = _mm256_insert_epi16(vHp, boundary[j], 0);
            vF = shift(vF);
            vF = _mm256_insert_epi16(vF, boundary[j+1]-open, 0);
            vFM = shift(vFM);
            vFL = shift(vFL);
            for (i=0; i<segLen; ++i) {
                __m256i case1not;
                __m256i case2not;
                __m256i case2;

                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm256_add_epi16(vHp, _mm256_load_si256(vP + i));
                vE = _mm256_load_si256(pvELoad + i);
                case1not = _mm256_or_si256(
                        _mm256_cmplt_epi16(vHp,vF),_mm256_cmplt_epi16(vHp,vE));
                case2not = _mm256_cmplt_epi16(vF,vE);
                case2 = _mm256_andnot_si256(case2not,case1not);

                vHM = _mm256_load_si256(pvHMStore + i);
                vHM = _mm256_andnot_si256(case2, vHM);
                vHM = _mm256_or_si256(vHM, _mm256_and_si256(case2, vFM));
                _mm256_store_si256(pvHMStore + i, vHM);
                _mm256_store_si256(pvEM + i, vHM);

                vHL = _mm256_load_si256(pvHLStore + i);
                vHL = _mm256_andnot_si256(case2, vHL);
                vHL = _mm256_or_si256(vHL, _mm256_and_si256(case2,
                            _mm256_add_epi16(vFL,vOne)));
                _mm256_store_si256(pvHLStore + i, vHL);
                _mm256_store_si256(pvEL + i, vHL);

                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi16(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
                arr_store_si256(result->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si256(result->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm256_sub_epi16(vH, vGapO);
                vF = _mm256_sub_epi16(vF, vGapE);
                if (! _mm256_movemask_epi8(_mm256_cmpgt_epi16(vF, vH))) goto end;
                vF = _mm256_max_epi16(vF, vH);
                vFM = vHM;
                vFL = vHL;
                vHp = _mm256_load_si256(pvHLoad + i);
            }
        }
end:
        {
        }
    }

    /* extract last value from the last column */
    {
        __m256i vH = _mm256_load_si256(pvHStore + offset);
        __m256i vHM = _mm256_load_si256(pvHMStore + offset);
        __m256i vHL = _mm256_load_si256(pvHLStore + offset);
        for (k=0; k<position; ++k) {
            vH = shift (vH);
            vHM = shift (vHM);
            vHL = shift (vHL);
        }
        score = (int16_t) _mm256_extract_epi16 (vH, 15);
        matches = (int16_t) _mm256_extract_epi16 (vHM, 15);
        length = (int16_t) _mm256_extract_epi16 (vHL, 15);
    }

    result->score = score;
    result->matches = matches;
    result->length = length;

    parasail_free(boundary);
    parasail_free(pvEL);
    parasail_free(pvEM);
    parasail_free(pvELoad);
    parasail_free(pvEStore);
    parasail_free(pvHLLoad);
    parasail_free(pvHLStore);
    parasail_free(pvHMLoad);
    parasail_free(pvHMStore);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);
    parasail_free(vProfileS);
    parasail_free(vProfile);

    return result;
}

