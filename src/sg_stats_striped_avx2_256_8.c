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

/* avx2 _mm256_packs_epi16 does not pack across 128-bit lanes, emulate it */
static inline __m256i pack16(__m256i a, __m256i b) {
    return _mm256_permute4x64_epi64(
            _mm256_packs_epi16(a, b),
            _MM_SHUFFLE(3,1,2,0));
}

/* avx2 does not have _mm256_cmplt_epi8, emulate it */
static inline __m256i _mm256_cmplt_epi8(__m256i a, __m256i b) {
    return _mm256_cmpgt_epi8(b,a);
}

/* avx2 does not have _mm256_cmplt_epi16, emulate it */
static inline __m256i _mm256_cmplt_epi16(__m256i a, __m256i b) {
    return _mm256_cmpgt_epi16(b,a);
}

#if HAVE_AVX2_MM256_INSERT_EPI8
#else
static inline __m256i _mm256_insert_epi8(__m256i a, int8_t b, int imm) {
    __m256i_8_t tmp;
    tmp.m = a;
    tmp.v[imm] = b;
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
        int32_t dlen)
{
    array[( 0*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  0);
    array[( 1*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  1);
    array[( 2*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  2);
    array[( 3*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  3);
    array[( 4*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  4);
    array[( 5*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  5);
    array[( 6*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  6);
    array[( 7*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  7);
    array[( 8*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  8);
    array[( 9*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH,  9);
    array[(10*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 10);
    array[(11*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 11);
    array[(12*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 12);
    array[(13*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 13);
    array[(14*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 14);
    array[(15*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 15);
    array[(16*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 16);
    array[(17*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 17);
    array[(18*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 18);
    array[(19*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 19);
    array[(20*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 20);
    array[(21*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 21);
    array[(22*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 22);
    array[(23*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 23);
    array[(24*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 24);
    array[(25*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 25);
    array[(26*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 26);
    array[(27*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 27);
    array[(28*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 28);
    array[(29*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 29);
    array[(30*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 30);
    array[(31*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8(vH, 31);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME sg_stats_table_striped_avx2_256_8
#else
#define FNAME sg_stats_striped_avx2_256_8
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
    __m256i vGapO = _mm256_set1_epi8(open);
    __m256i vGapE = _mm256_set1_epi8(gap);
    __m256i vNegInf = _mm256_set1_epi8(NEG_INF_8);
    __m256i vZero = _mm256_setzero_si256();
    __m256i vOne = _mm256_set1_epi8(1);
    int8_t score = NEG_INF_8;
    int8_t matches = NEG_INF_8;
    int8_t length = NEG_INF_8;
    __m256i vMaxH = vNegInf;
    __m256i vMaxHM = vNegInf;
    __m256i vMaxHL = vNegInf;
    __m256i vSaturationCheck = _mm256_setzero_si256();
    __m256i vNegLimit = _mm256_set1_epi8(INT8_MIN);
    __m256i vPosLimit = _mm256_set1_epi8(INT8_MAX);
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
                __m256i_8_t t;
                __m256i_8_t s;
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
            __m256i_8_t h;
            __m256i_8_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = 0;
                e.v[segNum] = -open;
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
        //vF = initialF;
        vF = vNegInf;
        vFM = vZero;
        vFL = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = shift(pvHStore[segLen - 1]);
        vHM = shift(pvHMStore[segLen - 1]);
        vHL = shift(pvHLStore[segLen - 1]);

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
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);

            /* calculate vM */
            vEM = _mm256_load_si256(pvEM + i);
            vHM = _mm256_andnot_si256(case1not,
                    _mm256_adds_epi8(vHM, _mm256_load_si256(vPS + i)));
            vHM = _mm256_or_si256(vHM, _mm256_and_si256(case2, vFM));
            vHM = _mm256_or_si256(vHM, _mm256_and_si256(case3, vEM));
            _mm256_store_si256(pvHMStore + i, vHM);

            /* calculate vL */
            vEL = _mm256_load_si256(pvEL + i);
            vHL = _mm256_andnot_si256(case1not, _mm256_adds_epi8(vHL, vOne));
            vHL = _mm256_or_si256(vHL, _mm256_and_si256(case2,
                        _mm256_adds_epi8(vFL, vOne)));
            vHL = _mm256_or_si256(vHL, _mm256_and_si256(case3,
                        _mm256_adds_epi8(vEL, vOne)));
            _mm256_store_si256(pvHLStore + i, vHL);

            /* check for saturation */
            {
                vSaturationCheck = _mm256_or_si256(vSaturationCheck,
                        _mm256_or_si256(
                            _mm256_or_si256(
                                _mm256_cmpeq_epi8(vH, vNegLimit),
                                _mm256_cmpeq_epi8(vH, vPosLimit)),
                            _mm256_or_si256(
                                _mm256_cmpeq_epi8(vHM, vPosLimit),
                                _mm256_cmpeq_epi8(vHL, vPosLimit))));
            }
#ifdef PARASAIL_TABLE
            arr_store_si256(result->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si256(result->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = _mm256_subs_epi8(vH, vGapO);
            vE = _mm256_subs_epi8(vE, vGapE);
            vE = _mm256_max_epi8(vE, vH);
            _mm256_store_si256(pvEStore + i, vE);
            _mm256_store_si256(pvEM + i, vHM);
            _mm256_store_si256(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm256_subs_epi8(vF, vGapE);
            vF = _mm256_max_epi8(vF, vH);
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
            vF = shift(vF);
            vF = _mm256_insert_epi8(vF, -open, 0);
            vFM = shift(vFM);
            vFL = shift(vFL);
            for (i=0; i<segLen; ++i) {
                __m256i case1not;
                __m256i case2not;
                __m256i case2;

                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm256_adds_epi8(vHp, _mm256_load_si256(vP + i));
                vE = _mm256_load_si256(pvELoad + i);
                case1not = _mm256_or_si256(
                        _mm256_cmplt_epi8(vHp,vF),_mm256_cmplt_epi8(vHp,vE));
                case2not = _mm256_cmplt_epi8(vF,vE);
                case2 = _mm256_andnot_si256(case2not,case1not);

                vHM = _mm256_load_si256(pvHMStore + i);
                vHM = _mm256_andnot_si256(case2, vHM);
                vHM = _mm256_or_si256(vHM, _mm256_and_si256(case2, vFM));
                _mm256_store_si256(pvHMStore + i, vHM);
                _mm256_store_si256(pvEM + i, vHM);

                vHL = _mm256_load_si256(pvHLStore + i);
                vHL = _mm256_andnot_si256(case2, vHL);
                vHL = _mm256_or_si256(vHL, _mm256_and_si256(case2,
                            _mm256_adds_epi8(vFL,vOne)));
                _mm256_store_si256(pvHLStore + i, vHL);
                _mm256_store_si256(pvEL + i, vHL);

                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi8(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
                /* check for saturation */
                {
                    vSaturationCheck = _mm256_or_si256(vSaturationCheck,
                            _mm256_or_si256(
                                _mm256_or_si256(
                                    _mm256_cmpeq_epi8(vH, vNegLimit),
                                    _mm256_cmpeq_epi8(vH, vPosLimit)),
                                _mm256_or_si256(
                                    _mm256_cmpeq_epi8(vHM, vPosLimit),
                                    _mm256_cmpeq_epi8(vHL, vPosLimit))));
                }
#ifdef PARASAIL_TABLE
                arr_store_si256(result->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si256(result->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm256_subs_epi8(vH, vGapO);
                vF = _mm256_subs_epi8(vF, vGapE);
                if (! _mm256_movemask_epi8(_mm256_cmpgt_epi8(vF, vH))) goto end;
                /*vF = _mm256_max_epi8(vF, vH);*/
                vFM = vHM;
                vFL = vHL;
                vHp = _mm256_load_si256(pvHLoad + i);
            }
        }
end:
        {
            /* extract vector containing last value from the column */
            __m256i cond_max;
            vH = _mm256_load_si256(pvHStore + offset);
            vHM = _mm256_load_si256(pvHMStore + offset);
            vHL = _mm256_load_si256(pvHLStore + offset);
            cond_max = _mm256_cmpgt_epi8(vH, vMaxH);
            vMaxH = _mm256_blendv_epi8(vMaxH, vH, cond_max);
            vMaxHM = _mm256_blendv_epi8(vMaxHM, vHM, cond_max);
            vMaxHL = _mm256_blendv_epi8(vMaxHL, vHL, cond_max);
        }
    }

    {
        /* extract last value from the column */
        int8_t tmp;
        for (k=0; k<position; ++k) {
            vMaxH  = shift (vMaxH);
            vMaxHM = shift (vMaxHM);
            vMaxHL = shift (vMaxHL);
        }
        /* max of last value in each column */
        tmp = (int8_t) _mm256_extract_epi8 (vMaxH, 31);
        if (tmp > score) {
            score = tmp;
            matches = (int8_t)_mm256_extract_epi8(vMaxHM, 31);
            length = (int8_t)_mm256_extract_epi8(vMaxHL, 31);
        }
    }

    /* max of last column */
    {
        vMaxH = vNegInf;
        vMaxHM = vNegInf;
        vMaxHL = vNegInf;

        for (i=0; i<segLen; ++i) {
            /* load the last stored values */
            __m256i vH = _mm256_load_si256(pvHStore + i);
            __m256i vHM = _mm256_load_si256(pvHMStore + i);
            __m256i vHL = _mm256_load_si256(pvHLStore + i);
            __m256i cond_max = _mm256_cmpgt_epi8(vH, vMaxH);
            vMaxH = _mm256_blendv_epi8(vMaxH, vH, cond_max);
            vMaxHM = _mm256_blendv_epi8(vMaxHM, vHM, cond_max);
            vMaxHL = _mm256_blendv_epi8(vMaxHL, vHL, cond_max);
        }

        /* max in vec */
        for (j=0; j<segWidth; ++j) {
            int8_t value = (int8_t) _mm256_extract_epi8(vMaxH, 31);
            if (value > score) {
                score = value;
                matches = (int8_t)_mm256_extract_epi8(vMaxHM, 31);
                length = (int8_t)_mm256_extract_epi8(vMaxHL, 31);
            }
            vMaxH = shift(vMaxH);
            vMaxHM = shift(vMaxHM);
            vMaxHL = shift(vMaxHL);
        }
    }

    if (_mm256_movemask_epi8(vSaturationCheck)) {
        result->saturated = 1;
        score = INT8_MAX;
        matches = 0;
        length = 0;
    }

    result->score = score;
    result->matches = matches;
    result->length = length;

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

