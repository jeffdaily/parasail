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

#define NEG_INF (INT32_MIN/(int32_t)(2))

#define _mm256_cmplt_epi32_rpl(a,b) _mm256_cmpgt_epi32(b,a)

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)


#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[(0*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 0);
    array[(1*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 1);
    array[(2*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 2);
    array[(3*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 3);
    array[(4*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 4);
    array[(5*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 5);
    array[(6*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 6);
    array[(7*seglen+t)*dlen + d] = (int32_t)_mm256_extract_epi32(vH, 7);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME sw_stats_table_striped_avx2_256_32
#else
#define FNAME sw_stats_striped_avx2_256_32
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
    const int32_t segWidth = 8; /* number of values in vector unit */
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
    __m256i* const restrict pvEM      = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvES      = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvEL      = parasail_memalign___m256i(32, segLen);
    __m256i vGapO = _mm256_set1_epi32(open);
    __m256i vGapE = _mm256_set1_epi32(gap);
    __m256i vZero = _mm256_setzero_si256();
    __m256i vOne = _mm256_set1_epi32(1);
    int32_t score = NEG_INF;
    int32_t matches = NEG_INF;
    int32_t similar = NEG_INF;
    int32_t length = NEG_INF;
    
    __m256i vMaxH = vZero;
    __m256i vMaxHM = vZero;
    __m256i vMaxHS = vZero;
    __m256i vMaxHL = vZero;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif

    parasail_memset___m256i(pvHMStore, vZero, segLen);
    parasail_memset___m256i(pvHSStore, vZero, segLen);
    parasail_memset___m256i(pvHLStore, vZero, segLen);

    /* Generate query profile.
     * Rearrange query sequence & calculate the weight of match/mismatch.
     * Don't alias. */
    {
        int32_t index = 0;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                __m256i_32_t p;
                __m256i_32_t m;
                __m256i_32_t s;
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
            __m256i_32_t h;
            __m256i_32_t e;
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
        vF = vZero;
        vFM = vZero;
        vFS = vZero;
        vFL = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = _mm256_slli_si256_rpl(pvHStore[segLen - 1], 4);
        vHM = _mm256_slli_si256_rpl(pvHMStore[segLen - 1], 4);
        vHS = _mm256_slli_si256_rpl(pvHSStore[segLen - 1], 4);
        vHL = _mm256_slli_si256_rpl(pvHLStore[segLen - 1], 4);

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

            vH = _mm256_add_epi32(vH, _mm256_load_si256(vP + i));
            vE = _mm256_load_si256(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            case1not = _mm256_or_si256(
                    _mm256_cmplt_epi32_rpl(vH,vF),_mm256_cmplt_epi32_rpl(vH,vE));
            case2not = _mm256_cmplt_epi32_rpl(vF,vE);
            case2 = _mm256_andnot_si256(case2not,case1not);
            case3 = _mm256_and_si256(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm256_max_epi32(vH, vE);
            vH = _mm256_max_epi32(vH, vF);
            vH = _mm256_max_epi32(vH, vZero);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);
            cond_zero = _mm256_cmpeq_epi32(vH, vZero);

            /* calculate vM */
            vEM = _mm256_load_si256(pvEM + i);
            vHM = _mm256_blendv_epi8(
                    _mm256_add_epi32(vHM, _mm256_load_si256(vPM + i)),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, vFM),
                        _mm256_and_si256(case3, vEM)),
                    case1not);
            vHM = _mm256_andnot_si256(cond_zero, vHM);
            _mm256_store_si256(pvHMStore + i, vHM);

            /* calculate vS */
            vES = _mm256_load_si256(pvES + i);
            vHS = _mm256_blendv_epi8(
                    _mm256_add_epi32(vHS, _mm256_load_si256(vPS + i)),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, vFS),
                        _mm256_and_si256(case3, vES)),
                    case1not);
            vHS = _mm256_andnot_si256(cond_zero, vHS);
            _mm256_store_si256(pvHSStore + i, vHS);

            /* calculate vL */
            vEL = _mm256_load_si256(pvEL + i);
            vHL = _mm256_blendv_epi8(
                    _mm256_add_epi32(vHL, vOne),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, _mm256_add_epi32(vFL, vOne)),
                        _mm256_and_si256(case3, _mm256_add_epi32(vEL, vOne))),
                    case1not);
            vHL = _mm256_andnot_si256(cond_zero, vHL);
            _mm256_store_si256(pvHLStore + i, vHL);
            
#ifdef PARASAIL_TABLE
            arr_store_si256(result->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si256(result->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si256(result->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif
            /* update max vector seen so far */
            {
                __m256i cond_max = _mm256_cmpgt_epi32(vH, vMaxH);
                vMaxH = _mm256_blendv_epi8(vMaxH, vH,  cond_max);
                vMaxHM = _mm256_blendv_epi8(vMaxHM, vHM, cond_max);
                vMaxHS = _mm256_blendv_epi8(vMaxHS, vHS, cond_max);
                vMaxHL = _mm256_blendv_epi8(vMaxHL, vHL, cond_max);
            }

            /* Update vE value. */
            vH = _mm256_sub_epi32(vH, vGapO);
            vE = _mm256_sub_epi32(vE, vGapE);
            vE = _mm256_max_epi32(vE, vH);
            _mm256_store_si256(pvEStore + i, vE);
            _mm256_store_si256(pvEM + i, vHM);
            _mm256_store_si256(pvES + i, vHS);
            _mm256_store_si256(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm256_sub_epi32(vF, vGapE);
            vF = _mm256_max_epi32(vF, vH);
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
            __m256i vHp = _mm256_slli_si256_rpl(pvHLoad[segLen - 1], 4);
            vF = _mm256_slli_si256_rpl(vF, 4);
            vFM = _mm256_slli_si256_rpl(vFM, 4);
            vFS = _mm256_slli_si256_rpl(vFS, 4);
            vFL = _mm256_slli_si256_rpl(vFL, 4);
            for (i=0; i<segLen; ++i) {
                __m256i case1not;
                __m256i case2not;
                __m256i case2;
                __m256i cond_zero;

                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm256_add_epi32(vHp, _mm256_load_si256(vP + i));
                vE = _mm256_load_si256(pvELoad + i);
                case1not = _mm256_or_si256(
                        _mm256_cmplt_epi32_rpl(vHp,vF),_mm256_cmplt_epi32_rpl(vHp,vE));
                case2not = _mm256_cmplt_epi32_rpl(vF,vE);
                case2 = _mm256_andnot_si256(case2not,case1not);

                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi32(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
                cond_zero = _mm256_cmpeq_epi32(vH, vZero);

                vHM = _mm256_load_si256(pvHMStore + i);
                vHM = _mm256_blendv_epi8(vHM, vFM, case2);
                vHM = _mm256_andnot_si256(cond_zero, vHM);
                _mm256_store_si256(pvHMStore + i, vHM);
                _mm256_store_si256(pvEM + i, vHM);

                vHS = _mm256_load_si256(pvHSStore + i);
                vHS = _mm256_blendv_epi8(vHS, vFS, case2);
                vHS = _mm256_andnot_si256(cond_zero, vHS);
                _mm256_store_si256(pvHSStore + i, vHS);
                _mm256_store_si256(pvES + i, vHS);

                vHL = _mm256_load_si256(pvHLStore + i);
                vHL = _mm256_blendv_epi8(vHL, _mm256_add_epi32(vFL,vOne), case2);
                vHL = _mm256_andnot_si256(cond_zero, vHL);
                _mm256_store_si256(pvHLStore + i, vHL);
                _mm256_store_si256(pvEL + i, vHL);

#ifdef PARASAIL_TABLE
                arr_store_si256(result->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si256(result->similar_table, vHS, i, segLen, j, s2Len);
                arr_store_si256(result->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm256_sub_epi32(vH, vGapO);
                vF = _mm256_sub_epi32(vF, vGapE);
                if (! _mm256_movemask_epi8(_mm256_cmpgt_epi32(vF, vH))) goto end;
                /*vF = _mm256_max_epi32(vF, vH);*/
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
        int32_t value = (int32_t) _mm256_extract_epi32(vMaxH, 7);
        if (value > score) {
            score = value;
            matches = (int32_t)_mm256_extract_epi32(vMaxHM, 7);
            similar = (int32_t)_mm256_extract_epi32(vMaxHS, 7);
            length = (int32_t)_mm256_extract_epi32(vMaxHL, 7);
        }
        vMaxH = _mm256_slli_si256_rpl(vMaxH, 4);
        vMaxHM = _mm256_slli_si256_rpl(vMaxHM, 4);
        vMaxHS = _mm256_slli_si256_rpl(vMaxHS, 4);
        vMaxHL = _mm256_slli_si256_rpl(vMaxHL, 4);
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


