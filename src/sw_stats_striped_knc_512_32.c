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
#include "parasail/memory.h"
#include "parasail/internal_knc.h"
#include "parasail/matrices/blosum_map.h"

#define NEG_INF_32 (INT32_MIN/(int32_t)(2))

static inline __m512i insert(__m512i a, int32_t b, int imm) {
    __m512i_32_t tmp;
    tmp.m = a;
    tmp.v[imm] = b;
    return tmp.m;
}

static inline int32_t extract(__m512i a, int imm) {
    __m512i_32_t tmp;
    tmp.m = a;
    return tmp.v[imm];
}

#ifdef PARASAIL_TABLE
static inline void arr_store_si512(
        int *array,
        __m512i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[( 0*seglen+t)*dlen + d] = (int32_t)extract(vH,  0);
    array[( 1*seglen+t)*dlen + d] = (int32_t)extract(vH,  1);
    array[( 2*seglen+t)*dlen + d] = (int32_t)extract(vH,  2);
    array[( 3*seglen+t)*dlen + d] = (int32_t)extract(vH,  3);
    array[( 4*seglen+t)*dlen + d] = (int32_t)extract(vH,  4);
    array[( 5*seglen+t)*dlen + d] = (int32_t)extract(vH,  5);
    array[( 6*seglen+t)*dlen + d] = (int32_t)extract(vH,  6);
    array[( 7*seglen+t)*dlen + d] = (int32_t)extract(vH,  7);
    array[( 8*seglen+t)*dlen + d] = (int32_t)extract(vH,  8);
    array[( 9*seglen+t)*dlen + d] = (int32_t)extract(vH,  9);
    array[(10*seglen+t)*dlen + d] = (int32_t)extract(vH, 10);
    array[(11*seglen+t)*dlen + d] = (int32_t)extract(vH, 11);
    array[(12*seglen+t)*dlen + d] = (int32_t)extract(vH, 12);
    array[(13*seglen+t)*dlen + d] = (int32_t)extract(vH, 13);
    array[(14*seglen+t)*dlen + d] = (int32_t)extract(vH, 14);
    array[(15*seglen+t)*dlen + d] = (int32_t)extract(vH, 15);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME sw_stats_table_striped_knc_512_32
#else
#define FNAME sw_stats_striped_knc_512_32
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
    const int32_t n = 24; /* number of amino acids in table */
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m512i* const restrict vProfile  = parasail_memalign_m512i(64, n * segLen);
    __m512i* const restrict vProfileS = parasail_memalign_m512i(64, n * segLen);
    __m512i* restrict pvHStore        = parasail_memalign_m512i(64, segLen);
    __m512i* restrict pvHLoad         = parasail_memalign_m512i(64, segLen);
    __m512i* restrict pvHMStore       = parasail_memalign_m512i(64, segLen);
    __m512i* restrict pvHMLoad        = parasail_memalign_m512i(64, segLen);
    __m512i* restrict pvHLStore       = parasail_memalign_m512i(64, segLen);
    __m512i* restrict pvHLLoad        = parasail_memalign_m512i(64, segLen);
    __m512i* restrict pvEStore        = parasail_memalign_m512i(64, segLen);
    __m512i* restrict pvELoad         = parasail_memalign_m512i(64, segLen);
    __m512i* const restrict pvEM      = parasail_memalign_m512i(64, segLen);
    __m512i* const restrict pvEL      = parasail_memalign_m512i(64, segLen);
    int32_t score = NEG_INF_32;
    int32_t matches = NEG_INF_32;
    int32_t length = NEG_INF_32;
    __m512i vSegLimit = _mm512_set1_epi32(segWidth-1);
    __m512i permute_idx = _mm512_set_16to16_pi(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15);
    __mmask16 permute_mask = _mm512_cmplt_epi32_mask(permute_idx, vSegLimit);
    __m512i vGapO = _mm512_set1_epi32(open);
    __m512i vGapE = _mm512_set1_epi32(gap);
    __m512i vZero = _mm512_set1_epi32(0);
    __m512i vOne = _mm512_set1_epi32(1);
    __m512i vNegInf = _mm512_set1_epi32(NEG_INF_32);
    __m512i vMaxH = vZero;
    __m512i vMaxHM = vZero;
    __m512i vMaxHL = vZero;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif

    parasail_memset_m512i(pvHMStore, vZero, segLen);
    parasail_memset_m512i(pvHLStore, vZero, segLen);

    /* Generate query profile.
     * Rearrange query sequence & calculate the weight of match/mismatch.
     * Don't alias. */
    {
        int32_t index = 0;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                __m512i_32_t t;
                __m512i_32_t s;
                j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = j >= s1Len ? 0 : matrix[k][parasail_blosum_map[(unsigned char)s1[j]]];
                    s.v[segNum] = j >= s1Len ? 0 : (k == parasail_blosum_map[(unsigned char)s1[j]]);
                    j += segLen;
                }
                _mm512_store_epi32(&vProfile[index], t.m);
                _mm512_store_epi32(&vProfileS[index], s.m);
                ++index;
            }
        }
    }

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m512i_32_t h;
            __m512i_32_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = 0;
                e.v[segNum] = NEG_INF_32;
            }
            _mm512_store_epi32(&pvHStore[index], h.m);
            _mm512_store_epi32(&pvEStore[index], e.m);
            ++index;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m512i vE;
        __m512i vEM;
        __m512i vEL;
        __m512i vF;
        __m512i vFM;
        __m512i vFL;
        __m512i vH;
        __m512i vHM;
        __m512i vHL;
        const __m512i* vP = NULL;
        const __m512i* vPS = NULL;
        __m512i* pv = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        vF = vZero;
        vFM = vZero;
        vFL = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = _mm512_mask_permutevar_epi32(vZero, permute_mask, permute_idx, _mm512_load_epi32(&pvHStore[segLen - 1]));
        vHM = _mm512_mask_permutevar_epi32(vZero, permute_mask, permute_idx, _mm512_load_epi32(&pvHMStore[segLen - 1]));
        vHL = _mm512_mask_permutevar_epi32(vZero, permute_mask, permute_idx, _mm512_load_epi32(&pvHLStore[segLen - 1]));

        /* Correct part of the vProfile */
        vP = vProfile + parasail_blosum_map[(unsigned char)s2[j]] * segLen;
        vPS = vProfileS + parasail_blosum_map[(unsigned char)s2[j]] * segLen;

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
            __mmask16 case1not;
            __mmask16 case2not;
            __mmask16 case2;
            __mmask16 case3;
            __mmask16 cond_zero;

            vH = _mm512_add_epi32(vH, _mm512_load_epi32(vP + i));
            vE = _mm512_load_epi32(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            case1not = _mm512_kor(
                    _mm512_cmplt_epi32_mask(vH,vF),
                    _mm512_cmplt_epi32_mask(vH,vE));
            case2not = _mm512_cmplt_epi32_mask(vF,vE);
            case2 = _mm512_kandn(case2not,case1not);
            case3 = _mm512_kand(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm512_max_epi32(vH, vE);
            vH = _mm512_max_epi32(vH, vF);
            vH = _mm512_max_epi32(vH, vZero);
            /* Save vH values. */
            _mm512_store_epi32(pvHStore + i, vH);
            cond_zero = _mm512_cmpeq_epi32_mask(vH, vZero);

            /* calculate vM */
            vEM = _mm512_load_epi32(pvEM + i);
            vHM = _mm512_mask_blend_epi32(case1not,
                    _mm512_add_epi32(vHM, _mm512_load_epi32(vPS + i)),
                    _mm512_mask_blend_epi32(case2, vEM, vFM));
            vHM = _mm512_mask_blend_epi32(cond_zero, vHM, vZero);
            _mm512_store_epi32(pvHMStore + i, vHM);

            /* calculate vL */
            vEL = _mm512_load_epi32(pvEL + i);
            vHL = _mm512_mask_blend_epi32(case1not,
                    _mm512_add_epi32(vHL, vOne),
                    _mm512_mask_blend_epi32(case2,
                        _mm512_add_epi32(vEL, vOne),
                        _mm512_add_epi32(vFL, vOne)));
            vHL = _mm512_mask_blend_epi32(cond_zero, vHL, vZero);
            _mm512_store_epi32(pvHLStore + i, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si512(result->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si512(result->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si512(result->score_table, vH, i, segLen, j, s2Len);
#endif

            /* update max vector seen so far */
            {
                __mmask16 cond_max = _mm512_cmpgt_epi32_mask(vH, vMaxH);
                vMaxH  = _mm512_mask_blend_epi32(cond_max, vMaxH, vH);
                vMaxHM = _mm512_mask_blend_epi32(cond_max, vMaxHM, vHM);
                vMaxHL = _mm512_mask_blend_epi32(cond_max, vMaxHL, vHL);
            }

            /* Update vE value. */
            vH = _mm512_sub_epi32(vH, vGapO);
            vE = _mm512_sub_epi32(vE, vGapE);
            vE = _mm512_max_epi32(vE, vH);
            _mm512_store_epi32(pvEStore + i, vE);
            _mm512_store_epi32(pvEM + i, vHM);
            _mm512_store_epi32(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm512_sub_epi32(vF, vGapE);
            vF = _mm512_max_epi32(vF, vH);
            vFM = vHM;
            vFL = vHL;

            /* Load the next vH. */
            vH = _mm512_load_epi32(pvHLoad + i);
            vHM = _mm512_load_epi32(pvHMLoad + i);
            vHL = _mm512_load_epi32(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            __m512i vHp = _mm512_mask_permutevar_epi32(
                    vZero, permute_mask, permute_idx,
                    _mm512_load_epi32(&pvHLoad[segLen - 1]));
            vF = _mm512_mask_permutevar_epi32(
                    vZero, permute_mask, permute_idx, vF);
            vFM = _mm512_mask_permutevar_epi32(
                    vZero, permute_mask, permute_idx, vFM);
            vFL = _mm512_mask_permutevar_epi32(
                    vZero, permute_mask, permute_idx, vFL);
            for (i=0; i<segLen; ++i) {
                __mmask16 case1not;
                __mmask16 case2not;
                __mmask16 case2;
                __mmask16 cond_zero;

                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm512_add_epi32(vHp, _mm512_load_epi32(vP + i));
                vE = _mm512_load_epi32(pvELoad + i);
                case1not = _mm512_kor(
                        _mm512_cmplt_epi32_mask(vHp,vF),
                        _mm512_cmplt_epi32_mask(vHp,vE));
                case2not = _mm512_cmplt_epi32_mask(vF,vE);
                case2 = _mm512_kandn(case2not,case1not);

                vH = _mm512_load_epi32(pvHStore + i);
                vH = _mm512_max_epi32(vH,vF);
                _mm512_store_epi32(pvHStore + i, vH);
                cond_zero = _mm512_cmpeq_epi32_mask(vH, vZero);

                vHM = _mm512_load_epi32(pvHMStore + i);
                vHM = _mm512_mask_blend_epi32(case2, vHM, vFM);
                vHM = _mm512_mask_blend_epi32(cond_zero, vHM, vZero);
                _mm512_store_epi32(pvHMStore + i, vHM);
                _mm512_store_epi32(pvEM + i, vHM);

                vHL = _mm512_load_epi32(pvHLStore + i);
                vHL = _mm512_mask_blend_epi32(case2, vHL, _mm512_add_epi32(vFL,vOne));
                vHL = _mm512_mask_blend_epi32(cond_zero, vHL, vZero);
                _mm512_store_epi32(pvHLStore + i, vHL);
                _mm512_store_epi32(pvEL + i, vHL);

#ifdef PARASAIL_TABLE
                arr_store_si512(result->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si512(result->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si512(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm512_sub_epi32(vH, vGapO);
                vF = _mm512_sub_epi32(vF, vGapE);
                if (! _mm512_mask2int(_mm512_cmpgt_epi32_mask(vF, vH))) goto end;
                /*vF = _mm512_max_epi32(vF, vH);*/
                vFM = vHM;
                vFL = vHL;
                vHp = _mm512_load_epi32(pvHLoad + i);
            }
        }
end:
        {
        }
    }

    /* max in vec */
    for (j=0; j<segWidth; ++j) {
        int32_t value = (int32_t) extract(vMaxH, 15);
        if (value > score) {
            score = value;
            matches = (int32_t) extract(vMaxHM, 15);
            length = (int32_t) extract(vMaxHL, 15);
        }
        vMaxH = _mm512_permutevar_epi32(permute_idx, vMaxH);
        vMaxHM = _mm512_permutevar_epi32(permute_idx, vMaxHM);
        vMaxHL = _mm512_permutevar_epi32(permute_idx, vMaxHL);
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

