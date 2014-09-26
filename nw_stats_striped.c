#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include <emmintrin.h>

#ifdef ALIGN_EXTRA
#include "align/align_striped_128_16_debug.h"
#else
#include "align/align_striped_128_16.h"
#endif
#include "blosum/blosum_map.h"


#if ALIGN_EXTRA
static inline void arr_store_si128(
        int *array,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[(0*seglen+t)*dlen + d] = _mm_extract_epi16(vH, 0);
    array[(1*seglen+t)*dlen + d] = _mm_extract_epi16(vH, 1);
    array[(2*seglen+t)*dlen + d] = _mm_extract_epi16(vH, 2);
    array[(3*seglen+t)*dlen + d] = _mm_extract_epi16(vH, 3);
    array[(4*seglen+t)*dlen + d] = _mm_extract_epi16(vH, 4);
    array[(5*seglen+t)*dlen + d] = _mm_extract_epi16(vH, 5);
    array[(6*seglen+t)*dlen + d] = _mm_extract_epi16(vH, 6);
    array[(7*seglen+t)*dlen + d] = _mm_extract_epi16(vH, 7);
}
#endif

#ifdef ALIGN_EXTRA
#define FNAME nw_stats_striped_128_16_debug
#else
#define FNAME nw_stats_striped_128_16
#endif

int FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict matches, int * const restrict length
#ifdef ALIGN_EXTRA
        , int * const restrict score_table
        , int * const restrict match_table
        , int * const restrict length_table
#endif
        )
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t nt = 0;
    int32_t segNum = 0;
    const int32_t n = 24; /* number of amino acids in table */
    int32_t segLen = (s1Len + 7) / 8;
    __m128i* vProfile = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    __m128i* vProfileS = (__m128i*)malloc(n * segLen * sizeof(__m128i));

    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = 7 - (s1Len - 1) / segLen;

    /* the max alignment score */
    int last_value;

    /* Define 16 byte 0 vector. */
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi16(1);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHMStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHMLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvEStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvELoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvEM = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvEL = (__m128i*) calloc(segLen, sizeof(__m128i));
    int16_t* boundary = (int16_t*) calloc(s2Len+1, sizeof(int16_t));

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    {
        int16_t *t = (int16_t*)vProfile;
        int16_t *s = (int16_t*)vProfileS;
        for (nt=0; nt<n; ++nt) {
            for (i=0; i<segLen; ++i) {
                int32_t j = i;
                for (segNum=0; segNum<8; ++segNum) {
                    *t++ = matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    *s++ = (nt == MAP_BLOSUM_[(unsigned char)s1[j]]);
                    j += segLen;
                }
            }
        }
    }

    /* initialize H and E */
    {
        int16_t *h = (int16_t*)pvHStore;
        int16_t *e = (int16_t*)pvEStore;
        for (i=0; i<segLen; ++i) {
            int32_t j = i;
            for (segNum=0; segNum<8; ++segNum) {
                *h = -open-gap*(segNum*segLen+i);
                *e = *h-open;
                j += segLen;
                ++h;
                ++e;
            }
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            boundary[i] = -open-gap*(i-1);
        }
    }

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    __m128i initialF = _mm_set_epi16(
            -open-open-7*segLen*gap,
            -open-open-6*segLen*gap,
            -open-open-5*segLen*gap,
            -open-open-4*segLen*gap,
            -open-open-3*segLen*gap,
            -open-open-2*segLen*gap,
            -open-open-1*segLen*gap,
            -open-open-0*segLen*gap);
    initialF = _mm_adds_epi16(initialF, vGapE);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vE;
        __m128i vEM;
        __m128i vEL;
        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        initialF = _mm_subs_epi16(initialF, vGapE);
        __m128i vF = initialF;
        __m128i vFM = vZero;
        __m128i vFL = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);
        __m128i vHM = _mm_slli_si128(pvHMStore[segLen - 1], 2);
        __m128i vHL = _mm_slli_si128(pvHLStore[segLen - 1], 2);

        /* insert upper boundary condition */
        vH = _mm_insert_epi16(vH, boundary[j], 0);

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;
        const __m128i* vPS = vProfileS + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
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
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            __m128i case1not = _mm_or_si128(
                    _mm_cmplt_epi16(vH,vF),_mm_cmplt_epi16(vH,vE));
            __m128i case2not = _mm_cmplt_epi16(vF,vE);
            __m128i case2 = _mm_andnot_si128(case2not,case1not);
            __m128i case3 = _mm_and_si128(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);

            /* calculate vM */
            vEM = _mm_load_si128(pvEM + i);
            vHM = _mm_andnot_si128(case1not,
                    _mm_add_epi16(vHM, _mm_load_si128(vPS + i)));
            vHM = _mm_or_si128(vHM, _mm_and_si128(case2, vFM));
            vHM = _mm_or_si128(vHM, _mm_and_si128(case3, vEM));
            _mm_store_si128(pvHMStore + i, vHM);

            /* calculate vL */
            vEL = _mm_load_si128(pvEL + i);
            vHL = _mm_andnot_si128(case1not, _mm_add_epi16(vHL, vOne));
            vHL = _mm_or_si128(vHL, _mm_and_si128(case2,
                        _mm_add_epi16(vFL, vOne)));
            vHL = _mm_or_si128(vHL, _mm_and_si128(case3,
                        _mm_add_epi16(vEL, vOne)));
            _mm_store_si128(pvHLStore + i, vHL);

#ifdef ALIGN_EXTRA
            arr_store_si128(match_table, vHM, i, segLen, j, s2Len);
            arr_store_si128(length_table, vHL, i, segLen, j, s2Len);
            arr_store_si128(score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = _mm_subs_epi16(vH, vGapO);
            vE = _mm_subs_epi16(vE, vGapE);
            vE = _mm_max_epi16(vE, vH);
            _mm_store_si128(pvEStore + i, vE);
            _mm_store_si128(pvEM + i, vHM);
            _mm_store_si128(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm_subs_epi16(vF, vGapE);
            vF = _mm_max_epi16(vF, vH);
            vFM = vHM;
            vFL = vHL;

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
            vHM = _mm_load_si128(pvHMLoad + i);
            vHL = _mm_load_si128(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<8; ++k) {
            __m128i vHp = _mm_slli_si128(pvHLoad[segLen - 1], 2);
            vHp = _mm_insert_epi16(vHp, boundary[j], 0);
            vF = _mm_slli_si128(vF, 2);
            vF = _mm_insert_epi16(vF, boundary[j+1]-open, 0);
            vFM = _mm_slli_si128(vFM, 2);
            vFL = _mm_slli_si128(vFL, 2);
            for (i=0; i<segLen; ++i) {
                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm_adds_epi16(vHp, _mm_load_si128(vP + i));
                vE = _mm_load_si128(pvELoad + i);
                __m128i case1not = _mm_or_si128(
                        _mm_cmplt_epi16(vHp,vF),_mm_cmplt_epi16(vHp,vE));
                __m128i case2not = _mm_cmplt_epi16(vF,vE);
                __m128i case2 = _mm_andnot_si128(case2not,case1not);

                vHM = _mm_load_si128(pvHMStore + i);
                vHM = _mm_andnot_si128(case2, vHM);
                vHM = _mm_or_si128(vHM, _mm_and_si128(case2, vFM));
                _mm_store_si128(pvHMStore + i, vHM);
                _mm_store_si128(pvEM + i, vHM);

                vHL = _mm_load_si128(pvHLStore + i);
                vHL = _mm_andnot_si128(case2, vHL);
                vHL = _mm_or_si128(vHL, _mm_and_si128(case2,
                            _mm_add_epi16(vFL,vOne)));
                _mm_store_si128(pvHLStore + i, vHL);
                _mm_store_si128(pvEL + i, vHL);

                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi16(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
#ifdef ALIGN_EXTRA
                arr_store_si128(match_table, vHM, i, segLen, j, s2Len);
                arr_store_si128(length_table, vHL, i, segLen, j, s2Len);
                arr_store_si128(score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm_subs_epi16(vH, vGapO);
                vF = _mm_subs_epi16(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) goto end;
                vF = _mm_max_epi16(vF, vH);
                vFM = vHM;
                vFL = vHL;
                vHp = _mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
        }
    }

    /* extract last value from the last column */
    __m128i vH = _mm_load_si128(pvHStore + offset);
    __m128i vHM = _mm_load_si128(pvHMStore + offset);
    __m128i vHL = _mm_load_si128(pvHLStore + offset);
    for (k=0; k<position; ++k) {
        vH = _mm_slli_si128 (vH, 2);
        vHM = _mm_slli_si128 (vHM, 2);
        vHL = _mm_slli_si128 (vHL, 2);
    }
    last_value = (int16_t) _mm_extract_epi16 (vH, 7);
    *matches = (int16_t) _mm_extract_epi16 (vHM, 7);
    *length = (int16_t) _mm_extract_epi16 (vHL, 7);

    free(vProfile);
    free(vProfileS);
    free(pvHStore);
    free(pvHLoad);
    free(pvHMStore);
    free(pvHMLoad);
    free(pvHLStore);
    free(pvHLLoad);
    free(pvEStore);
    free(pvELoad);
    free(pvEM);
    free(pvEL);
    free(boundary);

    return last_value;
}

