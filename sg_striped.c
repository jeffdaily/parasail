#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include <emmintrin.h>

#ifdef ALIGN_EXTRA
#include "align/align_debug.h"
#else
#include "align/align.h"
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
    array[(0*seglen+t)*dlen + d] = extract(vH, 0);
    array[(1*seglen+t)*dlen + d] = extract(vH, 1);
    array[(2*seglen+t)*dlen + d] = extract(vH, 2);
    array[(3*seglen+t)*dlen + d] = extract(vH, 3);
    array[(4*seglen+t)*dlen + d] = extract(vH, 4);
    array[(5*seglen+t)*dlen + d] = extract(vH, 5);
    array[(6*seglen+t)*dlen + d] = extract(vH, 6);
    array[(7*seglen+t)*dlen + d] = extract(vH, 7);
}
#endif

#ifdef ALIGN_EXTRA
#define FNAME sg_striped_debug
#else
#define FNAME sg_striped
#endif

int FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix
#ifdef ALIGN_EXTRA
        , int * const restrict score_table
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

    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = 7 - (s1Len - 1) / segLen;

    /* the max alignment score */
    int max = NEG_INF;

    /* Define 16 byte 0 vector. */
    __m128i vOne = _mm_set1_epi16(1);

    __m128i* pvHStore = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvHLoad = (__m128i*) calloc(segLen, sizeof(__m128i));
    __m128i* pvE = (__m128i*) calloc(segLen, sizeof(__m128i));

    /* Generate query profile rearrange query sequence & calculate the weight
     * of match/mismatch */
    {
        int16_t *t = (int16_t*)vProfile;
        for (nt=0; nt<n; ++nt) {
            for (i=0; i<segLen; ++i) {
                int32_t j = i;
                for (segNum=0; segNum<8; ++segNum) {
                    *t++ = matrix[nt*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    j += segLen;
                }
            }
        }
    }

    /* initialize E */
    {
        int16_t *e = (int16_t*)pvE;
        for (i=0; i<segLen; ++i) {
            for (segNum=0; segNum<8; ++segNum) {
                *e++ = -open;
            }
        }
    }

    /* 16 byte insertion begin vector */
    __m128i vGapO = _mm_set1_epi16(open);

    /* 16 byte insertion extension vector */
    __m128i vGapE = _mm_set1_epi16(gap);

    __m128i initialF = _mm_set_epi16(
            -open-7*segLen*gap,
            -open-6*segLen*gap,
            -open-5*segLen*gap,
            -open-4*segLen*gap,
            -open-3*segLen*gap,
            -open-2*segLen*gap,
            -open-1*segLen*gap,
            -open-0*segLen*gap);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vE;
        /* Initialize F value to 0.  Any errors to vH values will be corrected
         * in the Lazy_F loop.  */
        __m128i vF = initialF;

        /* load final segment of pvHStore and shift left by 2 bytes */
        __m128i vH = _mm_slli_si128(pvHStore[segLen - 1], 2);

        /* Correct part of the vProfile */
        const __m128i* vP = vProfile + MAP_BLOSUM_[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        __m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vH = _mm_adds_epi16(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi16(vH, vE);
            vH = _mm_max_epi16(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);
#ifdef ALIGN_EXTRA
            arr_store_si128(score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = _mm_subs_epi16(vH, vGapO);
            vE = _mm_subs_epi16(vE, vGapE);
            vE = _mm_max_epi16(vE, vH);
            _mm_store_si128(pvE + i, vE);

            /* Update vF value. */
            vF = _mm_subs_epi16(vF, vGapE);
            vF = _mm_max_epi16(vF, vH);

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<8; ++k) {
            vF = _mm_slli_si128(vF, 2);
            vF = _mm_insert_epi16(vF, -open, 0);
            for (i=0; i<segLen; ++i) {
                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi16(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
#ifdef ALIGN_EXTRA
                arr_store_si128(score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm_subs_epi16(vH, vGapO);
                vF = _mm_subs_epi16(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi16(vF, vH))) goto end;
                vF = _mm_max_epi16(vF, vH);
            }
        }
end:
        {
            /* extract last value from the column */
            vH = _mm_load_si128(pvHStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128 (vH, 2);
            }
            /* max of last value in each column */
            int16_t tmp = (signed short) _mm_extract_epi16 (vH, 7);
            if (tmp > max) {
                max = tmp;
            }
        }
    }

    /* max of last column */
    __m128i vNegInf = _mm_set1_epi16(NEG_INF);
    __m128i vMaxLastColH = vNegInf;
    __m128i vQIndex = _mm_set_epi16(
            7*segLen,
            6*segLen,
            5*segLen,
            4*segLen,
            3*segLen,
            2*segLen,
            1*segLen,
            0*segLen);
    __m128i vQLimit = _mm_set1_epi16(s1Len);

    for (i=0; i<segLen; ++i) {
        /* load the last stored values */
        __m128i vH = _mm_load_si128(pvHStore + i);
        /* mask off the values that were padded */
        __m128i vMask = _mm_cmplt_epi16(vQIndex, vQLimit);
        vH = _mm_or_si128(_mm_and_si128(vMask, vH),
                          _mm_andnot_si128(vMask, vNegInf));
        __m128i cond = _mm_cmplt_epi16(vH, vMaxLastColH);
        vMaxLastColH = _mm_or_si128(_mm_andnot_si128(cond, vH),
                                    _mm_and_si128(cond, vMaxLastColH));
        vQIndex = _mm_adds_epi16(vQIndex, vOne);
    }

    /* max in vec */
    for (j=0; j<8; ++j) {
        int16_t value = (signed short) _mm_extract_epi16(vMaxLastColH, j);
        if (value > max) {
            max = value;
        }
    }

    free(vProfile);
    free(pvHStore);
    free(pvHLoad);
    free(pvE);

    return max;
}

