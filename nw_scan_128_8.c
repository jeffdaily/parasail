#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include <emmintrin.h>
#include <smmintrin.h>

#ifdef ALIGN_EXTRA
#include "align/align_scan_128_8_debug.h"
#else
#include "align/align_scan_128_8.h"
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
    array[(0*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 0);
    array[(1*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 1);
    array[(2*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 2);
    array[(3*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 3);
    array[(4*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 4);
    array[(5*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 5);
    array[(6*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 6);
    array[(7*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 7);
    array[(8*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 8);
    array[(9*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 9);
    array[(10*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 10);
    array[(11*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 11);
    array[(12*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 12);
    array[(13*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 13);
    array[(14*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 14);
    array[(15*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 15);
}
#endif

#define PARALLEL_PREFIX_OP(vFt, gap, segLen)            \
{                                                       \
    union {                                             \
        __m128i m;                                      \
        int8_t v[16];                                   \
    } tmp;                                              \
    tmp.m = vFt;                                        \
    tmp.v[1]  = MAX(tmp.v[0] -segLen*gap, tmp.v[1]);    \
    tmp.v[2]  = MAX(tmp.v[1] -segLen*gap, tmp.v[2]);    \
    tmp.v[3]  = MAX(tmp.v[2] -segLen*gap, tmp.v[3]);    \
    tmp.v[4]  = MAX(tmp.v[3] -segLen*gap, tmp.v[4]);    \
    tmp.v[5]  = MAX(tmp.v[4] -segLen*gap, tmp.v[5]);    \
    tmp.v[6]  = MAX(tmp.v[5] -segLen*gap, tmp.v[6]);    \
    tmp.v[7]  = MAX(tmp.v[6] -segLen*gap, tmp.v[7]);    \
    tmp.v[8]  = MAX(tmp.v[7] -segLen*gap, tmp.v[8]);    \
    tmp.v[9]  = MAX(tmp.v[8] -segLen*gap, tmp.v[9]);    \
    tmp.v[10] = MAX(tmp.v[9] -segLen*gap, tmp.v[10]);   \
    tmp.v[11] = MAX(tmp.v[10]-segLen*gap, tmp.v[11]);   \
    tmp.v[12] = MAX(tmp.v[11]-segLen*gap, tmp.v[12]);   \
    tmp.v[13] = MAX(tmp.v[12]-segLen*gap, tmp.v[13]);   \
    tmp.v[14] = MAX(tmp.v[13]-segLen*gap, tmp.v[14]);   \
    tmp.v[15] = MAX(tmp.v[14]-segLen*gap, tmp.v[15]);   \
    vFt = tmp.m;                                        \
}

#ifdef ALIGN_EXTRA
#define FNAME nw_scan_128_8_debug
#else
#define FNAME nw_scan_128_8
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
    const int32_t n = 24; /* number of amino acids in table */
    const int32_t segWidth = 16; /* number of values in vector unit */
    int32_t segNum = 0;
    int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m128i* pvP = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    __m128i* pvE = (__m128i*)malloc(segLen * sizeof(__m128i));
    __m128i* pvHt = (__m128i*)malloc(segLen * sizeof(__m128i));
    __m128i* pvFt = (__m128i*)malloc(segLen * sizeof(__m128i));
    __m128i* pvH = (__m128i*)malloc(segLen * sizeof(__m128i));
    int8_t* boundary = (int8_t*)malloc((s2Len+1) * sizeof(int8_t));
    __m128i vGapO = _mm_set1_epi8(open);
    __m128i vGapE = _mm_set1_epi8(gap);
    __m128i vSaturationCheck = _mm_setzero_si128();
    __m128i vNegLimit = _mm_set1_epi8(INT8_MIN);
    __m128i vPosLimit = _mm_set1_epi8(INT8_MAX);
    int8_t score = 0;

    /* Generate query profile.
     * Rearrange query sequence & calculate the weight of match/mismatch */
    {
        int8_t *t = (int8_t*)pvP;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                int32_t j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    *t++ = matrix[k*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    j += segLen;
                }
            }
        }
    }

    /* initialize H and E */
    {
        int8_t *h = (int8_t*)pvH;
        int8_t *e = (int8_t*)pvE;
        for (i=0; i<segLen; ++i) {
            for (segNum=0; segNum<segWidth; ++segNum) {
                *h = -open-gap*(segNum*segLen+i);
                *e = NEG_INF_8;
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

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vE;
        __m128i vHt;
        __m128i vFt;
        __m128i vH;
        __m128i vHp;
        __m128i *pvW;
        __m128i vW;

#define LOOP_FUSION 1
#if LOOP_FUSION
        /* calculate E */
        /* calculate Ht */
        vHp = _mm_slli_si128(_mm_load_si128(pvH+(segLen-1)), 1);
        vHp = _mm_insert_epi8(vHp, boundary[j], 0);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            vE = _mm_max_epi8(
                    _mm_subs_epi8(vE, vGapE),
                    _mm_subs_epi8(vH, vGapO));
            vHt = _mm_max_epi8(
                    _mm_adds_epi8(vHp, vW),
                    vE);
            _mm_store_si128(pvE+i, vE);
            _mm_store_si128(pvHt+i, vHt);
            vHp = vH;
        }
#else
        /* calculate E */
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vE = _mm_max_epi8(
                    _mm_subs_epi8(vE, vGapE),
                    _mm_subs_epi8(vH, vGapO));
            _mm_store_si128(pvE+i, vE);
        }

        /* calculate Ht */
        vH = _mm_slli_si128(_mm_load_si128(pvH+(segLen-1)), 1);
        vH = _mm_insert_epi8(vH, boundary[j], 0);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            vHt = _mm_max_epi8(
                    _mm_adds_epi8(vH, vW),
                    vE);
            vH = _mm_load_si128(pvH+i);
            _mm_store_si128(pvHt+i, vHt);
        }
#endif

        /* calculate Ft */
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 1);
        vHt = _mm_insert_epi8(vHt, boundary[j+1], 0);
        vFt = _mm_set1_epi8(NEG_INF_8);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_max_epi8(
                    _mm_subs_epi8(vFt, vGapE),
                    vHt);
            vHt = _mm_load_si128(pvHt+i);
        }
        PARALLEL_PREFIX_OP(vFt, gap, segLen)
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 1);
        vHt = _mm_insert_epi8(vHt, boundary[j+1], 0);
        vFt = _mm_slli_si128(vFt, 1);
        vFt = _mm_insert_epi8(vFt, NEG_INF_8, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_max_epi8(
                    _mm_subs_epi8(vFt, vGapE),
                    vHt);
            vHt = _mm_load_si128(pvHt+i);
            _mm_store_si128(pvFt+i, vFt);
        }

        /* calculate H */
        for (i=0; i<segLen; ++i) {
            vHt = _mm_load_si128(pvHt+i);
            vFt = _mm_load_si128(pvFt+i);
            vH = _mm_max_epi8(
                    vHt,
                    _mm_subs_epi8(vFt, vGapO));
            _mm_store_si128(pvH+i, vH);
            /* check for saturation */
            {
                vSaturationCheck = _mm_or_si128(vSaturationCheck,
                        _mm_or_si128(
                            _mm_cmpeq_epi8(vH, vNegLimit),
                            _mm_cmpeq_epi8(vH, vPosLimit)));
            }
#ifdef ALIGN_EXTRA
            arr_store_si128(score_table, vH, i, segLen, j, s2Len);
#endif
        }
    }

    /* extract last value from the last column */
    {
        __m128i vH = _mm_load_si128(pvH + offset);
        for (k=0; k<position; ++k) {
            vH = _mm_slli_si128(vH, 1);
        }
        score = (int8_t) _mm_extract_epi8 (vH, 15);
    }

    if (_mm_movemask_epi8(vSaturationCheck)) {
        score = INT8_MAX;
    }

    free(pvP);
    free(pvE);
    free(pvHt);
    free(pvFt);
    free(pvH);
    free(boundary);

    return score;
}
