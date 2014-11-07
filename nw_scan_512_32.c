#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include <immintrin.h>
#include <zmmintrin.h>

#ifdef ALIGN_EXTRA
#include "align/align_scan_512_32_debug.h"
#else
#include "align/align_scan_512_32.h"
#endif
#include "blosum/blosum_map.h"


#if 1
#define SHIFT(v) _mm512_mask_permutevar_epi32(vZero, mShifter, vShifter, v)
#else
#define SHIFT(v) insert(_mm512_permutevar_epi32(vShifter, v), 0, 15)
#endif

static inline __m512i insert(__m512i v, int32_t value, int32_t position) {
    union {
        __m512i v;
        int32_t i[16];
    } u;
    u.v = v;
    u.i[position] = value;
    return u.v;
}

static inline int32_t extract(__m512i v, int32_t position) {
    union {
        __m512i v;
        int32_t i[16];
    } u;
    u.v = v;
    return u.i[position];
}

#ifdef ALIGN_EXTRA
static inline void arr_store_si512(
        int32_t *array,
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

#define PARALLEL_PREFIX_OP(vFt, gap, segLen)    \
{                                               \
    union {                                     \
        __m512i v;                              \
        int32_t i[16];                          \
    } u;                                        \
    u.v = vFt;                                  \
    u.i[ 1] = MAX(u.i[ 0]-segLen*gap, u.i[ 1]); \
    u.i[ 2] = MAX(u.i[ 1]-segLen*gap, u.i[ 2]); \
    u.i[ 3] = MAX(u.i[ 2]-segLen*gap, u.i[ 3]); \
    u.i[ 4] = MAX(u.i[ 3]-segLen*gap, u.i[ 4]); \
    u.i[ 5] = MAX(u.i[ 4]-segLen*gap, u.i[ 5]); \
    u.i[ 6] = MAX(u.i[ 5]-segLen*gap, u.i[ 6]); \
    u.i[ 7] = MAX(u.i[ 6]-segLen*gap, u.i[ 7]); \
    u.i[ 8] = MAX(u.i[ 7]-segLen*gap, u.i[ 8]); \
    u.i[ 9] = MAX(u.i[ 8]-segLen*gap, u.i[ 9]); \
    u.i[10] = MAX(u.i[ 9]-segLen*gap, u.i[10]); \
    u.i[11] = MAX(u.i[10]-segLen*gap, u.i[11]); \
    u.i[12] = MAX(u.i[11]-segLen*gap, u.i[12]); \
    u.i[13] = MAX(u.i[12]-segLen*gap, u.i[13]); \
    u.i[14] = MAX(u.i[13]-segLen*gap, u.i[14]); \
    u.i[15] = MAX(u.i[14]-segLen*gap, u.i[15]); \
    vFt = u.v;                                  \
}

#ifdef ALIGN_EXTRA
#define FNAME nw_scan_512_32_debug
#else
#define FNAME nw_scan_512_32
#endif

int32_t FNAME(
        const char * const restrict s1, const int32_t s1Len,
        const char * const restrict s2, const int32_t s2Len,
        const int32_t open, const int32_t gap,
        const int8_t * const restrict matrix
#ifdef ALIGN_EXTRA
        , int32_t * const restrict score_table
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
    __m512i* pvP = (__m512i*)_mm_malloc(n * segLen * sizeof(__m512i), 64);
    __m512i* pvE = (__m512i*)_mm_malloc(segLen * sizeof(__m512i), 64);
    __m512i* pvHt = (__m512i*)_mm_malloc(segLen * sizeof(__m512i), 64);
    __m512i* pvFt = (__m512i*)_mm_malloc(segLen * sizeof(__m512i), 64);
    __m512i* pvH = (__m512i*)_mm_malloc(segLen * sizeof(__m512i), 64);
    int32_t* boundary = (int32_t*)malloc((s2Len+1) * sizeof(int32_t));
    __m512i vGapO = _mm512_set1_epi32(open);
    __m512i vGapE = _mm512_set1_epi32(gap);
    int32_t score = 0;
    __m512i vZero = _mm512_set1_epi32(0);
    __m512i v15 = _mm512_set1_epi32(15);
    __m512i vShifter = {15,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    __mmask16 mShifter = _mm512_cmplt_epi32_mask(vShifter, v15); /* 65534 */

    /* Generate query profile.
     * Rearrange query sequence & calculate the weight of match/mismatch */
    {
        int32_t *t = (int32_t*)pvP;
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
        int32_t *h = (int32_t*)pvH;
        int32_t *e = (int32_t*)pvE;
        for (i=0; i<segLen; ++i) {
            for (segNum=0; segNum<segWidth; ++segNum) {
                *h = -open-gap*(segNum*segLen+i);
                *e = NEG_INF_32;
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
        __m512i vE;
        __m512i vHt;
        __m512i vFt;
        __m512i vH;
        __m512i vHp;
        __m512i *pvW;
        __m512i vW;

#define LOOP_FUSION 1
#if LOOP_FUSION
        /* calculate E */
        /* calculate Ht */
        vHp = _mm512_load_epi32(pvH+(segLen-1));
        vHp = SHIFT(vHp);
        vHp = insert(vHp, boundary[j], 0);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            vH = _mm512_load_epi32(pvH+i);
            vE = _mm512_load_epi32(pvE+i);
            vW = _mm512_load_epi32(pvW+i);
            vE = _mm512_max_epi32(
                    _mm512_sub_epi32(vE, vGapE),
                    _mm512_sub_epi32(vH, vGapO));
            vHt = _mm512_max_epi32(
                    _mm512_add_epi32(vHp, vW),
                    vE);
            _mm512_store_epi32(pvE+i, vE);
            _mm512_store_epi32(pvHt+i, vHt);
            vHp = vH;
        }
#else
        /* calculate E */
        for (i=0; i<segLen; ++i) {
            vH = _mm512_load_epi32(pvH+i);
            vE = _mm512_load_epi32(pvE+i);
            vE = _mm512_max_epi32(
                    _mm512_sub_epi32(vE, vGapE),
                    _mm512_sub_epi32(vH, vGapO));
            _mm512_store_epi32(pvE+i, vE);
        }

        /* calculate Ht */
        vH = _mm512_load_epi32(pvH+(segLen-1));
        vH = SHIFT(vH);
        vH = insert(vH, boundary[j], 0);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            vE = _mm512_load_epi32(pvE+i);
            vW = _mm512_load_epi32(pvW+i);
            vHt = _mm512_max_epi32(
                    _mm512_add_epi32(vH, vW),
                    vE);
            vH = _mm512_load_epi32(pvH+i);
            _mm512_store_epi32(pvHt+i, vHt);
        }
#endif

        /* calculate Ft */
        vHt = _mm512_load_epi32(pvHt+(segLen-1));
        vHt = SHIFT(vHt);
        vHt = insert(vHt, boundary[j+1], 0);
        vFt = _mm512_set1_epi32(NEG_INF_32);
        for (i=0; i<segLen; ++i) {
            vFt = _mm512_max_epi32(
                    _mm512_sub_epi32(vFt, vGapE),
                    vHt);
            vHt = _mm512_load_epi32(pvHt+i);
        }
        PARALLEL_PREFIX_OP(vFt, gap, segLen)
        vHt = _mm512_load_epi32(pvHt+(segLen-1));
        vHt = SHIFT(vHt);
        vHt = insert(vHt, boundary[j+1], 0);
        vFt = SHIFT(vFt);
        vFt = insert(vFt, NEG_INF_32, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm512_max_epi32(
                    _mm512_sub_epi32(vFt, vGapE),
                    vHt);
            vHt = _mm512_load_epi32(pvHt+i);
            _mm512_store_epi32(pvFt+i, vFt);
        }

        /* calculate H */
        for (i=0; i<segLen; ++i) {
            vHt = _mm512_load_epi32(pvHt+i);
            vFt = _mm512_load_epi32(pvFt+i);
            vH = _mm512_max_epi32(
                    vHt,
                    _mm512_sub_epi32(vFt, vGapO));
            _mm512_store_epi32(pvH+i, vH);
#ifdef ALIGN_EXTRA
            arr_store_si512(score_table, vH, i, segLen, j, s2Len);
#endif
        }
    }

    /* extract last value from the last column */
    {
        __m512i vH = _mm512_load_epi32(pvH + offset);
        for (k=0; k<position; ++k) {
            vH = SHIFT(vH);
        }
        score = (int32_t) extract (vH, 15);
    }

    _mm_free(pvP);
    _mm_free(pvE);
    _mm_free(pvHt);
    _mm_free(pvFt);
    _mm_free(pvH);
    free(boundary);

    return score;
}
