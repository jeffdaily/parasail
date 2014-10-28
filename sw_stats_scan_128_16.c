#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include <emmintrin.h>

#ifdef ALIGN_EXTRA
#include "align/align_scan_128_16_debug.h"
#else
#include "align/align_scan_128_16.h"
#endif
#include "blosum/blosum_map.h"


#ifdef ALIGN_EXTRA
static inline void arr_store_si128(
        int *array,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[(0*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 0);
    array[(1*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 1);
    array[(2*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 2);
    array[(3*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 3);
    array[(4*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 4);
    array[(5*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 5);
    array[(6*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 6);
    array[(7*seglen+t)*dlen + d] = (int16_t)_mm_extract_epi16(vH, 7);
}
#endif

#define PARALLEL_PREFIX_OP(vFt, gap, segLen)            \
{                                                       \
    union {                                             \
        __m128i m;                                      \
        int16_t v[8];                                   \
    } tmp;                                              \
    tmp.m = vFt;                                        \
    tmp.v[1] = MAX(tmp.v[0]-segLen*gap, tmp.v[1]);      \
    tmp.v[2] = MAX(tmp.v[1]-segLen*gap, tmp.v[2]);      \
    tmp.v[3] = MAX(tmp.v[2]-segLen*gap, tmp.v[3]);      \
    tmp.v[4] = MAX(tmp.v[3]-segLen*gap, tmp.v[4]);      \
    tmp.v[5] = MAX(tmp.v[4]-segLen*gap, tmp.v[5]);      \
    tmp.v[6] = MAX(tmp.v[5]-segLen*gap, tmp.v[6]);      \
    tmp.v[7] = MAX(tmp.v[6]-segLen*gap, tmp.v[7]);      \
    vFt = tmp.m;                                        \
}

#ifdef ALIGN_EXTRA
#define FNAME sw_stats_scan_128_16_debug
#else
#define FNAME sw_stats_scan_128_16
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
    const int32_t n = 24; /* number of amino acids in table */
    const int32_t segWidth = 8; /* number of values in vector unit */
    int32_t segNum = 0;
    int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    int32_t offset = (s1Len - 1) % segLen;
    int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m128i* pvP = (__m128i*)malloc(n * segLen * sizeof(__m128i));
    __m128i* pvPm= (__m128i*)malloc(n * segLen * sizeof(__m128i));
    __m128i* pvE = (__m128i*)calloc(segLen, sizeof(__m128i));
    __m128i* pvHt = (__m128i*)calloc(segLen, sizeof(__m128i));
    __m128i* pvFt = (__m128i*)calloc(segLen, sizeof(__m128i));
    __m128i* pvMt = (__m128i*)calloc(segLen, sizeof(__m128i));
    __m128i* pvLt = (__m128i*)calloc(segLen, sizeof(__m128i));
    __m128i* pvEx = (__m128i*)calloc(segLen, sizeof(__m128i));
    __m128i* pvH = (__m128i*)calloc(segLen, sizeof(__m128i));
    __m128i* pvM = (__m128i*)calloc(segLen, sizeof(__m128i));
    __m128i* pvL = (__m128i*)calloc(segLen, sizeof(__m128i));
    __m128i vGapO = _mm_set1_epi16(open);
    __m128i vGapE = _mm_set1_epi16(gap);
    __m128i vZero = _mm_setzero_si128();
    __m128i vOne = _mm_set1_epi16(1);
    int16_t score = 0;
    __m128i vMaxH = vZero;
    __m128i vMaxM = vZero;
    __m128i vMaxL = vZero;
    __m128i vQLimit = _mm_set1_epi16(s1Len);
    __m128i vQIndex_reset = _mm_set_epi16(
            7*segLen,
            6*segLen,
            5*segLen,
            4*segLen,
            3*segLen,
            2*segLen,
            1*segLen,
            0*segLen);

    /* Generate query profile and match profile.
     * Rearrange query sequence & calculate the weight of match/mismatch */
    {
        int16_t *t = (int16_t*)pvP;
        int16_t *s = (int16_t*)pvPm;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                int32_t j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    *t++ = matrix[k*n + MAP_BLOSUM_[(unsigned char)s1[j]]];
                    *s++ = (k == MAP_BLOSUM_[(unsigned char)s1[j]]);
                    j += segLen;
                }
            }
        }
    }

    /* initialize H and E */
    {
        int16_t *h = (int16_t*)pvH;
        int16_t *e = (int16_t*)pvE;
        for (i=0; i<segLen; ++i) {
            for (segNum=0; segNum<segWidth; ++segNum) {
                *h = 0;
                *e = NEG_INF_16;
                ++h;
                ++e;
            }
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vE;
        __m128i vHt;
        __m128i vFt;
        __m128i vH;
        __m128i *pvW;
        __m128i vW;
        __m128i *pvC;
        __m128i vC;
        __m128i vM;
        __m128i vMp;
        __m128i vMt;
        __m128i vL;
        __m128i vLp;
        __m128i vLt;
        __m128i vEx;
        __m128i cond_lmt;
        __m128i cond_max;
        __m128i cond_all;
        __m128i vQIndex;

        /* calculate E */
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vE = _mm_max_epi16(
                    _mm_sub_epi16(vE, vGapE),
                    _mm_sub_epi16(vH, vGapO));
            _mm_store_si128(pvE+i, vE);
        }

        /* calculate Ht */
        vH = _mm_slli_si128(_mm_load_si128(pvH+(segLen-1)), 2);
        vMp= _mm_slli_si128(_mm_load_si128(pvM+(segLen-1)), 2);
        vLp= _mm_slli_si128(_mm_load_si128(pvL+(segLen-1)), 2);
        vLp= _mm_add_epi16(vLp, vOne);
        pvW = pvP + MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        pvC = pvPm+ MAP_BLOSUM_[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            /* compute */
            vH = _mm_add_epi16(vH, vW);
            vHt = _mm_max_epi16(vH, vE);
            vHt = _mm_max_epi16(vHt, vZero);
            /* statistics */
            vC = _mm_load_si128(pvC+i);
            vMp = _mm_add_epi16(vMp, vC);
            vEx = _mm_cmpgt_epi16(vE, vH);
            vM = _mm_load_si128(pvM+i);
            vL = _mm_load_si128(pvL+i);
            vL = _mm_add_epi16(vL, vOne);
            vMt = _mm_and_si128(vEx, vM);
            vMt = _mm_or_si128(vMt, _mm_andnot_si128(vEx, vMp));
            vLt = _mm_and_si128(vEx, vL);
            vLt = _mm_or_si128(vLt, _mm_andnot_si128(vEx, vLp));
            cond_max = _mm_cmpeq_epi16(vHt, vZero);
            vEx = _mm_andnot_si128(cond_max, vEx);
            vMt = _mm_andnot_si128(cond_max, vMt);
            vLt = _mm_andnot_si128(cond_max, vLt);
            /* store results */
            _mm_store_si128(pvHt+i, vHt);
            _mm_store_si128(pvEx+i, vEx);
            _mm_store_si128(pvMt+i, vMt);
            _mm_store_si128(pvLt+i, vLt);
            /* prep for next iteration */
            vH = _mm_load_si128(pvH+i);
            vMp = vM;
            vLp = vL;
        }

        /* calculate Ft */
        vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 2);
        vFt = _mm_set1_epi16(NEG_INF_16);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_sub_epi16(vFt, vGapE);
            vFt = _mm_max_epi16(vFt, vHt);
            vHt = _mm_load_si128(pvHt+i);
        }
        PARALLEL_PREFIX_OP(vFt, gap, segLen)
            vHt = _mm_slli_si128(_mm_load_si128(pvHt+(segLen-1)), 2);
        vFt = _mm_slli_si128(vFt, 2);
        vFt = _mm_insert_epi16(vFt, NEG_INF_16, 0);
        for (i=0; i<segLen; ++i) {
            vFt = _mm_sub_epi16(vFt, vGapE);
            vFt = _mm_max_epi16(vFt, vHt);
            vHt = _mm_load_si128(pvHt+i);
            _mm_store_si128(pvFt+i, vFt);
        }

        /* calculate H,M,L */
        vMp = vZero;
        vLp = vOne;
        vC = _mm_cmpeq_epi16(vZero, vZero); /* check if prefix sum is needed */
        vC = _mm_srli_si128(vC, 2); /* zero out last value */
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vHt = _mm_load_si128(pvHt+i);
            vFt = _mm_load_si128(pvFt+i);
            /* compute */
            vFt = _mm_sub_epi16(vFt, vGapO);
            vH = _mm_max_epi16(vHt, vFt);
            /* statistics */
            vEx = _mm_load_si128(pvEx+i);
            vMt = _mm_load_si128(pvMt+i);
            vLt = _mm_load_si128(pvLt+i);
            vEx = _mm_or_si128(
                    _mm_and_si128(vEx, _mm_cmpeq_epi16(vHt, vFt)),
                    _mm_cmplt_epi16(vHt, vFt));
            vM = _mm_and_si128(vEx, vMp);
            vM = _mm_or_si128(vM, _mm_andnot_si128(vEx, vMt));
            vMp = vM;
            vL = _mm_and_si128(vEx, vLp);
            vL = _mm_or_si128(vL, _mm_andnot_si128(vEx, vLt));
            vLp = _mm_add_epi16(vL, vOne);
            vC = _mm_and_si128(vC, vEx);
            /* store results */
            _mm_store_si128(pvH+i, vH);
            _mm_store_si128(pvEx+i, vEx);
#ifdef ALIGN_EXTRA
            arr_store_si128(score_table, vH, i, segLen, j, s2Len);
#endif
        }
#if PREFIX_SUM_CHECK
        /* check if local prefix sum for L is needed */
        if (_mm_movemask_epi8(vC))
#endif
        {
            vLp = _mm_sub_epi16(vLp, vOne);
            {
                union {
                    __m128i m;
                    int16_t v[8];
                } uMp, uLp, uC;
                uC.m = vC;
                uMp.m = vMp;
                uMp.v[1] = uC.v[1] ? uMp.v[0] : uMp.v[1];
                uMp.v[2] = uC.v[2] ? uMp.v[1] : uMp.v[2];
                uMp.v[3] = uC.v[3] ? uMp.v[2] : uMp.v[3];
                uMp.v[4] = uC.v[4] ? uMp.v[3] : uMp.v[4];
                uMp.v[5] = uC.v[5] ? uMp.v[4] : uMp.v[5];
                uMp.v[6] = uC.v[6] ? uMp.v[5] : uMp.v[6];
                uMp.v[7] = uC.v[7] ? uMp.v[6] : uMp.v[7];
                vMp = uMp.m;
                uLp.m = vLp;
                uLp.v[1] = uC.v[1] ? uLp.v[1] + uLp.v[0] : uLp.v[1];
                uLp.v[2] = uC.v[2] ? uLp.v[2] + uLp.v[1] : uLp.v[2];
                uLp.v[3] = uC.v[3] ? uLp.v[3] + uLp.v[2] : uLp.v[3];
                uLp.v[4] = uC.v[4] ? uLp.v[4] + uLp.v[3] : uLp.v[4];
                uLp.v[5] = uC.v[5] ? uLp.v[5] + uLp.v[4] : uLp.v[5];
                uLp.v[6] = uC.v[6] ? uLp.v[6] + uLp.v[5] : uLp.v[6];
                uLp.v[7] = uC.v[7] ? uLp.v[7] + uLp.v[6] : uLp.v[7];
                vLp = uLp.m;
            }
            vLp = _mm_add_epi16(vLp, vOne);
        }
        /* final pass for M,L */
        vQIndex = vQIndex_reset;
        vMp = _mm_slli_si128(vMp, 2);
        vLp = _mm_slli_si128(vLp, 2);
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vH = _mm_load_si128(pvH+i);
            /* statistics */
            vEx = _mm_load_si128(pvEx+i);
            vMt = _mm_load_si128(pvMt+i);
            vLt = _mm_load_si128(pvLt+i);
            vM = _mm_and_si128(vEx, vMp);
            vM = _mm_or_si128(vM, _mm_andnot_si128(vEx, vMt));
            vMp = vM;
            vL = _mm_and_si128(vEx, vLp);
            vL = _mm_or_si128(vL, _mm_andnot_si128(vEx, vLt));
            vLp = _mm_add_epi16(vL, vOne);
            /* store results */
            _mm_store_si128(pvM+i, vM);
            _mm_store_si128(pvL+i, vL);
#ifdef ALIGN_EXTRA
            arr_store_si128(match_table, vM, i, segLen, j, s2Len);
            arr_store_si128(length_table, vL, i, segLen, j, s2Len);
#endif
            /* update max vector seen so far */
            {
                __m128i cond_lmt = _mm_cmplt_epi16(vQIndex, vQLimit);
                __m128i cond_max = _mm_cmpgt_epi16(vH, vMaxH);
                __m128i cond_all = _mm_and_si128(cond_max, cond_lmt);
                vMaxH = _mm_andnot_si128(cond_all, vMaxH);
                vMaxH = _mm_or_si128(vMaxH, _mm_and_si128(cond_all, vH));
                vMaxM = _mm_andnot_si128(cond_all, vMaxM);
                vMaxM = _mm_or_si128(vMaxM, _mm_and_si128(cond_all, vM));
                vMaxL = _mm_andnot_si128(cond_all, vMaxL);
                vMaxL = _mm_or_si128(vMaxL, _mm_and_si128(cond_all, vL));
                vQIndex = _mm_add_epi16(vQIndex, vOne);
            }
        }
    }

    /* max in vec */
    for (j=0; j<segWidth; ++j) {
        int16_t value = (int16_t) _mm_extract_epi16(vMaxH, 7);
        if (value > score) {
            *matches = (int16_t) _mm_extract_epi16(vMaxM, 7);
            *length = (int16_t) _mm_extract_epi16(vMaxL, 7);
            score = value;
        }
        vMaxH = _mm_slli_si128(vMaxH, 2);
        vMaxM = _mm_slli_si128(vMaxM, 2);
        vMaxL = _mm_slli_si128(vMaxL, 2);
    }

    free(pvP);
    free(pvPm);
    free(pvE);
    free(pvHt);
    free(pvFt);
    free(pvMt);
    free(pvLt);
    free(pvEx);
    free(pvH);
    free(pvM);
    free(pvL);

    return score;
}
