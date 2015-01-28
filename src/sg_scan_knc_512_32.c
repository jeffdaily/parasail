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
#include "parasail_internal_knc.h"
#include "blosum/blosum_map.h"

#define NEG_INF_32 (INT32_MIN/(int32_t)(2))
#define MAX(a,b) ((a)>(b)?(a):(b))

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


#ifdef PARASAIL_TABLE
#define FNAME sg_table_scan_knc_512_32
#else
#define FNAME sg_scan_knc_512_32
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    const int32_t n = 24; /* number of amino acids in table */
    const int32_t segWidth = 16; /* number of values in vector unit */
    int32_t segNum = 0;
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m512i* const restrict pvP  = parasail_memalign_m512i(64, n * segLen);
    __m512i* const restrict pvE  = parasail_memalign_m512i(64, segLen);
    __m512i* const restrict pvHt = parasail_memalign_m512i(64, segLen);
    __m512i* const restrict pvFt = parasail_memalign_m512i(64, segLen);
    __m512i* const restrict pvH  = parasail_memalign_m512i(64, segLen);
    __m512i vGapO = _mm512_set1_epi32(open);
    __m512i vGapE = _mm512_set1_epi32(gap);
    __m512i vNegInf = _mm512_set1_epi32(NEG_INF_32);
    __m512i vZero = _mm512_set1_epi32(0);
    __m512i vSegLimit = _mm512_set1_epi32(segWidth-1);
    int32_t score = NEG_INF_32;
    __m512i vMaxH = vNegInf;
    __m512i permute_idx = _mm512_set_16to16_pi(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15);
    __mmask16 permute_mask = _mm512_cmplt_epi32_mask(permute_idx, vSegLimit);
    __m512i segLenXgap_reset = _mm512_set_16to16_pi(
            NEG_INF_32, NEG_INF_32, NEG_INF_32, NEG_INF_32,
            NEG_INF_32, NEG_INF_32, NEG_INF_32, NEG_INF_32,
            NEG_INF_32, NEG_INF_32, NEG_INF_32, NEG_INF_32,
            NEG_INF_32, NEG_INF_32, NEG_INF_32, -segLen*gap);
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif

    /* Generate query profile.
     * Rearrange query sequence & calculate the weight of match/mismatch.
     * Don't alias. */
    {
        int32_t index = 0;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                __m512i_32_t t;
                j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = j >= s1Len ? 0 : matrix[k][MAP_BLOSUM_[(unsigned char)s1[j]]];
                    j += segLen;
                }
                _mm512_store_epi32(&pvP[index], t.m);
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
            _mm512_store_epi32(&pvH[index], h.m);
            _mm512_store_epi32(&pvE[index], e.m);
            ++index;
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

        /* calculate E */
        /* calculate Ht */
        vHp = _mm512_load_epi32(pvH+(segLen-1));
        vHp = _mm512_mask_permutevar_epi32(
                vZero, permute_mask, permute_idx, vHp);
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

        /* calculate Ft */
        vHt = _mm512_load_epi32(pvHt+(segLen-1));
        vHt = _mm512_mask_permutevar_epi32(
                vZero, permute_mask, permute_idx, vHt);
        vFt = vNegInf;
        for (i=0; i<segLen; ++i) {
            vFt = _mm512_sub_epi32(vFt, vGapE);
            vFt = _mm512_max_epi32(vFt, vHt);
            vHt = _mm512_load_epi32(pvHt+i);
        }
#if 1
        {
            __m512i_32_t tmp;
            tmp.m = vFt;
            tmp.v[ 1] = MAX(tmp.v[ 0]-segLen*gap, tmp.v[ 1]);
            tmp.v[ 2] = MAX(tmp.v[ 1]-segLen*gap, tmp.v[ 2]);
            tmp.v[ 3] = MAX(tmp.v[ 2]-segLen*gap, tmp.v[ 3]);
            tmp.v[ 4] = MAX(tmp.v[ 3]-segLen*gap, tmp.v[ 4]);
            tmp.v[ 5] = MAX(tmp.v[ 4]-segLen*gap, tmp.v[ 5]);
            tmp.v[ 6] = MAX(tmp.v[ 5]-segLen*gap, tmp.v[ 6]);
            tmp.v[ 7] = MAX(tmp.v[ 6]-segLen*gap, tmp.v[ 7]);
            tmp.v[ 8] = MAX(tmp.v[ 7]-segLen*gap, tmp.v[ 8]);
            tmp.v[ 9] = MAX(tmp.v[ 8]-segLen*gap, tmp.v[ 9]);
            tmp.v[10] = MAX(tmp.v[ 9]-segLen*gap, tmp.v[10]);
            tmp.v[11] = MAX(tmp.v[10]-segLen*gap, tmp.v[11]);
            tmp.v[12] = MAX(tmp.v[11]-segLen*gap, tmp.v[12]);
            tmp.v[13] = MAX(tmp.v[12]-segLen*gap, tmp.v[13]);
            tmp.v[14] = MAX(tmp.v[13]-segLen*gap, tmp.v[14]);
            tmp.v[15] = MAX(tmp.v[14]-segLen*gap, tmp.v[15]);
            vFt = tmp.m;
        }
#else
        {
            __m512i vFt_save = vFt;
            __m512i segLenXgap = segLenXgap_reset;
            for (i=0; i<segWidth-1; ++i) {
                __m512i vFtt = _mm512_mask_permutevar_epi32(
                        vNegInf, permute_mask, permute_idx, vFt);
                segLenXgap = _mm512_permutevar_epi32(permute_idx, segLenXgap);
                vFtt = _mm512_add_epi32(vFtt, segLenXgap);
                vFt = _mm512_max_epi32(vFt, vFtt);
            }
            //vFt = _mm512_blendv_epi8(vFt_save, vFt, insert_mask);
        }
#endif
        vHt = _mm512_load_epi32(pvHt+(segLen-1));
        vHt = _mm512_mask_permutevar_epi32(
                vZero, permute_mask, permute_idx, vHt);
        vFt = _mm512_mask_permutevar_epi32(
                vNegInf, permute_mask, permute_idx, vFt);
        for (i=0; i<segLen; ++i) {
            vFt = _mm512_sub_epi32(vFt, vGapE);
            vFt = _mm512_max_epi32(vFt, vHt);
            vHt = _mm512_load_epi32(pvHt+i);
            _mm512_store_epi32(pvFt+i, vFt);
        }

        /* calculate H */
        for (i=0; i<segLen; ++i) {
            vHt = _mm512_load_epi32(pvHt+i);
            vFt = _mm512_load_epi32(pvFt+i);
            vFt = _mm512_sub_epi32(vFt, vGapO);
            vH = _mm512_max_epi32(vHt, vFt);
            _mm512_store_epi32(pvH+i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si512(result->score_table, vH, i, segLen, j, s2Len);
#endif
        }

        /* extract vector containing last value from column */
        {
            vH = _mm512_load_epi32(pvH + offset);
            vMaxH = _mm512_max_epi32(vH, vMaxH);
        }
    }

    /* max last value from all columns */
    {
        int32_t value;
        for (k=0; k<position; ++k) {
            vMaxH = _mm512_permutevar_epi32(permute_idx, vMaxH);
        }
        value = (int32_t) extract(vMaxH, 15);
        if (value > score) {
            score = value;
        }
    }

    /* max of last column */
    {
        __m512i vOne = _mm512_set1_epi32(1);
        __m512i vQLimit = _mm512_set1_epi32(s1Len);
        __m512i vQIndex = _mm512_set_16to16_pi(
                segLen*15,
                segLen*14,
                segLen*13,
                segLen*12,
                segLen*11,
                segLen*10,
                segLen* 9,
                segLen* 8,
                segLen* 7,
                segLen* 6,
                segLen* 5,
                segLen* 4,
                segLen* 3,
                segLen* 2,
                segLen* 1,
                segLen* 0);
        vMaxH = vNegInf;

        for (i=0; i<segLen; ++i) {
            __m512i vH = _mm512_load_epi32(pvH + i);
            __mmask16 cond_lmt = _mm512_cmplt_epi32_mask(vQIndex, vQLimit);
            __mmask16 cond_max = _mm512_cmpgt_epi32_mask(vH, vMaxH);
            __mmask16 cond_all = _mm512_kand(cond_max, cond_lmt);
            vMaxH = _mm512_mask_blend_epi32(cond_all, vMaxH, vH);
            vQIndex = _mm512_add_epi32(vQIndex, vOne);
        }

#if 1
        /* max in vec */
        for (j=0; j<segWidth; ++j) {
            int32_t value = (int32_t) extract(vMaxH, 15);
            if (value > score) {
                score = value;
            }
            vMaxH = _mm512_permutevar_epi32(permute_idx, vMaxH);
        }
#else
        {
            int32_t value = (int32_t) _mm512_reduce_max_epi32(vMaxH);
            if (value > score) {
                score = value;
            }
        }
#endif
    }

    result->score = score;

    parasail_free(pvH);
    parasail_free(pvFt);
    parasail_free(pvHt);
    parasail_free(pvE);
    parasail_free(pvP);

    return result;
}

