/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include <immintrin.h>

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_knc.h"

#define NEG_INF (INT32_MIN/(int32_t)(2))

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

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m512i vH,
        int32_t t,
        int32_t seglen)
{
    col[ 0*seglen+t] = (int32_t)extract(vH,  0);
    col[ 1*seglen+t] = (int32_t)extract(vH,  1);
    col[ 2*seglen+t] = (int32_t)extract(vH,  2);
    col[ 3*seglen+t] = (int32_t)extract(vH,  3);
    col[ 4*seglen+t] = (int32_t)extract(vH,  4);
    col[ 5*seglen+t] = (int32_t)extract(vH,  5);
    col[ 6*seglen+t] = (int32_t)extract(vH,  6);
    col[ 7*seglen+t] = (int32_t)extract(vH,  7);
    col[ 8*seglen+t] = (int32_t)extract(vH,  8);
    col[ 9*seglen+t] = (int32_t)extract(vH,  9);
    col[10*seglen+t] = (int32_t)extract(vH, 10);
    col[11*seglen+t] = (int32_t)extract(vH, 11);
    col[12*seglen+t] = (int32_t)extract(vH, 12);
    col[13*seglen+t] = (int32_t)extract(vH, 13);
    col[14*seglen+t] = (int32_t)extract(vH, 14);
    col[15*seglen+t] = (int32_t)extract(vH, 15);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sg_table_scan_knc_512_32
#define PNAME parasail_sg_table_scan_profile_knc_512_32
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_rowcol_scan_knc_512_32
#define PNAME parasail_sg_rowcol_scan_profile_knc_512_32
#else
#define FNAME parasail_sg_scan_knc_512_32
#define PNAME parasail_sg_scan_profile_knc_512_32
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_knc_512_32(s1, s1Len, matrix);
    parasail_result_t *result = PNAME(profile, s2, s2Len, open, gap);
    parasail_profile_free(profile);
    return result;
}

parasail_result_t* PNAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m512i* const restrict pvP  = (__m512i*)profile->profile;
    __m512i* const restrict pvE  = parasail_memalign___m512i(64, segLen);
    __m512i* const restrict pvHt = parasail_memalign___m512i(64, segLen);
    __m512i* const restrict pvH  = parasail_memalign___m512i(64, segLen);
    __m512i vGapO = _mm512_set1_epi32(open);
    __m512i vGapE = _mm512_set1_epi32(gap);
    __m512i vNegInf = _mm512_set1_epi32(NEG_INF);
    int32_t score = NEG_INF;
    const int32_t segLenXgap = -segLen*gap;
    __m512i vSegLenXgap1 = _mm512_set1_epi32((segLen-1)*gap);
    __m512i vSegLenXgap = _mm512_set_16to16_pi(
            NEG_INF, segLenXgap, segLenXgap, segLenXgap,
            segLenXgap, segLenXgap, segLenXgap, segLenXgap,
            segLenXgap, segLenXgap, segLenXgap, segLenXgap,
            segLenXgap, segLenXgap, segLenXgap, segLenXgap);
    __m512i vSegLimit = _mm512_set1_epi32(segWidth-1);
    __m512i permute_idx = _mm512_set_16to16_pi(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15);
    __mmask16 permute_mask = _mm512_cmplt_epi32_mask(permute_idx, vSegLimit);
    __m512i vZero = _mm512_set1_epi32(0);
    __m512i vMaxH = vNegInf;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m512i_32_t h;
            __m512i_32_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = 0;
                e.v[segNum] = NEG_INF;
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
        /* calculate Ft */
        vHp = _mm512_load_epi32(pvH+(segLen-1));
        vHp = _mm512_mask_permutevar_epi32(
                vZero, permute_mask, permute_idx, vHp);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = vNegInf;
        vFt = vNegInf;
        for (i=0; i<segLen; ++i) {
            vH = _mm512_load_epi32(pvH+i);
            vE = _mm512_load_epi32(pvE+i);
            vW = _mm512_load_epi32(pvW+i);
            vE = _mm512_max_epi32(
                    _mm512_sub_epi32(vE, vGapE),
                    _mm512_sub_epi32(vH, vGapO));
            vFt = _mm512_sub_epi32(vFt, vGapE);
            vFt = _mm512_max_epi32(vFt, vHt);
            vHt = _mm512_max_epi32(
                    _mm512_add_epi32(vHp, vW),
                    vE);
            _mm512_store_epi32(pvE+i, vE);
            _mm512_store_epi32(pvHt+i, vHt);
            vHp = vH;
        }

        /* adjust Ft before local prefix scan */
        vHt = _mm512_mask_permutevar_epi32(
                vZero, permute_mask, permute_idx, vHt);
        vFt = _mm512_max_epi32(vFt,
                _mm512_sub_epi32(vHt, vSegLenXgap1));
        /* local prefix scan */
        for (i=0; i<segWidth-1; ++i) {
            __m512i vFtt = _mm512_mask_permutevar_epi32(
                    vNegInf, permute_mask, permute_idx, vFt);
            vFtt = _mm512_add_epi32(vFtt, vSegLenXgap);
            vFt = _mm512_max_epi32(vFt, vFtt);
        }
        vFt = _mm512_mask_permutevar_epi32(
                vNegInf, permute_mask, permute_idx, vFt);

        /* second Ft pass */
        /* calculate vH */
        for (i=0; i<segLen; ++i) {
            vFt = _mm512_sub_epi32(vFt, vGapE);
            vFt = _mm512_max_epi32(vFt, vHt);
            vHt = _mm512_load_epi32(pvHt+i);
            vH = _mm512_max_epi32(vHt, _mm512_sub_epi32(vFt, vGapO));
            _mm512_store_epi32(pvH+i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si512(result->score_table, vH, i, segLen, j, s2Len);
#endif
        }

        /* extract vector containing last value from column */
        {
            vH = _mm512_load_epi32(pvH + offset);
            vMaxH = _mm512_max_epi32(vH, vMaxH);
#ifdef PARASAIL_ROWCOL
            for (k=0; k<position; ++k) {
                vH = _mm512_permutevar_epi32(permute_idx, vH);
            }
            result->score_row[j] = (int32_t) extract (vH, 15);
#endif
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
        vMaxH = vNegInf;

        for (i=0; i<segLen; ++i) {
            __m512i vH = _mm512_load_epi32(pvH + i);
            vMaxH = _mm512_max_epi32(vH, vMaxH);
#ifdef PARASAIL_ROWCOL
            arr_store_col(result->score_col, vH, i, segLen);
#endif
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
    parasail_free(pvHt);
    parasail_free(pvE);

    return result;
}

