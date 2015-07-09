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
#define FNAME parasail_sg_stats_table_scan_knc_512_32
#define PNAME parasail_sg_stats_table_scan_profile_knc_512_32
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_stats_rowcol_scan_knc_512_32
#define PNAME parasail_sg_stats_rowcol_scan_profile_knc_512_32
#else
#define FNAME parasail_sg_stats_scan_knc_512_32
#define PNAME parasail_sg_stats_scan_profile_knc_512_32
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_knc_512_32(s1, s1Len, matrix);
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
    const int32_t segWidth = 16;
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m512i* const restrict pvP  = (__m512i*)profile->profile;
    __m512i* const restrict pvPm = (__m512i*)profile->profile_m;
    __m512i* const restrict pvPs = (__m512i*)profile->profile_s;
    __m512i* const restrict pvE  = parasail_memalign___m512i(64, segLen);
    __m512i* const restrict pvHt = parasail_memalign___m512i(64, segLen);
    __m512i* const restrict pvFt = parasail_memalign___m512i(64, segLen);
    __m512i* const restrict pvMt = parasail_memalign___m512i(64, segLen);
    __m512i* const restrict pvSt = parasail_memalign___m512i(64, segLen);
    __m512i* const restrict pvLt = parasail_memalign___m512i(64, segLen);
    int*     const restrict pvEx = parasail_memalign_int(64, segLen);
    __m512i* const restrict pvH  = parasail_memalign___m512i(64, segLen);
    __m512i* const restrict pvM  = parasail_memalign___m512i(64, segLen);
    __m512i* const restrict pvS  = parasail_memalign___m512i(64, segLen);
    __m512i* const restrict pvL  = parasail_memalign___m512i(64, segLen);
    __m512i vGapO = _mm512_set1_epi32(open);
    __m512i vGapE = _mm512_set1_epi32(gap);
    __m512i vZero = _mm512_set1_epi32(0);
    __m512i vOne = _mm512_set1_epi32(1);
    __m512i vNegInf = _mm512_set1_epi32(NEG_INF);
    int score = NEG_INF;
    int matches = 0;
    int similar = 0;
    int length = 0;
    __m512i vMaxH = vNegInf;
    __m512i vMaxM = vZero;
    __m512i vMaxS = vZero;
    __m512i vMaxL = vZero;
    __m512i permute_idx = _mm512_set_16to16_pi(14,13,12,11,10,9,8,7,6,5,4,3,2,1,0,15);
    __m512i vSegLimit = _mm512_set1_epi32(segWidth-1);
    __mmask16 permute_mask = _mm512_cmplt_epi32_mask(permute_idx, vSegLimit);
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif

    parasail_memset___m512i(pvM, vZero, segLen);
    parasail_memset___m512i(pvS, vZero, segLen);
    parasail_memset___m512i(pvL, vZero, segLen);

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
        __m512i *pvW;
        __m512i vW;
        __m512i *pvC;
        __m512i *pvD;
        __m512i vC;
        __m512i vD;
        __m512i vM;
        __m512i vMp;
        __m512i vMt;
        __m512i vS;
        __m512i vSp;
        __m512i vSt;
        __m512i vL;
        __m512i vLp;
        __m512i vLt;
        __mmask16 vEx;
        __mmask16 vCx;

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
        vH = _mm512_mask_permutevar_epi32(
                vZero, permute_mask, permute_idx, vH);
        vMp= _mm512_load_epi32(pvM+(segLen-1));
        vMp= _mm512_mask_permutevar_epi32(
                vZero, permute_mask, permute_idx, vMp);
        vSp= _mm512_load_epi32(pvS+(segLen-1));
        vSp= _mm512_mask_permutevar_epi32(
                vZero, permute_mask, permute_idx, vSp);
        vLp= _mm512_load_epi32(pvL+(segLen-1));
        vLp= _mm512_mask_permutevar_epi32(
                vZero, permute_mask, permute_idx, vLp);
        vLp= _mm512_add_epi32(vLp, vOne);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvC = pvPm+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        pvD = pvPs+ matrix->mapper[(unsigned char)s2[j]]*segLen;
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vE = _mm512_load_epi32(pvE+i);
            vW = _mm512_load_epi32(pvW+i);
            /* compute */
            vH = _mm512_add_epi32(vH, vW);
            vHt = _mm512_max_epi32(vH, vE);
            /* statistics */
            vC = _mm512_load_epi32(pvC+i);
            vD = _mm512_load_epi32(pvD+i);
            vMp = _mm512_add_epi32(vMp, vC);
            vSp = _mm512_add_epi32(vSp, vD);
            vEx = _mm512_cmpgt_epi32_mask(vE, vH);
            vM = _mm512_load_epi32(pvM+i);
            vS = _mm512_load_epi32(pvS+i);
            vL = _mm512_load_epi32(pvL+i);
            vL = _mm512_add_epi32(vL, vOne);
            vMt = _mm512_mask_blend_epi32(vEx, vMp, vM);
            vSt = _mm512_mask_blend_epi32(vEx, vSp, vS);
            vLt = _mm512_mask_blend_epi32(vEx, vLp, vL);
            /* store results */
            _mm512_store_epi32(pvHt+i, vHt);
            pvEx[i] = _mm512_mask2int(vEx);
            _mm512_store_epi32(pvMt+i, vMt);
            _mm512_store_epi32(pvSt+i, vSt);
            _mm512_store_epi32(pvLt+i, vLt);
            /* prep for next iteration */
            vH = _mm512_load_epi32(pvH+i);
            vMp = vM;
            vSp = vS;
            vLp = vL;
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

        /* calculate H,M,L */
        vMp = vZero;
        vSp = vZero;
        vLp = vOne;
        vCx = permute_mask;
        for (i=0; i<segLen; ++i) {
            /* load values we need */
            vHt = _mm512_load_epi32(pvHt+i);
            vFt = _mm512_load_epi32(pvFt+i);
            /* compute */
            vFt = _mm512_sub_epi32(vFt, vGapO);
            vH = _mm512_max_epi32(vHt, vFt);
            /* statistics */
            vEx = _mm512_int2mask(pvEx[i]);
            vMt = _mm512_load_epi32(pvMt+i);
            vSt = _mm512_load_epi32(pvSt+i);
            vLt = _mm512_load_epi32(pvLt+i);
            vEx = _mm512_kor(
                    _mm512_kand(vEx, _mm512_cmpeq_epi32_mask(vHt, vFt)),
                    _mm512_cmplt_epi32_mask(vHt, vFt));
            vM = _mm512_mask_blend_epi32(vEx, vMt, vMp);
            vS = _mm512_mask_blend_epi32(vEx, vSt, vSp);
            vL = _mm512_mask_blend_epi32(vEx, vLt, vLp);
            vMp = vM;
            vSp = vS;
            vLp = _mm512_add_epi32(vL, vOne);
            vCx = _mm512_kand(vCx, vEx);
            /* store results */
            _mm512_store_epi32(pvH+i, vH);
            pvEx[i] = _mm512_mask2int(vEx);
#ifdef PARASAIL_TABLE
            arr_store_si512(result->score_table, vH, i, segLen, j, s2Len);
#endif
        }
        {
            vLp = _mm512_sub_epi32(vLp, vOne);
            {
                __m512i_32_t uMp, uSp, uLp;
                int uC = _mm512_mask2int(vCx);
                uMp.m = vMp;
                uMp.v[ 1] = (uC&((int)1<< 1)) ? uMp.v[ 0] : uMp.v[ 1];
                uMp.v[ 2] = (uC&((int)1<< 2)) ? uMp.v[ 1] : uMp.v[ 2];
                uMp.v[ 3] = (uC&((int)1<< 3)) ? uMp.v[ 2] : uMp.v[ 3];
                uMp.v[ 4] = (uC&((int)1<< 4)) ? uMp.v[ 3] : uMp.v[ 4];
                uMp.v[ 5] = (uC&((int)1<< 5)) ? uMp.v[ 4] : uMp.v[ 5];
                uMp.v[ 6] = (uC&((int)1<< 6)) ? uMp.v[ 5] : uMp.v[ 6];
                uMp.v[ 7] = (uC&((int)1<< 7)) ? uMp.v[ 6] : uMp.v[ 7];
                uMp.v[ 8] = (uC&((int)1<< 8)) ? uMp.v[ 7] : uMp.v[ 8];
                uMp.v[ 9] = (uC&((int)1<< 9)) ? uMp.v[ 8] : uMp.v[ 9];
                uMp.v[10] = (uC&((int)1<<10)) ? uMp.v[ 9] : uMp.v[10];
                uMp.v[11] = (uC&((int)1<<11)) ? uMp.v[10] : uMp.v[11];
                uMp.v[12] = (uC&((int)1<<12)) ? uMp.v[11] : uMp.v[12];
                uMp.v[13] = (uC&((int)1<<13)) ? uMp.v[12] : uMp.v[13];
                uMp.v[14] = (uC&((int)1<<14)) ? uMp.v[13] : uMp.v[14];
                uMp.v[15] = (uC&((int)1<<15)) ? uMp.v[14] : uMp.v[15];
                vMp = uMp.m;
                uSp.m = vSp;
                uSp.v[ 1] = (uC&((int)1<< 1)) ? uSp.v[ 0] : uSp.v[ 1];
                uSp.v[ 2] = (uC&((int)1<< 2)) ? uSp.v[ 1] : uSp.v[ 2];
                uSp.v[ 3] = (uC&((int)1<< 3)) ? uSp.v[ 2] : uSp.v[ 3];
                uSp.v[ 4] = (uC&((int)1<< 4)) ? uSp.v[ 3] : uSp.v[ 4];
                uSp.v[ 5] = (uC&((int)1<< 5)) ? uSp.v[ 4] : uSp.v[ 5];
                uSp.v[ 6] = (uC&((int)1<< 6)) ? uSp.v[ 5] : uSp.v[ 6];
                uSp.v[ 7] = (uC&((int)1<< 7)) ? uSp.v[ 6] : uSp.v[ 7];
                uSp.v[ 8] = (uC&((int)1<< 8)) ? uSp.v[ 7] : uSp.v[ 8];
                uSp.v[ 9] = (uC&((int)1<< 9)) ? uSp.v[ 8] : uSp.v[ 9];
                uSp.v[10] = (uC&((int)1<<10)) ? uSp.v[ 9] : uSp.v[10];
                uSp.v[11] = (uC&((int)1<<11)) ? uSp.v[10] : uSp.v[11];
                uSp.v[12] = (uC&((int)1<<12)) ? uSp.v[11] : uSp.v[12];
                uSp.v[13] = (uC&((int)1<<13)) ? uSp.v[12] : uSp.v[13];
                uSp.v[14] = (uC&((int)1<<14)) ? uSp.v[13] : uSp.v[14];
                uSp.v[15] = (uC&((int)1<<15)) ? uSp.v[14] : uSp.v[15];
                vSp = uSp.m;
                uLp.m = vLp;
                uLp.v[ 1] = (uC&((int)1<< 1))?uLp.v[ 1] + uLp.v[ 0] : uLp.v[ 1];
                uLp.v[ 2] = (uC&((int)1<< 2))?uLp.v[ 2] + uLp.v[ 1] : uLp.v[ 2];
                uLp.v[ 3] = (uC&((int)1<< 3))?uLp.v[ 3] + uLp.v[ 2] : uLp.v[ 3];
                uLp.v[ 4] = (uC&((int)1<< 4))?uLp.v[ 4] + uLp.v[ 3] : uLp.v[ 4];
                uLp.v[ 5] = (uC&((int)1<< 5))?uLp.v[ 5] + uLp.v[ 4] : uLp.v[ 5];
                uLp.v[ 6] = (uC&((int)1<< 6))?uLp.v[ 6] + uLp.v[ 5] : uLp.v[ 6];
                uLp.v[ 7] = (uC&((int)1<< 7))?uLp.v[ 7] + uLp.v[ 6] : uLp.v[ 7];
                uLp.v[ 8] = (uC&((int)1<< 8))?uLp.v[ 8] + uLp.v[ 7] : uLp.v[ 8];
                uLp.v[ 9] = (uC&((int)1<< 9))?uLp.v[ 9] + uLp.v[ 8] : uLp.v[ 9];
                uLp.v[10] = (uC&((int)1<<10))?uLp.v[10] + uLp.v[ 9] : uLp.v[10];
                uLp.v[11] = (uC&((int)1<<11))?uLp.v[11] + uLp.v[10] : uLp.v[11];
                uLp.v[12] = (uC&((int)1<<12))?uLp.v[12] + uLp.v[11] : uLp.v[12];
                uLp.v[13] = (uC&((int)1<<13))?uLp.v[13] + uLp.v[12] : uLp.v[13];
                uLp.v[14] = (uC&((int)1<<14))?uLp.v[14] + uLp.v[13] : uLp.v[14];
                uLp.v[15] = (uC&((int)1<<15))?uLp.v[15] + uLp.v[14] : uLp.v[15];
                vLp = uLp.m;
            }
            vLp = _mm512_add_epi32(vLp, vOne);
        }
        /* final pass for M,L */
        vMp = _mm512_mask_permutevar_epi32(
                vZero, permute_mask, permute_idx, vMp);
        vSp = _mm512_mask_permutevar_epi32(
                vZero, permute_mask, permute_idx, vSp);
        vLp = _mm512_mask_permutevar_epi32(
                vZero, permute_mask, permute_idx, vLp);
        for (i=0; i<segLen; ++i) {
            /* statistics */
            vEx = _mm512_int2mask(pvEx[i]);
            vMt = _mm512_load_epi32(pvMt+i);
            vSt = _mm512_load_epi32(pvSt+i);
            vLt = _mm512_load_epi32(pvLt+i);
            vM = _mm512_mask_blend_epi32(vEx, vMt, vMp);
            vS = _mm512_mask_blend_epi32(vEx, vSt, vSp);
            vL = _mm512_mask_blend_epi32(vEx, vLt, vLp);
            vMp = vM;
            vSp = vS;
            vLp = _mm512_add_epi32(vL, vOne);
            /* store results */
            _mm512_store_epi32(pvM+i, vM);
            _mm512_store_epi32(pvS+i, vS);
            _mm512_store_epi32(pvL+i, vL);
#ifdef PARASAIL_TABLE
            arr_store_si512(result->matches_table, vM, i, segLen, j, s2Len);
            arr_store_si512(result->similar_table, vS, i, segLen, j, s2Len);
            arr_store_si512(result->length_table, vL, i, segLen, j, s2Len);
#endif
        }

        /* extract vector containing last value from column */
        {
            __mmask16 cond_max;
            vH = _mm512_load_epi32(pvH + offset);
            vM = _mm512_load_epi32(pvM + offset);
            vS = _mm512_load_epi32(pvS + offset);
            vL = _mm512_load_epi32(pvL + offset);
            cond_max = _mm512_cmpgt_epi32_mask(vH, vMaxH);
            vMaxH = _mm512_mask_blend_epi32(cond_max, vMaxH, vH);
            vMaxM = _mm512_mask_blend_epi32(cond_max, vMaxM, vM);
            vMaxS = _mm512_mask_blend_epi32(cond_max, vMaxS, vS);
            vMaxL = _mm512_mask_blend_epi32(cond_max, vMaxL, vL);
#ifdef PARASAIL_ROWCOL
            for (k=0; k<position; ++k) {
                vH = _mm512_permutevar_epi32(permute_idx, vH);
                vM = _mm512_permutevar_epi32(permute_idx, vM);
                vS = _mm512_permutevar_epi32(permute_idx, vS);
                vL = _mm512_permutevar_epi32(permute_idx, vL);
            }
            result->score_row[j] = (int32_t) extract (vH, 15);
            result->matches_row[j] = (int32_t) extract (vM, 15);
            result->similar_row[j] = (int32_t) extract (vS, 15);
            result->length_row[j] = (int32_t) extract (vL, 15);
#endif
        }
    }

    /* extract last value from column */
    {
        int32_t value;
        for (k=0; k<position; ++k) {
            vMaxH = _mm512_permutevar_epi32(permute_idx, vMaxH);
            vMaxM = _mm512_permutevar_epi32(permute_idx, vMaxM);
            vMaxS = _mm512_permutevar_epi32(permute_idx, vMaxS);
            vMaxL = _mm512_permutevar_epi32(permute_idx, vMaxL);
        }
        value = (int32_t) extract(vMaxH, 15);
        if (value > score) {
            score = value;
            matches = (int32_t)extract(vMaxM, 15);
            similar = (int32_t)extract(vMaxS, 15);
            length = (int32_t)extract(vMaxL, 15);
        }
    }

    /* max of last column */
    {
        vMaxH = vNegInf;
        vMaxM = vZero;
        vMaxS = vZero;
        vMaxL = vZero;

        for (i=0; i<segLen; ++i) {
            /* load the last stored values */
            __m512i vH = _mm512_load_epi32(pvH + i);
            __m512i vM = _mm512_load_epi32(pvM + i);
            __m512i vS = _mm512_load_epi32(pvS + i);
            __m512i vL = _mm512_load_epi32(pvL + i);
            __mmask16 cond_max = _mm512_cmpgt_epi32_mask(vH, vMaxH);
            vMaxH = _mm512_mask_blend_epi32(cond_max, vMaxH, vH);
            vMaxM = _mm512_mask_blend_epi32(cond_max, vMaxM, vM);
            vMaxS = _mm512_mask_blend_epi32(cond_max, vMaxS, vS);
            vMaxL = _mm512_mask_blend_epi32(cond_max, vMaxL, vL);
#ifdef PARASAIL_ROWCOL
            arr_store_col(result->score_col, vH, i, segLen);
            arr_store_col(result->matches_col, vM, i, segLen);
            arr_store_col(result->similar_col, vS, i, segLen);
            arr_store_col(result->length_col, vL, i, segLen);
#endif
        }

        /* max in vec */
        for (j=0; j<segWidth; ++j) {
            int32_t value = (int32_t) extract(vMaxH, 15);
            if (value > score) {
                score = value;
                matches = (int32_t)extract(vMaxM, 15);
                similar = (int32_t)extract(vMaxS, 15);
                length = (int32_t)extract(vMaxL, 15);
            }
            vMaxH = _mm512_permutevar_epi32(permute_idx, vMaxH);
            vMaxM = _mm512_permutevar_epi32(permute_idx, vMaxM);
            vMaxS = _mm512_permutevar_epi32(permute_idx, vMaxS);
            vMaxL = _mm512_permutevar_epi32(permute_idx, vMaxL);
        }
    }

    result->score = score;
    result->matches = matches;
    result->similar = similar;
    result->length = length;

    parasail_free(pvL);
    parasail_free(pvS);
    parasail_free(pvM);
    parasail_free(pvH);
    parasail_free(pvEx);
    parasail_free(pvLt);
    parasail_free(pvSt);
    parasail_free(pvMt);
    parasail_free(pvFt);
    parasail_free(pvHt);
    parasail_free(pvE);

    return result;
}

