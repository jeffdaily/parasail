/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>



#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_altivec.h"



static inline void arr_store(
        vec128i *array,
        vec128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    _mm_store_si128(array + (1LL*d*seglen+t), vH);
}

static inline vec128i arr_load(
        vec128i *array,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    return _mm_load_si128(array + (1LL*d*seglen+t));
}

#define FNAME parasail_nw_trace_scan_altivec_128_64
#define PNAME parasail_nw_trace_scan_profile_altivec_128_64

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    /* declare local variables */
    parasail_profile_t *profile = NULL;
    parasail_result_t *result = NULL;

    /* validate inputs */
    PARASAIL_CHECK_NULL(s2);
    PARASAIL_CHECK_GT0(s2Len);
    PARASAIL_CHECK_GE0(open);
    PARASAIL_CHECK_GE0(gap);
    PARASAIL_CHECK_NULL(matrix);
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        PARASAIL_CHECK_NULL(s1);
        PARASAIL_CHECK_GT0(s1Len);
    }

    /* initialize local variables */
    profile = parasail_profile_create_altivec_128_64(s1, s1Len, matrix);
    if (!profile) return NULL;
    result = PNAME(profile, s2, s2Len, open, gap);

    parasail_profile_free(profile);

    return result;
}

parasail_result_t* PNAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    /* declare local variables */
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t s1Len = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    const parasail_matrix_t *matrix = NULL;
    int32_t segWidth = 0;
    int32_t segLen = 0;
    int32_t offset = 0;
    int32_t position = 0;
    vec128i* restrict pvP = NULL;
    vec128i* restrict pvE = NULL;
    int64_t* restrict boundary = NULL;
    vec128i* restrict pvHt = NULL;
    vec128i* restrict pvH = NULL;
    vec128i* restrict pvGapper = NULL;
    vec128i vGapO;
    vec128i vGapE;
    int64_t NEG_LIMIT = 0;
    int64_t POS_LIMIT = 0;
    vec128i vZero;
    int64_t score = 0;
    vec128i vNegLimit;
    vec128i vPosLimit;
    vec128i vSaturationCheckMin;
    vec128i vSaturationCheckMax;
    vec128i vNegInfFront;
    vec128i vSegLenXgap;
    parasail_result_t *result = NULL;
    vec128i vTIns;
    vec128i vTDel;
    vec128i vTDiag;
    vec128i vTDiagE;
    vec128i vTInsE;
    vec128i vTDiagF;
    vec128i vTDelF;

    /* validate inputs */
    PARASAIL_CHECK_NULL(profile);
    PARASAIL_CHECK_NULL(profile->profile64.score);
    PARASAIL_CHECK_NULL(profile->matrix);
    PARASAIL_CHECK_GT0(profile->s1Len);
    PARASAIL_CHECK_NULL(s2);
    PARASAIL_CHECK_GT0(s2Len);
    PARASAIL_CHECK_GE0(open);
    PARASAIL_CHECK_GE0(gap);

    /* initialize stack variables */
    i = 0;
    j = 0;
    k = 0;
    s1Len = profile->s1Len;
    end_query = s1Len-1;
    end_ref = s2Len-1;
    matrix = profile->matrix;
    segWidth = 2; /* number of values in vector unit */
    segLen = (s1Len + segWidth - 1) / segWidth;
    offset = (s1Len - 1) % segLen;
    position = (segWidth - 1) - (s1Len - 1) / segLen;
    pvP = (vec128i*)profile->profile64.score;
    vGapO = _mm_set1_epi64(open);
    vGapE = _mm_set1_epi64(gap);
    NEG_LIMIT = (-open < matrix->min ? INT64_MIN + open : INT64_MIN - matrix->min) + 1;
    POS_LIMIT = INT64_MAX - matrix->max - 1;
    vZero = _mm_setzero_si128();
    score = NEG_LIMIT;
    vNegLimit = _mm_set1_epi64(NEG_LIMIT);
    vPosLimit = _mm_set1_epi64(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;
    vNegInfFront = vZero;
    vTIns  = _mm_set1_epi64(PARASAIL_INS);
    vTDel  = _mm_set1_epi64(PARASAIL_DEL);
    vTDiag = _mm_set1_epi64(PARASAIL_DIAG);
    vTDiagE = _mm_set1_epi64(PARASAIL_DIAG_E);
    vTInsE = _mm_set1_epi64(PARASAIL_INS_E);
    vTDiagF = _mm_set1_epi64(PARASAIL_DIAG_F);
    vTDelF = _mm_set1_epi64(PARASAIL_DEL_F);
    vNegInfFront = _mm_insert_epi64(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = _mm_add_epi64(vNegInfFront,
            _mm_slli_si128(_mm_set1_epi64(-segLen*gap), 8));

    /* initialize result */
    result = parasail_result_new_trace(segLen, s2Len, 16, sizeof(vec128i));
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_SCAN
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_2;

    /* initialize heap variables */
    pvE = parasail_memalign_vec128i(16, segLen);
    boundary = parasail_memalign_int64_t(16, s2Len+1);
    pvHt= parasail_memalign_vec128i(16, segLen);
    pvH = parasail_memalign_vec128i(16, segLen);
    pvGapper = parasail_memalign_vec128i(16, segLen);

    /* validate heap variables */
    if (!pvE) return NULL;
    if (!boundary) return NULL;
    if (!pvHt) return NULL;
    if (!pvH) return NULL;
    if (!pvGapper) return NULL;

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            int32_t segNum = 0;
            vec128i_64_t h;
            vec128i_64_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.v[segNum] = tmp < INT64_MIN ? INT64_MIN : tmp;
                tmp = tmp - open;
                e.v[segNum] = tmp < INT64_MIN ? INT64_MIN : tmp;
            }
            _mm_store_si128(&pvH[index], h.m);
            _mm_store_si128(&pvE[index], e.m);
            ++index;
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            int64_t tmp = -open-gap*(i-1);
            boundary[i] = tmp < INT64_MIN ? INT64_MIN : tmp;
        }
    }

    {
        vec128i vGapper = _mm_sub_epi64(vZero,vGapO);
        for (i=segLen-1; i>=0; --i) {
            _mm_store_si128(pvGapper+i, vGapper);
            vGapper = _mm_sub_epi64(vGapper, vGapE);
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        vec128i vE;
        vec128i vE_ext;
        vec128i vE_opn;
        vec128i vHt;
        vec128i vF;
        vec128i vF_ext;
        vec128i vF_opn;
        vec128i vH;
        vec128i vHp;
        vec128i *pvW;
        vec128i vW;
        vec128i case1;
        vec128i case2;
        vec128i vGapper;
        vec128i vT;
        vec128i vET;
        vec128i vFT;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = _mm_load_si128(pvH+(segLen-1));
        vHp = _mm_slli_si128(vHp, 8);
        vHp = _mm_insert_epi64(vHp, boundary[j], 0);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = _mm_sub_epi64(vNegLimit, pvGapper[0]);
        vF = vNegLimit;
        for (i=0; i<segLen; ++i) {
            vH = _mm_load_si128(pvH+i);
            vE = _mm_load_si128(pvE+i);
            vW = _mm_load_si128(pvW+i);
            vGapper = _mm_load_si128(pvGapper+i);
            vE_opn = _mm_sub_epi64(vH, vGapO);
            vE_ext = _mm_sub_epi64(vE, vGapE);
            case1 = _mm_cmpgt_epi64(vE_opn, vE_ext);
            vET = _mm_blendv_epi8(vTInsE, vTDiagE, case1);
            arr_store(result->trace->trace_table, vET, i, segLen, j);
            vE = _mm_max_epi64(vE_opn, vE_ext);
            vSaturationCheckMin = _mm_min_epi64(vSaturationCheckMin, vE);
            vGapper = _mm_add_epi64(vHt, vGapper);
            vF = _mm_max_epi64(vF, vGapper);
            vHp = _mm_add_epi64(vHp, vW);
            vHt = _mm_max_epi64(vE, vHp);
            _mm_store_si128(pvE+i, vE);
            _mm_store_si128(pvHt+i, vHt);
            _mm_store_si128(pvH+i, vHp);
            vHp = vH;
        }

        /* pseudo prefix scan on F and H */
        vHt = _mm_slli_si128(vHt, 8);
        vHt = _mm_insert_epi64(vHt, boundary[j+1], 0);
        vGapper = _mm_load_si128(pvGapper);
        vGapper = _mm_add_epi64(vHt, vGapper);
        vF = _mm_max_epi64(vF, vGapper);
        for (i=0; i<segWidth-2; ++i) {
            vec128i vFt = _mm_slli_si128(vF, 8);
            vFt = _mm_add_epi64(vFt, vSegLenXgap);
            vF = _mm_max_epi64(vF, vFt);
        }

        /* calculate final H */
        vF = _mm_slli_si128(vF, 8);
        vF = _mm_add_epi64(vF, vNegInfFront);
        vH = _mm_max_epi64(vF, vHt);
        for (i=0; i<segLen; ++i) {
            vET = arr_load(result->trace->trace_table, i, segLen, j);
            vHp = _mm_load_si128(pvH+i);
            vHt = _mm_load_si128(pvHt+i);
            vF_opn = _mm_sub_epi64(vH, vGapO);
            vF_ext = _mm_sub_epi64(vF, vGapE);
            vF = _mm_max_epi64(vF_opn, vF_ext);
            case1 = _mm_cmpgt_epi64(vF_opn, vF_ext);
            vFT = _mm_blendv_epi8(vTDelF, vTDiagF, case1);
            vH = _mm_max_epi64(vHt, vF);
            case1 = _mm_cmpeq_epi64(vH, vHp);
            case2 = _mm_cmpeq_epi64(vH, vF);
            vT = _mm_blendv_epi8(
                    _mm_blendv_epi8(vTIns, vTDel, case2),
                    vTDiag, case1);
            vT = _mm_or_si128(vT, vET);
            vT = _mm_or_si128(vT, vFT);
            arr_store(result->trace->trace_table, vT, i, segLen, j);
            _mm_store_si128(pvH+i, vH);
            vSaturationCheckMin = _mm_min_epi64(vSaturationCheckMin, vH);
            vSaturationCheckMin = _mm_min_epi64(vSaturationCheckMin, vF);
            vSaturationCheckMax = _mm_max_epi64(vSaturationCheckMax, vH);
        }
    }

    /* extract last value from the last column */
    {
        vec128i vH = _mm_load_si128(pvH + offset);
        for (k=0; k<position; ++k) {
            vH = _mm_slli_si128(vH, 8);
        }
        score = (int64_t) _mm_extract_epi64 (vH, 1);
    }

    if (_mm_movemask_epi8(_mm_or_si128(
            _mm_cmplt_epi64(vSaturationCheckMin, vNegLimit),
            _mm_cmpgt_epi64(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(pvGapper);
    parasail_free(pvH);
    parasail_free(pvHt);
    parasail_free(boundary);
    parasail_free(pvE);

    return result;
}


