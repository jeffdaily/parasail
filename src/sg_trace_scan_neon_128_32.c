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
#include "parasail/internal_neon.h"

#define SG_TRACE
#define SG_SUFFIX _scan_neon_128_32
#define SG_SUFFIX_PROF _scan_profile_neon_128_32
#include "sg_helper.h"



static inline void arr_store(
        simde__m128i *array,
        simde__m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    simde_mm_store_si128(array + (1LL*d*seglen+t), vH);
}

static inline simde__m128i arr_load(
        simde__m128i *array,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    return simde_mm_load_si128(array + (1LL*d*seglen+t));
}

#define FNAME parasail_sg_flags_trace_scan_neon_128_32
#define PNAME parasail_sg_flags_trace_scan_profile_neon_128_32

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix,
        int s1_beg, int s1_end, int s2_beg, int s2_end)
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
    profile = parasail_profile_create_neon_128_32(s1, s1Len, matrix);
    if (!profile) return NULL;
    result = PNAME(profile, s2, s2Len, open, gap, s1_beg, s1_end, s2_beg, s2_end);

    parasail_profile_free(profile);

    return result;
}

parasail_result_t* PNAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        int s1_beg, int s1_end, int s2_beg, int s2_end)
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
    simde__m128i* restrict pvP = NULL;
    simde__m128i* restrict pvE = NULL;
    int32_t* restrict boundary = NULL;
    simde__m128i* restrict pvHt = NULL;
    simde__m128i* restrict pvH = NULL;
    simde__m128i* restrict pvGapper = NULL;
    simde__m128i vGapO;
    simde__m128i vGapE;
    int32_t NEG_LIMIT = 0;
    int32_t POS_LIMIT = 0;
    simde__m128i vZero;
    int32_t score = 0;
    simde__m128i vNegLimit;
    simde__m128i vPosLimit;
    simde__m128i vSaturationCheckMin;
    simde__m128i vSaturationCheckMax;
    simde__m128i vMaxH;
    simde__m128i vPosMask;
    simde__m128i vNegInfFront;
    simde__m128i vSegLenXgap;
    parasail_result_t *result = NULL;
    simde__m128i vTIns;
    simde__m128i vTDel;
    simde__m128i vTDiag;
    simde__m128i vTDiagE;
    simde__m128i vTInsE;
    simde__m128i vTDiagF;
    simde__m128i vTDelF;

    /* validate inputs */
    PARASAIL_CHECK_NULL(profile);
    PARASAIL_CHECK_NULL(profile->profile32.score);
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
    segWidth = 4; /* number of values in vector unit */
    segLen = (s1Len + segWidth - 1) / segWidth;
    offset = (s1Len - 1) % segLen;
    position = (segWidth - 1) - (s1Len - 1) / segLen;
    pvP = (simde__m128i*)profile->profile32.score;
    vGapO = simde_mm_set1_epi32(open);
    vGapE = simde_mm_set1_epi32(gap);
    NEG_LIMIT = (-open < matrix->min ? INT32_MIN + open : INT32_MIN - matrix->min) + 1;
    POS_LIMIT = INT32_MAX - matrix->max - 1;
    vZero = simde_mm_setzero_si128();
    score = NEG_LIMIT;
    vNegLimit = simde_mm_set1_epi32(NEG_LIMIT);
    vPosLimit = simde_mm_set1_epi32(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;
    vMaxH = vNegLimit;
    vPosMask = simde_mm_cmpeq_epi32(simde_mm_set1_epi32(position),
            simde_mm_set_epi32(0,1,2,3));
    vNegInfFront = vZero;
    vNegInfFront = simde_mm_insert_epi32(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = simde_mm_add_epi32(vNegInfFront,
            simde_mm_slli_si128(simde_mm_set1_epi32(-segLen*gap), 4));
    vTIns  = simde_mm_set1_epi32(PARASAIL_INS);
    vTDel  = simde_mm_set1_epi32(PARASAIL_DEL);
    vTDiag = simde_mm_set1_epi32(PARASAIL_DIAG);
    vTDiagE= simde_mm_set1_epi32(PARASAIL_DIAG_E);
    vTInsE = simde_mm_set1_epi32(PARASAIL_INS_E);
    vTDiagF= simde_mm_set1_epi32(PARASAIL_DIAG_F);
    vTDelF = simde_mm_set1_epi32(PARASAIL_DEL_F);

    /* initialize result */
    result = parasail_result_new_trace(segLen, s2Len, 16, sizeof(simde__m128i));
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_SCAN
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_32 | PARASAIL_FLAG_LANES_4;
    result->flag |= s1_beg ? PARASAIL_FLAG_SG_S1_BEG : 0;
    result->flag |= s1_end ? PARASAIL_FLAG_SG_S1_END : 0;
    result->flag |= s2_beg ? PARASAIL_FLAG_SG_S2_BEG : 0;
    result->flag |= s2_end ? PARASAIL_FLAG_SG_S2_END : 0;

    /* initialize heap variables */
    pvE = parasail_memalign_simde__m128i(16, segLen);
    boundary = parasail_memalign_int32_t(16, s2Len+1);
    pvHt= parasail_memalign_simde__m128i(16, segLen);
    pvH = parasail_memalign_simde__m128i(16, segLen);
    pvGapper = parasail_memalign_simde__m128i(16, segLen);

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
            simde__m128i_private h_;
            simde__m128i_private e_;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = s1_beg ? 0 : (-open-gap*(segNum*segLen+i));
                h_.i32[segNum] = tmp < INT32_MIN ? INT32_MIN : tmp;
                tmp = tmp - open;
                e_.i32[segNum] = tmp < INT32_MIN ? INT32_MIN : tmp;
            }
            simde_mm_store_si128(&pvH[index], simde__m128i_from_private(h_));
            simde_mm_store_si128(&pvE[index], simde__m128i_from_private(e_));
            ++index;
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            int64_t tmp = s2_beg ? 0 : (-open-gap*(i-1));
            boundary[i] = tmp < INT32_MIN ? INT32_MIN : tmp;
        }
    }

    {
        simde__m128i vGapper = simde_mm_sub_epi32(vZero,vGapO);
        for (i=segLen-1; i>=0; --i) {
            simde_mm_store_si128(pvGapper+i, vGapper);
            vGapper = simde_mm_sub_epi32(vGapper, vGapE);
            /* long queries and/or large penalties will break the pseudo prefix scan */
            vSaturationCheckMin = simde_mm_min_epi32(vSaturationCheckMin, vGapper);
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        simde__m128i vE;
        simde__m128i vE_ext;
        simde__m128i vE_opn;
        simde__m128i vHt;
        simde__m128i vF;
        simde__m128i vF_ext;
        simde__m128i vF_opn;
        simde__m128i vH;
        simde__m128i vHp;
        simde__m128i *pvW;
        simde__m128i vW;
        simde__m128i case1;
        simde__m128i case2;
        simde__m128i vGapper;
        simde__m128i vT;
        simde__m128i vET;
        simde__m128i vFT;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = simde_mm_load_si128(pvH+(segLen-1));
        vHp = simde_mm_slli_si128(vHp, 4);
        vHp = simde_mm_insert_epi32(vHp, boundary[j], 0);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = simde_mm_sub_epi32(vNegLimit, pvGapper[0]);
        vF = vNegLimit;
        for (i=0; i<segLen; ++i) {
            vH = simde_mm_load_si128(pvH+i);
            vE = simde_mm_load_si128(pvE+i);
            vW = simde_mm_load_si128(pvW+i);
            vGapper = simde_mm_load_si128(pvGapper+i);
            vE_opn = simde_mm_sub_epi32(vH, vGapO);
            vE_ext = simde_mm_sub_epi32(vE, vGapE);
            case1 = simde_mm_cmpgt_epi32(vE_opn, vE_ext);
            vET = simde_mm_blendv_epi8(vTInsE, vTDiagE, case1);
            arr_store(result->trace->trace_table, vET, i, segLen, j);
            vE = simde_mm_max_epi32(vE_opn, vE_ext);
            vSaturationCheckMin = simde_mm_min_epi32(vSaturationCheckMin, vE);
            vGapper = simde_mm_add_epi32(vHt, vGapper);
            vF = simde_mm_max_epi32(vF, vGapper);
            vHp = simde_mm_add_epi32(vHp, vW);
            vHt = simde_mm_max_epi32(vE, vHp);
            simde_mm_store_si128(pvE+i, vE);
            simde_mm_store_si128(pvHt+i, vHt);
            simde_mm_store_si128(pvH+i, vHp);
            vHp = vH;
        }

        /* pseudo prefix scan on F and H */
        vHt = simde_mm_slli_si128(vHt, 4);
        vHt = simde_mm_insert_epi32(vHt, boundary[j+1], 0);
        vGapper = simde_mm_load_si128(pvGapper);
        vGapper = simde_mm_add_epi32(vHt, vGapper);
        vF = simde_mm_max_epi32(vF, vGapper);
        for (i=0; i<segWidth-2; ++i) {
            simde__m128i vFt = simde_mm_slli_si128(vF, 4);
            vFt = simde_mm_add_epi32(vFt, vSegLenXgap);
            vF = simde_mm_max_epi32(vF, vFt);
        }

        /* calculate final H */
        vF = simde_mm_slli_si128(vF, 4);
        vF = simde_mm_add_epi32(vF, vNegInfFront);
        vH = simde_mm_max_epi32(vF, vHt);
        for (i=0; i<segLen; ++i) {
            vET = arr_load(result->trace->trace_table, i, segLen, j);
            vHp = simde_mm_load_si128(pvH+i);
            vHt = simde_mm_load_si128(pvHt+i);
            vF_opn = simde_mm_sub_epi32(vH, vGapO);
            vF_ext = simde_mm_sub_epi32(vF, vGapE);
            vF = simde_mm_max_epi32(vF_opn, vF_ext);
            case1 = simde_mm_cmpgt_epi32(vF_opn, vF_ext);
            vFT = simde_mm_blendv_epi8(vTDelF, vTDiagF, case1);
            vH = simde_mm_max_epi32(vHt, vF);
            case1 = simde_mm_cmpeq_epi32(vH, vHp);
            case2 = simde_mm_cmpeq_epi32(vH, vF);
            vT = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vTIns, vTDel, case2),
                    vTDiag, case1);
            vT = simde_mm_or_si128(vT, vET);
            vT = simde_mm_or_si128(vT, vFT);
            arr_store(result->trace->trace_table, vT, i, segLen, j);
            simde_mm_store_si128(pvH+i, vH);
            vSaturationCheckMin = simde_mm_min_epi32(vSaturationCheckMin, vH);
            vSaturationCheckMin = simde_mm_min_epi32(vSaturationCheckMin, vF);
            vSaturationCheckMax = simde_mm_max_epi32(vSaturationCheckMax, vH);
        }

        /* extract vector containing last value from column */
        {
            simde__m128i vCompare;
            vH = simde_mm_load_si128(pvH + offset);
            vCompare = simde_mm_and_si128(vPosMask, simde_mm_cmpgt_epi32(vH, vMaxH));
            vMaxH = simde_mm_max_epi32(vH, vMaxH);
            if (simde_mm_movemask_epi8(vCompare)) {
                end_ref = j;
            }
        }
    } 

    /* max last value from all columns */
    if (s2_end)
    {
        for (k=0; k<position; ++k) {
            vMaxH = simde_mm_slli_si128(vMaxH, 4);
        }
        score = (int32_t) simde_mm_extract_epi32(vMaxH, 3);
        end_query = s1Len-1;
    }

    /* max of last column */
    if (s1_end)
    {
        /* Trace the alignment ending position on read. */
        int32_t *t = (int32_t*)pvH;
        int32_t column_len = segLen * segWidth;
        for (i = 0; i<column_len; ++i, ++t) {
            int32_t temp = i / segWidth + i % segWidth * segLen;
            if (temp >= s1Len) continue;
            if (*t > score) {
                score = *t;
                end_query = temp;
                end_ref = s2Len-1;
            }
            else if (*t == score && end_ref == s2Len-1 && temp < end_query) {
                end_query = temp;
            }
        }
    }

    if (!s1_end && !s2_end) {
        /* extract last value from the last column */
        {
            simde__m128i vH = simde_mm_load_si128(pvH + offset);
            for (k=0; k<position; ++k) {
                vH = simde_mm_slli_si128(vH, 4);
            }
            score = (int32_t) simde_mm_extract_epi32 (vH, 3);
            end_ref = s2Len - 1;
            end_query = s1Len - 1;
        }
    }

    if (simde_mm_movemask_epi8(simde_mm_or_si128(
            simde_mm_cmplt_epi32(vSaturationCheckMin, vNegLimit),
            simde_mm_cmpgt_epi32(vSaturationCheckMax, vPosLimit)))) {
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

SG_IMPL_ALL
SG_IMPL_PROF_ALL


