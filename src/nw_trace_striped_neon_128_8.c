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

#define SWAP(A,B) { simde__m128i* tmp = A; A = B; B = tmp; }

#define NEG_INF INT8_MIN


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

#define FNAME parasail_nw_trace_striped_neon_128_8
#define PNAME parasail_nw_trace_striped_profile_neon_128_8

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    /* declare local variables */
    parasail_profile_t *profile = NULL;
    parasail_result_t *result = NULL;

    /* validate inputs */
    PARASAIL_CHECK_NULL(s1);
    PARASAIL_CHECK_GT0(s1Len);
    PARASAIL_CHECK_NULL(s2);
    PARASAIL_CHECK_GT0(s2Len);
    PARASAIL_CHECK_GE0(open);
    PARASAIL_CHECK_GE0(gap);
    PARASAIL_CHECK_NULL(matrix);

    /* initialize local variables */
    profile = parasail_profile_create_neon_128_8(s1, s1Len, matrix);
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
    int s1Len = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    const parasail_matrix_t *matrix = NULL;
    int32_t segWidth = 0;
    int32_t segLen = 0;
    int32_t offset = 0;
    int32_t position = 0;
    simde__m128i* restrict vProfile = NULL;
    simde__m128i* restrict pvHStore = NULL;
    simde__m128i* restrict pvHLoad = NULL;
    simde__m128i* restrict pvE = NULL;
    simde__m128i* restrict pvEaStore = NULL;
    simde__m128i* restrict pvEaLoad = NULL;
    simde__m128i* restrict pvHT = NULL;
    int8_t* restrict boundary = NULL;
    simde__m128i vGapO;
    simde__m128i vGapE;
    simde__m128i vNegInf;
    int8_t score = 0;
    simde__m128i vNegLimit;
    simde__m128i vPosLimit;
    simde__m128i vSaturationCheckMin;
    simde__m128i vSaturationCheckMax;
    parasail_result_t *result = NULL;
    simde__m128i vTIns;
    simde__m128i vTDel;
    simde__m128i vTDiag;
    simde__m128i vTDiagE;
    simde__m128i vTInsE;
    simde__m128i vTDiagF;
    simde__m128i vTDelF;
    simde__m128i vTMask;
    simde__m128i vFTMask;
    
    /* validate inputs */
    PARASAIL_CHECK_NULL(profile);
    PARASAIL_CHECK_NULL(profile->profile8.score);
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
    segWidth = 16; /* number of values in vector unit */
    segLen = (s1Len + segWidth - 1) / segWidth;
    offset = (s1Len - 1) % segLen;
    position = (segWidth - 1) - (s1Len - 1) / segLen;
    vProfile = (simde__m128i*)profile->profile8.score;
    vGapO = simde_mm_set1_epi8(open);
    vGapE = simde_mm_set1_epi8(gap);
    vNegInf = simde_mm_set1_epi8(NEG_INF);
    score = NEG_INF;
    vTIns  = simde_mm_set1_epi8(PARASAIL_INS);
    vTDel  = simde_mm_set1_epi8(PARASAIL_DEL);
    vTDiag = simde_mm_set1_epi8(PARASAIL_DIAG);
    vTDiagE = simde_mm_set1_epi8(PARASAIL_DIAG_E);
    vTInsE = simde_mm_set1_epi8(PARASAIL_INS_E);
    vTDiagF = simde_mm_set1_epi8(PARASAIL_DIAG_F);
    vTDelF = simde_mm_set1_epi8(PARASAIL_DEL_F);
    vTMask = simde_mm_set1_epi8(PARASAIL_ZERO_MASK);
    vFTMask = simde_mm_set1_epi8(PARASAIL_F_MASK);
    vNegLimit = simde_mm_set1_epi8(INT8_MIN);
    vPosLimit = simde_mm_set1_epi8(INT8_MAX);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;

    /* initialize result */
    result = parasail_result_new_trace(segLen, s2Len, 16, sizeof(simde__m128i));
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_16;

    /* initialize heap variables */
    pvHStore = parasail_memalign_simde__m128i(16, segLen);
    pvHLoad = parasail_memalign_simde__m128i(16, segLen);
    pvE = parasail_memalign_simde__m128i(16, segLen);
    pvEaStore = parasail_memalign_simde__m128i(16, segLen);
    pvEaLoad = parasail_memalign_simde__m128i(16, segLen);
    pvHT = parasail_memalign_simde__m128i(16, segLen);
    boundary = parasail_memalign_int8_t(16, s2Len+1);

    /* validate heap variables */
    if (!pvHStore) return NULL;
    if (!pvHLoad) return NULL;
    if (!pvE) return NULL;
    if (!pvEaStore) return NULL;
    if (!pvEaLoad) return NULL;
    if (!pvHT) return NULL;
    if (!boundary) return NULL;

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            int32_t segNum = 0;
            simde__m128i h;
            simde__m128i e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.i8[segNum] = tmp < INT8_MIN ? INT8_MIN : tmp;
                tmp = tmp - open;
                e.i8[segNum] = tmp < INT8_MIN ? INT8_MIN : tmp;
            }
            simde_mm_store_si128(&pvHStore[index], h);
            simde_mm_store_si128(&pvE[index], e);
            simde_mm_store_si128(&pvEaStore[index], e);
            ++index;
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            int64_t tmp = -open-gap*(i-1);
            boundary[i] = tmp < INT8_MIN ? INT8_MIN : tmp;
        }
    }

    for (i=0; i<segLen; ++i) {
        arr_store(result->trace->trace_table, vTDiagE, i, segLen, 0);
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        simde__m128i vEF_opn;
        simde__m128i vE;
        simde__m128i vE_ext;
        simde__m128i vF;
        simde__m128i vF_ext;
        simde__m128i vFa;
        simde__m128i vFa_ext;
        simde__m128i vH;
        simde__m128i vH_dag;
        const simde__m128i* vP = NULL;

        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        vF = vNegInf;

        /* load final segment of pvHStore and shift left by 1 bytes */
        vH = simde_mm_load_si128(&pvHStore[segLen - 1]);
        vH = simde_mm_slli_si128(vH, 1);

        /* insert upper boundary condition */
        vH = simde_mm_insert_epi8(vH, boundary[j], 0);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        SWAP(pvHLoad, pvHStore)
        SWAP(pvEaLoad, pvEaStore)

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vE = simde_mm_load_si128(pvE + i);

            /* Get max from vH, vE and vF. */
            vH_dag = simde_mm_adds_epi8(vH, simde_mm_load_si128(vP + i));
            vH = simde_mm_max_epi8(vH_dag, vE);
            vH = simde_mm_max_epi8(vH, vF);
            /* Save vH values. */
            simde_mm_store_si128(pvHStore + i, vH);
            /* check for saturation */
            {
                vSaturationCheckMax = simde_mm_max_epi8(vSaturationCheckMax, vH);
                vSaturationCheckMin = simde_mm_min_epi8(vSaturationCheckMin, vH);
                vSaturationCheckMin = simde_mm_min_epi8(vSaturationCheckMin, vE);
                vSaturationCheckMin = simde_mm_min_epi8(vSaturationCheckMin, vF);
            }

            {
                simde__m128i vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                simde__m128i case1 = simde_mm_cmpeq_epi8(vH, vH_dag);
                simde__m128i case2 = simde_mm_cmpeq_epi8(vH, vF);
                simde__m128i vT = simde_mm_blendv_epi8(
                        simde_mm_blendv_epi8(vTIns, vTDel, case2),
                        vTDiag, case1);
                simde_mm_store_si128(pvHT + i, vT);
                vT = simde_mm_or_si128(vT, vTAll);
                arr_store(result->trace->trace_table, vT, i, segLen, j);
            }

            vEF_opn = simde_mm_subs_epi8(vH, vGapO);

            /* Update vE value. */
            vE_ext = simde_mm_subs_epi8(vE, vGapE);
            vE = simde_mm_max_epi8(vEF_opn, vE_ext);
            simde_mm_store_si128(pvE + i, vE);
            {
                simde__m128i vEa = simde_mm_load_si128(pvEaLoad + i);
                simde__m128i vEa_ext = simde_mm_subs_epi8(vEa, vGapE);
                vEa = simde_mm_max_epi8(vEF_opn, vEa_ext);
                simde_mm_store_si128(pvEaStore + i, vEa);
                if (j+1<s2Len) {
                    simde__m128i cond = simde_mm_cmpgt_epi8(vEF_opn, vEa_ext);
                    simde__m128i vT = simde_mm_blendv_epi8(vTInsE, vTDiagE, cond);
                    arr_store(result->trace->trace_table, vT, i, segLen, j+1);
                }
            }

            /* Update vF value. */
            vF_ext = simde_mm_subs_epi8(vF, vGapE);
            vF = simde_mm_max_epi8(vEF_opn, vF_ext);
            if (i+1<segLen) {
                simde__m128i vTAll = arr_load(result->trace->trace_table, i+1, segLen, j);
                simde__m128i cond = simde_mm_cmpgt_epi8(vEF_opn, vF_ext);
                simde__m128i vT = simde_mm_blendv_epi8(vTDelF, vTDiagF, cond);
                vT = simde_mm_or_si128(vT, vTAll);
                arr_store(result->trace->trace_table, vT, i+1, segLen, j);
            }

            /* Load the next vH. */
            vH = simde_mm_load_si128(pvHLoad + i);
        }


        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        vFa_ext = vF_ext;
        vFa = vF;
        for (k=0; k<segWidth; ++k) {
            int64_t tmp = boundary[j+1]-open;
            int8_t tmp2 = tmp < INT8_MIN ? INT8_MIN : tmp;
            simde__m128i vHp = simde_mm_load_si128(&pvHLoad[segLen - 1]);
            vHp = simde_mm_slli_si128(vHp, 1);
            vHp = simde_mm_insert_epi8(vHp, boundary[j], 0);
            vEF_opn = simde_mm_slli_si128(vEF_opn, 1);
            vEF_opn = simde_mm_insert_epi8(vEF_opn, tmp2, 0);
            vF_ext = simde_mm_slli_si128(vF_ext, 1);
            vF_ext = simde_mm_insert_epi8(vF_ext, NEG_INF, 0);
            vF = simde_mm_slli_si128(vF, 1);
            vF = simde_mm_insert_epi8(vF, tmp2, 0);
            vFa_ext = simde_mm_slli_si128(vFa_ext, 1);
            vFa_ext = simde_mm_insert_epi8(vFa_ext, NEG_INF, 0);
            vFa = simde_mm_slli_si128(vFa, 1);
            vFa = simde_mm_insert_epi8(vFa, tmp2, 0);
            for (i=0; i<segLen; ++i) {
                vH = simde_mm_load_si128(pvHStore + i);
                vH = simde_mm_max_epi8(vH,vF);
                simde_mm_store_si128(pvHStore + i, vH);
                /* check for saturation */
            {
                vSaturationCheckMax = simde_mm_max_epi8(vSaturationCheckMax, vH);
                vSaturationCheckMin = simde_mm_min_epi8(vSaturationCheckMin, vH);
                vSaturationCheckMin = simde_mm_min_epi8(vSaturationCheckMin, vE);
                vSaturationCheckMin = simde_mm_min_epi8(vSaturationCheckMin, vF);
            }
                {
                    simde__m128i vTAll;
                    simde__m128i vT;
                    simde__m128i case1;
                    simde__m128i case2;
                    simde__m128i cond;
                    vHp = simde_mm_adds_epi8(vHp, simde_mm_load_si128(vP + i));
                    case1 = simde_mm_cmpeq_epi8(vH, vHp);
                    case2 = simde_mm_cmpeq_epi8(vH, vF);
                    cond = simde_mm_andnot_si128(case1,case2);
                    vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                    vT = simde_mm_load_si128(pvHT + i);
                    vT = simde_mm_blendv_epi8(vT, vTDel, cond);
                    simde_mm_store_si128(pvHT + i, vT);
                    vTAll = simde_mm_and_si128(vTAll, vTMask);
                    vTAll = simde_mm_or_si128(vTAll, vT);
                    arr_store(result->trace->trace_table, vTAll, i, segLen, j);
                }
                /* Update vF value. */
                {
                    simde__m128i vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                    simde__m128i cond = simde_mm_cmpgt_epi8(vEF_opn, vFa_ext);
                    simde__m128i vT = simde_mm_blendv_epi8(vTDelF, vTDiagF, cond);
                    vTAll = simde_mm_and_si128(vTAll, vFTMask);
                    vTAll = simde_mm_or_si128(vTAll, vT);
                    arr_store(result->trace->trace_table, vTAll, i, segLen, j);
                }
                vEF_opn = simde_mm_subs_epi8(vH, vGapO);
                vF_ext = simde_mm_subs_epi8(vF, vGapE);
                {
                    simde__m128i vEa = simde_mm_load_si128(pvEaLoad + i);
                    simde__m128i vEa_ext = simde_mm_subs_epi8(vEa, vGapE);
                    vEa = simde_mm_max_epi8(vEF_opn, vEa_ext);
                    simde_mm_store_si128(pvEaStore + i, vEa);
                    if (j+1<s2Len) {
                        simde__m128i cond = simde_mm_cmpgt_epi8(vEF_opn, vEa_ext);
                        simde__m128i vT = simde_mm_blendv_epi8(vTInsE, vTDiagE, cond);
                        arr_store(result->trace->trace_table, vT, i, segLen, j+1);
                    }
                }
                if (! simde_mm_movemask_epi8(
                            simde_mm_or_si128(
                                simde_mm_cmpgt_epi8(vF_ext, vEF_opn),
                                simde_mm_cmpeq_epi8(vF_ext, vEF_opn))))
                    goto end;
                /*vF = simde_mm_max_epi8(vEF_opn, vF_ext);*/
                vF = vF_ext;
                vFa_ext = simde_mm_subs_epi8(vFa, vGapE);
                vFa = simde_mm_max_epi8(vEF_opn, vFa_ext);
                vHp = simde_mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
        }
    }

    /* extract last value from the last column */
    {
        simde__m128i vH = simde_mm_load_si128(pvHStore + offset);
        for (k=0; k<position; ++k) {
            vH = simde_mm_slli_si128 (vH, 1);
        }
        score = (int8_t) simde_mm_extract_epi8 (vH, 15);
    }

    if (simde_mm_movemask_epi8(simde_mm_or_si128(
            simde_mm_cmpeq_epi8(vSaturationCheckMin, vNegLimit),
            simde_mm_cmpeq_epi8(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = INT8_MAX;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(boundary);
    parasail_free(pvHT);
    parasail_free(pvEaLoad);
    parasail_free(pvEaStore);
    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}


