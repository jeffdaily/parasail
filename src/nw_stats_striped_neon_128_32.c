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



#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        simde__m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[1LL*(0*seglen+t)*dlen + d] = (int32_t)simde_mm_extract_epi32(vH, 0);
    array[1LL*(1*seglen+t)*dlen + d] = (int32_t)simde_mm_extract_epi32(vH, 1);
    array[1LL*(2*seglen+t)*dlen + d] = (int32_t)simde_mm_extract_epi32(vH, 2);
    array[1LL*(3*seglen+t)*dlen + d] = (int32_t)simde_mm_extract_epi32(vH, 3);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        simde__m128i vH,
        int32_t t,
        int32_t seglen)
{
    col[0*seglen+t] = (int32_t)simde_mm_extract_epi32(vH, 0);
    col[1*seglen+t] = (int32_t)simde_mm_extract_epi32(vH, 1);
    col[2*seglen+t] = (int32_t)simde_mm_extract_epi32(vH, 2);
    col[3*seglen+t] = (int32_t)simde_mm_extract_epi32(vH, 3);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_nw_stats_table_striped_neon_128_32
#define PNAME parasail_nw_stats_table_striped_profile_neon_128_32
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_stats_rowcol_striped_neon_128_32
#define PNAME parasail_nw_stats_rowcol_striped_profile_neon_128_32
#else
#define FNAME parasail_nw_stats_striped_neon_128_32
#define PNAME parasail_nw_stats_striped_profile_neon_128_32
#endif
#endif

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
    if (matrix->type == PARASAIL_MATRIX_TYPE_PSSM) {
        PARASAIL_CHECK_NULL_PSSM_STATS(s1);
    }
    else {
        PARASAIL_CHECK_NULL(s1);
        PARASAIL_CHECK_GT0(s1Len);
    }

    /* initialize local variables */
    profile = parasail_profile_create_stats_neon_128_32(s1, s1Len, matrix);
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
    simde__m128i* restrict vProfile = NULL;
    simde__m128i* restrict vProfileM = NULL;
    simde__m128i* restrict vProfileS = NULL;
    simde__m128i* restrict pvHStore = NULL;
    simde__m128i* restrict pvHLoad = NULL;
    simde__m128i* restrict pvHMStore = NULL;
    simde__m128i* restrict pvHMLoad = NULL;
    simde__m128i* restrict pvHSStore = NULL;
    simde__m128i* restrict pvHSLoad = NULL;
    simde__m128i* restrict pvHLStore = NULL;
    simde__m128i* restrict pvHLLoad = NULL;
    simde__m128i* restrict pvE = NULL;
    simde__m128i* restrict pvEM = NULL;
    simde__m128i* restrict pvES = NULL;
    simde__m128i* restrict pvEL = NULL;
    int32_t* restrict boundary = NULL;
    simde__m128i vGapO;
    simde__m128i vGapE;
    int32_t NEG_LIMIT = 0;
    int32_t POS_LIMIT = 0;
    simde__m128i vZero;
    simde__m128i vOne;
    int32_t score = 0;
    int32_t matches = 0;
    int32_t similar = 0;
    int32_t length = 0;
    simde__m128i vNegLimit;
    simde__m128i vPosLimit;
    simde__m128i vSaturationCheckMin;
    simde__m128i vSaturationCheckMax;
    parasail_result_t *result = NULL;

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
    vProfile  = (simde__m128i*)profile->profile32.score;
    vProfileM = (simde__m128i*)profile->profile32.matches;
    vProfileS = (simde__m128i*)profile->profile32.similar;
    vGapO = simde_mm_set1_epi32(open);
    vGapE = simde_mm_set1_epi32(gap);
    NEG_LIMIT = (-open < matrix->min ? INT32_MIN + open : INT32_MIN - matrix->min) + 1;
    POS_LIMIT = INT32_MAX - matrix->max - 1;
    vZero = simde_mm_setzero_si128();
    vOne = simde_mm_set1_epi32(1);
    score = NEG_LIMIT;
    matches = 0;
    similar = 0;
    length = 0;
    vNegLimit = simde_mm_set1_epi32(NEG_LIMIT);
    vPosLimit = simde_mm_set1_epi32(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;

    /* initialize result */
#ifdef PARASAIL_TABLE
    result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    result = parasail_result_new_stats();
#endif
#endif
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_32 | PARASAIL_FLAG_LANES_4;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    /* initialize heap variables */
    pvHStore  = parasail_memalign_simde__m128i(16, segLen);
    pvHLoad   = parasail_memalign_simde__m128i(16, segLen);
    pvHMStore = parasail_memalign_simde__m128i(16, segLen);
    pvHMLoad  = parasail_memalign_simde__m128i(16, segLen);
    pvHSStore = parasail_memalign_simde__m128i(16, segLen);
    pvHSLoad  = parasail_memalign_simde__m128i(16, segLen);
    pvHLStore = parasail_memalign_simde__m128i(16, segLen);
    pvHLLoad  = parasail_memalign_simde__m128i(16, segLen);
    pvE       = parasail_memalign_simde__m128i(16, segLen);
    pvEM      = parasail_memalign_simde__m128i(16, segLen);
    pvES      = parasail_memalign_simde__m128i(16, segLen);
    pvEL      = parasail_memalign_simde__m128i(16, segLen);
    boundary  = parasail_memalign_int32_t(16, s2Len+1);

    /* validate heap variables */
    if (!pvHStore) return NULL;
    if (!pvHLoad) return NULL;
    if (!pvHMStore) return NULL;
    if (!pvHMLoad) return NULL;
    if (!pvHSStore) return NULL;
    if (!pvHSLoad) return NULL;
    if (!pvHLStore) return NULL;
    if (!pvHLLoad) return NULL;
    if (!pvE) return NULL;
    if (!pvEM) return NULL;
    if (!pvES) return NULL;
    if (!pvEL) return NULL;
    if (!boundary) return NULL;

    parasail_memset_simde__m128i(pvHMStore, vZero, segLen);
    parasail_memset_simde__m128i(pvHSStore, vZero, segLen);
    parasail_memset_simde__m128i(pvHLStore, vZero, segLen);
    parasail_memset_simde__m128i(pvEM, vZero, segLen);
    parasail_memset_simde__m128i(pvES, vZero, segLen);
    parasail_memset_simde__m128i(pvEL, vOne, segLen);

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            int32_t segNum = 0;
            simde__m128i_private h_;
            simde__m128i_private e_;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h_.i32[segNum] = tmp < INT32_MIN ? INT32_MIN : tmp;
                tmp = tmp - open;
                e_.i32[segNum] = tmp < INT32_MIN ? INT32_MIN : tmp;
            }
            simde_mm_store_si128(&pvHStore[index], simde__m128i_from_private(h_));
            simde_mm_store_si128(&pvE[index], simde__m128i_from_private(e_));
            ++index;
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            int64_t tmp = -open-gap*(i-1);
            boundary[i] = tmp < INT32_MIN ? INT32_MIN : tmp;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        simde__m128i vEF_opn;
        simde__m128i vE;
        simde__m128i vE_ext;
        simde__m128i vEM;
        simde__m128i vES;
        simde__m128i vEL;
        simde__m128i vF;
        simde__m128i vF_ext;
        simde__m128i vFM;
        simde__m128i vFS;
        simde__m128i vFL;
        simde__m128i vH;
        simde__m128i vH_dag;
        simde__m128i vHM;
        simde__m128i vHS;
        simde__m128i vHL;
        const simde__m128i* vP = NULL;
        const simde__m128i* vPM = NULL;
        const simde__m128i* vPS = NULL;

        /* Initialize F value to neg inf.  Any errors to vH values will
         * be corrected in the Lazy_F loop.  */
        vF = vNegLimit;
        vFM = vZero;
        vFS = vZero;
        vFL = vOne;

        /* load final segment of pvHStore and shift left by 4 bytes */
        vH = simde_mm_load_si128(&pvHStore[segLen - 1]);
        vHM = simde_mm_load_si128(&pvHMStore[segLen - 1]);
        vHS = simde_mm_load_si128(&pvHSStore[segLen - 1]);
        vHL = simde_mm_load_si128(&pvHLStore[segLen - 1]);
        vH = simde_mm_slli_si128(vH, 4);
        vHM = simde_mm_slli_si128(vHM, 4);
        vHS = simde_mm_slli_si128(vHS, 4);
        vHL = simde_mm_slli_si128(vHL, 4);

        /* insert upper boundary condition */
        vH = simde_mm_insert_epi32(vH, boundary[j], 0);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;
        vPM = vProfileM + matrix->mapper[(unsigned char)s2[j]] * segLen;
        vPS = vProfileS + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        SWAP(pvHLoad,  pvHStore)
        SWAP(pvHMLoad, pvHMStore)
        SWAP(pvHSLoad, pvHSStore)
        SWAP(pvHLLoad, pvHLStore)

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            simde__m128i case1;
            simde__m128i case2;

            vE = simde_mm_load_si128(pvE+ i);
            vEM = simde_mm_load_si128(pvEM+ i);
            vES = simde_mm_load_si128(pvES+ i);
            vEL = simde_mm_load_si128(pvEL+ i);

            /* Get max from vH, vE and vF. */
            vH_dag = simde_mm_add_epi32(vH, simde_mm_load_si128(vP + i));
            vH = simde_mm_max_epi32(vH_dag, vE);
            vH = simde_mm_max_epi32(vH, vF);
            /* Save vH values. */
            simde_mm_store_si128(pvHStore + i, vH);

            case1 = simde_mm_cmpeq_epi32(vH, vH_dag);
            case2 = simde_mm_cmpeq_epi32(vH, vF);

            /* calculate vM */
            vHM = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vEM, vFM, case2),
                    simde_mm_add_epi32(vHM, simde_mm_load_si128(vPM + i)),
                    case1);
            simde_mm_store_si128(pvHMStore + i, vHM);

            /* calculate vS */
            vHS = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vES, vFS, case2),
                    simde_mm_add_epi32(vHS, simde_mm_load_si128(vPS + i)),
                    case1);
            simde_mm_store_si128(pvHSStore + i, vHS);

            /* calculate vL */
            vHL = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vEL, vFL, case2),
                    simde_mm_add_epi32(vHL, vOne),
                    case1);
            simde_mm_store_si128(pvHLStore + i, vHL);

            vSaturationCheckMin = simde_mm_min_epi32(vSaturationCheckMin, vH);
            vSaturationCheckMax = simde_mm_max_epi32(vSaturationCheckMax, vH);
            vSaturationCheckMax = simde_mm_max_epi32(vSaturationCheckMax, vHM);
            vSaturationCheckMax = simde_mm_max_epi32(vSaturationCheckMax, vHS);
            vSaturationCheckMax = simde_mm_max_epi32(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si128(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
#endif
            vEF_opn = simde_mm_sub_epi32(vH, vGapO);

            /* Update vE value. */
            vE_ext = simde_mm_sub_epi32(vE, vGapE);
            vE = simde_mm_max_epi32(vEF_opn, vE_ext);
            case1 = simde_mm_cmpgt_epi32(vEF_opn, vE_ext);
            vEM = simde_mm_blendv_epi8(vEM, vHM, case1);
            vES = simde_mm_blendv_epi8(vES, vHS, case1);
            vEL = simde_mm_blendv_epi8(
                    simde_mm_add_epi32(vEL, vOne),
                    simde_mm_add_epi32(vHL, vOne),
                    case1);
            simde_mm_store_si128(pvE + i, vE);
            simde_mm_store_si128(pvEM + i, vEM);
            simde_mm_store_si128(pvES + i, vES);
            simde_mm_store_si128(pvEL + i, vEL);

            /* Update vF value. */
            vF_ext = simde_mm_sub_epi32(vF, vGapE);
            vF = simde_mm_max_epi32(vEF_opn, vF_ext);
            case1 = simde_mm_cmpgt_epi32(vEF_opn, vF_ext);
            vFM = simde_mm_blendv_epi8(vFM, vHM, case1);
            vFS = simde_mm_blendv_epi8(vFS, vHS, case1);
            vFL = simde_mm_blendv_epi8(
                    simde_mm_add_epi32(vFL, vOne),
                    simde_mm_add_epi32(vHL, vOne),
                    case1);

            /* Load the next vH. */
            vH = simde_mm_load_si128(pvHLoad + i);
            vHM = simde_mm_load_si128(pvHMLoad + i);
            vHS = simde_mm_load_si128(pvHSLoad + i);
            vHL = simde_mm_load_si128(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            simde__m128i vHp = simde_mm_load_si128(&pvHLoad[segLen - 1]);
            int64_t tmp = boundary[j+1]-open;
            int32_t tmp2 = tmp < INT32_MIN ? INT32_MIN : tmp;
            vHp = simde_mm_slli_si128(vHp, 4);
            vF = simde_mm_slli_si128(vF, 4);
            vFM = simde_mm_slli_si128(vFM, 4);
            vFS = simde_mm_slli_si128(vFS, 4);
            vFL = simde_mm_slli_si128(vFL, 4);
            vHp = simde_mm_insert_epi32(vHp, boundary[j], 0);
            vF = simde_mm_insert_epi32(vF, tmp2, 0);
            vFL = simde_mm_insert_epi32(vFL, 1, 0);
            for (i=0; i<segLen; ++i) {
                simde__m128i case1;
                simde__m128i case2;
                simde__m128i cond;

                vHp = simde_mm_add_epi32(vHp, simde_mm_load_si128(vP + i));
                vH = simde_mm_load_si128(pvHStore + i);
                vH = simde_mm_max_epi32(vH,vF);
                simde_mm_store_si128(pvHStore + i, vH);
                case1 = simde_mm_cmpeq_epi32(vH, vHp);
                case2 = simde_mm_cmpeq_epi32(vH, vF);
                cond = simde_mm_andnot_si128(case1, case2);

                /* calculate vM */
                vHM = simde_mm_load_si128(pvHMStore + i);
                vHM = simde_mm_blendv_epi8(vHM, vFM, cond);
                simde_mm_store_si128(pvHMStore + i, vHM);

                /* calculate vS */
                vHS = simde_mm_load_si128(pvHSStore + i);
                vHS = simde_mm_blendv_epi8(vHS, vFS, cond);
                simde_mm_store_si128(pvHSStore + i, vHS);

                /* calculate vL */
                vHL = simde_mm_load_si128(pvHLStore + i);
                vHL = simde_mm_blendv_epi8(vHL, vFL, cond);
                simde_mm_store_si128(pvHLStore + i, vHL);

                vSaturationCheckMin = simde_mm_min_epi32(vSaturationCheckMin, vH);
                vSaturationCheckMax = simde_mm_max_epi32(vSaturationCheckMax, vH);
                vSaturationCheckMax = simde_mm_max_epi32(vSaturationCheckMax, vHM);
                vSaturationCheckMax = simde_mm_max_epi32(vSaturationCheckMax, vHS);
                vSaturationCheckMax = simde_mm_max_epi32(vSaturationCheckMax, vHL);
#ifdef PARASAIL_TABLE
                arr_store_si128(result->stats->tables->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si128(result->stats->tables->similar_table, vHS, i, segLen, j, s2Len);
                arr_store_si128(result->stats->tables->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si128(result->stats->tables->score_table, vH, i, segLen, j, s2Len);
#endif
                /* Update vF value. */
                vEF_opn = simde_mm_sub_epi32(vH, vGapO);
                vF_ext = simde_mm_sub_epi32(vF, vGapE);
                if (! simde_mm_movemask_epi8(
                            simde_mm_or_si128(
                                simde_mm_cmpgt_epi32(vF_ext, vEF_opn),
                                simde_mm_cmpeq_epi32(vF_ext, vEF_opn))))
                    goto end;
                /*vF = simde_mm_max_epi32(vEF_opn, vF_ext);*/
                vF = vF_ext;
                cond = simde_mm_cmpgt_epi32(vEF_opn, vF_ext);
                vFM = simde_mm_blendv_epi8(vFM, vHM, cond);
                vFS = simde_mm_blendv_epi8(vFS, vHS, cond);
                vFL = simde_mm_blendv_epi8(
                        simde_mm_add_epi32(vFL, vOne),
                        simde_mm_add_epi32(vHL, vOne),
                        cond);
                vHp = simde_mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = simde_mm_load_si128(pvHStore + offset);
            vHM = simde_mm_load_si128(pvHMStore + offset);
            vHS = simde_mm_load_si128(pvHSStore + offset);
            vHL = simde_mm_load_si128(pvHLStore + offset);
            for (k=0; k<position; ++k) {
                vH = simde_mm_slli_si128 (vH, 4);
                vHM = simde_mm_slli_si128 (vHM, 4);
                vHS = simde_mm_slli_si128 (vHS, 4);
                vHL = simde_mm_slli_si128 (vHL, 4);
            }
            result->stats->rowcols->score_row[j] = (int32_t) simde_mm_extract_epi32 (vH, 3);
            result->stats->rowcols->matches_row[j] = (int32_t) simde_mm_extract_epi32 (vHM, 3);
            result->stats->rowcols->similar_row[j] = (int32_t) simde_mm_extract_epi32 (vHS, 3);
            result->stats->rowcols->length_row[j] = (int32_t) simde_mm_extract_epi32 (vHL, 3);
        }
#endif
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        simde__m128i vH = simde_mm_load_si128(pvHStore+i);
        simde__m128i vHM = simde_mm_load_si128(pvHMStore+i);
        simde__m128i vHS = simde_mm_load_si128(pvHSStore+i);
        simde__m128i vHL = simde_mm_load_si128(pvHLStore+i);
        arr_store_col(result->stats->rowcols->score_col, vH, i, segLen);
        arr_store_col(result->stats->rowcols->matches_col, vHM, i, segLen);
        arr_store_col(result->stats->rowcols->similar_col, vHS, i, segLen);
        arr_store_col(result->stats->rowcols->length_col, vHL, i, segLen);
    }
#endif

    /* extract last value from the last column */
    {
        simde__m128i vH = simde_mm_load_si128(pvHStore + offset);
        simde__m128i vHM = simde_mm_load_si128(pvHMStore + offset);
        simde__m128i vHS = simde_mm_load_si128(pvHSStore + offset);
        simde__m128i vHL = simde_mm_load_si128(pvHLStore + offset);
        for (k=0; k<position; ++k) {
            vH = simde_mm_slli_si128 (vH, 4);
            vHM = simde_mm_slli_si128 (vHM, 4);
            vHS = simde_mm_slli_si128 (vHS, 4);
            vHL = simde_mm_slli_si128 (vHL, 4);
        }
        score = (int32_t) simde_mm_extract_epi32 (vH, 3);
        matches = (int32_t) simde_mm_extract_epi32 (vHM, 3);
        similar = (int32_t) simde_mm_extract_epi32 (vHS, 3);
        length = (int32_t) simde_mm_extract_epi32 (vHL, 3);
    }

    if (simde_mm_movemask_epi8(simde_mm_or_si128(
            simde_mm_cmplt_epi32(vSaturationCheckMin, vNegLimit),
            simde_mm_cmpgt_epi32(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = 0;
        matches = 0;
        similar = 0;
        length = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->stats->matches = matches;
    result->stats->similar = similar;
    result->stats->length = length;

    parasail_free(boundary);
    parasail_free(pvEL);
    parasail_free(pvES);
    parasail_free(pvEM);
    parasail_free(pvE);
    parasail_free(pvHLLoad);
    parasail_free(pvHLStore);
    parasail_free(pvHSLoad);
    parasail_free(pvHSStore);
    parasail_free(pvHMLoad);
    parasail_free(pvHMStore);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}


