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



#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        simde__m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[1LL*(0*seglen+t)*dlen + d] = (int64_t)simde_mm_extract_epi64(vH, 0);
    array[1LL*(1*seglen+t)*dlen + d] = (int64_t)simde_mm_extract_epi64(vH, 1);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        simde__m128i vH,
        int32_t t,
        int32_t seglen)
{
    col[0*seglen+t] = (int64_t)simde_mm_extract_epi64(vH, 0);
    col[1*seglen+t] = (int64_t)simde_mm_extract_epi64(vH, 1);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_nw_table_striped_neon_128_64
#define PNAME parasail_nw_table_striped_profile_neon_128_64
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_rowcol_striped_neon_128_64
#define PNAME parasail_nw_rowcol_striped_profile_neon_128_64
#else
#define FNAME parasail_nw_striped_neon_128_64
#define PNAME parasail_nw_striped_profile_neon_128_64
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
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        PARASAIL_CHECK_NULL(s1);
        PARASAIL_CHECK_GT0(s1Len);
    }

    /* initialize local variables */
    profile = parasail_profile_create_neon_128_64(s1, s1Len, matrix);
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
    int64_t* restrict boundary = NULL;
    simde__m128i vGapO;
    simde__m128i vGapE;
    int64_t NEG_LIMIT = 0;
    int64_t POS_LIMIT = 0;
    int64_t score = 0;
    simde__m128i vNegLimit;
    simde__m128i vPosLimit;
    simde__m128i vSaturationCheckMin;
    simde__m128i vSaturationCheckMax;
    parasail_result_t *result = NULL;

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
    vProfile = (simde__m128i*)profile->profile64.score;
    vGapO = simde_mm_set1_epi64x(open);
    vGapE = simde_mm_set1_epi64x(gap);
    NEG_LIMIT = (-open < matrix->min ? INT64_MIN + open : INT64_MIN - matrix->min) + 1;
    POS_LIMIT = INT64_MAX - matrix->max - 1;
    score = NEG_LIMIT;
    vNegLimit = simde_mm_set1_epi64x(NEG_LIMIT);
    vPosLimit = simde_mm_set1_epi64x(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;

    /* initialize result */
#ifdef PARASAIL_TABLE
    result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    result = parasail_result_new_rowcol1(segLen*segWidth, s2Len);
#else
    result = parasail_result_new();
#endif
#endif
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_2;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    /* initialize heap variables */
    pvHStore = parasail_memalign_simde__m128i(16, segLen);
    pvHLoad  = parasail_memalign_simde__m128i(16, segLen);
    pvE      = parasail_memalign_simde__m128i(16, segLen);
    boundary = parasail_memalign_int64_t(16, s2Len+1);

    /* validate heap variables */
    if (!pvHStore) return NULL;
    if (!pvHLoad) return NULL;
    if (!pvE) return NULL;
    if (!boundary) return NULL;

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            int32_t segNum = 0;
            simde__m128i_private h_;
            simde__m128i_private e_;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h_.i64[segNum] = tmp < INT64_MIN ? INT64_MIN : tmp;
                tmp = tmp - open;
                e_.i64[segNum] = tmp < INT64_MIN ? INT64_MIN : tmp;
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
            boundary[i] = tmp < INT64_MIN ? INT64_MIN : tmp;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        simde__m128i vE;
        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        simde__m128i vF = vNegLimit;

        /* load final segment of pvHStore and shift left by 2 bytes */
        simde__m128i vH = simde_mm_slli_si128(pvHStore[segLen - 1], 8);

        /* Correct part of the vProfile */
        const simde__m128i* vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        simde__m128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* insert upper boundary condition */
        vH = simde_mm_insert_epi64(vH, boundary[j], 0);

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vH = simde_mm_add_epi64(vH, simde_mm_load_si128(vP + i));
            vE = simde_mm_load_si128(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = simde_mm_max_epi64(vH, vE);
            vH = simde_mm_max_epi64(vH, vF);
            /* Save vH values. */
            simde_mm_store_si128(pvHStore + i, vH);
            vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vH);
            vSaturationCheckMin = simde_mm_min_epi64(vSaturationCheckMin, vH);
            vSaturationCheckMin = simde_mm_min_epi64(vSaturationCheckMin, vE);
            vSaturationCheckMin = simde_mm_min_epi64(vSaturationCheckMin, vF);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->tables->score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = simde_mm_sub_epi64(vH, vGapO);
            vE = simde_mm_sub_epi64(vE, vGapE);
            vE = simde_mm_max_epi64(vE, vH);
            simde_mm_store_si128(pvE + i, vE);

            /* Update vF value. */
            vF = simde_mm_sub_epi64(vF, vGapE);
            vF = simde_mm_max_epi64(vF, vH);

            /* Load the next vH. */
            vH = simde_mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            int64_t tmp = boundary[j+1]-open;
            int64_t tmp2 = tmp < INT64_MIN ? INT64_MIN : tmp;
            vF = simde_mm_slli_si128(vF, 8);
            vF = simde_mm_insert_epi64(vF, tmp2, 0);
            for (i=0; i<segLen; ++i) {
                vH = simde_mm_load_si128(pvHStore + i);
                vH = simde_mm_max_epi64(vH,vF);
                simde_mm_store_si128(pvHStore + i, vH);
                vSaturationCheckMin = simde_mm_min_epi64(vSaturationCheckMin, vH);
                vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vH);
#ifdef PARASAIL_TABLE
                arr_store_si128(result->tables->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = simde_mm_sub_epi64(vH, vGapO);
                vF = simde_mm_sub_epi64(vF, vGapE);
                if (! simde_mm_movemask_epi8(simde_mm_cmpgt_epi64(vF, vH))) goto end;
                /*vF = simde_mm_max_epi64(vF, vH);*/
            }
        }
end:
        {
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = simde_mm_load_si128(pvHStore + offset);
            for (k=0; k<position; ++k) {
                vH = simde_mm_slli_si128(vH, 8);
            }
            result->rowcols->score_row[j] = (int64_t) simde_mm_extract_epi64 (vH, 1);
        }
#endif
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        simde__m128i vH = simde_mm_load_si128(pvHStore+i);
        arr_store_col(result->rowcols->score_col, vH, i, segLen);
    }
#endif

    /* extract last value from the last column */
    {
        simde__m128i vH = simde_mm_load_si128(pvHStore + offset);
        for (k=0; k<position; ++k) {
            vH = simde_mm_slli_si128 (vH, 8);
        }
        score = (int64_t) simde_mm_extract_epi64 (vH, 1);
    }

    if (simde_mm_movemask_epi8(simde_mm_or_si128(
            simde_mm_cmplt_epi64(vSaturationCheckMin, vNegLimit),
            simde_mm_cmpgt_epi64(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(boundary);
    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

