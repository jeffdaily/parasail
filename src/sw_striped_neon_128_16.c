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
        int32_t dlen,
        int32_t bias)
{
    array[1LL*(0*seglen+t)*dlen + d] = (int16_t)simde_mm_extract_epi16(vH, 0) - bias;
    array[1LL*(1*seglen+t)*dlen + d] = (int16_t)simde_mm_extract_epi16(vH, 1) - bias;
    array[1LL*(2*seglen+t)*dlen + d] = (int16_t)simde_mm_extract_epi16(vH, 2) - bias;
    array[1LL*(3*seglen+t)*dlen + d] = (int16_t)simde_mm_extract_epi16(vH, 3) - bias;
    array[1LL*(4*seglen+t)*dlen + d] = (int16_t)simde_mm_extract_epi16(vH, 4) - bias;
    array[1LL*(5*seglen+t)*dlen + d] = (int16_t)simde_mm_extract_epi16(vH, 5) - bias;
    array[1LL*(6*seglen+t)*dlen + d] = (int16_t)simde_mm_extract_epi16(vH, 6) - bias;
    array[1LL*(7*seglen+t)*dlen + d] = (int16_t)simde_mm_extract_epi16(vH, 7) - bias;
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        simde__m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t bias)
{
    col[0*seglen+t] = (int16_t)simde_mm_extract_epi16(vH, 0) - bias;
    col[1*seglen+t] = (int16_t)simde_mm_extract_epi16(vH, 1) - bias;
    col[2*seglen+t] = (int16_t)simde_mm_extract_epi16(vH, 2) - bias;
    col[3*seglen+t] = (int16_t)simde_mm_extract_epi16(vH, 3) - bias;
    col[4*seglen+t] = (int16_t)simde_mm_extract_epi16(vH, 4) - bias;
    col[5*seglen+t] = (int16_t)simde_mm_extract_epi16(vH, 5) - bias;
    col[6*seglen+t] = (int16_t)simde_mm_extract_epi16(vH, 6) - bias;
    col[7*seglen+t] = (int16_t)simde_mm_extract_epi16(vH, 7) - bias;
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_table_striped_neon_128_16
#define PNAME parasail_sw_table_striped_profile_neon_128_16
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_rowcol_striped_neon_128_16
#define PNAME parasail_sw_rowcol_striped_profile_neon_128_16
#else
#define FNAME parasail_sw_striped_neon_128_16
#define PNAME parasail_sw_striped_profile_neon_128_16
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
    profile = parasail_profile_create_neon_128_16(s1, s1Len, matrix);
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
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int32_t s1Len = 0;
    const parasail_matrix_t *matrix = NULL;
    int32_t segWidth = 0;
    int32_t segLen = 0;
#ifdef PARASAIL_ROWCOL
    int32_t offset = 0;
    int32_t position = 0;
#endif
    simde__m128i* restrict vProfile = NULL;
    simde__m128i* restrict pvHStore = NULL;
    simde__m128i* restrict pvHLoad = NULL;
    simde__m128i* restrict pvHMax = NULL;
    simde__m128i* restrict pvE = NULL;
    simde__m128i vGapO;
    simde__m128i vGapE;
    simde__m128i vZero;
    int16_t bias = 0;
    int16_t score = 0;
    simde__m128i vBias;
    simde__m128i vMaxH;
    simde__m128i vMaxHUnit;
    int16_t maxp = 0;
    simde__m128i insert_mask;
    /*int16_t stop = 0;*/
    parasail_result_t *result = NULL;

    /* validate inputs */
    PARASAIL_CHECK_NULL(profile);
    PARASAIL_CHECK_NULL(profile->profile16.score);
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
    end_query = 0;
    end_ref = 0;
    s1Len = profile->s1Len;
    matrix = profile->matrix;
    segWidth = 8; /* number of values in vector unit */
    segLen = (s1Len + segWidth - 1) / segWidth;
#ifdef PARASAIL_ROWCOL
    offset = (s1Len - 1) % segLen;
    position = (segWidth - 1) - (s1Len - 1) / segLen;
#endif
    vProfile = (simde__m128i*)profile->profile16.score;
    vGapO = simde_mm_set1_epi16(open);
    vGapE = simde_mm_set1_epi16(gap);
    vZero = simde_mm_set1_epi16(0);
    bias = INT16_MIN;
    score = bias;
    vBias = simde_mm_set1_epi16(bias);
    vMaxH = vBias;
    vMaxHUnit = vBias;
    maxp = INT16_MAX - (int16_t)(matrix->max+1);
    insert_mask = simde_mm_cmpgt_epi16(
            simde_mm_set_epi16(0,0,0,0,0,0,0,1),
            vZero);
    /*stop = profile->stop == INT32_MAX ?  INT16_MAX : (int16_t)profile->stop-bias;*/

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
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_BITS_16 | PARASAIL_FLAG_LANES_8;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    /* initialize heap variables */
    pvHStore = parasail_memalign_simde__m128i(16, segLen);
    pvHLoad = parasail_memalign_simde__m128i(16, segLen);
    pvHMax = parasail_memalign_simde__m128i(16, segLen);
    pvE = parasail_memalign_simde__m128i(16, segLen);

    /* validate heap variables */
    if (!pvHStore) return NULL;
    if (!pvHLoad) return NULL;
    if (!pvHMax) return NULL;
    if (!pvE) return NULL;

    /* initialize H and E */
    parasail_memset_simde__m128i(pvHStore, vBias, segLen);
    parasail_memset_simde__m128i(pvE, vBias, segLen);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        simde__m128i vE;
        simde__m128i vF;
        simde__m128i vH;
        const simde__m128i* vP = NULL;
        simde__m128i* pv = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        vF = vBias;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = simde_mm_slli_si128(pvHStore[segLen - 1], 2);
        vH = simde_mm_blendv_epi8(vH, vBias, insert_mask);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        if (end_ref == j-2) {
            /* Swap in the max buffer. */
            pv = pvHMax;
            pvHMax = pvHLoad;
            pvHLoad = pvHStore;
            pvHStore = pv;
        }
        else {
            /* Swap the 2 H buffers. */
            pv = pvHLoad;
            pvHLoad = pvHStore;
            pvHStore = pv;
        }

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vH = simde_mm_adds_epi16(vH, simde_mm_load_si128(vP + i));
            vE = simde_mm_load_si128(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = simde_mm_max_epi16(vH, vE);
            vH = simde_mm_max_epi16(vH, vF);
            /* Save vH values. */
            simde_mm_store_si128(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->tables->score_table, vH, i, segLen, j, s2Len, bias);
#endif
            vMaxH = simde_mm_max_epi16(vH, vMaxH);

            /* Update vE value. */
            vH = simde_mm_subs_epi16(vH, vGapO);
            vE = simde_mm_subs_epi16(vE, vGapE);
            vE = simde_mm_max_epi16(vE, vH);
            simde_mm_store_si128(pvE + i, vE);

            /* Update vF value. */
            vF = simde_mm_subs_epi16(vF, vGapE);
            vF = simde_mm_max_epi16(vF, vH);

            /* Load the next vH. */
            vH = simde_mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            vF = simde_mm_slli_si128(vF, 2);
            vF = simde_mm_blendv_epi8(vF, vBias, insert_mask);
            for (i=0; i<segLen; ++i) {
                vH = simde_mm_load_si128(pvHStore + i);
                vH = simde_mm_max_epi16(vH,vF);
                simde_mm_store_si128(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
                arr_store_si128(result->tables->score_table, vH, i, segLen, j, s2Len, bias);
#endif
                vMaxH = simde_mm_max_epi16(vH, vMaxH);
                vH = simde_mm_subs_epi16(vH, vGapO);
                vF = simde_mm_subs_epi16(vF, vGapE);
                if (! simde_mm_movemask_epi8(simde_mm_cmpgt_epi16(vF, vH))) goto end;
                /*vF = simde_mm_max_epi16(vF, vH);*/
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
                vH = simde_mm_slli_si128(vH, 2);
            }
            result->rowcols->score_row[j] = (int16_t) simde_mm_extract_epi16 (vH, 7) - bias;
        }
#endif

        {
            simde__m128i vCompare = simde_mm_cmpgt_epi16(vMaxH, vMaxHUnit);
            if (simde_mm_movemask_epi8(vCompare)) {
                score = simde_mm_hmax_epi16(vMaxH);
                /* if score has potential to overflow, abort early */
                if (score > maxp) {
                    result->flag |= PARASAIL_FLAG_SATURATED;
                    break;
                }
                vMaxHUnit = simde_mm_set1_epi16(score);
                end_ref = j;
            }
        }

        /*if (score == stop) break;*/
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        simde__m128i vH = simde_mm_load_si128(pvHStore+i);
        arr_store_col(result->rowcols->score_col, vH, i, segLen, bias);
    }
#endif

    if (score == INT16_MAX) {
        result->flag |= PARASAIL_FLAG_SATURATED;
    }

    if (parasail_result_is_saturated(result)) {
        score = INT16_MAX;
        end_query = 0;
        end_ref = 0;
    }
    else {
        if (end_ref == j-1) {
            /* end_ref was the last store column */
            simde__m128i *pv = pvHMax;
            pvHMax = pvHStore;
            pvHStore = pv;
        }
        else if (end_ref == j-2) {
            /* end_ref was the last load column */
            simde__m128i *pv = pvHMax;
            pvHMax = pvHLoad;
            pvHLoad = pv;
        }
        /* Trace the alignment ending position on read. */
        {
            int16_t *t = (int16_t*)pvHMax;
            int32_t column_len = segLen * segWidth;
            end_query = s1Len - 1;
            for (i = 0; i<column_len; ++i, ++t) {
                if (*t == score) {
                    int32_t temp = i / segWidth + i % segWidth * segLen;
                    if (temp < end_query) {
                        end_query = temp;
                    }
                }
            }
        }
    }

    result->score = score - bias;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(pvE);
    parasail_free(pvHMax);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

