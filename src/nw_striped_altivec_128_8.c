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



#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        vec128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[1LL*( 0*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  0);
    array[1LL*( 1*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  1);
    array[1LL*( 2*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  2);
    array[1LL*( 3*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  3);
    array[1LL*( 4*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  4);
    array[1LL*( 5*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  5);
    array[1LL*( 6*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  6);
    array[1LL*( 7*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  7);
    array[1LL*( 8*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  8);
    array[1LL*( 9*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH,  9);
    array[1LL*(10*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 10);
    array[1LL*(11*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 11);
    array[1LL*(12*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 12);
    array[1LL*(13*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 13);
    array[1LL*(14*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 14);
    array[1LL*(15*seglen+t)*dlen + d] = (int8_t)_mm_extract_epi8(vH, 15);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        vec128i vH,
        int32_t t,
        int32_t seglen)
{
    col[ 0*seglen+t] = (int8_t)_mm_extract_epi8(vH,  0);
    col[ 1*seglen+t] = (int8_t)_mm_extract_epi8(vH,  1);
    col[ 2*seglen+t] = (int8_t)_mm_extract_epi8(vH,  2);
    col[ 3*seglen+t] = (int8_t)_mm_extract_epi8(vH,  3);
    col[ 4*seglen+t] = (int8_t)_mm_extract_epi8(vH,  4);
    col[ 5*seglen+t] = (int8_t)_mm_extract_epi8(vH,  5);
    col[ 6*seglen+t] = (int8_t)_mm_extract_epi8(vH,  6);
    col[ 7*seglen+t] = (int8_t)_mm_extract_epi8(vH,  7);
    col[ 8*seglen+t] = (int8_t)_mm_extract_epi8(vH,  8);
    col[ 9*seglen+t] = (int8_t)_mm_extract_epi8(vH,  9);
    col[10*seglen+t] = (int8_t)_mm_extract_epi8(vH, 10);
    col[11*seglen+t] = (int8_t)_mm_extract_epi8(vH, 11);
    col[12*seglen+t] = (int8_t)_mm_extract_epi8(vH, 12);
    col[13*seglen+t] = (int8_t)_mm_extract_epi8(vH, 13);
    col[14*seglen+t] = (int8_t)_mm_extract_epi8(vH, 14);
    col[15*seglen+t] = (int8_t)_mm_extract_epi8(vH, 15);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_nw_table_striped_altivec_128_8
#define PNAME parasail_nw_table_striped_profile_altivec_128_8
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_rowcol_striped_altivec_128_8
#define PNAME parasail_nw_rowcol_striped_profile_altivec_128_8
#else
#define FNAME parasail_nw_striped_altivec_128_8
#define PNAME parasail_nw_striped_profile_altivec_128_8
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
    profile = parasail_profile_create_altivec_128_8(s1, s1Len, matrix);
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
    vec128i* restrict vProfile = NULL;
    vec128i* restrict pvHStore = NULL;
    vec128i* restrict pvHLoad = NULL;
    vec128i* restrict pvE = NULL;
    int8_t* restrict boundary = NULL;
    vec128i vGapO;
    vec128i vGapE;
    int8_t NEG_LIMIT = 0;
    int8_t POS_LIMIT = 0;
    int8_t score = 0;
    vec128i vNegLimit;
    vec128i vPosLimit;
    vec128i vSaturationCheckMin;
    vec128i vSaturationCheckMax;
    parasail_result_t *result = NULL;

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
    vProfile = (vec128i*)profile->profile8.score;
    vGapO = _mm_set1_epi8(open);
    vGapE = _mm_set1_epi8(gap);
    NEG_LIMIT = (-open < matrix->min ? INT8_MIN + open : INT8_MIN - matrix->min) + 1;
    POS_LIMIT = INT8_MAX - matrix->max - 1;
    score = NEG_LIMIT;
    vNegLimit = _mm_set1_epi8(NEG_LIMIT);
    vPosLimit = _mm_set1_epi8(POS_LIMIT);
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
        | PARASAIL_FLAG_BITS_8 | PARASAIL_FLAG_LANES_16;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    /* initialize heap variables */
    pvHStore = parasail_memalign_vec128i(16, segLen);
    pvHLoad  = parasail_memalign_vec128i(16, segLen);
    pvE      = parasail_memalign_vec128i(16, segLen);
    boundary = parasail_memalign_int8_t(16, s2Len+1);

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
            vec128i_8_t h;
            vec128i_8_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.v[segNum] = tmp < INT8_MIN ? INT8_MIN : tmp;
                tmp = tmp - open;
                e.v[segNum] = tmp < INT8_MIN ? INT8_MIN : tmp;
            }
            _mm_store_si128(&pvHStore[index], h.m);
            _mm_store_si128(&pvE[index], e.m);
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

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        vec128i vE;
        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        vec128i vF = vNegLimit;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vec128i vH = _mm_slli_si128(pvHStore[segLen - 1], 1);

        /* Correct part of the vProfile */
        const vec128i* vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        vec128i* pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* insert upper boundary condition */
        vH = _mm_insert_epi8(vH, boundary[j], 0);

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vH = _mm_adds_epi8(vH, _mm_load_si128(vP + i));
            vE = _mm_load_si128(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = _mm_max_epi8(vH, vE);
            vH = _mm_max_epi8(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);
            vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vH);
            vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vH);
            vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vE);
            vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vF);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->tables->score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = _mm_subs_epi8(vH, vGapO);
            vE = _mm_subs_epi8(vE, vGapE);
            vE = _mm_max_epi8(vE, vH);
            _mm_store_si128(pvE + i, vE);

            /* Update vF value. */
            vF = _mm_subs_epi8(vF, vGapE);
            vF = _mm_max_epi8(vF, vH);

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            int64_t tmp = boundary[j+1]-open;
            int8_t tmp2 = tmp < INT8_MIN ? INT8_MIN : tmp;
            vF = _mm_slli_si128(vF, 1);
            vF = _mm_insert_epi8(vF, tmp2, 0);
            for (i=0; i<segLen; ++i) {
                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi8(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
                vSaturationCheckMin = _mm_min_epi8(vSaturationCheckMin, vH);
                vSaturationCheckMax = _mm_max_epi8(vSaturationCheckMax, vH);
#ifdef PARASAIL_TABLE
                arr_store_si128(result->tables->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm_subs_epi8(vH, vGapO);
                vF = _mm_subs_epi8(vF, vGapE);
                if (! _mm_movemask_epi8(_mm_cmpgt_epi8(vF, vH))) goto end;
                /*vF = _mm_max_epi8(vF, vH);*/
            }
        }
end:
        {
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = _mm_load_si128(pvHStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128(vH, 1);
            }
            result->rowcols->score_row[j] = (int8_t) _mm_extract_epi8 (vH, 15);
        }
#endif
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        vec128i vH = _mm_load_si128(pvHStore+i);
        arr_store_col(result->rowcols->score_col, vH, i, segLen);
    }
#endif

    /* extract last value from the last column */
    {
        vec128i vH = _mm_load_si128(pvHStore + offset);
        for (k=0; k<position; ++k) {
            vH = _mm_slli_si128 (vH, 1);
        }
        score = (int8_t) _mm_extract_epi8 (vH, 15);
    }

    if (_mm_movemask_epi8(_mm_or_si128(
            _mm_cmplt_epi8(vSaturationCheckMin, vNegLimit),
            _mm_cmpgt_epi8(vSaturationCheckMax, vPosLimit)))) {
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

