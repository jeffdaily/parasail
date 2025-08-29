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
#define FNAME parasail_nw_table_scan_neon_128_64
#define PNAME parasail_nw_table_scan_profile_neon_128_64
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_rowcol_scan_neon_128_64
#define PNAME parasail_nw_rowcol_scan_profile_neon_128_64
#else
#define FNAME parasail_nw_scan_neon_128_64
#define PNAME parasail_nw_scan_profile_neon_128_64
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
    int64_t* restrict boundary = NULL;
    simde__m128i* restrict pvHt = NULL;
    simde__m128i* restrict pvH = NULL;
    simde__m128i* restrict pvGapper = NULL;
    simde__m128i vGapO;
    simde__m128i vGapE;
    int64_t NEG_LIMIT = 0;
    int64_t POS_LIMIT = 0;
    simde__m128i vZero;
    int64_t score = 0;
    simde__m128i vNegLimit;
    simde__m128i vPosLimit;
    simde__m128i vSaturationCheckMin;
    simde__m128i vSaturationCheckMax;
    simde__m128i vNegInfFront;
    simde__m128i vSegLenXgap;
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
    pvP = (simde__m128i*)profile->profile64.score;
    vGapO = simde_mm_set1_epi64x(open);
    vGapE = simde_mm_set1_epi64x(gap);
    NEG_LIMIT = (-open < matrix->min ? INT64_MIN + open : INT64_MIN - matrix->min) + 1;
    POS_LIMIT = INT64_MAX - matrix->max - 1;
    vZero = simde_mm_setzero_si128();
    score = NEG_LIMIT;
    vNegLimit = simde_mm_set1_epi64x(NEG_LIMIT);
    vPosLimit = simde_mm_set1_epi64x(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;
    vNegInfFront = vZero;
    vNegInfFront = simde_mm_insert_epi64(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = simde_mm_add_epi64(vNegInfFront,
            simde_mm_slli_si128(simde_mm_set1_epi64x(-segLen*gap), 8));

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
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_SCAN
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_2;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    /* initialize heap variables */
    pvE = parasail_memalign_simde__m128i(16, segLen);
    boundary = parasail_memalign_int64_t(16, s2Len+1);
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
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h_.i64[segNum] = tmp < INT64_MIN ? INT64_MIN : tmp;
                tmp = tmp - open;
                e_.i64[segNum] = tmp < INT64_MIN ? INT64_MIN : tmp;
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
            int64_t tmp = -open-gap*(i-1);
            boundary[i] = tmp < INT64_MIN ? INT64_MIN : tmp;
        }
    }

    {
        simde__m128i vGapper = simde_mm_sub_epi64(vZero,vGapO);
        for (i=segLen-1; i>=0; --i) {
            simde_mm_store_si128(pvGapper+i, vGapper);
            vGapper = simde_mm_sub_epi64(vGapper, vGapE);
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        simde__m128i vE;
        simde__m128i vHt;
        simde__m128i vF;
        simde__m128i vH;
        simde__m128i vHp;
        simde__m128i *pvW;
        simde__m128i vW;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = simde_mm_load_si128(pvH+(segLen-1));
        vHp = simde_mm_slli_si128(vHp, 8);
        vHp = simde_mm_insert_epi64(vHp, boundary[j], 0);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = simde_mm_sub_epi64(vNegLimit, pvGapper[0]);
        vF = vNegLimit;
        for (i=0; i<segLen; ++i) {
            vH = simde_mm_load_si128(pvH+i);
            vE = simde_mm_load_si128(pvE+i);
            vW = simde_mm_load_si128(pvW+i);
            vE = simde_mm_max_epi64(
                    simde_mm_sub_epi64(vE, vGapE),
                    simde_mm_sub_epi64(vH, vGapO));
            vHp = simde_mm_add_epi64(vHp, vW);
            vF = simde_mm_max_epi64(vF, simde_mm_add_epi64(vHt, pvGapper[i]));
            vHt = simde_mm_max_epi64(vE, vHp);
            simde_mm_store_si128(pvE+i, vE);
            simde_mm_store_si128(pvHt+i, vHt);
            vHp = vH;
        }

        /* pseudo prefix scan on F and H */
        vHt = simde_mm_slli_si128(vHt, 8);
        vHt = simde_mm_insert_epi64(vHt, boundary[j+1], 0);
        vF = simde_mm_max_epi64(vF, simde_mm_add_epi64(vHt, pvGapper[0]));
        for (i=0; i<segWidth-2; ++i) {
            simde__m128i vFt = simde_mm_slli_si128(vF, 8);
            vFt = simde_mm_add_epi64(vFt, vSegLenXgap);
            vF = simde_mm_max_epi64(vF, vFt);
        }

        /* calculate final H */
        vF = simde_mm_slli_si128(vF, 8);
        vF = simde_mm_add_epi64(vF, vNegInfFront);
        vH = simde_mm_max_epi64(vHt, vF);
        for (i=0; i<segLen; ++i) {
            vHt = simde_mm_load_si128(pvHt+i);
            vF = simde_mm_max_epi64(
                    simde_mm_sub_epi64(vF, vGapE),
                    simde_mm_sub_epi64(vH, vGapO));
            vH = simde_mm_max_epi64(vHt, vF);
            simde_mm_store_si128(pvH+i, vH);
            vSaturationCheckMin = simde_mm_min_epi64(vSaturationCheckMin, vH);
            vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vH);
#ifdef PARASAIL_TABLE
            arr_store_si128(result->tables->score_table, vH, i, segLen, j, s2Len);
#endif
        } 


#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = simde_mm_load_si128(pvH + offset);
            for (k=0; k<position; ++k) {
                vH = simde_mm_slli_si128(vH, 8);
            }
            result->rowcols->score_row[j] = (int64_t) simde_mm_extract_epi64 (vH, 1);
        }
#endif
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        simde__m128i vH = simde_mm_load_si128(pvH+i);
        arr_store_col(result->rowcols->score_col, vH, i, segLen);
    }
#endif

    /* extract last value from the last column */
    {
        simde__m128i vH = simde_mm_load_si128(pvH + offset);
        for (k=0; k<position; ++k) {
            vH = simde_mm_slli_si128(vH, 8);
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

    parasail_free(pvGapper);
    parasail_free(pvH);
    parasail_free(pvHt);
    parasail_free(boundary);
    parasail_free(pvE);

    return result;
}


