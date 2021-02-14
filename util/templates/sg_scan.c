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

%(HEADER)s

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_%(ISA)s.h"

#define SG_SUFFIX %(SUFFIX)s
#define SG_SUFFIX_PROF %(SUFFIX_PROF)s
#include "sg_helper.h"

%(FIXES)s

#ifdef PARASAIL_TABLE
static inline void arr_store_si%(BITS)s(
        int *array,
        %(VTYPE)s vH,
        %(INDEX)s t,
        %(INDEX)s seglen,
        %(INDEX)s d,
        %(INDEX)s dlen)
{
%(PRINTER)s
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        %(VTYPE)s vH,
        %(INDEX)s t,
        %(INDEX)s seglen)
{
%(PRINTER_ROWCOL)s
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME %(NAME_TABLE)s
#define PNAME %(PNAME_TABLE)s
#else
#ifdef PARASAIL_ROWCOL
#define FNAME %(NAME_ROWCOL)s
#define PNAME %(PNAME_ROWCOL)s
#else
#define FNAME %(NAME)s
#define PNAME %(PNAME)s
#endif
#endif

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
    profile = parasail_profile_create_%(ISA)s_%(BITS)s_%(WIDTH)s(s1, s1Len, matrix);
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
    %(INDEX)s i = 0;
    %(INDEX)s j = 0;
    %(INDEX)s k = 0;
    %(INDEX)s end_query = 0;
    %(INDEX)s end_ref = 0;
    %(INDEX)s s1Len = 0;
    const parasail_matrix_t *matrix = NULL;
    %(INDEX)s segWidth = 0;
    %(INDEX)s segLen = 0;
    %(INDEX)s offset = 0;
    %(INDEX)s position = 0;
    %(VTYPE)s* restrict pvP = NULL;
    %(VTYPE)s* restrict pvE = NULL;
    %(INT)s* restrict boundary = NULL;
    %(VTYPE)s* restrict pvHt = NULL;
    %(VTYPE)s* restrict pvH = NULL;
    %(VTYPE)s* restrict pvGapper = NULL;
    %(VTYPE)s vGapO;
    %(VTYPE)s vGapE;
    %(INT)s NEG_LIMIT = 0;
    %(INT)s POS_LIMIT = 0;
    %(VTYPE)s vZero;
    %(INT)s score = 0;
    %(VTYPE)s vNegLimit;
    %(VTYPE)s vPosLimit;
    %(VTYPE)s vSaturationCheckMin;
    %(VTYPE)s vSaturationCheckMax;
    %(VTYPE)s vMaxH;
    %(VTYPE)s vPosMask;
    %(VTYPE)s vNegInfFront;
    %(VTYPE)s vSegLenXgap;
    parasail_result_t *result = NULL;

    /* validate inputs */
    PARASAIL_CHECK_NULL(profile);
    PARASAIL_CHECK_NULL(profile->profile%(WIDTH)s.score);
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
    segWidth = %(LANES)s; /* number of values in vector unit */
    segLen = (s1Len + segWidth - 1) / segWidth;
    offset = (s1Len - 1) %% segLen;
    position = (segWidth - 1) - (s1Len - 1) / segLen;
    pvP = (%(VTYPE)s*)profile->profile%(WIDTH)s.score;
    vGapO = %(VSET1)s(open);
    vGapE = %(VSET1)s(gap);
    NEG_LIMIT = (-open < matrix->min ? INT%(WIDTH)s_MIN + open : INT%(WIDTH)s_MIN - matrix->min) + 1;
    POS_LIMIT = INT%(WIDTH)s_MAX - matrix->max - 1;
    vZero = %(VSET0)s();
    score = NEG_LIMIT;
    vNegLimit = %(VSET1)s(NEG_LIMIT);
    vPosLimit = %(VSET1)s(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;
    vMaxH = vNegLimit;
    vPosMask = %(VCMPEQ)s(%(VSET1)s(position),
            %(VSET)s(%(POSITION_MASK)s));
    vNegInfFront = vZero;
    vNegInfFront = %(VINSERT)s(vNegInfFront, NEG_LIMIT, 0);
    vSegLenXgap = %(VADD)s(vNegInfFront,
            %(VSHIFT)s(%(VSET1)s(-segLen*gap), %(BYTES)s));

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
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_SCAN
        | PARASAIL_FLAG_BITS_%(WIDTH)s | PARASAIL_FLAG_LANES_%(LANES)s;
    result->flag |= s1_beg ? PARASAIL_FLAG_SG_S1_BEG : 0;
    result->flag |= s1_end ? PARASAIL_FLAG_SG_S1_END : 0;
    result->flag |= s2_beg ? PARASAIL_FLAG_SG_S2_BEG : 0;
    result->flag |= s2_end ? PARASAIL_FLAG_SG_S2_END : 0;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    /* initialize heap variables */
    pvE = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    boundary = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+1);
    pvHt= parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    pvH = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    pvGapper = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);

    /* validate heap variables */
    if (!pvE) return NULL;
    if (!boundary) return NULL;
    if (!pvHt) return NULL;
    if (!pvH) return NULL;
    if (!pvGapper) return NULL;

%(INIT_H_AND_E)s

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            int64_t tmp = s2_beg ? 0 : (-open-gap*(i-1));
            boundary[i] = tmp < INT%(WIDTH)s_MIN ? INT%(WIDTH)s_MIN : tmp;
        }
    }

    {
        %(VTYPE)s vGapper = %(VSUB)s(vZero,vGapO);
        for (i=segLen-1; i>=0; --i) {
            %(VSTORE)s(pvGapper+i, vGapper);
            vGapper = %(VSUB)s(vGapper, vGapE);
            /* long queries and/or large penalties will break the pseudo prefix scan */
            vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vGapper);
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        %(VTYPE)s vE;
        %(VTYPE)s vHt;
        %(VTYPE)s vF;
        %(VTYPE)s vH;
        %(VTYPE)s vHp;
        %(VTYPE)s *pvW;
        %(VTYPE)s vW;

        /* calculate E */
        /* calculate Ht */
        /* calculate F and H first pass */
        vHp = %(VLOAD)s(pvH+(segLen-1));
        vHp = %(VSHIFT)s(vHp, %(BYTES)s);
        vHp = %(VINSERT)s(vHp, boundary[j], 0);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = %(VSUB)s(vNegLimit, pvGapper[0]);
        vF = vNegLimit;
        for (i=0; i<segLen; ++i) {
            vH = %(VLOAD)s(pvH+i);
            vE = %(VLOAD)s(pvE+i);
            vW = %(VLOAD)s(pvW+i);
            vE = %(VMAX)s(
                    %(VSUB)s(vE, vGapE),
                    %(VSUB)s(vH, vGapO));
            vHp = %(VADD)s(vHp, vW);
            vF = %(VMAX)s(vF, %(VADD)s(vHt, pvGapper[i]));
            vHt = %(VMAX)s(vE, vHp);
            %(VSTORE)s(pvE+i, vE);
            %(VSTORE)s(pvHt+i, vHt);
            vHp = vH;
        }

        /* pseudo prefix scan on F and H */
        vHt = %(VSHIFT)s(vHt, %(BYTES)s);
        vHt = %(VINSERT)s(vHt, boundary[j+1], 0);
        vF = %(VMAX)s(vF, %(VADD)s(vHt, pvGapper[0]));
        for (i=0; i<segWidth-2; ++i) {
            %(VTYPE)s vFt = %(VSHIFT)s(vF, %(BYTES)s);
            vFt = %(VADD)s(vFt, vSegLenXgap);
            vF = %(VMAX)s(vF, vFt);
        }

        /* calculate final H */
        vF = %(VSHIFT)s(vF, %(BYTES)s);
        vF = %(VADD)s(vF, vNegInfFront);
        vH = %(VMAX)s(vHt, vF);
        for (i=0; i<segLen; ++i) {
            vHt = %(VLOAD)s(pvHt+i);
            vF = %(VMAX)s(
                    %(VSUB)s(vF, vGapE),
                    %(VSUB)s(vH, vGapO));
            vH = %(VMAX)s(vHt, vF);
            %(VSTORE)s(pvH+i, vH);
            vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vH);
            vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vF);
            vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vH);
#ifdef PARASAIL_TABLE
            arr_store_si%(BITS)s(result->tables->score_table, vH, i, segLen, j, s2Len);
#endif
        } 

        /* extract vector containing last value from column */
        {
            %(VTYPE)s vCompare;
            vH = %(VLOAD)s(pvH + offset);
            vCompare = %(VAND)s(vPosMask, %(VCMPGT)s(vH, vMaxH));
            vMaxH = %(VMAX)s(vH, vMaxH);
            if (%(VMOVEMASK)s(vCompare)) {
                end_ref = j;
            }
#ifdef PARASAIL_ROWCOL
            for (k=0; k<position; ++k) {
                vH = %(VSHIFT)s(vH, %(BYTES)s);
            }
            result->rowcols->score_row[j] = (%(INT)s) %(VEXTRACT)s (vH, %(LAST_POS)s);
#endif
        }
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        %(VTYPE)s vH = %(VLOAD)s(pvH + i);
        arr_store_col(result->rowcols->score_col, vH, i, segLen);
    }
#endif

    /* max last value from all columns */
    if (s2_end)
    {
        for (k=0; k<position; ++k) {
            vMaxH = %(VSHIFT)s(vMaxH, %(BYTES)s);
        }
        score = (%(INT)s) %(VEXTRACT)s(vMaxH, %(LAST_POS)s);
        end_query = s1Len-1;
    }

    /* max of last column */
    if (s1_end)
    {
        /* Trace the alignment ending position on read. */
        %(INT)s *t = (%(INT)s*)pvH;
        %(INDEX)s column_len = segLen * segWidth;
        for (i = 0; i<column_len; ++i, ++t) {
            %(INDEX)s temp = i / segWidth + i %% segWidth * segLen;
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
            %(VTYPE)s vH = %(VLOAD)s(pvH + offset);
            for (k=0; k<position; ++k) {
                vH = %(VSHIFT)s(vH, %(BYTES)s);
            }
            score = (%(INT)s) %(VEXTRACT)s (vH, %(LAST_POS)s);
            end_ref = s2Len - 1;
            end_query = s1Len - 1;
        }
    }

    if (%(VMOVEMASK)s(%(VOR)s(
            %(VCMPLT)s(vSaturationCheckMin, vNegLimit),
            %(VCMPGT)s(vSaturationCheckMax, vPosLimit)))) {
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

