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
#include <string.h>

%(HEADER)s

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_%(ISA)s.h"

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
    profile = parasail_profile_create_%(ISA)s_%(BITS)s_%(WIDTH)s(s1, s1Len, matrix);
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
    %(INDEX)s i = 0;
    %(INDEX)s j = 0;
    %(INDEX)s end_query = 0;
    %(INDEX)s end_ref = 0;
    %(INDEX)s s1Len = 0;
    const parasail_matrix_t *matrix = NULL;
    %(INDEX)s segWidth = 0;
    %(INDEX)s segLen = 0;
#ifdef PARASAIL_ROWCOL
    %(INDEX)s offset = 0;
    %(INDEX)s position = 0;
#endif
    %(VTYPE)s* restrict pvP = NULL;
    %(VTYPE)s* restrict pvE = NULL;
    %(VTYPE)s* restrict pvHt = NULL;
    %(VTYPE)s* restrict pvH = NULL;
    %(VTYPE)s* restrict pvHMax = NULL;
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
    %(VTYPE)s vMaxHUnit;
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
    end_query = 0;
    end_ref = 0;
    s1Len = profile->s1Len;
    matrix = profile->matrix;
    segWidth = %(LANES)s; /* number of values in vector unit */
    segLen = (s1Len + segWidth - 1) / segWidth;
#ifdef PARASAIL_ROWCOL
    offset = (s1Len - 1) %% segLen;
    position = (segWidth - 1) - (s1Len - 1) / segLen;
#endif
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
    vMaxHUnit = vNegLimit;
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
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_SCAN
        | PARASAIL_FLAG_BITS_%(WIDTH)s | PARASAIL_FLAG_LANES_%(LANES)s;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    /* initialize heap variables */
    pvE = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    pvHt= parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    pvH = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    pvHMax = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    pvGapper = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);

    /* validate heap variables */
    if (!pvE) return NULL;
    if (!pvHt) return NULL;
    if (!pvH) return NULL;
    if (!pvHMax) return NULL;
    if (!pvGapper) return NULL;

    parasail_memset_%(VTYPE)s(pvH, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvE, vNegLimit, segLen);
    {
        %(VTYPE)s vGapper = %(VSUB)s(vZero,vGapO);
        for (i=segLen-1; i>=0; --i) {
            %(VSTORE)s(pvGapper+i, vGapper);
            vGapper = %(VSUB)s(vGapper, vGapE);
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
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = vZero;
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
        vH = %(VMAX)s(vH, vZero);
        for (i=0; i<segLen; ++i) {
            vHt = %(VLOAD)s(pvHt+i);
            vF = %(VMAX)s(
                    %(VSUB)s(vF, vGapE),
                    %(VSUB)s(vH, vGapO));
            vH = %(VMAX)s(vHt, vF);
            vH = %(VMAX)s(vH, vZero);
            %(VSTORE)s(pvH+i, vH);
            vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vH);
            vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vH);
#ifdef PARASAIL_TABLE
            arr_store_si%(BITS)s(result->tables->score_table, vH, i, segLen, j, s2Len);
#endif
            vMaxH = %(VMAX)s(vH, vMaxH);
        } 

        {
            %(VTYPE)s vCompare = %(VCMPGT)s(vMaxH, vMaxHUnit);
            if (%(VMOVEMASK)s(vCompare)) {
                score = %(VHMAX)s(vMaxH);
                vMaxHUnit = %(VSET1)s(score);
                end_ref = j;
                (void)memcpy(pvHMax, pvH, sizeof(%(VTYPE)s)*segLen);
            }
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            %(INDEX)s k = 0;
            %(VTYPE)s vH = %(VLOAD)s(pvH + offset);
            for (k=0; k<position; ++k) {
                vH = %(VSHIFT)s(vH, %(BYTES)s);
            }
            result->rowcols->score_row[j] = (%(INT)s) %(VEXTRACT)s (vH, %(LAST_POS)s);
        }
#endif
    }

    /* Trace the alignment ending position on read. */
    {
        %(INT)s *t = (%(INT)s*)pvHMax;
        %(INDEX)s column_len = segLen * segWidth;
        end_query = s1Len;
        for (i = 0; i<column_len; ++i, ++t) {
            if (*t == score) {
                %(INDEX)s temp = i / segWidth + i %% segWidth * segLen;
                if (temp < end_query) {
                    end_query = temp;
                }
            }
        }
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        %(VTYPE)s vH = %(VLOAD)s(pvH+i);
        arr_store_col(result->rowcols->score_col, vH, i, segLen);
    }
#endif

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
    parasail_free(pvHMax);
    parasail_free(pvH);
    parasail_free(pvHt);
    parasail_free(pvE);

    return result;
}
