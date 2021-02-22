/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>

%(HEADER)s

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_%(ISA)s.h"

%(FIXES)s

#ifdef PARASAIL_TABLE
static inline void arr_store_si%(BITS)s(
        int *array,
        %(VTYPE)s vWH,
        %(INDEX)s i,
        %(INDEX)s s1Len,
        %(INDEX)s j,
        %(INDEX)s s2Len)
{
%(PRINTER)s
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_rowcol(
        int *row,
        int *col,
        %(VTYPE)s vWH,
        %(INDEX)s i,
        %(INDEX)s s1Len,
        %(INDEX)s j,
        %(INDEX)s s2Len)
{
%(PRINTER_ROWCOL)s
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME %(NAME_TABLE)s
#else
#ifdef PARASAIL_ROWCOL
#define FNAME %(NAME_ROWCOL)s
#else
#define FNAME %(NAME)s
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int _s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    /* declare local variables */
    %(INDEX)s N = 0;
    %(INDEX)s PAD = 0;
    %(INDEX)s PAD2 = 0;
    %(INDEX)s s1Len = 0;
    %(INDEX)s s1Len_PAD = 0;
    %(INDEX)s s2Len_PAD = 0;
    %(INT)s * restrict s1 = NULL;
    %(INT)s * restrict s2B = NULL;
    %(INT)s * restrict _H_pr = NULL;
    %(INT)s * restrict _F_pr = NULL;
    %(INT)s * restrict s2 = NULL;
    %(INT)s * restrict H_pr = NULL;
    %(INT)s * restrict F_pr = NULL;
    parasail_result_t *result = NULL;
    %(INDEX)s i = 0;
    %(INDEX)s j = 0;
    %(INDEX)s end_query = 0;
    %(INDEX)s end_ref = 0;
    %(INT)s NEG_LIMIT = 0;
    %(INT)s POS_LIMIT = 0;
    %(INT)s score = 0;
    %(VTYPE)s vNegLimit;
    %(VTYPE)s vPosLimit;
    %(VTYPE)s vSaturationCheckMin;
    %(VTYPE)s vSaturationCheckMax;
    %(VTYPE)s vNegInf;
    %(VTYPE)s vOpen;
    %(VTYPE)s vGap;
    %(VTYPE)s vOne;
    %(VTYPE)s vN;
    %(VTYPE)s vGapN;
    %(VTYPE)s vNegOne;
    %(VTYPE)s vI;
    %(VTYPE)s vJreset;
    %(VTYPE)s vMax;
    %(VTYPE)s vILimit;
    %(VTYPE)s vILimit1;
    %(VTYPE)s vJLimit;
    %(VTYPE)s vJLimit1;
    %(VTYPE)s vIBoundary;

    /* validate inputs */
    PARASAIL_CHECK_NULL(_s2);
    PARASAIL_CHECK_GT0(s2Len);
    PARASAIL_CHECK_GE0(open);
    PARASAIL_CHECK_GE0(gap);
    PARASAIL_CHECK_NULL(matrix);
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        PARASAIL_CHECK_NULL(_s1);
        PARASAIL_CHECK_GT0(_s1Len);
    }

    /* initialize stack variables */
    N = %(LANES)s; /* number of values in vector */
    PAD = N-1;
    PAD2 = PAD*2;
    s1Len = matrix->type == PARASAIL_MATRIX_TYPE_SQUARE ? _s1Len : matrix->length;
    s1Len_PAD = s1Len+PAD;
    s2Len_PAD = s2Len+PAD;
    i = 0;
    j = 0;
    end_query = s1Len-1;
    end_ref = s2Len-1;
    NEG_LIMIT = (-open < matrix->min ? INT%(WIDTH)s_MIN + open : INT%(WIDTH)s_MIN - matrix->min) + 1;
    POS_LIMIT = INT%(WIDTH)s_MAX - matrix->max - 1;
    score = NEG_LIMIT;
    vNegLimit = %(VSET1)s(NEG_LIMIT);
    vPosLimit = %(VSET1)s(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;
    vNegInf = %(VSET1)s(NEG_LIMIT);
    vOpen = %(VSET1)s(open);
    vGap  = %(VSET1)s(gap);
    vOne = %(VSET1)s(1);
    vN = %(VSET1)s(N);
    vGapN = %(VSET1)s(gap*N);
    vNegOne = %(VSET1)s(-1);
    vI = %(VSET)s(%(DIAG_I)s);
    vJreset = %(VSET)s(%(DIAG_J)s);
    vMax = vNegInf;
    vILimit = %(VSET1)s(s1Len);
    vILimit1 = %(VSUB)s(vILimit, vOne);
    vJLimit = %(VSET1)s(s2Len);
    vJLimit1 = %(VSUB)s(vJLimit, vOne);
    vIBoundary = %(VSET)s(
            %(DIAG_IBoundary)s);

    /* initialize result */
#ifdef PARASAIL_TABLE
    result = parasail_result_new_table1(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    result = parasail_result_new_rowcol1(s1Len, s2Len);
#else
    result = parasail_result_new();
#endif
#endif
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_NW | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_BITS_%(WIDTH)s | PARASAIL_FLAG_LANES_%(LANES)s;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    /* initialize heap variables */
    s2B= parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    _H_pr = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    _F_pr = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    s2 = s2B+PAD; /* will allow later for negative indices */
    H_pr = _H_pr+PAD;
    F_pr = _F_pr+PAD;

    /* validate heap variables */
    if (!s2B) return NULL;
    if (!_H_pr) return NULL;
    if (!_F_pr) return NULL;

    /* convert _s1 from char to int in range 0-23 */
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        s1 = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s1Len+PAD);
        if (!s1) return NULL;
        for (i=0; i<s1Len; ++i) {
            s1[i] = matrix->mapper[(unsigned char)_s1[i]];
        }
        /* pad back of s1 with dummy values */
        for (i=s1Len; i<s1Len_PAD; ++i) {
            s1[i] = 0; /* point to first matrix row because we don't care */
        }
    }

    /* convert _s2 from char to int in range 0-23 */
    for (j=0; j<s2Len; ++j) {
        s2[j] = matrix->mapper[(unsigned char)_s2[j]];
    }
    /* pad front of s2 with dummy values */
    for (j=-PAD; j<0; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }
    /* pad back of s2 with dummy values */
    for (j=s2Len; j<s2Len_PAD; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }

    /* set initial values for stored row */
    for (j=0; j<s2Len; ++j) {
        H_pr[j] = -open - j*gap;
        F_pr[j] = NEG_LIMIT;
    }
    /* pad front of stored row values */
    for (j=-PAD; j<0; ++j) {
        H_pr[j] = NEG_LIMIT;
        F_pr[j] = NEG_LIMIT;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        H_pr[j] = NEG_LIMIT;
        F_pr[j] = NEG_LIMIT;
    }
    H_pr[-1] = 0; /* upper left corner */

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        %(VTYPE)s vNH = vNegInf;
        %(VTYPE)s vWH = vNegInf;
        %(VTYPE)s vE = vNegInf;
        %(VTYPE)s vF = vNegInf;
        %(VTYPE)s vJ = vJreset;
        %(DIAG_MATROW_DECL)s
        vNH = %(VRSHIFT)s(vNH, %(BYTES)s);
        vNH = %(VINSERT)s(vNH, H_pr[-1], %(LAST_POS)s);
        vWH = %(VRSHIFT)s(vWH, %(BYTES)s);
        vWH = %(VINSERT)s(vWH, -open - i*gap, %(LAST_POS)s);
        H_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            %(VTYPE)s vMat;
            %(VTYPE)s vNWH = vNH;
            vNH = %(VRSHIFT)s(vWH, %(BYTES)s);
            vNH = %(VINSERT)s(vNH, H_pr[j], %(LAST_POS)s);
            vF = %(VRSHIFT)s(vF, %(BYTES)s);
            vF = %(VINSERT)s(vF, F_pr[j], %(LAST_POS)s);
            vF = %(VMAX)s(
                    %(VSUB)s(vNH, vOpen),
                    %(VSUB)s(vF, vGap));
            vE = %(VMAX)s(
                    %(VSUB)s(vWH, vOpen),
                    %(VSUB)s(vE, vGap));
            vMat = %(VSET)s(
                    %(DIAG_MATROW_USE)s
                    );
            vNWH = %(VADD)s(vNWH, vMat);
            vWH = %(VMAX)s(vNWH, vE);
            vWH = %(VMAX)s(vWH, vF);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                %(VTYPE)s cond = %(VCMPEQ)s(vJ,vNegOne);
                vWH = %(VBLEND)s(vWH, vIBoundary, cond);
                vF = %(VBLEND)s(vF, vNegInf, cond);
                vE = %(VBLEND)s(vE, vNegInf, cond);
            }
            /* cannot start checking sat until after J clears boundary */
            if (j > PAD) {
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vWH);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWH);
            }
#ifdef PARASAIL_TABLE
            arr_store_si%(BITS)s(result->tables->score_table, vWH, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->rowcols->score_row, result->rowcols->score_col, vWH, i, s1Len, j, s2Len);
#endif
            H_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vWH,0);
            F_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vF,0);
            /* as minor diagonal vector passes across table, extract
               last table value at the i,j bound */
            {
                %(VTYPE)s cond_valid_I = %(VCMPEQ)s(vI, vILimit1);
                %(VTYPE)s cond_valid_J = %(VCMPEQ)s(vJ, vJLimit1);
                %(VTYPE)s cond_all = %(VAND)s(cond_valid_I, cond_valid_J);
                vMax = %(VBLEND)s(vMax, vWH, cond_all);
            }
            vJ = %(VADD)s(vJ, vOne);
        }
        vI = %(VADD)s(vI, vN);
        vIBoundary = %(VSUB)s(vIBoundary, vGapN);
    }

    /* max in vMax */
    for (i=0; i<N; ++i) {
        %(INT)s value;
        value = (%(INT)s) %(VEXTRACT)s(vMax, %(LAST_POS)s);
        if (value > score) {
            score = value;
        }
        vMax = %(VSHIFT)s(vMax, %(BYTES)s);
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

    parasail_free(_F_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        parasail_free(s1);
    }

    return result;
}

