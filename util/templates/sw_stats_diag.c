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
    %(INT)s * restrict _HM_pr = NULL;
    %(INT)s * restrict _HS_pr = NULL;
    %(INT)s * restrict _HL_pr = NULL;
    %(INT)s * restrict _F_pr = NULL;
    %(INT)s * restrict _FM_pr = NULL;
    %(INT)s * restrict _FS_pr = NULL;
    %(INT)s * restrict _FL_pr = NULL;
    %(INT)s * restrict s2 = NULL;
    %(INT)s * restrict H_pr = NULL;
    %(INT)s * restrict HM_pr = NULL;
    %(INT)s * restrict HS_pr = NULL;
    %(INT)s * restrict HL_pr = NULL;
    %(INT)s * restrict F_pr = NULL;
    %(INT)s * restrict FM_pr = NULL;
    %(INT)s * restrict FS_pr = NULL;
    %(INT)s * restrict FL_pr = NULL;
    parasail_result_t *result = NULL;
    %(INDEX)s i = 0;
    %(INDEX)s j = 0;
    %(INDEX)s end_query = 0;
    %(INDEX)s end_ref = 0;
    %(INT)s NEG_LIMIT = 0;
    %(INT)s POS_LIMIT = 0;
    %(INT)s score = 0;
    %(INT)s matches = 0;
    %(INT)s similar = 0;
    %(INT)s length = 0;
    %(VTYPE)s vNegLimit;
    %(VTYPE)s vPosLimit;
    %(VTYPE)s vSaturationCheckMin;
    %(VTYPE)s vSaturationCheckMax;
    %(VTYPE)s vNegInf;
    %(VTYPE)s vNegInf0;
    %(VTYPE)s vOpen;
    %(VTYPE)s vGap;
    %(VTYPE)s vZero;
    %(VTYPE)s vOne;
    %(VTYPE)s vN;
    %(VTYPE)s vNegOne;
    %(VTYPE)s vI;
    %(VTYPE)s vJreset;
    %(VTYPE)s vMaxH;
    %(VTYPE)s vMaxM;
    %(VTYPE)s vMaxS;
    %(VTYPE)s vMaxL;
    %(VTYPE)s vEndI;
    %(VTYPE)s vEndJ;
    %(VTYPE)s vILimit;
    %(VTYPE)s vJLimit;

    /* validate inputs */
    PARASAIL_CHECK_NULL(_s2);
    PARASAIL_CHECK_GT0(s2Len);
    PARASAIL_CHECK_GE0(open);
    PARASAIL_CHECK_GE0(gap);
    PARASAIL_CHECK_NULL(matrix);
    if (matrix->type == PARASAIL_MATRIX_TYPE_PSSM) {
        PARASAIL_CHECK_NULL_PSSM_STATS(_s1);
    }
    else {
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
    end_query = 0;
    end_ref = 0;
    NEG_LIMIT = (-open < matrix->min ? INT%(WIDTH)s_MIN + open : INT%(WIDTH)s_MIN - matrix->min) + 1;
    POS_LIMIT = INT%(WIDTH)s_MAX - matrix->max - 1;
    score = NEG_LIMIT;
    matches = NEG_LIMIT;
    similar = NEG_LIMIT;
    length = NEG_LIMIT;
    vNegLimit = %(VSET1)s(NEG_LIMIT);
    vPosLimit = %(VSET1)s(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;
    vNegInf = %(VSET1)s(NEG_LIMIT);
    vNegInf0 = %(VRSHIFT)s(vNegInf, %(BYTES)s); /* shift in a 0 */
    vOpen = %(VSET1)s(open);
    vGap  = %(VSET1)s(gap);
    vZero = %(VSET1)s(0);
    vOne = %(VSET1)s(1);
    vN = %(VSET1)s(N);
    vNegOne = %(VSET1)s(-1);
    vI = %(VSET)s(%(DIAG_I)s);
    vJreset = %(VSET)s(%(DIAG_J)s);
    vMaxH = vNegInf;
    vMaxM = vNegInf;
    vMaxS = vNegInf;
    vMaxL = vNegInf;
    vEndI = vNegInf;
    vEndJ = vNegInf;
    vILimit = %(VSET1)s(s1Len);
    vJLimit = %(VSET1)s(s2Len);

    /* initialize result */
#ifdef PARASAIL_TABLE
    result = parasail_result_new_table3(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    result = parasail_result_new_rowcol3(s1Len, s2Len);
#else
    result = parasail_result_new_stats();
#endif
#endif
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_%(WIDTH)s | PARASAIL_FLAG_LANES_%(LANES)s;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    /* initialize heap variables */
    s1     = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s1Len+PAD);
    s2B    = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    _H_pr  = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    _HM_pr = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    _HS_pr = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    _HL_pr = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    _F_pr  = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    _FM_pr = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    _FS_pr = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    _FL_pr = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    s2 = s2B+PAD; /* will allow later for negative indices */
    H_pr = _H_pr+PAD;
    HM_pr = _HM_pr+PAD;
    HS_pr = _HS_pr+PAD;
    HL_pr = _HL_pr+PAD;
    F_pr = _F_pr+PAD;
    FM_pr = _FM_pr+PAD;
    FS_pr = _FS_pr+PAD;
    FL_pr = _FL_pr+PAD;

    /* validate heap variables */
    if (!s1) return NULL;
    if (!s2B) return NULL;
    if (!_H_pr) return NULL;
    if (!_HM_pr) return NULL;
    if (!_HS_pr) return NULL;
    if (!_HL_pr) return NULL;
    if (!_F_pr) return NULL;
    if (!_FM_pr) return NULL;
    if (!_FS_pr) return NULL;
    if (!_FL_pr) return NULL;

    /* convert _s1 from char to int in range 0-23 */
    for (i=0; i<s1Len; ++i) {
        s1[i] = matrix->mapper[(unsigned char)_s1[i]];
    }
    /* pad back of s1 with dummy values */
    for (i=s1Len; i<s1Len_PAD; ++i) {
        s1[i] = 0; /* point to first matrix row because we don't care */
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
        H_pr[j] = 0;
        HM_pr[j] = 0;
        HS_pr[j] = 0;
        HL_pr[j] = 0;
        F_pr[j] = NEG_LIMIT;
        FM_pr[j] = 0;
        FS_pr[j] = 0;
        FL_pr[j] = 0;
    }
    /* pad front of stored row values */
    for (j=-PAD; j<0; ++j) {
        H_pr[j] = NEG_LIMIT;
        HM_pr[j] = 0;
        HS_pr[j] = 0;
        HL_pr[j] = 0;
        F_pr[j] = NEG_LIMIT;
        FM_pr[j] = 0;
        FS_pr[j] = 0;
        FL_pr[j] = 0;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        H_pr[j] = NEG_LIMIT;
        HM_pr[j] = 0;
        HS_pr[j] = 0;
        HL_pr[j] = 0;
        F_pr[j] = NEG_LIMIT;
        FM_pr[j] = 0;
        FS_pr[j] = 0;
        FL_pr[j] = 0;
    }
    H_pr[-1] = 0; /* upper left corner */

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        %(VTYPE)s case1 = vZero;
        %(VTYPE)s case2 = vZero;
        %(VTYPE)s case0 = vZero;
        %(VTYPE)s vNH = vNegInf0;
        %(VTYPE)s vNM = vZero;
        %(VTYPE)s vNS = vZero;
        %(VTYPE)s vNL = vZero;
        %(VTYPE)s vWH = vNegInf0;
        %(VTYPE)s vWM = vZero;
        %(VTYPE)s vWS = vZero;
        %(VTYPE)s vWL = vZero;
        %(VTYPE)s vE = vNegInf;
        %(VTYPE)s vE_opn = vNegInf;
        %(VTYPE)s vE_ext = vNegInf;
        %(VTYPE)s vEM = vZero;
        %(VTYPE)s vES = vZero;
        %(VTYPE)s vEL = vZero;
        %(VTYPE)s vF = vNegInf;
        %(VTYPE)s vF_opn = vNegInf;
        %(VTYPE)s vF_ext = vNegInf;
        %(VTYPE)s vFM = vZero;
        %(VTYPE)s vFS = vZero;
        %(VTYPE)s vFL = vZero;
        %(VTYPE)s vJ = vJreset;
        %(VTYPE)s vs1 = %(VSET)s(
                %(DIAG_VS1)s);
        %(VTYPE)s vs2 = vNegInf;
        %(DIAG_MATROW_DECL)s
        %(VTYPE)s vIltLimit = %(VCMPLT)s(vI, vILimit);
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            %(VTYPE)s vMat;
            %(VTYPE)s vNWH = vNH;
            %(VTYPE)s vNWM = vNM;
            %(VTYPE)s vNWS = vNS;
            %(VTYPE)s vNWL = vNL;
            vNH = %(VRSHIFT)s(vWH, %(BYTES)s);
            vNH = %(VINSERT)s(vNH, H_pr[j], %(LAST_POS)s);
            vNM = %(VRSHIFT)s(vWM, %(BYTES)s);
            vNM = %(VINSERT)s(vNM, HM_pr[j], %(LAST_POS)s);
            vNS = %(VRSHIFT)s(vWS, %(BYTES)s);
            vNS = %(VINSERT)s(vNS, HS_pr[j], %(LAST_POS)s);
            vNL = %(VRSHIFT)s(vWL, %(BYTES)s);
            vNL = %(VINSERT)s(vNL, HL_pr[j], %(LAST_POS)s);
            vF = %(VRSHIFT)s(vF, %(BYTES)s);
            vF = %(VINSERT)s(vF, F_pr[j], %(LAST_POS)s);
            vFM = %(VRSHIFT)s(vFM, %(BYTES)s);
            vFM = %(VINSERT)s(vFM, FM_pr[j], %(LAST_POS)s);
            vFS = %(VRSHIFT)s(vFS, %(BYTES)s);
            vFS = %(VINSERT)s(vFS, FS_pr[j], %(LAST_POS)s);
            vFL = %(VRSHIFT)s(vFL, %(BYTES)s);
            vFL = %(VINSERT)s(vFL, FL_pr[j], %(LAST_POS)s);
            vF_opn = %(VSUB)s(vNH, vOpen);
            vF_ext = %(VSUB)s(vF, vGap);
            vF = %(VMAX)s(vF_opn, vF_ext);
            case1 = %(VCMPGT)s(vF_opn, vF_ext);
            vFM = %(VBLEND)s(vFM, vNM, case1);
            vFS = %(VBLEND)s(vFS, vNS, case1);
            vFL = %(VBLEND)s(vFL, vNL, case1);
            vFL = %(VADD)s(vFL, vOne);
            vE_opn = %(VSUB)s(vWH, vOpen);
            vE_ext = %(VSUB)s(vE, vGap);
            vE = %(VMAX)s(vE_opn, vE_ext);
            case1 = %(VCMPGT)s(vE_opn, vE_ext);
            vEM = %(VBLEND)s(vEM, vWM, case1);
            vES = %(VBLEND)s(vES, vWS, case1);
            vEL = %(VBLEND)s(vEL, vWL, case1);
            vEL = %(VADD)s(vEL, vOne);
            vs2 = %(VRSHIFT)s(vs2, %(BYTES)s);
            vs2 = %(VINSERT)s(vs2, s2[j], %(LAST_POS)s);
            vMat = %(VSET)s(
                    %(DIAG_MATROW_USE)s
                    );
            vNWH = %(VADD)s(vNWH, vMat);
            vWH = %(VMAX)s(vNWH, vE);
            vWH = %(VMAX)s(vWH, vF);
            vWH = %(VMAX)s(vWH, vZero);
            case1 = %(VCMPEQ)s(vWH, vNWH);
            case2 = %(VCMPEQ)s(vWH, vF);
            case0 = %(VCMPEQ)s(vWH, vZero);
            vWM = %(VBLEND)s(
                    %(VBLEND)s(vEM, vFM, case2),
                    %(VADD)s(vNWM,
                        %(VAND)s(
                            %(VCMPEQ)s(vs1,vs2),
                            vOne)),
                    case1);
            vWM = %(VBLEND)s(vWM, vZero, case0);
            vWS = %(VBLEND)s(
                    %(VBLEND)s(vES, vFS, case2),
                    %(VADD)s(vNWS,
                        %(VAND)s(
                            %(VCMPGT)s(vMat,vZero),
                            vOne)),
                    case1);
            vWS = %(VBLEND)s(vWS, vZero, case0);
            vWL = %(VBLEND)s(
                    %(VBLEND)s(vEL, vFL, case2),
                    %(VADD)s(vNWL, vOne), case1);
            vWL = %(VBLEND)s(vWL, vZero, case0);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                %(VTYPE)s cond = %(VCMPEQ)s(vJ,vNegOne);
                vWH = %(VANDNOT)s(cond, vWH);
                vWM = %(VANDNOT)s(cond, vWM);
                vWS = %(VANDNOT)s(cond, vWS);
                vWL = %(VANDNOT)s(cond, vWL);
                vE = %(VBLEND)s(vE, vNegInf, cond);
                vEM = %(VANDNOT)s(cond, vEM);
                vES = %(VANDNOT)s(cond, vES);
                vEL = %(VANDNOT)s(cond, vEL);
            }
            /* cannot start checking sat until after J clears boundary */
            if (j > PAD) {
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vWH);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWH);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWM);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWS);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWL);
            }
#ifdef PARASAIL_TABLE
            arr_store_si%(BITS)s(result->stats->tables->score_table, vWH, i, s1Len, j, s2Len);
            arr_store_si%(BITS)s(result->stats->tables->matches_table, vWM, i, s1Len, j, s2Len);
            arr_store_si%(BITS)s(result->stats->tables->similar_table, vWS, i, s1Len, j, s2Len);
            arr_store_si%(BITS)s(result->stats->tables->length_table, vWL, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->stats->rowcols->score_row,   result->stats->rowcols->score_col, vWH, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->matches_row, result->stats->rowcols->matches_col, vWM, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->similar_row, result->stats->rowcols->similar_col, vWS, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->length_row,  result->stats->rowcols->length_col, vWL, i, s1Len, j, s2Len);
#endif
            H_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vWH,0);
            HM_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vWM,0);
            HS_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vWS,0);
            HL_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vWL,0);
            F_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vF,0);
            FM_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vFM,0);
            FS_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vFS,0);
            FL_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vFL,0);
            /* as minor diagonal vector passes across table, extract
             * max values within the i,j bounds */
            {
                %(VTYPE)s cond_valid_J = %(VAND)s(
                        %(VCMPGT)s(vJ, vNegOne),
                        %(VCMPLT)s(vJ, vJLimit));
                %(VTYPE)s cond_valid_IJ = %(VAND)s(cond_valid_J, vIltLimit);
                %(VTYPE)s cond_eq = %(VCMPEQ)s(vWH, vMaxH);
                %(VTYPE)s cond_max = %(VCMPGT)s(vWH, vMaxH);
                %(VTYPE)s cond_all = %(VAND)s(cond_max, cond_valid_IJ);
                %(VTYPE)s cond_Jlt = %(VCMPLT)s(vJ, vEndJ);
                vMaxH = %(VBLEND)s(vMaxH, vWH, cond_all);
                vMaxM = %(VBLEND)s(vMaxM, vWM, cond_all);
                vMaxS = %(VBLEND)s(vMaxS, vWS, cond_all);
                vMaxL = %(VBLEND)s(vMaxL, vWL, cond_all);
                vEndI = %(VBLEND)s(vEndI, vI, cond_all);
                vEndJ = %(VBLEND)s(vEndJ, vJ, cond_all);
                cond_all = %(VAND)s(cond_Jlt, cond_eq);
                cond_all = %(VAND)s(cond_all, cond_valid_IJ);
                vMaxM = %(VBLEND)s(vMaxM, vWM, cond_all);
                vMaxS = %(VBLEND)s(vMaxS, vWS, cond_all);
                vMaxL = %(VBLEND)s(vMaxL, vWL, cond_all);
                vEndI = %(VBLEND)s(vEndI, vI, cond_all);
                vEndJ = %(VBLEND)s(vEndJ, vJ, cond_all);
            }
            vJ = %(VADD)s(vJ, vOne);
        }
        vI = %(VADD)s(vI, vN);
    }

    /* alignment ending position */
    {
        %(INT)s *t = (%(INT)s*)&vMaxH;
        %(INT)s *m = (%(INT)s*)&vMaxM;
        %(INT)s *s = (%(INT)s*)&vMaxS;
        %(INT)s *l = (%(INT)s*)&vMaxL;
        %(INT)s *i = (%(INT)s*)&vEndI;
        %(INT)s *j = (%(INT)s*)&vEndJ;
        %(INDEX)s k;
        for (k=0; k<N; ++k, ++t, ++m, ++s, ++l, ++i, ++j) {
            if (*t > score) {
                score = *t;
                matches = *m;
                similar = *s;
                length = *l;
                end_query = *i;
                end_ref = *j;
            }
            else if (*t == score) {
                if (*j < end_ref) {
                    matches = *m;
                    similar = *s;
                    length = *l;
                    end_query = *i;
                    end_ref = *j;
                }
                else if (*j == end_ref && *i < end_query) {
                    matches = *m;
                    similar = *s;
                    length = *l;
                    end_query = *i;
                    end_ref = *j;
                }
            }
        }
    }

    if (%(VMOVEMASK)s(%(VOR)s(
            %(VCMPLT)s(vSaturationCheckMin, vNegLimit),
            %(VCMPGT)s(vSaturationCheckMax, vPosLimit)))) {
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

    parasail_free(_FL_pr);
    parasail_free(_FS_pr);
    parasail_free(_FM_pr);
    parasail_free(_F_pr);
    parasail_free(_HL_pr);
    parasail_free(_HS_pr);
    parasail_free(_HM_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}

