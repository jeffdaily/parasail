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

#define SG_TRACE
#define SG_SUFFIX %(SUFFIX)s
#include "sg_helper.h"

%(FIXES)s

static inline void arr_store_si%(BITS)s(
        int8_t *array,
        %(VTYPE)s vWH,
        %(INDEX)s i,
        %(INDEX)s s1Len,
        %(INDEX)s j,
        %(INDEX)s s2Len)
{
%(PRINTER_TRACE)s
}

#define FNAME %(NAME_TRACE)s

parasail_result_t* FNAME(
        const char * const restrict _s1, const int _s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix,
        int s1_beg, int s1_end, int s2_beg, int s2_end)
{
    /* declare local variables */
    %(INDEX)s N = 0;
    %(INDEX)s PAD = 0;
    %(INDEX)s PAD2 = 0;
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
    %(INDEX)s s1Len = 0;
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
    %(VTYPE)s vNegInf0;
    %(VTYPE)s vOpen;
    %(VTYPE)s vGap;
    %(VTYPE)s vOne16;
    %(VTYPE)s vN16;
    %(VTYPE)s vNegOne16;
    %(VTYPE)s vILo16;
    %(VTYPE)s vIHi16;
    %(VTYPE)s vJresetLo16;
    %(VTYPE)s vJresetHi16;
    %(VTYPE)s vMaxH;
    %(VTYPE)s vEndILo;
    %(VTYPE)s vEndIHi;
    %(VTYPE)s vEndJLo;
    %(VTYPE)s vEndJHi;
    %(VTYPE)s vILimit16;
    %(VTYPE)s vILimit116;
    %(VTYPE)s vJLimit16;
    %(VTYPE)s vJLimit116;
    %(VTYPE)s vTDiag;
    %(VTYPE)s vTIns;
    %(VTYPE)s vTDel;
    %(VTYPE)s vTDiagE;
    %(VTYPE)s vTInsE;
    %(VTYPE)s vTDiagF;
    %(VTYPE)s vTDelF;

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
    end_query = 0;
    end_ref = 0;
    NEG_LIMIT = (-open < matrix->min ? INT%(WIDTH)s_MIN + open : INT%(WIDTH)s_MIN - matrix->min) + 1;
    POS_LIMIT = INT%(WIDTH)s_MAX - matrix->max - 1;
    score = NEG_LIMIT;
    vNegLimit = %(VSET1)s(NEG_LIMIT);
    vPosLimit = %(VSET1)s(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;
    vNegInf = %(VSET1)s(NEG_LIMIT);
    vNegInf0 = %(VRSHIFT)s(vNegInf, %(BYTES)s); /* shift in a 0 */
    vOpen = %(VSET1)s(open);
    vGap  = %(VSET1)s(gap);
    vOne16 = %(VSET1x16)s(1);
    vN16 = %(VSET1x16)s(N);
    vNegOne16 = %(VSET1x16)s(-1);
    vILo16 = %(VSETx16)s(%(DIAG_ILO)s);
    vIHi16 = %(VSETx16)s(%(DIAG_IHI)s);
    vJresetLo16 = %(VSETx16)s(%(DIAG_JLO)s);
    vJresetHi16 = %(VSETx16)s(%(DIAG_JHI)s);
    vMaxH = vNegInf;
    vEndILo = vNegInf;
    vEndIHi = vNegInf;
    vEndJLo = vNegInf;
    vEndJHi = vNegInf;
    vILimit16 = %(VSET1x16)s(s1Len);
    vILimit116 = %(VSUBx16)s(vILimit16, vOne16);
    vJLimit16 = %(VSET1x16)s(s2Len);
    vJLimit116 = %(VSUBx16)s(vJLimit16, vOne16);
    vTDiag = %(VSET1)s(PARASAIL_DIAG);
    vTIns = %(VSET1)s(PARASAIL_INS);
    vTDel = %(VSET1)s(PARASAIL_DEL);
    vTDiagE = %(VSET1)s(PARASAIL_DIAG_E);
    vTInsE = %(VSET1)s(PARASAIL_INS_E);
    vTDiagF = %(VSET1)s(PARASAIL_DIAG_F);
    vTDelF = %(VSET1)s(PARASAIL_DEL_F);

    /* initialize result */
    result = parasail_result_new_trace(s1Len, s2Len, %(ALIGNMENT)s, sizeof(int8_t));
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_%(WIDTH)s | PARASAIL_FLAG_LANES_%(LANES)s;

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
        H_pr[j] = 0;
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

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        %(VTYPE)s vNH = vNegInf0;
        %(VTYPE)s vWH = vNegInf0;
        %(VTYPE)s vE = vNegInf;
        %(VTYPE)s vE_opn = vNegInf;
        %(VTYPE)s vE_ext = vNegInf;
        %(VTYPE)s vF = vNegInf;
        %(VTYPE)s vF_opn = vNegInf;
        %(VTYPE)s vF_ext = vNegInf;
        %(VTYPE)s vJLo16 = vJresetLo16;
        %(VTYPE)s vJHi16 = vJresetHi16;
        %(DIAG_MATROW_DECL)s
        %(VTYPE)s vIltLimit = %(VPACKS)s(
                %(VCMPLTx16)s(vILo16, vILimit16),
                %(VCMPLTx16)s(vIHi16, vILimit16));
        %(VTYPE)s vIeqLimit1 = %(VPACKS)s(
                %(VCMPEQx16)s(vILo16, vILimit116),
                %(VCMPEQx16)s(vIHi16, vILimit116));
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            %(VTYPE)s vMat;
            %(VTYPE)s vNWH = vNH;
            vNH = %(VRSHIFT)s(vWH, %(BYTES)s);
            vNH = %(VINSERT)s(vNH, H_pr[j], %(LAST_POS)s);
            vF = %(VRSHIFT)s(vF, %(BYTES)s);
            vF = %(VINSERT)s(vF, F_pr[j], %(LAST_POS)s);
            vF_opn = %(VSUB)s(vNH, vOpen);
            vF_ext = %(VSUB)s(vF, vGap);
            vF = %(VMAX)s(vF_opn, vF_ext);
            vE_opn = %(VSUB)s(vWH, vOpen);
            vE_ext = %(VSUB)s(vE, vGap);
            vE = %(VMAX)s(vE_opn, vE_ext);
            vMat = %(VSET)s(
                    %(DIAG_MATROW_USE)s
                    );
            vNWH = %(VADD)s(vNWH, vMat);
            vWH = %(VMAX)s(vNWH, vE);
            vWH = %(VMAX)s(vWH, vF);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                %(VTYPE)s cond = %(VPACKS)s(
                        %(VCMPEQx16)s(vJLo16,vNegOne16),
                        %(VCMPEQx16)s(vJHi16,vNegOne16));
                vWH = %(VANDNOT)s(cond, vWH);
                vF = %(VBLEND)s(vF, vNegInf, cond);
                vE = %(VBLEND)s(vE, vNegInf, cond);
            }
            /* cannot start checking sat until after J clears boundary */
            if (j > PAD) {
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vWH);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWH);
            }
            /* trace table */
            {
                %(VTYPE)s case1 = %(VCMPEQ)s(vWH, vNWH);
                %(VTYPE)s case2 = %(VCMPEQ)s(vWH, vF);
                %(VTYPE)s vT = %(VBLEND)s(
                        %(VBLEND)s(vTIns, vTDel, case2),
                        vTDiag,
                        case1);
                %(VTYPE)s condE = %(VCMPGT)s(vE_opn, vE_ext);
                %(VTYPE)s condF = %(VCMPGT)s(vF_opn, vF_ext);
                %(VTYPE)s vET = %(VBLEND)s(vTInsE, vTDiagE, condE);
                %(VTYPE)s vFT = %(VBLEND)s(vTDelF, vTDiagF, condF);
                vT = %(VOR)s(vT, vET);
                vT = %(VOR)s(vT, vFT);
                arr_store_si%(BITS)s(result->trace->trace_table, vT, i, s1Len, j, s2Len);
            }
            H_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vWH,0);
            F_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vF,0);
            /* as minor diagonal vector passes across the i or j limit
             * boundary, extract the last value of the column or row */
            {
                %(VTYPE)s vJeqLimit1 = %(VPACKS)s(
                        %(VCMPEQx16)s(vJLo16, vJLimit116),
                        %(VCMPEQx16)s(vJHi16, vJLimit116));
                %(VTYPE)s vJgtNegOne = %(VPACKS)s(
                        %(VCMPGTx16)s(vJLo16, vNegOne16),
                        %(VCMPGTx16)s(vJHi16, vNegOne16));
                %(VTYPE)s vJltLimit = %(VPACKS)s(
                        %(VCMPLTx16)s(vJLo16, vJLimit16),
                        %(VCMPLTx16)s(vJHi16, vJLimit16));
                %(VTYPE)s cond_j = %(VAND)s(vIltLimit, vJeqLimit1);
                %(VTYPE)s cond_i = %(VAND)s(vIeqLimit1,
                        %(VAND)s(vJgtNegOne, vJltLimit));
                %(VTYPE)s cond_valid_IJ = %(VOR)s(cond_i, cond_j);
                %(VTYPE)s cond_eq = %(VCMPEQ)s(vWH, vMaxH);
                %(VTYPE)s cond_max = %(VCMPGT)s(vWH, vMaxH);
                %(VTYPE)s cond_all = %(VAND)s(cond_max, cond_valid_IJ);
                %(VTYPE)s cond_Jlt = %(VPACKS)s(
                        %(VCMPLTx16)s(vJLo16, vEndJLo),
                        %(VCMPLTx16)s(vJHi16, vEndJHi));
                %(VTYPE)s cond_lo = %(VUNPACKLO)s(cond_all, cond_all);
                %(VTYPE)s cond_hi = %(VUNPACKHI)s(cond_all, cond_all);
                vMaxH = %(VBLEND)s(vMaxH, vWH, cond_all);
                vEndILo = %(VBLEND)s(vEndILo, vILo16, cond_lo);
                vEndIHi = %(VBLEND)s(vEndIHi, vIHi16, cond_hi);
                vEndJLo = %(VBLEND)s(vEndJLo, vJLo16, cond_lo);
                vEndJHi = %(VBLEND)s(vEndJHi, vJHi16, cond_hi);
                cond_all = %(VAND)s(cond_Jlt, cond_eq);
                cond_all = %(VAND)s(cond_all, cond_valid_IJ);
                cond_lo = %(VUNPACKLO)s(cond_all, cond_all);
                cond_hi = %(VUNPACKHI)s(cond_all, cond_all);
                vEndILo = %(VBLEND)s(vEndILo, vILo16, cond_lo);
                vEndIHi = %(VBLEND)s(vEndIHi, vIHi16, cond_hi);
                vEndJLo = %(VBLEND)s(vEndJLo, vJLo16, cond_lo);
                vEndJHi = %(VBLEND)s(vEndJHi, vJHi16, cond_hi);
            }
            vJLo16 = %(VADDx16)s(vJLo16, vOne16);
            vJHi16 = %(VADDx16)s(vJHi16, vOne16);
        }
        vILo16 = %(VADDx16)s(vILo16, vN16);
        vIHi16 = %(VADDx16)s(vIHi16, vN16);
    }

    /* alignment ending position */
    {
        %(INT)s *t = (%(INT)s*)&vMaxH;
        int16_t *ilo = (int16_t*)&vEndILo;
        int16_t *jlo = (int16_t*)&vEndJLo;
        int16_t *ihi = (int16_t*)&vEndIHi;
        int16_t *jhi = (int16_t*)&vEndJHi;
        %(INDEX)s k;
        for (k=0; k<N/2; ++k, ++t, ++ilo, ++jlo) {
            if (*t > score) {
                score = *t;
                end_query = *ilo;
                end_ref = *jlo;
            }
            else if (*t == score) {
                if (*jlo < end_ref) {
                    end_query = *ilo;
                    end_ref = *jlo;
                }
                else if (*jlo == end_ref && *ilo < end_query) {
                    end_query = *ilo;
                    end_ref = *jlo;
                }
            }
        }
        for (k=N/2; k<N; ++k, ++t, ++ihi, ++jhi) {
            if (*t > score) {
                score = *t;
                end_query = *ihi;
                end_ref = *jhi;
            }
            else if (*t == score) {
                if (*jhi < end_ref) {
                    end_query = *ihi;
                    end_ref = *jhi;
                }
                else if (*jhi == end_ref && *ihi < end_query) {
                    end_query = *ihi;
                    end_ref = *jhi;
                }
            }
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

    parasail_free(_F_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        parasail_free(s1);
    }

    return result;
}

SG_IMPL_ALL

