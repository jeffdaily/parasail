/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>

%(HEADER)s

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_%(ISA)s.h"

#define NEG_INF %(NEG_INF)s
%(FIXES)s

#ifdef PARASAIL_TABLE
static inline void arr_store_si%(BITS)s(
        int *array,
        %(VTYPE)s vWscore,
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
        %(VTYPE)s vWscore,
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
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    const %(INDEX)s N = %(LANES)s; /* number of values in vector */
    const %(INDEX)s PAD = N-1;
    const %(INDEX)s PAD2 = PAD*2;
    const %(INDEX)s s1Len_PAD = s1Len+PAD;
    const %(INDEX)s s2Len_PAD = s2Len+PAD;
    %(INT)s * const restrict s1 = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s1Len+PAD);
    %(INT)s * const restrict s2B= parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    %(INT)s * const restrict _tbl_pr = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    %(INT)s * const restrict _del_pr = parasail_memalign_%(INT)s(%(ALIGNMENT)s, s2Len+PAD2);
    %(INT)s * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    %(INT)s * const restrict tbl_pr = _tbl_pr+PAD;
    %(INT)s * const restrict del_pr = _del_pr+PAD;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif
    %(INDEX)s i = 0;
    %(INDEX)s j = 0;
    %(INDEX)s end_query = 0;
    %(INDEX)s end_ref = 0;
    %(INT)s score = NEG_INF;
    %(VTYPE)s vNegInf = %(VSET1)s(NEG_INF);
    %(VTYPE)s vNegInf0 = %(VRSHIFT)s(vNegInf, %(BYTES)s); /* shift in a 0 */
    %(VTYPE)s vOpen = %(VSET1)s(open);
    %(VTYPE)s vGap  = %(VSET1)s(gap);
    %(VTYPE)s vZero = %(VSET1)s(0);
    %(VTYPE)s vOne16 = %(VSET1x16)s(1);
    %(VTYPE)s vNegOne16 = %(VSET1x16)s(-1);
    %(VTYPE)s vN16 = %(VSET1x16)s(N);
    %(VTYPE)s vILo16 = %(VSETx16)s(%(DIAG_ILO)s);
    %(VTYPE)s vIHi16 = %(VSETx16)s(%(DIAG_IHI)s);
    %(VTYPE)s vJresetLo16 = %(VSETx16)s(%(DIAG_JLO)s);
    %(VTYPE)s vJresetHi16 = %(VSETx16)s(%(DIAG_JHI)s);
    %(VTYPE)s vMax = vNegInf;
    %(VTYPE)s vEndILo = vNegInf;
    %(VTYPE)s vEndIHi = vNegInf;
    %(VTYPE)s vEndJLo = vNegInf;
    %(VTYPE)s vEndJHi = vNegInf;
    %(VTYPE)s vILimit16 = %(VSET1x16)s(s1Len);
    %(VTYPE)s vJLimit16 = %(VSET1x16)s(s2Len);
    %(SATURATION_CHECK_INIT)s

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
        tbl_pr[j] = 0;
        del_pr[j] = NEG_INF;
    }
    /* pad front of stored row values */
    for (j=-PAD; j<0; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
    }

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        %(VTYPE)s vNscore = vNegInf0;
        %(VTYPE)s vWscore = vNegInf0;
        %(VTYPE)s vIns = vNegInf;
        %(VTYPE)s vDel = vNegInf;
        %(VTYPE)s vJLo16 = vJresetLo16;
        %(VTYPE)s vJHi16 = vJresetHi16;
        %(DIAG_MATROW_DECL)s
        %(VTYPE)s vIltLimit = %(VPACKS)s(
                    %(VCMPLTx16)s(vILo16, vILimit16),
                    %(VCMPLTx16)s(vIHi16, vILimit16));
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            %(VTYPE)s vMat;
            %(VTYPE)s vNWscore = vNscore;
            vNscore = %(VRSHIFT)s(vWscore, %(BYTES)s);
            vNscore = %(VINSERT)s(vNscore, tbl_pr[j], %(LAST_POS)s);
            vDel = %(VRSHIFT)s(vDel, %(BYTES)s);
            vDel = %(VINSERT)s(vDel, del_pr[j], %(LAST_POS)s);
            vDel = %(VMAX)s(
                    %(VSUB)s(vNscore, vOpen),
                    %(VSUB)s(vDel, vGap));
            vIns = %(VMAX)s(
                    %(VSUB)s(vWscore, vOpen),
                    %(VSUB)s(vIns, vGap));
            vMat = %(VSET)s(
                    %(DIAG_MATROW_USE)s
                    );
            vNWscore = %(VADD)s(vNWscore, vMat);
            vWscore = %(VMAX)s(vNWscore, vIns);
            vWscore = %(VMAX)s(vWscore, vDel);
            vWscore = %(VMAX)s(vWscore, vZero);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                %(VTYPE)s cond = %(VPACKS)s(
                        %(VCMPEQx16)s(vJLo16,vNegOne16),
                        %(VCMPEQx16)s(vJHi16,vNegOne16));
                vWscore = %(VANDNOT)s(cond, vWscore);
                vDel = %(VBLEND)s(vDel, vNegInf, cond);
                vIns = %(VBLEND)s(vIns, vNegInf, cond);
            }
            %(SATURATION_CHECK_MID)s
#ifdef PARASAIL_TABLE
            arr_store_si%(BITS)s(result->score_table, vWscore, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->score_row, result->score_col, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vWscore,0);
            del_pr[j-%(LAST_POS)s] = (%(INT)s)%(VEXTRACT)s(vDel,0);
            /* as minor diagonal vector passes across table, extract
             * max values within the i,j bounds */
            {
                %(VTYPE)s cond_valid_J = %(VAND)s(
                        %(VPACKS)s(
                            %(VCMPGTx16)s(vJLo16, vNegOne16),
                            %(VCMPGTx16)s(vJHi16, vNegOne16)),
                        %(VPACKS)s(
                            %(VCMPLTx16)s(vJLo16, vJLimit16),
                            %(VCMPLTx16)s(vJHi16, vJLimit16)));
                %(VTYPE)s cond_valid_IJ = %(VAND)s(cond_valid_J, vIltLimit);
                %(VTYPE)s cond_eq = %(VCMPEQ)s(vWscore, vMax);
                %(VTYPE)s cond_max = %(VCMPGT)s(vWscore, vMax);
                %(VTYPE)s cond_all = %(VAND)s(cond_max, cond_valid_IJ);
                %(VTYPE)s cond_Jlt = %(VPACKS)s(
                        %(VCMPLTx16)s(vJLo16, vEndJLo),
                        %(VCMPLTx16)s(vJHi16, vEndJHi));
                %(VTYPE)s cond_lo = %(VUNPACKLO)s(cond_all, cond_all);
                %(VTYPE)s cond_hi = %(VUNPACKHI)s(cond_all, cond_all);
                vMax = %(VBLEND)s(vMax, vWscore, cond_all);
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
        %(INT)s *t = (%(INT)s*)&vMax;
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

    %(SATURATION_CHECK_FINAL)s

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(_del_pr);
    parasail_free(_tbl_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}

