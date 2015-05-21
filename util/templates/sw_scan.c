/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
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

#define NEG_INF %(NEG_INF)s
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

#ifdef PARASAIL_TABLE
#define FNAME %(NAME_TABLE)s
#else
#define FNAME %(NAME)s
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    %(INDEX)s i = 0;
    %(INDEX)s j = 0;
    %(INDEX)s k = 0;
    %(INDEX)s segNum = 0;
    const %(INDEX)s n = 24; /* number of amino acids in table */
    const %(INDEX)s segWidth = %(LANES)s; /* number of values in vector unit */
    const %(INDEX)s segLen = (s1Len + segWidth - 1) / segWidth;
    %(VTYPE)s* const restrict pvP = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, n * segLen);
    %(VTYPE)s* const restrict pvE = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvHt= parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvH = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s vGapO = %(VSET1)s(open);
    %(VTYPE)s vGapE = %(VSET1)s(gap);
    %(VTYPE)s vNegInf = %(VSET1)s(NEG_INF);
    %(VTYPE)s vZero = %(VSET0)s();
    %(INT)s score = NEG_INF;
    %(VTYPE)s vMaxH = vNegInf;
    const %(INT)s segLenXgap = -segLen*gap;
    %(VTYPE)s insert_mask = %(VCMPEQ)s(vZero,
            %(VSET)s(%(SCAN_INSERT_MASK)s));
    %(VTYPE)s vSegLenXgap1 = %(VSET1)s((segLen-1)*gap);
    %(VTYPE)s vSegLenXgap = %(VBLEND)s(vNegInf,
            %(VSET1)s(segLenXgap),
            insert_mask);
    %(SATURATION_CHECK_INIT)s
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif

    /* Generate query profile.
     * Rearrange query sequence & calculate the weight of match/mismatch.
     * Don't alias. */
    {
        %(INDEX)s index = 0;
        for (k=0; k<n; ++k) {
            for (i=0; i<segLen; ++i) {
                %(VTYPE)s_%(WIDTH)s_t t;
                j = i;
                for (segNum=0; segNum<segWidth; ++segNum) {
                    t.v[segNum] = j >= s1Len ? 0 : matrix->matrix[matrix->size*k+matrix->mapper[(unsigned char)s1[j]]];
                    j += segLen;
                }
                %(VSTORE)s(&pvP[index], t.m);
                ++index;
            }
        }
    }

    /* initialize H and E */
    {
        %(INDEX)s index = 0;
        for (i=0; i<segLen; ++i) {
            %(VTYPE)s_%(WIDTH)s_t h;
            %(VTYPE)s_%(WIDTH)s_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = 0;
                e.v[segNum] = NEG_INF;
            }
            %(VSTORE)s(&pvH[index], h.m);
            %(VSTORE)s(&pvE[index], e.m);
            ++index;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        %(VTYPE)s vE;
        %(VTYPE)s vHt;
        %(VTYPE)s vFt;
        %(VTYPE)s vH;
        %(VTYPE)s vHp;
        %(VTYPE)s *pvW;
        %(VTYPE)s vW;

        /* calculate E */
        /* calculate Ht */
        /* calculate Ft first pass */
        vHp = %(VLOAD)s(pvH+(segLen-1));
        vHp = %(VSHIFT)s(vHp, %(BYTES)s);
        pvW = pvP + matrix->mapper[(unsigned char)s2[j]]*segLen;
        vHt = vNegInf;
        vFt = vNegInf;
        for (i=0; i<segLen; ++i) {
            vH = %(VLOAD)s(pvH+i);
            vE = %(VLOAD)s(pvE+i);
            vW = %(VLOAD)s(pvW+i);
            vE = %(VMAX)s(
                    %(VSUB)s(vE, vGapE),
                    %(VSUB)s(vH, vGapO));
            vFt = %(VSUB)s(vFt, vGapE);
            vFt = %(VMAX)s(vFt, vHt);
            vHt = %(VMAX)s(
                    %(VADD)s(vHp, vW),
                    vE);
            %(VSTORE)s(pvE+i, vE);
            %(VSTORE)s(pvHt+i, vHt);
            vHp = vH;
        }

        /* adjust Ft before local prefix scan */
        vHt = %(VSHIFT)s(vHt, %(BYTES)s);
        vFt = %(VMAX)s(vFt,
                %(VSUB)s(vHt, vSegLenXgap1));
        /* local prefix scan */
        vFt = %(VBLEND)s(vNegInf, vFt, insert_mask);
        for (i=0; i<segWidth-1; ++i) {
                %(VTYPE)s vFtt = %(VROTATE)s(vFt, %(BYTES)s);
                vFtt = %(VADD)s(vFtt, vSegLenXgap);
                vFt = %(VMAX)s(vFt, vFtt);
        }
        vFt = %(VROTATE)s(vFt, %(BYTES)s);

        /* second Ft pass */
        /* calculate vH */
        for (i=0; i<segLen; ++i) {
            vFt = %(VSUB)s(vFt, vGapE);
            vFt = %(VMAX)s(vFt, vHt);
            vHt = %(VLOAD)s(pvHt+i);
            vH = %(VMAX)s(vHt, %(VSUB)s(vFt, vGapO));
            vH = %(VMAX)s(vH, vZero);
            %(VSTORE)s(pvH+i, vH);
            %(SATURATION_CHECK_MID)s
#ifdef PARASAIL_TABLE
            arr_store_si%(BITS)s(result->score_table, vH, i, segLen, j, s2Len);
#endif
            vMaxH = %(VMAX)s(vH, vMaxH);
        }
    }

    /* max in vec */
    for (j=0; j<segWidth; ++j) {
        %(INT)s value = (%(INT)s) %(VEXTRACT)s(vMaxH, %(LAST_POS)s);
        if (value > score) {
            score = value;
        }
        vMaxH = %(VSHIFT)s(vMaxH, %(BYTES)s);
    }

    %(SATURATION_CHECK_FINAL)s

    result->score = score;

    parasail_free(pvH);
    parasail_free(pvHt);
    parasail_free(pvE);
    parasail_free(pvP);

    return result;
}
