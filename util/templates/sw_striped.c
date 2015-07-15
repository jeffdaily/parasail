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
#include <string.h>

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
    parasail_profile_t *profile = parasail_profile_create_%(ISA)s_%(BITS)s_%(WIDTH)s(s1, s1Len, matrix);
    parasail_result_t *result = PNAME(profile, s2, s2Len, open, gap);
    parasail_profile_free(profile);
    return result;
}

parasail_result_t* PNAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    %(INDEX)s i = 0;
    %(INDEX)s j = 0;
    %(INDEX)s k = 0;
    %(INDEX)s end_query = 0;
    %(INDEX)s end_ref = 0;
    %(INDEX)s segNum = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const %(INDEX)s segWidth = %(LANES)s; /* number of values in vector unit */
    const %(INDEX)s segLen = (s1Len + segWidth - 1) / segWidth;
    %(VTYPE)s* const restrict vProfile = (%(VTYPE)s*)profile->profile%(WIDTH)s.score;
    %(VTYPE)s* restrict pvHStore = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHLoad =  parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvE = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvHMax = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s vGapO = %(VSET1)s(open);
    %(VTYPE)s vGapE = %(VSET1)s(gap);
    %(VTYPE)s vZero = %(VSET0)s();
    %(VTYPE)s vNegInf = %(VSET1)s(NEG_INF);
    %(INT)s score = NEG_INF;
    %(VTYPE)s vMaxH = vNegInf;
    %(VTYPE)s vMaxHUnit = vNegInf;
    %(SATURATION_CHECK_INIT)s
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(segLen*segWidth, s2Len);
    const %(INDEX)s offset = (s1Len - 1) %% segLen;
    const %(INDEX)s position = (segWidth - 1) - (s1Len - 1) / segLen;
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif

    /* initialize H and E */
    {
        %(INDEX)s index = 0;
        for (i=0; i<segLen; ++i) {
            %(VTYPE)s_%(WIDTH)s_t h;
            %(VTYPE)s_%(WIDTH)s_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                h.v[segNum] = 0;
                e.v[segNum] = -open;
            }
            %(VSTORE)s(&pvHStore[index], h.m);
            %(VSTORE)s(&pvE[index], e.m);
            ++index;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        %(VTYPE)s vE;
        %(VTYPE)s vF;
        %(VTYPE)s vH;
        const %(VTYPE)s* vP = NULL;
        %(VTYPE)s* pv = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be
         * corrected in the Lazy_F loop.  */
        vF = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = %(VSHIFT)s(pvHStore[segLen - 1], %(BYTES)s);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vH = %(VADD)s(vH, %(VLOAD)s(vP + i));
            vE = %(VLOAD)s(pvE + i);

            /* Get max from vH, vE and vF. */
            vH = %(VMAX)s(vH, vE);
            vH = %(VMAX)s(vH, vF);
            vH = %(VMAX)s(vH, vZero);
            /* Save vH values. */
            %(VSTORE)s(pvHStore + i, vH);
            %(SATURATION_CHECK_MID)s
#ifdef PARASAIL_TABLE
            arr_store_si%(BITS)s(result->score_table, vH, i, segLen, j, s2Len);
#endif
            vMaxH = %(VMAX)s(vH, vMaxH);

            /* Update vE value. */
            vH = %(VSUB)s(vH, vGapO);
            vE = %(VSUB)s(vE, vGapE);
            vE = %(VMAX)s(vE, vH);
            %(VSTORE)s(pvE + i, vE);

            /* Update vF value. */
            vF = %(VSUB)s(vF, vGapE);
            vF = %(VMAX)s(vF, vH);

            /* Load the next vH. */
            vH = %(VLOAD)s(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            vF = %(VSHIFT)s(vF, %(BYTES)s);
            for (i=0; i<segLen; ++i) {
                vH = %(VLOAD)s(pvHStore + i);
                vH = %(VMAX)s(vH,vF);
                %(VSTORE)s(pvHStore + i, vH);
                %(SATURATION_CHECK_MID)s
#ifdef PARASAIL_TABLE
                arr_store_si%(BITS)s(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vMaxH = %(VMAX)s(vH, vMaxH);
                vH = %(VSUB)s(vH, vGapO);
                vF = %(VSUB)s(vF, vGapE);
                if (! %(VMOVEMASK)s(%(VCMPGT)s(vF, vH))) goto end;
                /*vF = %(VMAX)s(vF, vH);*/
            }
        }
end:
        {
        }

        {
            %(VTYPE)s vCompare = %(VCMPGT)s(vMaxH, vMaxHUnit);
            if (%(VMOVEMASK)s(vCompare)) {
                score = %(VHMAX)s(vMaxH);
                vMaxHUnit = %(VSET1)s(score);
                end_ref = j;
                (void)memcpy(pvHMax, pvHStore, sizeof(%(VTYPE)s)*segLen);
            }
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = %(VLOAD)s(pvHStore + offset);
            for (k=0; k<position; ++k) {
                vH = %(VSHIFT)s(vH, %(BYTES)s);
            }
            result->score_row[j] = (%(INT)s) %(VEXTRACT)s (vH, %(LAST_POS)s);
        }
#endif
    }

    /* Trace the alignment ending position on read. */
    {
        %(INT)s *t = (%(INT)s*)pvHMax;
        %(INDEX)s column_len = segLen * segWidth;
        end_query = s1Len - 1;
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
        %(VTYPE)s vH = %(VLOAD)s(pvHStore+i);
        arr_store_col(result->score_col, vH, i, segLen);
    }
#endif

    %(SATURATION_CHECK_FINAL)s

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(pvHMax);
    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}
