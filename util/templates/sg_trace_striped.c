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

#define SWAP(A,B) { %(VTYPE)s* tmp = A; A = B; B = tmp; }

#define NEG_INF %(NEG_INF)s
%(FIXES)s

static inline void arr_store(
        %(VTYPE)s *array,
        %(VTYPE)s vH,
        %(INDEX)s t,
        %(INDEX)s seglen,
        %(INDEX)s d)
{
    %(VSTORE)s(array + (d*seglen+t), vH);
}

#define FNAME %(NAME_TRACE)s
#define PNAME %(PNAME_TRACE)s

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
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const %(INDEX)s segWidth = %(LANES)s; /* number of values in vector unit */
    const %(INDEX)s segLen = (s1Len + segWidth - 1) / segWidth;
    const %(INDEX)s offset = (s1Len - 1) %% segLen;
    const %(INDEX)s position = (segWidth - 1) - (s1Len - 1) / segLen;
    %(VTYPE)s* const restrict vProfile = (%(VTYPE)s*)profile->profile%(WIDTH)s.score;
    %(VTYPE)s* restrict pvHStore = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHLoad = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvE = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvEaStore = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvEaLoad = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvHT = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s vGapO = %(VSET1)s(open);
    %(VTYPE)s vGapE = %(VSET1)s(gap);
    %(VTYPE)s vNegInf = %(VSET1)s(NEG_INF);
    %(INT)s score = NEG_INF;
    %(VTYPE)s vMaxH = vNegInf;
    %(VTYPE)s vPosMask = %(VCMPEQ)s(%(VSET1)s(position),
            %(VSET)s(%(POSITION_MASK)s));
    %(SATURATION_CHECK_INIT)s
    parasail_result_t *result = parasail_result_new_trace(segLen*segWidth, s2Len, %(BYTES)s);
    %(VTYPE)s vTIns  = %(VSET1)s(PARASAIL_INS);
    %(VTYPE)s vTDel  = %(VSET1)s(PARASAIL_DEL);
    %(VTYPE)s vTDiag = %(VSET1)s(PARASAIL_DIAG);

    /* initialize H and E */
    parasail_memset_%(VTYPE)s(pvHStore, %(VSET1)s(0), segLen);
    parasail_memset_%(VTYPE)s(pvE, %(VSET1)s(-open), segLen);
    parasail_memset_%(VTYPE)s(pvEaStore, %(VSET1)s(-open), segLen);

    for (i=0; i<segLen; ++i) {
        arr_store(result->trace_ins_table, vTDiag, i, segLen, 0);
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        %(VTYPE)s vEF_opn;
        %(VTYPE)s vE;
        %(VTYPE)s vE_ext;
        %(VTYPE)s vF;
        %(VTYPE)s vF_ext;
        %(VTYPE)s vFa;
        %(VTYPE)s vFa_ext;
        %(VTYPE)s vH;
        %(VTYPE)s vH_dag;
        const %(VTYPE)s* vP = NULL;

        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop. */
        vF = vNegInf;

        /* load final segment of pvHStore and shift left by %(BYTES)s bytes */
        vH = %(VLOAD)s(&pvHStore[segLen - 1]);
        vH = %(VSHIFT)s(vH, %(BYTES)s);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        SWAP(pvHLoad, pvHStore)
        SWAP(pvEaLoad, pvEaStore)

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vE = %(VLOAD)s(pvE + i);

            /* Get max from vH, vE and vF. */
            vH_dag = %(VADD)s(vH, %(VLOAD)s(vP + i));
            vH = %(VMAX)s(vH_dag, vE);
            vH = %(VMAX)s(vH, vF);
            /* Save vH values. */
            %(VSTORE)s(pvHStore + i, vH);
            %(SATURATION_CHECK_MID)s

            {
                %(VTYPE)s case1 = %(VCMPEQ)s(vH, vH_dag);
                %(VTYPE)s case2 = %(VCMPEQ)s(vH, vF);
                %(VTYPE)s vT = %(VBLEND)s(
                        %(VBLEND)s(vTIns, vTDel, case2),
                        vTDiag, case1);
                %(VSTORE)s(pvHT + i, vT);
                arr_store(result->trace_table, vT, i, segLen, j);
            }

            vEF_opn = %(VSUB)s(vH, vGapO);

            /* Update vE value. */
            vE_ext = %(VSUB)s(vE, vGapE);
            vE = %(VMAX)s(vEF_opn, vE_ext);
            %(VSTORE)s(pvE + i, vE);
            {
                %(VTYPE)s vEa = %(VLOAD)s(pvEaLoad + i);
                %(VTYPE)s vEa_ext = %(VSUB)s(vEa, vGapE);
                vEa = %(VMAX)s(vEF_opn, vEa_ext);
                %(VSTORE)s(pvEaStore + i, vEa);
                if (j+1<s2Len) {
                    %(VTYPE)s cond = %(VCMPGT)s(vEF_opn, vEa_ext);
                    %(VTYPE)s vT = %(VBLEND)s(vTIns, vTDiag, cond);
                    arr_store(result->trace_ins_table, vT, i, segLen, j+1);
                }
            }

            /* Update vF value. */
            vF_ext = %(VSUB)s(vF, vGapE);
            vF = %(VMAX)s(vEF_opn, vF_ext);
            {
                %(VTYPE)s cond = %(VCMPGT)s(vEF_opn, vF_ext);
                %(VTYPE)s vT = %(VBLEND)s(vTDel, vTDiag, cond);
                if (i+1<segLen) {
                    arr_store(result->trace_del_table, vT, i+1, segLen, j);
                }
            }

            /* Load the next vH. */
            vH = %(VLOAD)s(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        vFa_ext = vF_ext;
        vFa = vF;
        for (k=0; k<segWidth; ++k) {
            %(VTYPE)s vHp = %(VLOAD)s(&pvHLoad[segLen - 1]);
            vHp = %(VSHIFT)s(vHp, %(BYTES)s);
            vEF_opn = %(VSHIFT)s(vEF_opn, %(BYTES)s);
            vEF_opn = %(VINSERT)s(vEF_opn, -open, 0);
            vF_ext = %(VSHIFT)s(vF_ext, %(BYTES)s);
            vF_ext = %(VINSERT)s(vF_ext, NEG_INF, 0);
            vF = %(VSHIFT)s(vF, %(BYTES)s);
            vF = %(VINSERT)s(vF, -open, 0);
            vFa_ext = %(VSHIFT)s(vFa_ext, %(BYTES)s);
            vFa_ext = %(VINSERT)s(vFa_ext, NEG_INF, 0);
            vFa = %(VSHIFT)s(vFa, %(BYTES)s);
            vFa = %(VINSERT)s(vFa, -open, 0);
            for (i=0; i<segLen; ++i) {
                vH = %(VLOAD)s(pvHStore + i);
                vH = %(VMAX)s(vH,vF);
                %(VSTORE)s(pvHStore + i, vH);
                %(SATURATION_CHECK_MID)s
                {
                    %(VTYPE)s vT;
                    %(VTYPE)s case1;
                    %(VTYPE)s case2;
                    %(VTYPE)s cond;
                    vHp = %(VADD)s(vHp, %(VLOAD)s(vP + i));
                    case1 = %(VCMPEQ)s(vH, vHp);
                    case2 = %(VCMPEQ)s(vH, vF);
                    cond = %(VANDNOT)s(case1,case2);
                    vT = %(VLOAD)s(pvHT + i);
                    vT = %(VBLEND)s(vT, vTDel, cond);
                    %(VSTORE)s(pvHT + i, vT);
                    arr_store(result->trace_table, vT, i, segLen, j);
                }
                /* Update vF value. */
                {
                    %(VTYPE)s cond = %(VCMPGT)s(vEF_opn, vFa_ext);
                    %(VTYPE)s vT = %(VBLEND)s(vTDel, vTDiag, cond);
                    arr_store(result->trace_del_table, vT, i, segLen, j);
                }
                vEF_opn = %(VSUB)s(vH, vGapO);
                vF_ext = %(VSUB)s(vF, vGapE);
                {
                    %(VTYPE)s vT;
                    %(VTYPE)s cond;
                    %(VTYPE)s vEa = %(VLOAD)s(pvEaLoad + i);
                    %(VTYPE)s vEa_ext = %(VSUB)s(vEa, vGapE);
                    vEa = %(VMAX)s(vEF_opn, vEa_ext);
                    %(VSTORE)s(pvEaStore + i, vEa);
                    cond = %(VCMPGT)s(vEF_opn, vEa_ext);
                    vT = %(VBLEND)s(vTIns, vTDiag, cond);
                    if (j+1<s2Len) {
                        arr_store(result->trace_ins_table, vT, i, segLen, j+1);
                    }
                }
                if (! %(VMOVEMASK)s(
                            %(VOR)s(
                                %(VCMPGT)s(vF_ext, vEF_opn),
                                %(VCMPEQ)s(vF_ext, vEF_opn))))
                    goto end;
                /*vF = %(VMAX)s(vEF_opn, vF_ext);*/
                vF = vF_ext;
                vFa_ext = %(VSUB)s(vFa, vGapE);
                vFa = %(VMAX)s(vEF_opn, vFa_ext);
                vHp = %(VLOAD)s(pvHLoad + i);
            }
        }
end:
        {
            /* extract vector containing last value from the column */
            %(VTYPE)s vCompare;
            vH = %(VLOAD)s(pvHStore + offset);
            vCompare = %(VAND)s(vPosMask, %(VCMPGT)s(vH, vMaxH));
            vMaxH = %(VMAX)s(vH, vMaxH);
            if (%(VMOVEMASK)s(vCompare)) {
                end_ref = j;
                end_query = s1Len - 1;
            }
        }
    }

    /* max last value from all columns */
    {
        for (k=0; k<position; ++k) {
            vMaxH = %(VSHIFT)s(vMaxH, %(BYTES)s);
        }
        score = (%(INT)s) %(VEXTRACT)s(vMaxH, %(LAST_POS)s);
    }

    /* max of last column */
    {
        %(INT)s score_last;
        vMaxH = vNegInf;

        for (i=0; i<segLen; ++i) {
            %(VTYPE)s vH = %(VLOAD)s(pvHStore + i);
            vMaxH = %(VMAX)s(vH, vMaxH);
        }

        /* max in vec */
        score_last = %(VHMAX)s(vMaxH);
        if (score_last > score) {
            score = score_last;
            end_ref = s2Len - 1;
            end_query = s1Len;
            /* Trace the alignment ending position on read. */
            {
                %(INT)s *t = (%(INT)s*)pvHStore;
                %(INDEX)s column_len = segLen * segWidth;
                for (i = 0; i<column_len; ++i, ++t) {
                    if (*t == score) {
                        %(INDEX)s temp = i / segWidth + i %% segWidth * segLen;
                        if (temp < end_query) {
                            end_query = temp;
                        }
                    }
                }
            }
        }
    }

    %(SATURATION_CHECK_FINAL)s

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->flag = PARASAIL_FLAG_SG | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_%(WIDTH)s | PARASAIL_FLAG_LANES_%(LANES)s;

    parasail_free(pvHT);
    parasail_free(pvEaLoad);
    parasail_free(pvEaStore);
    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

