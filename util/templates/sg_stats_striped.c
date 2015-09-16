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

#define FASTSTATS

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
#define INAME PNAME
#define STATIC
#else
#ifdef PARASAIL_ROWCOL
#define FNAME %(NAME_ROWCOL)s
#define PNAME %(PNAME_ROWCOL)s
#define INAME PNAME
#define STATIC
#else
#define FNAME %(NAME)s
#ifdef FASTSTATS
#define PNAME %(PNAME)s_internal
#define INAME %(PNAME)s
#define STATIC static
#else
#define PNAME %(PNAME)s
#define INAME PNAME
#define STATIC
#endif
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_%(ISA)s_%(BITS)s_%(WIDTH)s(s1, s1Len, matrix);
    parasail_result_t *result = INAME(profile, s2, s2Len, open, gap);
    parasail_profile_free(profile);
    return result;
}

STATIC parasail_result_t* PNAME(
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
    %(VTYPE)s* const restrict vProfile  = (%(VTYPE)s*)profile->profile%(WIDTH)s.score;
    %(VTYPE)s* const restrict vProfileM = (%(VTYPE)s*)profile->profile%(WIDTH)s.matches;
    %(VTYPE)s* const restrict vProfileS = (%(VTYPE)s*)profile->profile%(WIDTH)s.similar;
    %(VTYPE)s* restrict pvHStore        = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHLoad         = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHMStore       = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHMLoad        = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHSStore       = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHSLoad        = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHLStore       = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvHLLoad        = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvEStore        = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* restrict pvELoad         = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvEM      = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvES      = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s* const restrict pvEL      = parasail_memalign_%(VTYPE)s(%(ALIGNMENT)s, segLen);
    %(VTYPE)s vGapO = %(VSET1)s(open);
    %(VTYPE)s vGapE = %(VSET1)s(gap);
    %(VTYPE)s vNegInf = %(VSET1)s(NEG_INF);
    %(VTYPE)s vZero = %(VSET0)s();
    %(VTYPE)s vOne = %(VSET1)s(1);
    %(INT)s score = NEG_INF;
    %(INT)s matches = NEG_INF;
    %(INT)s similar = NEG_INF;
    %(INT)s length = NEG_INF;
    %(STATS_SATURATION_CHECK_INIT)s
    %(VTYPE)s vMaxH = vNegInf;
    %(VTYPE)s vMaxHM = vNegInf;
    %(VTYPE)s vMaxHS = vNegInf;
    %(VTYPE)s vMaxHL = vNegInf;
    %(VTYPE)s vPosMask = %(VCMPEQ)s(%(VSET1)s(position),
            %(VSET)s(%(POSITION_MASK)s));
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif

    parasail_memset_%(VTYPE)s(pvHMStore, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvHSStore, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvHLStore, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvEM, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvES, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvEL, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvHStore, vZero, segLen);
    parasail_memset_%(VTYPE)s(pvEStore, %(VSET1)s(-open), segLen);

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        %(VTYPE)s vE;
        %(VTYPE)s vEM;
        %(VTYPE)s vES;
        %(VTYPE)s vEL;
        %(VTYPE)s vF;
        %(VTYPE)s vFM;
        %(VTYPE)s vFS;
        %(VTYPE)s vFL;
        %(VTYPE)s vH;
        %(VTYPE)s vHM;
        %(VTYPE)s vHS;
        %(VTYPE)s vHL;
        const %(VTYPE)s* vP = NULL;
        const %(VTYPE)s* vPM = NULL;
        const %(VTYPE)s* vPS = NULL;
        %(VTYPE)s* pv = NULL;

        /* Initialize F value to neg inf.  Any errors to vH values will
         * be corrected in the Lazy_F loop.  */
        vF = vNegInf;
        vFM = vZero;
        vFS = vZero;
        vFL = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = %(VSHIFT)s(pvHStore[segLen - 1], %(BYTES)s);
        vHM = %(VSHIFT)s(pvHMStore[segLen - 1], %(BYTES)s);
        vHS = %(VSHIFT)s(pvHSStore[segLen - 1], %(BYTES)s);
        vHL = %(VSHIFT)s(pvHLStore[segLen - 1], %(BYTES)s);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;
        vPM = vProfileM + matrix->mapper[(unsigned char)s2[j]] * segLen;
        vPS = vProfileS + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        pv = pvHLoad;
        pvHLoad = pvHStore;
        pvHStore = pv;
        pv = pvHMLoad;
        pvHMLoad = pvHMStore;
        pvHMStore = pv;
        pv = pvHSLoad;
        pvHSLoad = pvHSStore;
        pvHSStore = pv;
        pv = pvHLLoad;
        pvHLLoad = pvHLStore;
        pvHLStore = pv;
        pv = pvELoad;
        pvELoad = pvEStore;
        pvEStore = pv;

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            %(VTYPE)s case1not;
            %(VTYPE)s case2not;
            %(VTYPE)s case2;
            %(VTYPE)s case3;

            vH = %(VADD)s(vH, %(VLOAD)s(vP + i));
            vE = %(VLOAD)s(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            case1not = %(VOR)s(
                    %(VCMPLT)s(vH,vF),%(VCMPLT)s(vH,vE));
            case2not = %(VCMPLT)s(vF,vE);
            case2 = %(VANDNOT)s(case2not,case1not);
            case3 = %(VAND)s(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = %(VMAX)s(vH, vE);
            vH = %(VMAX)s(vH, vF);
            /* Save vH values. */
            %(VSTORE)s(pvHStore + i, vH);

            /* calculate vM */
            vEM = %(VLOAD)s(pvEM + i);
            vHM = %(VBLEND)s(
                    %(VADD)s(vHM, %(VLOAD)s(vPM + i)),
                    %(VOR)s(
                        %(VAND)s(case2, vFM),
                        %(VAND)s(case3, vEM)),
                    case1not);
            %(VSTORE)s(pvHMStore + i, vHM);

            /* calculate vS */
            vES = %(VLOAD)s(pvES + i);
            vHS = %(VBLEND)s(
                    %(VADD)s(vHS, %(VLOAD)s(vPS + i)),
                    %(VOR)s(
                        %(VAND)s(case2, vFS),
                        %(VAND)s(case3, vES)),
                    case1not);
            %(VSTORE)s(pvHSStore + i, vHS);

            /* calculate vL */
            vEL = %(VLOAD)s(pvEL + i);
            vHL = %(VBLEND)s(
                    %(VADD)s(vHL, vOne),
                    %(VOR)s(
                        %(VAND)s(case2, %(VADD)s(vFL, vOne)),
                        %(VAND)s(case3, %(VADD)s(vEL, vOne))),
                    case1not);
            %(VSTORE)s(pvHLStore + i, vHL);
            %(STATS_SATURATION_CHECK_MID)s
#ifdef PARASAIL_TABLE
            arr_store_si%(BITS)s(result->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si%(BITS)s(result->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si%(BITS)s(result->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si%(BITS)s(result->score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = %(VSUB)s(vH, vGapO);
            vE = %(VSUB)s(vE, vGapE);
            vE = %(VMAX)s(vE, vH);
            %(VSTORE)s(pvEStore + i, vE);
            %(VSTORE)s(pvEM + i, vHM);
            %(VSTORE)s(pvES + i, vHS);
            %(VSTORE)s(pvEL + i, vHL);

            /* Update vF value. */
            vF = %(VSUB)s(vF, vGapE);
            vF = %(VMAX)s(vF, vH);
            vFM = vHM;
            vFS = vHS;
            vFL = vHL;

            /* Load the next vH. */
            vH = %(VLOAD)s(pvHLoad + i);
            vHM = %(VLOAD)s(pvHMLoad + i);
            vHS = %(VLOAD)s(pvHSLoad + i);
            vHL = %(VLOAD)s(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            %(VTYPE)s vHp = %(VSHIFT)s(pvHLoad[segLen - 1], %(BYTES)s);
            vF = %(VSHIFT)s(vF, %(BYTES)s);
            vF = %(VINSERT)s(vF, -open, 0);
            vFM = %(VSHIFT)s(vFM, %(BYTES)s);
            vFS = %(VSHIFT)s(vFS, %(BYTES)s);
            vFL = %(VSHIFT)s(vFL, %(BYTES)s);
            for (i=0; i<segLen; ++i) {
                %(VTYPE)s case1not;
                %(VTYPE)s case2not;
                %(VTYPE)s case2;

                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = %(VADD)s(vHp, %(VLOAD)s(vP + i));
                vE = %(VLOAD)s(pvELoad + i);
                case1not = %(VOR)s(
                        %(VCMPLT)s(vHp,vF),%(VCMPLT)s(vHp,vE));
                case2not = %(VCMPLT)s(vF,vE);
                case2 = %(VANDNOT)s(case2not,case1not);

                vHM = %(VLOAD)s(pvHMStore + i);
                vHM = %(VBLEND)s(vHM, vFM, case2);
                %(VSTORE)s(pvHMStore + i, vHM);
                %(VSTORE)s(pvEM + i, vHM);

                vHS = %(VLOAD)s(pvHSStore + i);
                vHS = %(VBLEND)s(vHS, vFS, case2);
                %(VSTORE)s(pvHSStore + i, vHS);
                %(VSTORE)s(pvES + i, vHS);

                vHL = %(VLOAD)s(pvHLStore + i);
                vHL = %(VBLEND)s(vHL, %(VADD)s(vFL,vOne), case2);
                %(VSTORE)s(pvHLStore + i, vHL);
                %(VSTORE)s(pvEL + i, vHL);

                vH = %(VLOAD)s(pvHStore + i);
                vH = %(VMAX)s(vH,vF);
                %(VSTORE)s(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
                arr_store_si%(BITS)s(result->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si%(BITS)s(result->similar_table, vHS, i, segLen, j, s2Len);
                arr_store_si%(BITS)s(result->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si%(BITS)s(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = %(VSUB)s(vH, vGapO);
                vF = %(VSUB)s(vF, vGapE);
                if (! %(VMOVEMASK)s(%(VCMPGT)s(vF, vH))) goto end;
                /*vF = %(VMAX)s(vF, vH);*/
                vFM = vHM;
                vFS = vHS;
                vFL = vHL;
                vHp = %(VLOAD)s(pvHLoad + i);
            }
        }
end:
        {
            /* extract vector containing last value from the column */
            %(VTYPE)s cond_max;
            vH = %(VLOAD)s(pvHStore + offset);
            vHM = %(VLOAD)s(pvHMStore + offset);
            vHS = %(VLOAD)s(pvHSStore + offset);
            vHL = %(VLOAD)s(pvHLStore + offset);
            cond_max = %(VCMPGT)s(vH, vMaxH);
            vMaxH = %(VBLEND)s(vMaxH, vH, cond_max);
            vMaxHM = %(VBLEND)s(vMaxHM, vHM, cond_max);
            vMaxHS = %(VBLEND)s(vMaxHS, vHS, cond_max);
            vMaxHL = %(VBLEND)s(vMaxHL, vHL, cond_max);
            if (%(VMOVEMASK)s(%(VAND)s(vPosMask, cond_max))) {
                end_ref = j;
                end_query = s1Len - 1;
            }
#ifdef PARASAIL_ROWCOL
            for (k=0; k<position; ++k) {
                vH = %(VSHIFT)s(vH, %(BYTES)s);
                vHM = %(VSHIFT)s(vHM, %(BYTES)s);
                vHS = %(VSHIFT)s(vHS, %(BYTES)s);
                vHL = %(VSHIFT)s(vHL, %(BYTES)s);
            }
            result->score_row[j] = (%(INT)s) %(VEXTRACT)s (vH, %(LAST_POS)s);
            result->matches_row[j] = (%(INT)s) %(VEXTRACT)s (vHM, %(LAST_POS)s);
            result->similar_row[j] = (%(INT)s) %(VEXTRACT)s (vHS, %(LAST_POS)s);
            result->length_row[j] = (%(INT)s) %(VEXTRACT)s (vHL, %(LAST_POS)s);
#endif
        }
    }

    {
        /* extract last value from the column */
        for (k=0; k<position; ++k) {
            vMaxH  = %(VSHIFT)s (vMaxH, %(BYTES)s);
            vMaxHM = %(VSHIFT)s (vMaxHM, %(BYTES)s);
            vMaxHS = %(VSHIFT)s (vMaxHS, %(BYTES)s);
            vMaxHL = %(VSHIFT)s (vMaxHL, %(BYTES)s);
        }
        score = (%(INT)s) %(VEXTRACT)s (vMaxH, %(LAST_POS)s);
        matches = (%(INT)s)%(VEXTRACT)s(vMaxHM, %(LAST_POS)s);
        similar = (%(INT)s)%(VEXTRACT)s(vMaxHS, %(LAST_POS)s);
        length = (%(INT)s)%(VEXTRACT)s(vMaxHL, %(LAST_POS)s);
    }

    /* max of last column */
    if (INT32_MAX == profile->stop || 0 == profile->stop)
    {
        %(INT)s score_last;
        vMaxH = vNegInf;

        if (0 == profile->stop) {
            /* ignore last row contributions */
            score = NEG_INF;
            matches = NEG_INF;
            similar = NEG_INF;
            length = NEG_INF;
            end_query = s1Len;
            end_ref = s2Len - 1;
        }

        for (i=0; i<segLen; ++i) {
            /* load the last stored values */
            %(VTYPE)s vH = %(VLOAD)s(pvHStore + i);
#ifdef PARASAIL_ROWCOL
            %(VTYPE)s vHM = %(VLOAD)s(pvHMStore + i);
            %(VTYPE)s vHS = %(VLOAD)s(pvHSStore + i);
            %(VTYPE)s vHL = %(VLOAD)s(pvHLStore + i);
            arr_store_col(result->score_col, vH, i, segLen);
            arr_store_col(result->matches_col, vHM, i, segLen);
            arr_store_col(result->similar_col, vHS, i, segLen);
            arr_store_col(result->length_col, vHL, i, segLen);
#endif
            vMaxH = %(VMAX)s(vH, vMaxH);
        }

        /* max in vec */
        score_last = %(VHMAX)s(vMaxH);
        if (score_last > score) {
            end_query = s1Len;
            end_ref = s2Len - 1;
            /* Trace the alignment ending position on read. */
            {
                %(INT)s *t = (%(INT)s*)pvHStore;
                %(INT)s *m = (%(INT)s*)pvHMStore;
                %(INT)s *s = (%(INT)s*)pvHSStore;
                %(INT)s *l = (%(INT)s*)pvHLStore;
                %(INDEX)s column_len = segLen * segWidth;
                for (i = 0; i<column_len; ++i, ++t, ++m, ++s, ++l) {
                    %(INDEX)s temp = i / segWidth + i %% segWidth * segLen;
                    if (temp < s1Len) {
                        if (*t > score || (*t == score && temp < end_query)) {
                            score = *t;
                            end_query = temp;
                            matches = *m;
                            similar = *s;
                            length = *l;
                        }
                    }
                }
            }
        }
    }

    %(STATS_SATURATION_CHECK_FINAL)s

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->matches = matches;
    result->similar = similar;
    result->length = length;

    parasail_free(pvEL);
    parasail_free(pvES);
    parasail_free(pvEM);
    parasail_free(pvELoad);
    parasail_free(pvEStore);
    parasail_free(pvHLLoad);
    parasail_free(pvHLStore);
    parasail_free(pvHSLoad);
    parasail_free(pvHSStore);
    parasail_free(pvHMLoad);
    parasail_free(pvHMStore);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

#ifdef FASTSTATS
#ifdef PARASAIL_TABLE
#else
#ifdef PARASAIL_ROWCOL
#else
#include <assert.h>
parasail_result_t* INAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    const char *s1 = profile->s1;
    const parasail_matrix_t *matrix = profile->matrix;

    /* find the end loc first with the faster implementation */
    parasail_result_t *result = %(PNAME_BASE)s(profile, s2, s2Len, open, gap);
    if (!result->saturated) {
        int s1Len_new = 0;
        int s2Len_new = 0;
        parasail_result_t *result_final = NULL;

        /* using the end loc, call the original stats function */
        s1Len_new = result->end_query+1;
        s2Len_new = result->end_ref+1;

        if (s1Len_new == profile->s1Len) {
            /* special 'stop' value tells stats function not to
             * consider last column results */
            int stop_save = profile->stop;
            ((parasail_profile_t*)profile)->stop = 1;
            result_final = PNAME(
                    profile, s2, s2Len_new, open, gap);
            ((parasail_profile_t*)profile)->stop = stop_save;
        }
        else {
            parasail_profile_t *profile_final = NULL;
            profile_final = parasail_profile_create_stats_%(ISA)s_%(BITS)s_%(WIDTH)s(
                    s1, s1Len_new, matrix);
            /* special 'stop' value tells stats function not to
             * consider last row results */
            profile_final->stop = 0;
            result_final = PNAME(
                    profile_final, s2, s2Len_new, open, gap);

            parasail_profile_free(profile_final);
        }

        parasail_result_free(result);

        /* correct the end locations before returning */
        result_final->end_query = s1Len_new-1;
        result_final->end_ref = s2Len_new-1;
        return result_final;
    }
    else {
        return result;
    }
}
#endif
#endif
#endif

