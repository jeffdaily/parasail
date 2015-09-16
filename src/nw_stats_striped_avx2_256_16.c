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

#include <immintrin.h>

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_avx.h"

#define NEG_INF (INT16_MIN/(int16_t)(2))

#define _mm256_cmplt_epi16_rpl(a,b) _mm256_cmpgt_epi16(b,a)

#if HAVE_AVX2_MM256_INSERT_EPI16
#define _mm256_insert_epi16_rpl _mm256_insert_epi16
#else
static inline __m256i _mm256_insert_epi16_rpl(__m256i a, int16_t i, int imm) {
    __m256i_16_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI16
#define _mm256_extract_epi16_rpl _mm256_extract_epi16
#else
static inline int16_t _mm256_extract_epi16_rpl(__m256i a, int imm) {
    __m256i_16_t A;
    A.m = a;
    return A.v[imm];
}
#endif

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)


#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[( 0*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  0);
    array[( 1*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  1);
    array[( 2*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  2);
    array[( 3*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  3);
    array[( 4*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  4);
    array[( 5*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  5);
    array[( 6*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  6);
    array[( 7*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  7);
    array[( 8*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  8);
    array[( 9*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH,  9);
    array[(10*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 10);
    array[(11*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 11);
    array[(12*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 12);
    array[(13*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 13);
    array[(14*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 14);
    array[(15*seglen+t)*dlen + d] = (int16_t)_mm256_extract_epi16_rpl(vH, 15);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m256i vH,
        int32_t t,
        int32_t seglen)
{
    col[ 0*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  0);
    col[ 1*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  1);
    col[ 2*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  2);
    col[ 3*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  3);
    col[ 4*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  4);
    col[ 5*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  5);
    col[ 6*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  6);
    col[ 7*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  7);
    col[ 8*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  8);
    col[ 9*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH,  9);
    col[10*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 10);
    col[11*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 11);
    col[12*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 12);
    col[13*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 13);
    col[14*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 14);
    col[15*seglen+t] = (int16_t)_mm256_extract_epi16_rpl(vH, 15);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_nw_stats_table_striped_avx2_256_16
#define PNAME parasail_nw_stats_table_striped_profile_avx2_256_16
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_stats_rowcol_striped_avx2_256_16
#define PNAME parasail_nw_stats_rowcol_striped_profile_avx2_256_16
#else
#define FNAME parasail_nw_stats_striped_avx2_256_16
#define PNAME parasail_nw_stats_striped_profile_avx2_256_16
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_stats_avx_256_16(s1, s1Len, matrix);
    parasail_result_t *result = PNAME(profile, s2, s2Len, open, gap);
    parasail_profile_free(profile);
    return result;
}

parasail_result_t* PNAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int32_t segNum = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m256i* const restrict vProfile  = (__m256i*)profile->profile16.score;
    __m256i* const restrict vProfileM = (__m256i*)profile->profile16.matches;
    __m256i* const restrict vProfileS = (__m256i*)profile->profile16.similar;
    __m256i* restrict pvHStore        = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLoad         = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHMStore       = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHMLoad        = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHSStore       = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHSLoad        = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLStore       = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvHLLoad        = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvEStore        = parasail_memalign___m256i(32, segLen);
    __m256i* restrict pvELoad         = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvEM      = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvES      = parasail_memalign___m256i(32, segLen);
    __m256i* const restrict pvEL      = parasail_memalign___m256i(32, segLen);
    int16_t* const restrict boundary  = parasail_memalign_int16_t(32, s2Len+1);
    __m256i vGapO = _mm256_set1_epi16(open);
    __m256i vGapE = _mm256_set1_epi16(gap);
    __m256i vNegInf = _mm256_set1_epi16(NEG_INF);
    __m256i vZero = _mm256_setzero_si256();
    __m256i vOne = _mm256_set1_epi16(1);
    int16_t score = NEG_INF;
    int16_t matches = NEG_INF;
    int16_t similar = NEG_INF;
    int16_t length = NEG_INF;
    
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(segLen*segWidth, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif

    parasail_memset___m256i(pvHMStore, vZero, segLen);
    parasail_memset___m256i(pvHSStore, vZero, segLen);
    parasail_memset___m256i(pvHLStore, vZero, segLen);
    parasail_memset___m256i(pvEM, vZero, segLen);
    parasail_memset___m256i(pvES, vZero, segLen);
    parasail_memset___m256i(pvEL, vZero, segLen);

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            __m256i_16_t h;
            __m256i_16_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = -open-gap*(segNum*segLen+i);
                h.v[segNum] = tmp < INT16_MIN ? INT16_MIN : tmp;
                tmp = tmp - open;
                e.v[segNum] = tmp < INT16_MIN ? INT16_MIN : tmp;
            }
            _mm256_store_si256(&pvHStore[index], h.m);
            _mm256_store_si256(&pvEStore[index], e.m);
            ++index;
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            int64_t tmp = -open-gap*(i-1);
            boundary[i] = tmp < INT16_MIN ? INT16_MIN : tmp;
        }
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m256i vE;
        __m256i vEM;
        __m256i vES;
        __m256i vEL;
        __m256i vF;
        __m256i vFM;
        __m256i vFS;
        __m256i vFL;
        __m256i vH;
        __m256i vHM;
        __m256i vHS;
        __m256i vHL;
        const __m256i* vP = NULL;
        const __m256i* vPM = NULL;
        const __m256i* vPS = NULL;
        __m256i* pv = NULL;

        /* Initialize F value to neg inf.  Any errors to vH values will
         * be corrected in the Lazy_F loop.  */
        vF = vNegInf;
        vFM = vZero;
        vFS = vZero;
        vFL = vZero;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = _mm256_slli_si256_rpl(pvHStore[segLen - 1], 2);
        vHM = _mm256_slli_si256_rpl(pvHMStore[segLen - 1], 2);
        vHS = _mm256_slli_si256_rpl(pvHSStore[segLen - 1], 2);
        vHL = _mm256_slli_si256_rpl(pvHLStore[segLen - 1], 2);

        /* insert upper boundary condition */
        vH = _mm256_insert_epi16_rpl(vH, boundary[j], 0);

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
            __m256i case1not;
            __m256i case2not;
            __m256i case2;
            __m256i case3;

            vH = _mm256_add_epi16(vH, _mm256_load_si256(vP + i));
            vE = _mm256_load_si256(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            case1not = _mm256_or_si256(
                    _mm256_cmplt_epi16_rpl(vH,vF),_mm256_cmplt_epi16_rpl(vH,vE));
            case2not = _mm256_cmplt_epi16_rpl(vF,vE);
            case2 = _mm256_andnot_si256(case2not,case1not);
            case3 = _mm256_and_si256(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm256_max_epi16(vH, vE);
            vH = _mm256_max_epi16(vH, vF);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);

            /* calculate vM */
            vEM = _mm256_load_si256(pvEM + i);
            vHM = _mm256_blendv_epi8(
                    _mm256_add_epi16(vHM, _mm256_load_si256(vPM + i)),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, vFM),
                        _mm256_and_si256(case3, vEM)),
                    case1not);
            _mm256_store_si256(pvHMStore + i, vHM);

            /* calculate vS */
            vES = _mm256_load_si256(pvES + i);
            vHS = _mm256_blendv_epi8(
                    _mm256_add_epi16(vHS, _mm256_load_si256(vPS + i)),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, vFS),
                        _mm256_and_si256(case3, vES)),
                    case1not);
            _mm256_store_si256(pvHSStore + i, vHS);

            /* calculate vL */
            vEL = _mm256_load_si256(pvEL + i);
            vHL = _mm256_blendv_epi8(
                    _mm256_add_epi16(vHL, vOne),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, _mm256_add_epi16(vFL, vOne)),
                        _mm256_and_si256(case3, _mm256_add_epi16(vEL, vOne))),
                    case1not);
            _mm256_store_si256(pvHLStore + i, vHL);
            
#ifdef PARASAIL_TABLE
            arr_store_si256(result->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si256(result->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si256(result->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = _mm256_sub_epi16(vH, vGapO);
            vE = _mm256_sub_epi16(vE, vGapE);
            vE = _mm256_max_epi16(vE, vH);
            _mm256_store_si256(pvEStore + i, vE);
            _mm256_store_si256(pvEM + i, vHM);
            _mm256_store_si256(pvES + i, vHS);
            _mm256_store_si256(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm256_sub_epi16(vF, vGapE);
            vF = _mm256_max_epi16(vF, vH);
            vFM = vHM;
            vFS = vHS;
            vFL = vHL;

            /* Load the next vH. */
            vH = _mm256_load_si256(pvHLoad + i);
            vHM = _mm256_load_si256(pvHMLoad + i);
            vHS = _mm256_load_si256(pvHSLoad + i);
            vHL = _mm256_load_si256(pvHLLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
            __m256i vHp = _mm256_slli_si256_rpl(pvHLoad[segLen - 1], 2);
            int64_t tmp = boundary[j+1]-open;
            int16_t tmp2 = tmp < INT16_MIN ? INT16_MIN : tmp;
            vF = _mm256_slli_si256_rpl(vF, 2);
            vF = _mm256_insert_epi16_rpl(vF, tmp2, 0);
            vFM = _mm256_slli_si256_rpl(vFM, 2);
            vFS = _mm256_slli_si256_rpl(vFS, 2);
            vFL = _mm256_slli_si256_rpl(vFL, 2);
            for (i=0; i<segLen; ++i) {
                __m256i case1not;
                __m256i case2not;
                __m256i case2;

                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm256_add_epi16(vHp, _mm256_load_si256(vP + i));
                vE = _mm256_load_si256(pvELoad + i);
                case1not = _mm256_or_si256(
                        _mm256_cmplt_epi16_rpl(vHp,vF),_mm256_cmplt_epi16_rpl(vHp,vE));
                case2not = _mm256_cmplt_epi16_rpl(vF,vE);
                case2 = _mm256_andnot_si256(case2not,case1not);

                vHM = _mm256_load_si256(pvHMStore + i);
                vHM = _mm256_blendv_epi8(vHM, vFM, case2);
                _mm256_store_si256(pvHMStore + i, vHM);
                _mm256_store_si256(pvEM + i, vHM);

                vHS = _mm256_load_si256(pvHSStore + i);
                vHS = _mm256_blendv_epi8(vHS, vFS, case2);
                _mm256_store_si256(pvHSStore + i, vHS);
                _mm256_store_si256(pvES + i, vHS);

                vHL = _mm256_load_si256(pvHLStore + i);
                vHL = _mm256_blendv_epi8(vHL, _mm256_add_epi16(vFL,vOne), case2);
                _mm256_store_si256(pvHLStore + i, vHL);
                _mm256_store_si256(pvEL + i, vHL);

                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi16(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
                
#ifdef PARASAIL_TABLE
                arr_store_si256(result->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si256(result->similar_table, vHS, i, segLen, j, s2Len);
                arr_store_si256(result->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm256_sub_epi16(vH, vGapO);
                vF = _mm256_sub_epi16(vF, vGapE);
                if (! _mm256_movemask_epi8(_mm256_cmpgt_epi16(vF, vH))) goto end;
                /*vF = _mm256_max_epi16(vF, vH);*/
                vFM = vHM;
                vFS = vHS;
                vFL = vHL;
                vHp = _mm256_load_si256(pvHLoad + i);
            }
        }
end:
        {
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = _mm256_load_si256(pvHStore + offset);
            vHM = _mm256_load_si256(pvHMStore + offset);
            vHS = _mm256_load_si256(pvHSStore + offset);
            vHL = _mm256_load_si256(pvHLStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm256_slli_si256_rpl (vH, 2);
                vHM = _mm256_slli_si256_rpl (vHM, 2);
                vHS = _mm256_slli_si256_rpl (vHS, 2);
                vHL = _mm256_slli_si256_rpl (vHL, 2);
            }
            result->score_row[j] = (int16_t) _mm256_extract_epi16_rpl (vH, 15);
            result->matches_row[j] = (int16_t) _mm256_extract_epi16_rpl (vHM, 15);
            result->similar_row[j] = (int16_t) _mm256_extract_epi16_rpl (vHS, 15);
            result->length_row[j] = (int16_t) _mm256_extract_epi16_rpl (vHL, 15);
        }
#endif
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        __m256i vH = _mm256_load_si256(pvHStore+i);
        __m256i vHM = _mm256_load_si256(pvHMStore+i);
        __m256i vHS = _mm256_load_si256(pvHSStore+i);
        __m256i vHL = _mm256_load_si256(pvHLStore+i);
        arr_store_col(result->score_col, vH, i, segLen);
        arr_store_col(result->matches_col, vHM, i, segLen);
        arr_store_col(result->similar_col, vHS, i, segLen);
        arr_store_col(result->length_col, vHL, i, segLen);
    }
#endif

    /* extract last value from the last column */
    {
        __m256i vH = _mm256_load_si256(pvHStore + offset);
        __m256i vHM = _mm256_load_si256(pvHMStore + offset);
        __m256i vHS = _mm256_load_si256(pvHSStore + offset);
        __m256i vHL = _mm256_load_si256(pvHLStore + offset);
        for (k=0; k<position; ++k) {
            vH = _mm256_slli_si256_rpl (vH, 2);
            vHM = _mm256_slli_si256_rpl (vHM, 2);
            vHS = _mm256_slli_si256_rpl (vHS, 2);
            vHL = _mm256_slli_si256_rpl (vHL, 2);
        }
        score = (int16_t) _mm256_extract_epi16_rpl (vH, 15);
        matches = (int16_t) _mm256_extract_epi16_rpl (vHM, 15);
        similar = (int16_t) _mm256_extract_epi16_rpl (vHS, 15);
        length = (int16_t) _mm256_extract_epi16_rpl (vHL, 15);
    }

    

    result->score = score;
    result->matches = matches;
    result->similar = similar;
    result->length = length;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(boundary);
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


