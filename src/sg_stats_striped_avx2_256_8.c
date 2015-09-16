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

#define FASTSTATS

#define NEG_INF INT8_MIN

#if HAVE_AVX2_MM256_INSERT_EPI8
#define _mm256_insert_epi8_rpl _mm256_insert_epi8
#else
static inline __m256i _mm256_insert_epi8_rpl(__m256i a, int8_t i, int imm) {
    __m256i_8_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI8
#define _mm256_extract_epi8_rpl _mm256_extract_epi8
#else
static inline int8_t _mm256_extract_epi8_rpl(__m256i a, int imm) {
    __m256i_8_t A;
    A.m = a;
    return A.v[imm];
}
#endif

#define _mm256_cmplt_epi8_rpl(a,b) _mm256_cmpgt_epi8(b,a)

#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)

static inline int8_t _mm256_hmax_epi8_rpl(__m256i a) {
    a = _mm256_max_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 8));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 4));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 2));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 1));
    return _mm256_extract_epi8_rpl(a, 31);
}


#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[( 0*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  0);
    array[( 1*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  1);
    array[( 2*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  2);
    array[( 3*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  3);
    array[( 4*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  4);
    array[( 5*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  5);
    array[( 6*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  6);
    array[( 7*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  7);
    array[( 8*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  8);
    array[( 9*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH,  9);
    array[(10*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 10);
    array[(11*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 11);
    array[(12*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 12);
    array[(13*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 13);
    array[(14*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 14);
    array[(15*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 15);
    array[(16*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 16);
    array[(17*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 17);
    array[(18*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 18);
    array[(19*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 19);
    array[(20*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 20);
    array[(21*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 21);
    array[(22*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 22);
    array[(23*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 23);
    array[(24*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 24);
    array[(25*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 25);
    array[(26*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 26);
    array[(27*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 27);
    array[(28*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 28);
    array[(29*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 29);
    array[(30*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 30);
    array[(31*seglen+t)*dlen + d] = (int8_t)_mm256_extract_epi8_rpl(vH, 31);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m256i vH,
        int32_t t,
        int32_t seglen)
{
    col[ 0*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  0);
    col[ 1*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  1);
    col[ 2*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  2);
    col[ 3*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  3);
    col[ 4*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  4);
    col[ 5*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  5);
    col[ 6*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  6);
    col[ 7*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  7);
    col[ 8*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  8);
    col[ 9*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH,  9);
    col[10*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 10);
    col[11*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 11);
    col[12*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 12);
    col[13*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 13);
    col[14*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 14);
    col[15*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 15);
    col[16*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 16);
    col[17*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 17);
    col[18*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 18);
    col[19*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 19);
    col[20*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 20);
    col[21*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 21);
    col[22*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 22);
    col[23*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 23);
    col[24*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 24);
    col[25*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 25);
    col[26*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 26);
    col[27*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 27);
    col[28*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 28);
    col[29*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 29);
    col[30*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 30);
    col[31*seglen+t] = (int8_t)_mm256_extract_epi8_rpl(vH, 31);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sg_stats_table_striped_avx2_256_8
#define PNAME parasail_sg_stats_table_striped_profile_avx2_256_8
#define INAME PNAME
#define STATIC
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_stats_rowcol_striped_avx2_256_8
#define PNAME parasail_sg_stats_rowcol_striped_profile_avx2_256_8
#define INAME PNAME
#define STATIC
#else
#define FNAME parasail_sg_stats_striped_avx2_256_8
#ifdef FASTSTATS
#define PNAME parasail_sg_stats_striped_profile_avx2_256_8_internal
#define INAME parasail_sg_stats_striped_profile_avx2_256_8
#define STATIC static
#else
#define PNAME parasail_sg_stats_striped_profile_avx2_256_8
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
    parasail_profile_t *profile = parasail_profile_create_stats_avx_256_8(s1, s1Len, matrix);
    parasail_result_t *result = INAME(profile, s2, s2Len, open, gap);
    parasail_profile_free(profile);
    return result;
}

STATIC parasail_result_t* PNAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const int32_t segWidth = 32; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
    __m256i* const restrict vProfile  = (__m256i*)profile->profile8.score;
    __m256i* const restrict vProfileM = (__m256i*)profile->profile8.matches;
    __m256i* const restrict vProfileS = (__m256i*)profile->profile8.similar;
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
    __m256i vGapO = _mm256_set1_epi8(open);
    __m256i vGapE = _mm256_set1_epi8(gap);
    __m256i vNegInf = _mm256_set1_epi8(NEG_INF);
    __m256i vZero = _mm256_setzero_si256();
    __m256i vOne = _mm256_set1_epi8(1);
    int8_t score = NEG_INF;
    int8_t matches = NEG_INF;
    int8_t similar = NEG_INF;
    int8_t length = NEG_INF;
    __m256i vNegLimit = _mm256_set1_epi8(INT8_MIN);
    __m256i vPosLimit = _mm256_set1_epi8(INT8_MAX);
    __m256i vSaturationCheckMin = vPosLimit;
    __m256i vSaturationCheckMax = vNegLimit;
    __m256i vMaxH = vNegInf;
    __m256i vMaxHM = vNegInf;
    __m256i vMaxHS = vNegInf;
    __m256i vMaxHL = vNegInf;
    __m256i vPosMask = _mm256_cmpeq_epi8(_mm256_set1_epi8(position),
            _mm256_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31));
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
    parasail_memset___m256i(pvHStore, vZero, segLen);
    parasail_memset___m256i(pvEStore, _mm256_set1_epi8(-open), segLen);

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
        vH = _mm256_slli_si256_rpl(pvHStore[segLen - 1], 1);
        vHM = _mm256_slli_si256_rpl(pvHMStore[segLen - 1], 1);
        vHS = _mm256_slli_si256_rpl(pvHSStore[segLen - 1], 1);
        vHL = _mm256_slli_si256_rpl(pvHLStore[segLen - 1], 1);

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

            vH = _mm256_adds_epi8(vH, _mm256_load_si256(vP + i));
            vE = _mm256_load_si256(pvELoad + i);

            /* determine which direction of length and match to
             * propagate, before vH is finished calculating */
            case1not = _mm256_or_si256(
                    _mm256_cmplt_epi8_rpl(vH,vF),_mm256_cmplt_epi8_rpl(vH,vE));
            case2not = _mm256_cmplt_epi8_rpl(vF,vE);
            case2 = _mm256_andnot_si256(case2not,case1not);
            case3 = _mm256_and_si256(case1not,case2not);

            /* Get max from vH, vE and vF. */
            vH = _mm256_max_epi8(vH, vE);
            vH = _mm256_max_epi8(vH, vF);
            /* Save vH values. */
            _mm256_store_si256(pvHStore + i, vH);

            /* calculate vM */
            vEM = _mm256_load_si256(pvEM + i);
            vHM = _mm256_blendv_epi8(
                    _mm256_adds_epi8(vHM, _mm256_load_si256(vPM + i)),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, vFM),
                        _mm256_and_si256(case3, vEM)),
                    case1not);
            _mm256_store_si256(pvHMStore + i, vHM);

            /* calculate vS */
            vES = _mm256_load_si256(pvES + i);
            vHS = _mm256_blendv_epi8(
                    _mm256_adds_epi8(vHS, _mm256_load_si256(vPS + i)),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, vFS),
                        _mm256_and_si256(case3, vES)),
                    case1not);
            _mm256_store_si256(pvHSStore + i, vHS);

            /* calculate vL */
            vEL = _mm256_load_si256(pvEL + i);
            vHL = _mm256_blendv_epi8(
                    _mm256_adds_epi8(vHL, vOne),
                    _mm256_or_si256(
                        _mm256_and_si256(case2, _mm256_adds_epi8(vFL, vOne)),
                        _mm256_and_si256(case3, _mm256_adds_epi8(vEL, vOne))),
                    case1not);
            _mm256_store_si256(pvHLStore + i, vHL);
            /* check for saturation */
            {
                vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vH);
                vSaturationCheckMin = _mm256_min_epi8(vSaturationCheckMin, vH);
                vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vHM);
                vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vHS);
                vSaturationCheckMax = _mm256_max_epi8(vSaturationCheckMax, vHL);
            }
#ifdef PARASAIL_TABLE
            arr_store_si256(result->matches_table, vHM, i, segLen, j, s2Len);
            arr_store_si256(result->similar_table, vHS, i, segLen, j, s2Len);
            arr_store_si256(result->length_table, vHL, i, segLen, j, s2Len);
            arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif

            /* Update vE value. */
            vH = _mm256_subs_epi8(vH, vGapO);
            vE = _mm256_subs_epi8(vE, vGapE);
            vE = _mm256_max_epi8(vE, vH);
            _mm256_store_si256(pvEStore + i, vE);
            _mm256_store_si256(pvEM + i, vHM);
            _mm256_store_si256(pvES + i, vHS);
            _mm256_store_si256(pvEL + i, vHL);

            /* Update vF value. */
            vF = _mm256_subs_epi8(vF, vGapE);
            vF = _mm256_max_epi8(vF, vH);
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
            __m256i vHp = _mm256_slli_si256_rpl(pvHLoad[segLen - 1], 1);
            vF = _mm256_slli_si256_rpl(vF, 1);
            vF = _mm256_insert_epi8_rpl(vF, -open, 0);
            vFM = _mm256_slli_si256_rpl(vFM, 1);
            vFS = _mm256_slli_si256_rpl(vFS, 1);
            vFL = _mm256_slli_si256_rpl(vFL, 1);
            for (i=0; i<segLen; ++i) {
                __m256i case1not;
                __m256i case2not;
                __m256i case2;

                /* need to know where match and length come from so
                 * recompute the cases as in the main loop */
                vHp = _mm256_adds_epi8(vHp, _mm256_load_si256(vP + i));
                vE = _mm256_load_si256(pvELoad + i);
                case1not = _mm256_or_si256(
                        _mm256_cmplt_epi8_rpl(vHp,vF),_mm256_cmplt_epi8_rpl(vHp,vE));
                case2not = _mm256_cmplt_epi8_rpl(vF,vE);
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
                vHL = _mm256_blendv_epi8(vHL, _mm256_adds_epi8(vFL,vOne), case2);
                _mm256_store_si256(pvHLStore + i, vHL);
                _mm256_store_si256(pvEL + i, vHL);

                vH = _mm256_load_si256(pvHStore + i);
                vH = _mm256_max_epi8(vH,vF);
                _mm256_store_si256(pvHStore + i, vH);
#ifdef PARASAIL_TABLE
                arr_store_si256(result->matches_table, vHM, i, segLen, j, s2Len);
                arr_store_si256(result->similar_table, vHS, i, segLen, j, s2Len);
                arr_store_si256(result->length_table, vHL, i, segLen, j, s2Len);
                arr_store_si256(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vH = _mm256_subs_epi8(vH, vGapO);
                vF = _mm256_subs_epi8(vF, vGapE);
                if (! _mm256_movemask_epi8(_mm256_cmpgt_epi8(vF, vH))) goto end;
                /*vF = _mm256_max_epi8(vF, vH);*/
                vFM = vHM;
                vFS = vHS;
                vFL = vHL;
                vHp = _mm256_load_si256(pvHLoad + i);
            }
        }
end:
        {
            /* extract vector containing last value from the column */
            __m256i cond_max;
            vH = _mm256_load_si256(pvHStore + offset);
            vHM = _mm256_load_si256(pvHMStore + offset);
            vHS = _mm256_load_si256(pvHSStore + offset);
            vHL = _mm256_load_si256(pvHLStore + offset);
            cond_max = _mm256_cmpgt_epi8(vH, vMaxH);
            vMaxH = _mm256_blendv_epi8(vMaxH, vH, cond_max);
            vMaxHM = _mm256_blendv_epi8(vMaxHM, vHM, cond_max);
            vMaxHS = _mm256_blendv_epi8(vMaxHS, vHS, cond_max);
            vMaxHL = _mm256_blendv_epi8(vMaxHL, vHL, cond_max);
            if (_mm256_movemask_epi8(_mm256_and_si256(vPosMask, cond_max))) {
                end_ref = j;
                end_query = s1Len - 1;
            }
#ifdef PARASAIL_ROWCOL
            for (k=0; k<position; ++k) {
                vH = _mm256_slli_si256_rpl(vH, 1);
                vHM = _mm256_slli_si256_rpl(vHM, 1);
                vHS = _mm256_slli_si256_rpl(vHS, 1);
                vHL = _mm256_slli_si256_rpl(vHL, 1);
            }
            result->score_row[j] = (int8_t) _mm256_extract_epi8_rpl (vH, 31);
            result->matches_row[j] = (int8_t) _mm256_extract_epi8_rpl (vHM, 31);
            result->similar_row[j] = (int8_t) _mm256_extract_epi8_rpl (vHS, 31);
            result->length_row[j] = (int8_t) _mm256_extract_epi8_rpl (vHL, 31);
#endif
        }
    }

    {
        /* extract last value from the column */
        for (k=0; k<position; ++k) {
            vMaxH  = _mm256_slli_si256_rpl (vMaxH, 1);
            vMaxHM = _mm256_slli_si256_rpl (vMaxHM, 1);
            vMaxHS = _mm256_slli_si256_rpl (vMaxHS, 1);
            vMaxHL = _mm256_slli_si256_rpl (vMaxHL, 1);
        }
        score = (int8_t) _mm256_extract_epi8_rpl (vMaxH, 31);
        matches = (int8_t)_mm256_extract_epi8_rpl(vMaxHM, 31);
        similar = (int8_t)_mm256_extract_epi8_rpl(vMaxHS, 31);
        length = (int8_t)_mm256_extract_epi8_rpl(vMaxHL, 31);
    }

    /* max of last column */
    if (INT32_MAX == profile->stop || 0 == profile->stop)
    {
        int8_t score_last;
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
            __m256i vH = _mm256_load_si256(pvHStore + i);
#ifdef PARASAIL_ROWCOL
            __m256i vHM = _mm256_load_si256(pvHMStore + i);
            __m256i vHS = _mm256_load_si256(pvHSStore + i);
            __m256i vHL = _mm256_load_si256(pvHLStore + i);
            arr_store_col(result->score_col, vH, i, segLen);
            arr_store_col(result->matches_col, vHM, i, segLen);
            arr_store_col(result->similar_col, vHS, i, segLen);
            arr_store_col(result->length_col, vHL, i, segLen);
#endif
            vMaxH = _mm256_max_epi8(vH, vMaxH);
        }

        /* max in vec */
        score_last = _mm256_hmax_epi8_rpl(vMaxH);
        if (score_last > score) {
            end_query = s1Len;
            end_ref = s2Len - 1;
            /* Trace the alignment ending position on read. */
            {
                int8_t *t = (int8_t*)pvHStore;
                int8_t *m = (int8_t*)pvHMStore;
                int8_t *s = (int8_t*)pvHSStore;
                int8_t *l = (int8_t*)pvHLStore;
                int32_t column_len = segLen * segWidth;
                for (i = 0; i<column_len; ++i, ++t, ++m, ++s, ++l) {
                    int32_t temp = i / segWidth + i % segWidth * segLen;
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

    if (_mm256_movemask_epi8(_mm256_or_si256(
            _mm256_cmpeq_epi8(vSaturationCheckMin, vNegLimit),
            _mm256_cmpeq_epi8(vSaturationCheckMax, vPosLimit)))) {
        result->saturated = 1;
        score = INT8_MAX;
        matches = 0;
        similar = 0;
        length = 0;
        end_query = 0;
        end_ref = 0;
    }

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
    parasail_result_t *result = parasail_sg_striped_profile_avx2_256_8(profile, s2, s2Len, open, gap);
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
            profile_final = parasail_profile_create_stats_avx_256_8(
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


