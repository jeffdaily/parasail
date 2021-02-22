/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>



#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_neon.h"



#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        simde__m128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[1LL*(i+0)*s2Len + (j-0)] = (int64_t)simde_mm_extract_epi64(vWH, 1);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[1LL*(i+1)*s2Len + (j-1)] = (int64_t)simde_mm_extract_epi64(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_rowcol(
        int *row,
        int *col,
        simde__m128i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (i+0 == s1Len-1 && 0 <= j-0 && j-0 < s2Len) {
        row[j-0] = (int64_t)simde_mm_extract_epi64(vWH, 1);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int64_t)simde_mm_extract_epi64(vWH, 1);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int64_t)simde_mm_extract_epi64(vWH, 0);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int64_t)simde_mm_extract_epi64(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_stats_table_diag_neon_128_64
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_stats_rowcol_diag_neon_128_64
#else
#define FNAME parasail_sw_stats_diag_neon_128_64
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int _s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    /* declare local variables */
    int32_t N = 0;
    int32_t PAD = 0;
    int32_t PAD2 = 0;
    int32_t s1Len = 0;
    int32_t s1Len_PAD = 0;
    int32_t s2Len_PAD = 0;
    int64_t * restrict s1 = NULL;
    int64_t * restrict s2B = NULL;
    int64_t * restrict _H_pr = NULL;
    int64_t * restrict _HM_pr = NULL;
    int64_t * restrict _HS_pr = NULL;
    int64_t * restrict _HL_pr = NULL;
    int64_t * restrict _F_pr = NULL;
    int64_t * restrict _FM_pr = NULL;
    int64_t * restrict _FS_pr = NULL;
    int64_t * restrict _FL_pr = NULL;
    int64_t * restrict s2 = NULL;
    int64_t * restrict H_pr = NULL;
    int64_t * restrict HM_pr = NULL;
    int64_t * restrict HS_pr = NULL;
    int64_t * restrict HL_pr = NULL;
    int64_t * restrict F_pr = NULL;
    int64_t * restrict FM_pr = NULL;
    int64_t * restrict FS_pr = NULL;
    int64_t * restrict FL_pr = NULL;
    parasail_result_t *result = NULL;
    int32_t i = 0;
    int32_t j = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int64_t NEG_LIMIT = 0;
    int64_t POS_LIMIT = 0;
    int64_t score = 0;
    int64_t matches = 0;
    int64_t similar = 0;
    int64_t length = 0;
    simde__m128i vNegLimit;
    simde__m128i vPosLimit;
    simde__m128i vSaturationCheckMin;
    simde__m128i vSaturationCheckMax;
    simde__m128i vNegInf;
    simde__m128i vNegInf0;
    simde__m128i vOpen;
    simde__m128i vGap;
    simde__m128i vZero;
    simde__m128i vOne;
    simde__m128i vN;
    simde__m128i vNegOne;
    simde__m128i vI;
    simde__m128i vJreset;
    simde__m128i vMaxH;
    simde__m128i vMaxM;
    simde__m128i vMaxS;
    simde__m128i vMaxL;
    simde__m128i vEndI;
    simde__m128i vEndJ;
    simde__m128i vILimit;
    simde__m128i vJLimit;

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
    N = 2; /* number of values in vector */
    PAD = N-1;
    PAD2 = PAD*2;
    s1Len = matrix->type == PARASAIL_MATRIX_TYPE_SQUARE ? _s1Len : matrix->length;
    s1Len_PAD = s1Len+PAD;
    s2Len_PAD = s2Len+PAD;
    i = 0;
    j = 0;
    end_query = 0;
    end_ref = 0;
    NEG_LIMIT = (-open < matrix->min ? INT64_MIN + open : INT64_MIN - matrix->min) + 1;
    POS_LIMIT = INT64_MAX - matrix->max - 1;
    score = NEG_LIMIT;
    matches = NEG_LIMIT;
    similar = NEG_LIMIT;
    length = NEG_LIMIT;
    vNegLimit = simde_mm_set1_epi64x(NEG_LIMIT);
    vPosLimit = simde_mm_set1_epi64x(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;
    vNegInf = simde_mm_set1_epi64x(NEG_LIMIT);
    vNegInf0 = simde_mm_srli_si128(vNegInf, 8); /* shift in a 0 */
    vOpen = simde_mm_set1_epi64x(open);
    vGap  = simde_mm_set1_epi64x(gap);
    vZero = simde_mm_set1_epi64x(0);
    vOne = simde_mm_set1_epi64x(1);
    vN = simde_mm_set1_epi64x(N);
    vNegOne = simde_mm_set1_epi64x(-1);
    vI = simde_mm_set_epi64x(0,1);
    vJreset = simde_mm_set_epi64x(0,-1);
    vMaxH = vNegInf;
    vMaxM = vNegInf;
    vMaxS = vNegInf;
    vMaxL = vNegInf;
    vEndI = vNegInf;
    vEndJ = vNegInf;
    vILimit = simde_mm_set1_epi64x(s1Len);
    vJLimit = simde_mm_set1_epi64x(s2Len);

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
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_2;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    /* initialize heap variables */
    s1     = parasail_memalign_int64_t(16, s1Len+PAD);
    s2B    = parasail_memalign_int64_t(16, s2Len+PAD2);
    _H_pr  = parasail_memalign_int64_t(16, s2Len+PAD2);
    _HM_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    _HS_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    _HL_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    _F_pr  = parasail_memalign_int64_t(16, s2Len+PAD2);
    _FM_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    _FS_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    _FL_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
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
        simde__m128i case1 = vZero;
        simde__m128i case2 = vZero;
        simde__m128i case0 = vZero;
        simde__m128i vNH = vNegInf0;
        simde__m128i vNM = vZero;
        simde__m128i vNS = vZero;
        simde__m128i vNL = vZero;
        simde__m128i vWH = vNegInf0;
        simde__m128i vWM = vZero;
        simde__m128i vWS = vZero;
        simde__m128i vWL = vZero;
        simde__m128i vE = vNegInf;
        simde__m128i vE_opn = vNegInf;
        simde__m128i vE_ext = vNegInf;
        simde__m128i vEM = vZero;
        simde__m128i vES = vZero;
        simde__m128i vEL = vZero;
        simde__m128i vF = vNegInf;
        simde__m128i vF_opn = vNegInf;
        simde__m128i vF_ext = vNegInf;
        simde__m128i vFM = vZero;
        simde__m128i vFS = vZero;
        simde__m128i vFL = vZero;
        simde__m128i vJ = vJreset;
        simde__m128i vs1 = simde_mm_set_epi64x(
                s1[i+0],
                s1[i+1]);
        simde__m128i vs2 = vNegInf;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+0] : ((i+0 >= s1Len) ? s1Len-1 : i+0))];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+1] : ((i+1 >= s1Len) ? s1Len-1 : i+1))];
        simde__m128i vIltLimit = simde_mm_cmplt_epi64(vI, vILimit);
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            simde__m128i vMat;
            simde__m128i vNWH = vNH;
            simde__m128i vNWM = vNM;
            simde__m128i vNWS = vNS;
            simde__m128i vNWL = vNL;
            vNH = simde_mm_srli_si128(vWH, 8);
            vNH = simde_mm_insert_epi64(vNH, H_pr[j], 1);
            vNM = simde_mm_srli_si128(vWM, 8);
            vNM = simde_mm_insert_epi64(vNM, HM_pr[j], 1);
            vNS = simde_mm_srli_si128(vWS, 8);
            vNS = simde_mm_insert_epi64(vNS, HS_pr[j], 1);
            vNL = simde_mm_srli_si128(vWL, 8);
            vNL = simde_mm_insert_epi64(vNL, HL_pr[j], 1);
            vF = simde_mm_srli_si128(vF, 8);
            vF = simde_mm_insert_epi64(vF, F_pr[j], 1);
            vFM = simde_mm_srli_si128(vFM, 8);
            vFM = simde_mm_insert_epi64(vFM, FM_pr[j], 1);
            vFS = simde_mm_srli_si128(vFS, 8);
            vFS = simde_mm_insert_epi64(vFS, FS_pr[j], 1);
            vFL = simde_mm_srli_si128(vFL, 8);
            vFL = simde_mm_insert_epi64(vFL, FL_pr[j], 1);
            vF_opn = simde_mm_sub_epi64(vNH, vOpen);
            vF_ext = simde_mm_sub_epi64(vF, vGap);
            vF = simde_mm_max_epi64(vF_opn, vF_ext);
            case1 = simde_mm_cmpgt_epi64(vF_opn, vF_ext);
            vFM = simde_mm_blendv_epi8(vFM, vNM, case1);
            vFS = simde_mm_blendv_epi8(vFS, vNS, case1);
            vFL = simde_mm_blendv_epi8(vFL, vNL, case1);
            vFL = simde_mm_add_epi64(vFL, vOne);
            vE_opn = simde_mm_sub_epi64(vWH, vOpen);
            vE_ext = simde_mm_sub_epi64(vE, vGap);
            vE = simde_mm_max_epi64(vE_opn, vE_ext);
            case1 = simde_mm_cmpgt_epi64(vE_opn, vE_ext);
            vEM = simde_mm_blendv_epi8(vEM, vWM, case1);
            vES = simde_mm_blendv_epi8(vES, vWS, case1);
            vEL = simde_mm_blendv_epi8(vEL, vWL, case1);
            vEL = simde_mm_add_epi64(vEL, vOne);
            vs2 = simde_mm_srli_si128(vs2, 8);
            vs2 = simde_mm_insert_epi64(vs2, s2[j], 1);
            vMat = simde_mm_set_epi64x(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]]
                    );
            vNWH = simde_mm_add_epi64(vNWH, vMat);
            vWH = simde_mm_max_epi64(vNWH, vE);
            vWH = simde_mm_max_epi64(vWH, vF);
            vWH = simde_mm_max_epi64(vWH, vZero);
            case1 = simde_mm_cmpeq_epi64(vWH, vNWH);
            case2 = simde_mm_cmpeq_epi64(vWH, vF);
            case0 = simde_mm_cmpeq_epi64(vWH, vZero);
            vWM = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vEM, vFM, case2),
                    simde_mm_add_epi64(vNWM,
                        simde_mm_and_si128(
                            simde_mm_cmpeq_epi64(vs1,vs2),
                            vOne)),
                    case1);
            vWM = simde_mm_blendv_epi8(vWM, vZero, case0);
            vWS = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vES, vFS, case2),
                    simde_mm_add_epi64(vNWS,
                        simde_mm_and_si128(
                            simde_mm_cmpgt_epi64(vMat,vZero),
                            vOne)),
                    case1);
            vWS = simde_mm_blendv_epi8(vWS, vZero, case0);
            vWL = simde_mm_blendv_epi8(
                    simde_mm_blendv_epi8(vEL, vFL, case2),
                    simde_mm_add_epi64(vNWL, vOne), case1);
            vWL = simde_mm_blendv_epi8(vWL, vZero, case0);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                simde__m128i cond = simde_mm_cmpeq_epi64(vJ,vNegOne);
                vWH = simde_mm_andnot_si128(cond, vWH);
                vWM = simde_mm_andnot_si128(cond, vWM);
                vWS = simde_mm_andnot_si128(cond, vWS);
                vWL = simde_mm_andnot_si128(cond, vWL);
                vE = simde_mm_blendv_epi8(vE, vNegInf, cond);
                vEM = simde_mm_andnot_si128(cond, vEM);
                vES = simde_mm_andnot_si128(cond, vES);
                vEL = simde_mm_andnot_si128(cond, vEL);
            }
            /* cannot start checking sat until after J clears boundary */
            if (j > PAD) {
                vSaturationCheckMin = simde_mm_min_epi64(vSaturationCheckMin, vWH);
                vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vWH);
                vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vWM);
                vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vWS);
                vSaturationCheckMax = simde_mm_max_epi64(vSaturationCheckMax, vWL);
            }
#ifdef PARASAIL_TABLE
            arr_store_si128(result->stats->tables->score_table, vWH, i, s1Len, j, s2Len);
            arr_store_si128(result->stats->tables->matches_table, vWM, i, s1Len, j, s2Len);
            arr_store_si128(result->stats->tables->similar_table, vWS, i, s1Len, j, s2Len);
            arr_store_si128(result->stats->tables->length_table, vWL, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->stats->rowcols->score_row,   result->stats->rowcols->score_col, vWH, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->matches_row, result->stats->rowcols->matches_col, vWM, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->similar_row, result->stats->rowcols->similar_col, vWS, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->length_row,  result->stats->rowcols->length_col, vWL, i, s1Len, j, s2Len);
#endif
            H_pr[j-1] = (int64_t)simde_mm_extract_epi64(vWH,0);
            HM_pr[j-1] = (int64_t)simde_mm_extract_epi64(vWM,0);
            HS_pr[j-1] = (int64_t)simde_mm_extract_epi64(vWS,0);
            HL_pr[j-1] = (int64_t)simde_mm_extract_epi64(vWL,0);
            F_pr[j-1] = (int64_t)simde_mm_extract_epi64(vF,0);
            FM_pr[j-1] = (int64_t)simde_mm_extract_epi64(vFM,0);
            FS_pr[j-1] = (int64_t)simde_mm_extract_epi64(vFS,0);
            FL_pr[j-1] = (int64_t)simde_mm_extract_epi64(vFL,0);
            /* as minor diagonal vector passes across table, extract
             * max values within the i,j bounds */
            {
                simde__m128i cond_valid_J = simde_mm_and_si128(
                        simde_mm_cmpgt_epi64(vJ, vNegOne),
                        simde_mm_cmplt_epi64(vJ, vJLimit));
                simde__m128i cond_valid_IJ = simde_mm_and_si128(cond_valid_J, vIltLimit);
                simde__m128i cond_eq = simde_mm_cmpeq_epi64(vWH, vMaxH);
                simde__m128i cond_max = simde_mm_cmpgt_epi64(vWH, vMaxH);
                simde__m128i cond_all = simde_mm_and_si128(cond_max, cond_valid_IJ);
                simde__m128i cond_Jlt = simde_mm_cmplt_epi64(vJ, vEndJ);
                vMaxH = simde_mm_blendv_epi8(vMaxH, vWH, cond_all);
                vMaxM = simde_mm_blendv_epi8(vMaxM, vWM, cond_all);
                vMaxS = simde_mm_blendv_epi8(vMaxS, vWS, cond_all);
                vMaxL = simde_mm_blendv_epi8(vMaxL, vWL, cond_all);
                vEndI = simde_mm_blendv_epi8(vEndI, vI, cond_all);
                vEndJ = simde_mm_blendv_epi8(vEndJ, vJ, cond_all);
                cond_all = simde_mm_and_si128(cond_Jlt, cond_eq);
                cond_all = simde_mm_and_si128(cond_all, cond_valid_IJ);
                vMaxM = simde_mm_blendv_epi8(vMaxM, vWM, cond_all);
                vMaxS = simde_mm_blendv_epi8(vMaxS, vWS, cond_all);
                vMaxL = simde_mm_blendv_epi8(vMaxL, vWL, cond_all);
                vEndI = simde_mm_blendv_epi8(vEndI, vI, cond_all);
                vEndJ = simde_mm_blendv_epi8(vEndJ, vJ, cond_all);
            }
            vJ = simde_mm_add_epi64(vJ, vOne);
        }
        vI = simde_mm_add_epi64(vI, vN);
    }

    /* alignment ending position */
    {
        int64_t *t = (int64_t*)&vMaxH;
        int64_t *m = (int64_t*)&vMaxM;
        int64_t *s = (int64_t*)&vMaxS;
        int64_t *l = (int64_t*)&vMaxL;
        int64_t *i = (int64_t*)&vEndI;
        int64_t *j = (int64_t*)&vEndJ;
        int32_t k;
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

    if (simde_mm_movemask_epi8(simde_mm_or_si128(
            simde_mm_cmplt_epi64(vSaturationCheckMin, vNegLimit),
            simde_mm_cmpgt_epi64(vSaturationCheckMax, vPosLimit)))) {
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


