/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>

#include <emmintrin.h>
#include <smmintrin.h>

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_sse.h"

#define NEG_INF (INT16_MIN/(int16_t)(2))


#ifdef PARASAIL_TABLE
static inline void arr_store_si128(
        int *array,
        __m128i vWscore,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[(i+0)*s2Len + (j-0)] = (int16_t)_mm_extract_epi16(vWscore, 7);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[(i+1)*s2Len + (j-1)] = (int16_t)_mm_extract_epi16(vWscore, 6);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[(i+2)*s2Len + (j-2)] = (int16_t)_mm_extract_epi16(vWscore, 5);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[(i+3)*s2Len + (j-3)] = (int16_t)_mm_extract_epi16(vWscore, 4);
    }
    if (0 <= i+4 && i+4 < s1Len && 0 <= j-4 && j-4 < s2Len) {
        array[(i+4)*s2Len + (j-4)] = (int16_t)_mm_extract_epi16(vWscore, 3);
    }
    if (0 <= i+5 && i+5 < s1Len && 0 <= j-5 && j-5 < s2Len) {
        array[(i+5)*s2Len + (j-5)] = (int16_t)_mm_extract_epi16(vWscore, 2);
    }
    if (0 <= i+6 && i+6 < s1Len && 0 <= j-6 && j-6 < s2Len) {
        array[(i+6)*s2Len + (j-6)] = (int16_t)_mm_extract_epi16(vWscore, 1);
    }
    if (0 <= i+7 && i+7 < s1Len && 0 <= j-7 && j-7 < s2Len) {
        array[(i+7)*s2Len + (j-7)] = (int16_t)_mm_extract_epi16(vWscore, 0);
    }
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_rowcol(
        int *row,
        int *col,
        __m128i vWscore,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (i+0 == s1Len-1 && 0 <= j-0 && j-0 < s2Len) {
        row[j-0] = (int16_t)_mm_extract_epi16(vWscore, 7);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int16_t)_mm_extract_epi16(vWscore, 7);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int16_t)_mm_extract_epi16(vWscore, 6);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int16_t)_mm_extract_epi16(vWscore, 6);
    }
    if (i+2 == s1Len-1 && 0 <= j-2 && j-2 < s2Len) {
        row[j-2] = (int16_t)_mm_extract_epi16(vWscore, 5);
    }
    if (j-2 == s2Len-1 && 0 <= i+2 && i+2 < s1Len) {
        col[(i+2)] = (int16_t)_mm_extract_epi16(vWscore, 5);
    }
    if (i+3 == s1Len-1 && 0 <= j-3 && j-3 < s2Len) {
        row[j-3] = (int16_t)_mm_extract_epi16(vWscore, 4);
    }
    if (j-3 == s2Len-1 && 0 <= i+3 && i+3 < s1Len) {
        col[(i+3)] = (int16_t)_mm_extract_epi16(vWscore, 4);
    }
    if (i+4 == s1Len-1 && 0 <= j-4 && j-4 < s2Len) {
        row[j-4] = (int16_t)_mm_extract_epi16(vWscore, 3);
    }
    if (j-4 == s2Len-1 && 0 <= i+4 && i+4 < s1Len) {
        col[(i+4)] = (int16_t)_mm_extract_epi16(vWscore, 3);
    }
    if (i+5 == s1Len-1 && 0 <= j-5 && j-5 < s2Len) {
        row[j-5] = (int16_t)_mm_extract_epi16(vWscore, 2);
    }
    if (j-5 == s2Len-1 && 0 <= i+5 && i+5 < s1Len) {
        col[(i+5)] = (int16_t)_mm_extract_epi16(vWscore, 2);
    }
    if (i+6 == s1Len-1 && 0 <= j-6 && j-6 < s2Len) {
        row[j-6] = (int16_t)_mm_extract_epi16(vWscore, 1);
    }
    if (j-6 == s2Len-1 && 0 <= i+6 && i+6 < s1Len) {
        col[(i+6)] = (int16_t)_mm_extract_epi16(vWscore, 1);
    }
    if (i+7 == s1Len-1 && 0 <= j-7 && j-7 < s2Len) {
        row[j-7] = (int16_t)_mm_extract_epi16(vWscore, 0);
    }
    if (j-7 == s2Len-1 && 0 <= i+7 && i+7 < s1Len) {
        col[(i+7)] = (int16_t)_mm_extract_epi16(vWscore, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sg_stats_table_diag_sse41_128_16
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_stats_rowcol_diag_sse41_128_16
#else
#define FNAME parasail_sg_stats_diag_sse41_128_16
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    const int32_t N = 8; /* number of values in vector */
    const int32_t PAD = N-1;
    const int32_t PAD2 = PAD*2;
    const int32_t s1Len_PAD = s1Len+PAD;
    const int32_t s2Len_PAD = s2Len+PAD;
    int16_t * const restrict s1      = parasail_memalign_int16_t(16, s1Len+PAD);
    int16_t * const restrict s2B     = parasail_memalign_int16_t(16, s2Len+PAD2);
    int16_t * const restrict _tbl_pr = parasail_memalign_int16_t(16, s2Len+PAD2);
    int16_t * const restrict _del_pr = parasail_memalign_int16_t(16, s2Len+PAD2);
    int16_t * const restrict _mch_pr = parasail_memalign_int16_t(16, s2Len+PAD2);
    int16_t * const restrict _sim_pr = parasail_memalign_int16_t(16, s2Len+PAD2);
    int16_t * const restrict _len_pr = parasail_memalign_int16_t(16, s2Len+PAD2);
    int16_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int16_t * const restrict tbl_pr = _tbl_pr+PAD;
    int16_t * const restrict del_pr = _del_pr+PAD;
    int16_t * const restrict mch_pr = _mch_pr+PAD;
    int16_t * const restrict sim_pr = _sim_pr+PAD;
    int16_t * const restrict len_pr = _len_pr+PAD;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif
    int32_t i = 0;
    int32_t j = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int16_t score = NEG_INF;
    int16_t matches = NEG_INF;
    int16_t similar = NEG_INF;
    int16_t length = NEG_INF;
    
    __m128i vNegInf = _mm_set1_epi16(NEG_INF);
    __m128i vNegInf0 = _mm_srli_si128(vNegInf, 2); /* shift in a 0 */
    __m128i vOpen = _mm_set1_epi16(open);
    __m128i vGap  = _mm_set1_epi16(gap);
    __m128i vZero = _mm_set1_epi16(0);
    __m128i vOne = _mm_set1_epi16(1);
    __m128i vN = _mm_set1_epi16(N);
    __m128i vNegOne = _mm_set1_epi16(-1);
    __m128i vI = _mm_set_epi16(0,1,2,3,4,5,6,7);
    __m128i vJreset = _mm_set_epi16(0,-1,-2,-3,-4,-5,-6,-7);
    __m128i vMaxScore = vNegInf;
    __m128i vMaxMatch = vNegInf;
    __m128i vMaxSimilar = vNegInf;
    __m128i vMaxLength = vNegInf;
    __m128i vEndI = vNegInf;
    __m128i vEndJ = vNegInf;
    __m128i vILimit = _mm_set1_epi16(s1Len);
    __m128i vILimit1 = _mm_sub_epi16(vILimit, vOne);
    __m128i vJLimit = _mm_set1_epi16(s2Len);
    __m128i vJLimit1 = _mm_sub_epi16(vJLimit, vOne);

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
        mch_pr[j] = 0;
        sim_pr[j] = 0;
        len_pr[j] = 0;
    }
    /* pad front of stored row values */
    for (j=-PAD; j<0; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        sim_pr[j] = 0;
        len_pr[j] = 0;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        sim_pr[j] = 0;
        len_pr[j] = 0;
    }
    tbl_pr[-1] = 0; /* upper left corner */

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        __m128i vNscore = vNegInf0;
        __m128i vNmatch = vZero;
        __m128i vNsimilar = vZero;
        __m128i vNlength = vZero;
        __m128i vWscore = vNegInf0;
        __m128i vWmatch = vZero;
        __m128i vWsimilar = vZero;
        __m128i vWlength = vZero;
        __m128i vIns = vNegInf;
        __m128i vDel = vNegInf;
        __m128i vJ = vJreset;
        __m128i vs1 = _mm_set_epi16(
                s1[i+0],
                s1[i+1],
                s1[i+2],
                s1[i+3],
                s1[i+4],
                s1[i+5],
                s1[i+6],
                s1[i+7]);
        __m128i vs2 = vNegInf;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size*s1[i+2]];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size*s1[i+3]];
        const int * const restrict matrow4 = &matrix->matrix[matrix->size*s1[i+4]];
        const int * const restrict matrow5 = &matrix->matrix[matrix->size*s1[i+5]];
        const int * const restrict matrow6 = &matrix->matrix[matrix->size*s1[i+6]];
        const int * const restrict matrow7 = &matrix->matrix[matrix->size*s1[i+7]];
        __m128i vIltLimit = _mm_cmplt_epi16(vI, vILimit);
        __m128i vIeqLimit1 = _mm_cmpeq_epi16(vI, vILimit1);
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            __m128i vNWmatch = vNmatch;
            __m128i vNWsimilar = vNsimilar;
            __m128i vNWlength = vNlength;
            vNscore = _mm_srli_si128(vWscore, 2);
            vNscore = _mm_insert_epi16(vNscore, tbl_pr[j], 7);
            vNmatch = _mm_srli_si128(vWmatch, 2);
            vNmatch = _mm_insert_epi16(vNmatch, mch_pr[j], 7);
            vNsimilar = _mm_srli_si128(vWsimilar, 2);
            vNsimilar = _mm_insert_epi16(vNsimilar, sim_pr[j], 7);
            vNlength = _mm_srli_si128(vWlength, 2);
            vNlength = _mm_insert_epi16(vNlength, len_pr[j], 7);
            vDel = _mm_srli_si128(vDel, 2);
            vDel = _mm_insert_epi16(vDel, del_pr[j], 7);
            vDel = _mm_max_epi16(
                    _mm_sub_epi16(vNscore, vOpen),
                    _mm_sub_epi16(vDel, vGap));
            vIns = _mm_max_epi16(
                    _mm_sub_epi16(vWscore, vOpen),
                    _mm_sub_epi16(vIns, vGap));
            vs2 = _mm_srli_si128(vs2, 2);
            vs2 = _mm_insert_epi16(vs2, s2[j], 7);
            vMat = _mm_set_epi16(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]],
                    matrow4[s2[j-4]],
                    matrow5[s2[j-5]],
                    matrow6[s2[j-6]],
                    matrow7[s2[j-7]]
                    );
            vNWscore = _mm_add_epi16(vNWscore, vMat);
            vWscore = _mm_max_epi16(vNWscore, vIns);
            vWscore = _mm_max_epi16(vWscore, vDel);
            /* conditional block */
            {
                __m128i case1not;
                __m128i case2not;
                __m128i case2;
                __m128i case3;
                __m128i vCmatch;
                __m128i vCsimilar;
                __m128i vClength;
                case1not = _mm_or_si128(
                        _mm_cmplt_epi16(vNWscore,vDel),
                        _mm_cmplt_epi16(vNWscore,vIns));
                case2not = _mm_cmplt_epi16(vDel,vIns);
                case2 = _mm_andnot_si128(case2not,case1not);
                case3 = _mm_and_si128(case1not,case2not);
                vCmatch = _mm_andnot_si128(case1not,
                        _mm_add_epi16(vNWmatch, _mm_and_si128(
                                _mm_cmpeq_epi16(vs1,vs2),vOne)));
                vCmatch = _mm_or_si128(vCmatch, _mm_and_si128(case2, vNmatch));
                vCmatch = _mm_or_si128(vCmatch, _mm_and_si128(case3, vWmatch));
                vCsimilar = _mm_andnot_si128(case1not,
                        _mm_add_epi16(vNWsimilar, _mm_and_si128(
                                _mm_cmpgt_epi16(vMat,vZero),vOne)));
                vCsimilar = _mm_or_si128(vCsimilar, _mm_and_si128(case2, vNsimilar));
                vCsimilar = _mm_or_si128(vCsimilar, _mm_and_si128(case3, vWsimilar));
                vClength= _mm_andnot_si128(case1not,
                        _mm_add_epi16(vNWlength, vOne));
                vClength= _mm_or_si128(vClength,_mm_and_si128(case2,
                            _mm_add_epi16(vNlength, vOne)));
                vClength= _mm_or_si128(vClength,_mm_and_si128(case3,
                            _mm_add_epi16(vWlength, vOne)));
                vWmatch = vCmatch;
                vWsimilar = vCsimilar;
                vWlength = vClength;
            }

            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m128i cond = _mm_cmpeq_epi16(vJ,vNegOne);
                vWscore = _mm_andnot_si128(cond, vWscore);
                vWmatch = _mm_andnot_si128(cond, vWmatch);
                vWsimilar = _mm_andnot_si128(cond, vWsimilar);
                vWlength = _mm_andnot_si128(cond, vWlength);
                vDel = _mm_blendv_epi8(vDel, vNegInf, cond);
                vIns = _mm_blendv_epi8(vIns, vNegInf, cond);
            }
            
#ifdef PARASAIL_TABLE
            arr_store_si128(result->score_table, vWscore, i, s1Len, j, s2Len);
            arr_store_si128(result->matches_table, vWmatch, i, s1Len, j, s2Len);
            arr_store_si128(result->similar_table, vWsimilar, i, s1Len, j, s2Len);
            arr_store_si128(result->length_table, vWlength, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->score_row, result->score_col, vWscore, i, s1Len, j, s2Len);
            arr_store_rowcol(result->matches_row, result->matches_col, vWmatch, i, s1Len, j, s2Len);
            arr_store_rowcol(result->similar_row, result->similar_col, vWsimilar, i, s1Len, j, s2Len);
            arr_store_rowcol(result->length_row, result->length_col, vWlength, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-7] = (int16_t)_mm_extract_epi16(vWscore,0);
            mch_pr[j-7] = (int16_t)_mm_extract_epi16(vWmatch,0);
            sim_pr[j-7] = (int16_t)_mm_extract_epi16(vWsimilar,0);
            len_pr[j-7] = (int16_t)_mm_extract_epi16(vWlength,0);
            del_pr[j-7] = (int16_t)_mm_extract_epi16(vDel,0);
            /* as minor diagonal vector passes across the i or j limit
             * boundary, extract the last value of the column or row */
            {
                __m128i vJeqLimit1 = _mm_cmpeq_epi16(vJ, vJLimit1);
                __m128i vJgtNegOne = _mm_cmpgt_epi16(vJ, vNegOne);
                __m128i vJltLimit = _mm_cmplt_epi16(vJ, vJLimit);
                __m128i cond_j = _mm_and_si128(vIltLimit, vJeqLimit1);
                __m128i cond_i = _mm_and_si128(vIeqLimit1,
                        _mm_and_si128(vJgtNegOne, vJltLimit));
                __m128i cond_ij = _mm_or_si128(cond_i, cond_j);
                __m128i cond_max = _mm_cmpgt_epi16(vWscore, vMaxScore);
                __m128i cond_eq = _mm_cmpeq_epi16(vWscore, vMaxScore);
                __m128i cond_all = _mm_and_si128(cond_max, cond_ij);
                __m128i cond_Jlt = _mm_cmplt_epi16(vJ, vEndJ);
                vMaxScore = _mm_blendv_epi8(vMaxScore, vWscore, cond_all);
                vMaxMatch = _mm_blendv_epi8(vMaxMatch, vWmatch, cond_all);
                vMaxSimilar = _mm_blendv_epi8(vMaxSimilar, vWsimilar, cond_all);
                vMaxLength = _mm_blendv_epi8(vMaxLength, vWlength, cond_all);
                vEndI = _mm_blendv_epi8(vEndI, vI, cond_all);
                vEndJ = _mm_blendv_epi8(vEndJ, vJ, cond_all);
                cond_all = _mm_and_si128(cond_Jlt, cond_eq);
                cond_all = _mm_and_si128(cond_all, cond_ij);
                vMaxMatch = _mm_blendv_epi8(vMaxMatch, vWmatch, cond_all);
                vMaxSimilar = _mm_blendv_epi8(vMaxSimilar, vWsimilar, cond_all);
                vMaxLength = _mm_blendv_epi8(vMaxLength, vWlength, cond_all);
                vEndI = _mm_blendv_epi8(vEndI, vI, cond_all);
                vEndJ = _mm_blendv_epi8(vEndJ, vJ, cond_all);
            }
            vJ = _mm_add_epi16(vJ, vOne);
        }
        vI = _mm_add_epi16(vI, vN);
    }

    /* alignment ending position */
    {
        int16_t *t = (int16_t*)&vMaxScore;
        int16_t *m = (int16_t*)&vMaxMatch;
        int16_t *s = (int16_t*)&vMaxSimilar;
        int16_t *l = (int16_t*)&vMaxLength;
        int16_t *i = (int16_t*)&vEndI;
        int16_t *j = (int16_t*)&vEndJ;
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

    

    result->score = score;
    result->matches = matches;
    result->similar = similar;
    result->length = length;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(_len_pr);
    parasail_free(_sim_pr);
    parasail_free(_mch_pr);
    parasail_free(_del_pr);
    parasail_free(_tbl_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}


