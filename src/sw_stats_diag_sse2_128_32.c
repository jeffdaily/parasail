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

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_sse.h"

#define NEG_INF (INT32_MIN/(int32_t)(2))

static inline __m128i _mm_blendv_epi8_rpl(__m128i a, __m128i b, __m128i mask) {
    a = _mm_andnot_si128(mask, a);
    a = _mm_or_si128(a, _mm_and_si128(mask, b));
    return a;
}

static inline __m128i _mm_insert_epi32_rpl(__m128i a, int32_t i, const int imm) {
    __m128i_32_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}

static inline __m128i _mm_max_epi32_rpl(__m128i a, __m128i b) {
    __m128i mask = _mm_cmpgt_epi32(a, b);
    a = _mm_and_si128(a, mask);
    b = _mm_andnot_si128(mask, b);
    return _mm_or_si128(a, b);
}

static inline int32_t _mm_extract_epi32_rpl(__m128i a, const int imm) {
    __m128i_32_t A;
    A.m = a;
    return A.v[imm];
}


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
        array[(i+0)*s2Len + (j-0)] = (int32_t)_mm_extract_epi32_rpl(vWscore, 3);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[(i+1)*s2Len + (j-1)] = (int32_t)_mm_extract_epi32_rpl(vWscore, 2);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[(i+2)*s2Len + (j-2)] = (int32_t)_mm_extract_epi32_rpl(vWscore, 1);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[(i+3)*s2Len + (j-3)] = (int32_t)_mm_extract_epi32_rpl(vWscore, 0);
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
        row[j-0] = (int32_t)_mm_extract_epi32_rpl(vWscore, 3);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int32_t)_mm_extract_epi32_rpl(vWscore, 3);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int32_t)_mm_extract_epi32_rpl(vWscore, 2);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int32_t)_mm_extract_epi32_rpl(vWscore, 2);
    }
    if (i+2 == s1Len-1 && 0 <= j-2 && j-2 < s2Len) {
        row[j-2] = (int32_t)_mm_extract_epi32_rpl(vWscore, 1);
    }
    if (j-2 == s2Len-1 && 0 <= i+2 && i+2 < s1Len) {
        col[(i+2)] = (int32_t)_mm_extract_epi32_rpl(vWscore, 1);
    }
    if (i+3 == s1Len-1 && 0 <= j-3 && j-3 < s2Len) {
        row[j-3] = (int32_t)_mm_extract_epi32_rpl(vWscore, 0);
    }
    if (j-3 == s2Len-1 && 0 <= i+3 && i+3 < s1Len) {
        col[(i+3)] = (int32_t)_mm_extract_epi32_rpl(vWscore, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_stats_table_diag_sse2_128_32
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_stats_rowcol_diag_sse2_128_32
#else
#define FNAME parasail_sw_stats_diag_sse2_128_32
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    const int32_t N = 4; /* number of values in vector */
    const int32_t PAD = N-1;
    const int32_t PAD2 = PAD*2;
    const int32_t s1Len_PAD = s1Len+PAD;
    const int32_t s2Len_PAD = s2Len+PAD;
    int32_t * const restrict s1      = parasail_memalign_int32_t(16, s1Len+PAD);
    int32_t * const restrict s2B     = parasail_memalign_int32_t(16, s2Len+PAD2);
    int32_t * const restrict _tbl_pr = parasail_memalign_int32_t(16, s2Len+PAD2);
    int32_t * const restrict _del_pr = parasail_memalign_int32_t(16, s2Len+PAD2);
    int32_t * const restrict _mch_pr = parasail_memalign_int32_t(16, s2Len+PAD2);
    int32_t * const restrict _sim_pr = parasail_memalign_int32_t(16, s2Len+PAD2);
    int32_t * const restrict _len_pr = parasail_memalign_int32_t(16, s2Len+PAD2);
    int32_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int32_t * const restrict tbl_pr = _tbl_pr+PAD;
    int32_t * const restrict del_pr = _del_pr+PAD;
    int32_t * const restrict mch_pr = _mch_pr+PAD;
    int32_t * const restrict sim_pr = _sim_pr+PAD;
    int32_t * const restrict len_pr = _len_pr+PAD;
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
    int32_t score = NEG_INF;
    int32_t matches = NEG_INF;
    int32_t similar = NEG_INF;
    int32_t length = NEG_INF;
    
    __m128i vNegInf = _mm_set1_epi32(NEG_INF);
    __m128i vNegInf0 = _mm_srli_si128(vNegInf, 4); /* shift in a 0 */
    __m128i vOpen = _mm_set1_epi32(open);
    __m128i vGap  = _mm_set1_epi32(gap);
    __m128i vZero = _mm_set1_epi32(0);
    __m128i vOne = _mm_set1_epi32(1);
    __m128i vN = _mm_set1_epi32(N);
    __m128i vNegOne = _mm_set1_epi32(-1);
    __m128i vI = _mm_set_epi32(0,1,2,3);
    __m128i vJreset = _mm_set_epi32(0,-1,-2,-3);
    __m128i vMaxScore = vNegInf;
    __m128i vMaxMatch = vNegInf;
    __m128i vMaxSimilar = vNegInf;
    __m128i vMaxLength = vNegInf;
    __m128i vEndI = vNegInf;
    __m128i vEndJ = vNegInf;
    __m128i vILimit = _mm_set1_epi32(s1Len);
    __m128i vJLimit = _mm_set1_epi32(s2Len);

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
        __m128i vs1 = _mm_set_epi32(
                s1[i+0],
                s1[i+1],
                s1[i+2],
                s1[i+3]);
        __m128i vs2 = vNegInf;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size*s1[i+2]];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size*s1[i+3]];
        __m128i vIltLimit = _mm_cmplt_epi32(vI, vILimit);
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            __m128i vNWmatch = vNmatch;
            __m128i vNWsimilar = vNsimilar;
            __m128i vNWlength = vNlength;
            vNscore = _mm_srli_si128(vWscore, 4);
            vNscore = _mm_insert_epi32_rpl(vNscore, tbl_pr[j], 3);
            vNmatch = _mm_srli_si128(vWmatch, 4);
            vNmatch = _mm_insert_epi32_rpl(vNmatch, mch_pr[j], 3);
            vNsimilar = _mm_srli_si128(vWsimilar, 4);
            vNsimilar = _mm_insert_epi32_rpl(vNsimilar, sim_pr[j], 3);
            vNlength = _mm_srli_si128(vWlength, 4);
            vNlength = _mm_insert_epi32_rpl(vNlength, len_pr[j], 3);
            vDel = _mm_srli_si128(vDel, 4);
            vDel = _mm_insert_epi32_rpl(vDel, del_pr[j], 3);
            vDel = _mm_max_epi32_rpl(
                    _mm_sub_epi32(vNscore, vOpen),
                    _mm_sub_epi32(vDel, vGap));
            vIns = _mm_max_epi32_rpl(
                    _mm_sub_epi32(vWscore, vOpen),
                    _mm_sub_epi32(vIns, vGap));
            vs2 = _mm_srli_si128(vs2, 4);
            vs2 = _mm_insert_epi32_rpl(vs2, s2[j], 3);
            vMat = _mm_set_epi32(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]]
                    );
            vNWscore = _mm_add_epi32(vNWscore, vMat);
            vWscore = _mm_max_epi32_rpl(vNWscore, vIns);
            vWscore = _mm_max_epi32_rpl(vWscore, vDel);
            vWscore = _mm_max_epi32_rpl(vWscore, vZero);
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
                        _mm_cmplt_epi32(vNWscore,vDel),
                        _mm_cmplt_epi32(vNWscore,vIns));
                case2not = _mm_cmplt_epi32(vDel,vIns);
                case2 = _mm_andnot_si128(case2not,case1not);
                case3 = _mm_and_si128(case1not,case2not);
                vCmatch = _mm_andnot_si128(case1not,
                        _mm_add_epi32(vNWmatch, _mm_and_si128(
                                _mm_cmpeq_epi32(vs1,vs2),vOne)));
                vCmatch = _mm_or_si128(vCmatch, _mm_and_si128(case2, vNmatch));
                vCmatch = _mm_or_si128(vCmatch, _mm_and_si128(case3, vWmatch));
                vCsimilar = _mm_andnot_si128(case1not,
                        _mm_add_epi32(vNWsimilar, _mm_and_si128(
                                _mm_cmpgt_epi32(vMat,vZero),vOne)));
                vCsimilar = _mm_or_si128(vCsimilar, _mm_and_si128(case2, vNsimilar));
                vCsimilar = _mm_or_si128(vCsimilar, _mm_and_si128(case3, vWsimilar));
                vClength= _mm_andnot_si128(case1not,
                        _mm_add_epi32(vNWlength, vOne));
                vClength= _mm_or_si128(vClength,_mm_and_si128(case2,
                            _mm_add_epi32(vNlength, vOne)));
                vClength= _mm_or_si128(vClength,_mm_and_si128(case3,
                            _mm_add_epi32(vWlength, vOne)));
                vWmatch = vCmatch;
                vWsimilar = vCsimilar;
                vWlength = vClength;
            }
            {
                __m128i cond_zero = _mm_cmpeq_epi32(vWscore, vZero);
                vWmatch = _mm_andnot_si128(cond_zero, vWmatch);
                vWsimilar = _mm_andnot_si128(cond_zero, vWsimilar);
                vWlength = _mm_andnot_si128(cond_zero, vWlength);
            }

            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m128i cond = _mm_cmpeq_epi32(vJ,vNegOne);
                vWscore = _mm_andnot_si128(cond, vWscore);
                vWmatch = _mm_andnot_si128(cond, vWmatch);
                vWsimilar = _mm_andnot_si128(cond, vWsimilar);
                vWlength = _mm_andnot_si128(cond, vWlength);
                vDel = _mm_blendv_epi8_rpl(vDel, vNegInf, cond);
                vIns = _mm_blendv_epi8_rpl(vIns, vNegInf, cond);
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
            tbl_pr[j-3] = (int32_t)_mm_extract_epi32_rpl(vWscore,0);
            mch_pr[j-3] = (int32_t)_mm_extract_epi32_rpl(vWmatch,0);
            sim_pr[j-3] = (int32_t)_mm_extract_epi32_rpl(vWsimilar,0);
            len_pr[j-3] = (int32_t)_mm_extract_epi32_rpl(vWlength,0);
            del_pr[j-3] = (int32_t)_mm_extract_epi32_rpl(vDel,0);
            /* as minor diagonal vector passes across table, extract
             * max values within the i,j bounds */
            {
                __m128i cond_valid_J = _mm_and_si128(
                        _mm_cmpgt_epi32(vJ, vNegOne),
                        _mm_cmplt_epi32(vJ, vJLimit));
                __m128i cond_valid_IJ = _mm_and_si128(cond_valid_J, vIltLimit);
                __m128i cond_eq = _mm_cmpeq_epi32(vWscore, vMaxScore);
                __m128i cond_max = _mm_cmpgt_epi32(vWscore, vMaxScore);
                __m128i cond_all = _mm_and_si128(cond_max, cond_valid_IJ);
                __m128i cond_Jlt = _mm_cmplt_epi32(vJ, vEndJ);
                vMaxScore = _mm_blendv_epi8_rpl(vMaxScore, vWscore, cond_all);
                vMaxMatch = _mm_blendv_epi8_rpl(vMaxMatch, vWmatch, cond_all);
                vMaxSimilar = _mm_blendv_epi8_rpl(vMaxSimilar, vWsimilar, cond_all);
                vMaxLength = _mm_blendv_epi8_rpl(vMaxLength, vWlength, cond_all);
                vEndI = _mm_blendv_epi8_rpl(vEndI, vI, cond_all);
                vEndJ = _mm_blendv_epi8_rpl(vEndJ, vJ, cond_all);
                cond_all = _mm_and_si128(cond_Jlt, cond_eq);
                cond_all = _mm_and_si128(cond_all, cond_valid_IJ);
                vMaxMatch = _mm_blendv_epi8_rpl(vMaxMatch, vWmatch, cond_all);
                vMaxSimilar = _mm_blendv_epi8_rpl(vMaxSimilar, vWsimilar, cond_all);
                vMaxLength = _mm_blendv_epi8_rpl(vMaxLength, vWlength, cond_all);
                vEndI = _mm_blendv_epi8_rpl(vEndI, vI, cond_all);
                vEndJ = _mm_blendv_epi8_rpl(vEndJ, vJ, cond_all);
            }
            vJ = _mm_add_epi32(vJ, vOne);
        }
        vI = _mm_add_epi32(vI, vN);
    }

    /* alignment ending position */
    {
        int32_t *t = (int32_t*)&vMaxScore;
        int32_t *m = (int32_t*)&vMaxMatch;
        int32_t *s = (int32_t*)&vMaxSimilar;
        int32_t *l = (int32_t*)&vMaxLength;
        int32_t *i = (int32_t*)&vEndI;
        int32_t *j = (int32_t*)&vEndJ;
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


