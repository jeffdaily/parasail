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

#define NEG_INF (INT64_MIN/(int64_t)(2))

static inline __m128i _mm_cmpgt_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? 0xFFFFFFFFFFFFFFFF : 0;
    A.v[1] = (A.v[1]>B.v[1]) ? 0xFFFFFFFFFFFFFFFF : 0;
    return A.m;
}

static inline __m128i _mm_max_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]>B.v[1]) ? A.v[1] : B.v[1];
    return A.m;
}

static inline __m128i _mm_cmplt_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]<B.v[0]) ? 0xFFFFFFFFFFFFFFFF : 0;
    A.v[1] = (A.v[1]<B.v[1]) ? 0xFFFFFFFFFFFFFFFF : 0;
    return A.m;
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
        array[(i+0)*s2Len + (j-0)] = (int64_t)_mm_extract_epi64(vWscore, 1);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[(i+1)*s2Len + (j-1)] = (int64_t)_mm_extract_epi64(vWscore, 0);
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
        row[j-0] = (int64_t)_mm_extract_epi64(vWscore, 1);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int64_t)_mm_extract_epi64(vWscore, 1);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int64_t)_mm_extract_epi64(vWscore, 0);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int64_t)_mm_extract_epi64(vWscore, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_nw_stats_table_diag_sse41_128_64
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_nw_stats_rowcol_diag_sse41_128_64
#else
#define FNAME parasail_nw_stats_diag_sse41_128_64
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    const int32_t N = 2; /* number of values in vector */
    const int32_t PAD = N-1;
    const int32_t PAD2 = PAD*2;
    const int32_t s1Len_PAD = s1Len+PAD;
    const int32_t s2Len_PAD = s2Len+PAD;
    int64_t * const restrict s1      = parasail_memalign_int64_t(16, s1Len+PAD);
    int64_t * const restrict s2B     = parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict _tbl_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict _del_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict _mch_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict _sim_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict _len_pr = parasail_memalign_int64_t(16, s2Len+PAD2);
    int64_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int64_t * const restrict tbl_pr = _tbl_pr+PAD;
    int64_t * const restrict del_pr = _del_pr+PAD;
    int64_t * const restrict mch_pr = _mch_pr+PAD;
    int64_t * const restrict sim_pr = _sim_pr+PAD;
    int64_t * const restrict len_pr = _len_pr+PAD;
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
    int64_t score = NEG_INF;
    int64_t matches = NEG_INF;
    int64_t similar = NEG_INF;
    int64_t length = NEG_INF;
    
    __m128i vNegInf = _mm_set1_epi64x(NEG_INF);
    __m128i vOpen = _mm_set1_epi64x(open);
    __m128i vGap  = _mm_set1_epi64x(gap);
    __m128i vZero = _mm_set1_epi64x(0);
    __m128i vOne = _mm_set1_epi64x(1);
    __m128i vN = _mm_set1_epi64x(N);
    __m128i vGapN = _mm_set1_epi64x(gap*N);
    __m128i vNegOne = _mm_set1_epi64x(-1);
    __m128i vI = _mm_set_epi64x(0,1);
    __m128i vJreset = _mm_set_epi64x(0,-1);
    __m128i vMaxScore = vNegInf;
    __m128i vMaxMatch = vNegInf;
    __m128i vMaxSimilar = vNegInf;
    __m128i vMaxLength = vNegInf;
    __m128i vILimit = _mm_set1_epi64x(s1Len);
    __m128i vILimit1 = _mm_sub_epi64(vILimit, vOne);
    __m128i vJLimit = _mm_set1_epi64x(s2Len);
    __m128i vJLimit1 = _mm_sub_epi64(vJLimit, vOne);
    __m128i vIBoundary = _mm_set_epi64x(
            -open-0*gap,
            -open-1*gap);

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
        tbl_pr[j] = -open - j*gap;
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
        __m128i vNscore = vNegInf;
        __m128i vNmatch = vZero;
        __m128i vNsimilar = vZero;
        __m128i vNlength = vZero;
        __m128i vWscore = vNegInf;
        __m128i vWmatch = vZero;
        __m128i vWsimilar = vZero;
        __m128i vWlength = vZero;
        __m128i vIns = vNegInf;
        __m128i vDel = vNegInf;
        __m128i vJ = vJreset;
        __m128i vs1 = _mm_set_epi64x(
                s1[i+0],
                s1[i+1]);
        __m128i vs2 = vNegInf;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        vNscore = _mm_srli_si128(vNscore, 8);
        vNscore = _mm_insert_epi64(vNscore, tbl_pr[-1], 1);
        vNmatch = _mm_srli_si128(vNmatch, 8);
        vNmatch = _mm_insert_epi64(vNmatch, 0, 1);
        vNsimilar = _mm_srli_si128(vNsimilar, 8);
        vNsimilar = _mm_insert_epi64(vNsimilar, 0, 1);
        vNlength = _mm_srli_si128(vNlength, 8);
        vNlength = _mm_insert_epi64(vNlength, 0, 1);
        vWscore = _mm_srli_si128(vWscore, 8);
        vWscore = _mm_insert_epi64(vWscore, -open - i*gap, 1);
        vWmatch = _mm_srli_si128(vWmatch, 8);
        vWmatch = _mm_insert_epi64(vWmatch, 0, 1);
        vWsimilar = _mm_srli_si128(vWsimilar, 8);
        vWsimilar = _mm_insert_epi64(vWsimilar, 0, 1);
        vWlength = _mm_srli_si128(vWlength, 8);
        vWlength = _mm_insert_epi64(vWlength, 0, 1);
        tbl_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            __m128i vNWmatch = vNmatch;
            __m128i vNWsimilar = vNsimilar;
            __m128i vNWlength = vNlength;
            vNscore = _mm_srli_si128(vWscore, 8);
            vNscore = _mm_insert_epi64(vNscore, tbl_pr[j], 1);
            vNmatch = _mm_srli_si128(vWmatch, 8);
            vNmatch = _mm_insert_epi64(vNmatch, mch_pr[j], 1);
            vNsimilar = _mm_srli_si128(vWsimilar, 8);
            vNsimilar = _mm_insert_epi64(vNsimilar, sim_pr[j], 1);
            vNlength = _mm_srli_si128(vWlength, 8);
            vNlength = _mm_insert_epi64(vNlength, len_pr[j], 1);
            vDel = _mm_srli_si128(vDel, 8);
            vDel = _mm_insert_epi64(vDel, del_pr[j], 1);
            vDel = _mm_max_epi64_rpl(
                    _mm_sub_epi64(vNscore, vOpen),
                    _mm_sub_epi64(vDel, vGap));
            vIns = _mm_max_epi64_rpl(
                    _mm_sub_epi64(vWscore, vOpen),
                    _mm_sub_epi64(vIns, vGap));
            vs2 = _mm_srli_si128(vs2, 8);
            vs2 = _mm_insert_epi64(vs2, s2[j], 1);
            vMat = _mm_set_epi64x(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]]
                    );
            vNWscore = _mm_add_epi64(vNWscore, vMat);
            vWscore = _mm_max_epi64_rpl(vNWscore, vIns);
            vWscore = _mm_max_epi64_rpl(vWscore, vDel);
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
                        _mm_cmplt_epi64_rpl(vNWscore,vDel),
                        _mm_cmplt_epi64_rpl(vNWscore,vIns));
                case2not = _mm_cmplt_epi64_rpl(vDel,vIns);
                case2 = _mm_andnot_si128(case2not,case1not);
                case3 = _mm_and_si128(case1not,case2not);
                vCmatch = _mm_andnot_si128(case1not,
                        _mm_add_epi64(vNWmatch, _mm_and_si128(
                                _mm_cmpeq_epi64(vs1,vs2),vOne)));
                vCmatch = _mm_or_si128(vCmatch, _mm_and_si128(case2, vNmatch));
                vCmatch = _mm_or_si128(vCmatch, _mm_and_si128(case3, vWmatch));
                vCsimilar = _mm_andnot_si128(case1not,
                        _mm_add_epi64(vNWsimilar, _mm_and_si128(
                                _mm_cmpgt_epi64_rpl(vMat,vZero),vOne)));
                vCsimilar = _mm_or_si128(vCsimilar, _mm_and_si128(case2, vNsimilar));
                vCsimilar = _mm_or_si128(vCsimilar, _mm_and_si128(case3, vWsimilar));
                vClength= _mm_andnot_si128(case1not,
                        _mm_add_epi64(vNWlength, vOne));
                vClength= _mm_or_si128(vClength,_mm_and_si128(case2,
                            _mm_add_epi64(vNlength, vOne)));
                vClength= _mm_or_si128(vClength,_mm_and_si128(case3,
                            _mm_add_epi64(vWlength, vOne)));
                vWmatch = vCmatch;
                vWsimilar = vCsimilar;
                vWlength = vClength;
            }

            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m128i cond = _mm_cmpeq_epi64(vJ,vNegOne);
                vWscore = _mm_blendv_epi8(vWscore, vIBoundary, cond);
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
            tbl_pr[j-1] = (int64_t)_mm_extract_epi64(vWscore,0);
            mch_pr[j-1] = (int64_t)_mm_extract_epi64(vWmatch,0);
            sim_pr[j-1] = (int64_t)_mm_extract_epi64(vWsimilar,0);
            len_pr[j-1] = (int64_t)_mm_extract_epi64(vWlength,0);
            del_pr[j-1] = (int64_t)_mm_extract_epi64(vDel,0);
            /* as minor diagonal vector passes across table, extract
               last table value at the i,j bound */
            {
                __m128i cond_valid_I = _mm_cmpeq_epi64(vI, vILimit1);
                __m128i cond_valid_J = _mm_cmpeq_epi64(vJ, vJLimit1);
                __m128i cond_all = _mm_and_si128(cond_valid_I, cond_valid_J);
                vMaxScore = _mm_blendv_epi8(vMaxScore, vWscore, cond_all);
                vMaxMatch = _mm_blendv_epi8(vMaxMatch, vWmatch, cond_all);
                vMaxSimilar = _mm_blendv_epi8(vMaxSimilar, vWsimilar, cond_all);
                vMaxLength = _mm_blendv_epi8(vMaxLength, vWlength, cond_all);
            }
            vJ = _mm_add_epi64(vJ, vOne);
        }
        vI = _mm_add_epi64(vI, vN);
        vIBoundary = _mm_sub_epi64(vIBoundary, vGapN);
    }

    /* max in vMaxScore */
    for (i=0; i<N; ++i) {
        int64_t value;
        value = (int64_t) _mm_extract_epi64(vMaxScore, 1);
        if (value > score) {
            score = value;
            matches = (int64_t) _mm_extract_epi64(vMaxMatch, 1);
            similar = (int64_t) _mm_extract_epi64(vMaxSimilar, 1);
            length= (int64_t) _mm_extract_epi64(vMaxLength, 1);
        }
        vMaxScore = _mm_slli_si128(vMaxScore, 8);
        vMaxMatch = _mm_slli_si128(vMaxMatch, 8);
        vMaxSimilar = _mm_slli_si128(vMaxSimilar, 8);
        vMaxLength = _mm_slli_si128(vMaxLength, 8);
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


