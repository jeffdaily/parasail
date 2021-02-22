/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>

#include <immintrin.h>

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_avx.h"


#if HAVE_AVX2_MM256_INSERT_EPI64
#define _mm256_insert_epi64_rpl _mm256_insert_epi64
#else
static inline __m256i _mm256_insert_epi64_rpl(__m256i a, int64_t i, int imm) {
    __m256i_64_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif

#if HAVE_AVX2_MM256_SET1_EPI64X
#define _mm256_set1_epi64x_rpl _mm256_set1_epi64x
#else
static inline __m256i _mm256_set1_epi64x_rpl(int64_t i) {
    __m256i_64_t A;
    A.v[0] = i;
    A.v[1] = i;
    A.v[2] = i;
    A.v[3] = i;
    return A.m;
}
#endif

static inline __m256i _mm256_max_epi64_rpl(__m256i a, __m256i b) {
    __m256i_64_t A;
    __m256i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]>B.v[1]) ? A.v[1] : B.v[1];
    A.v[2] = (A.v[2]>B.v[2]) ? A.v[2] : B.v[2];
    A.v[3] = (A.v[3]>B.v[3]) ? A.v[3] : B.v[3];
    return A.m;
}

#if HAVE_AVX2_MM256_SET_EPI64X
#define _mm256_set_epi64x_rpl _mm256_set_epi64x
#else
static inline __m256i _mm256_set_epi64x_rpl(int64_t e3, int64_t e2, int64_t e1, int64_t e0) {
    __m256i_64_t A;
    A.v[0] = e0;
    A.v[1] = e1;
    A.v[2] = e2;
    A.v[3] = e3;
    return A.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI64
#define _mm256_extract_epi64_rpl _mm256_extract_epi64
#else
static inline int64_t _mm256_extract_epi64_rpl(__m256i a, int imm) {
    __m256i_64_t A;
    A.m = a;
    return A.v[imm];
}
#endif

static inline __m256i _mm256_min_epi64_rpl(__m256i a, __m256i b) {
    __m256i_64_t A;
    __m256i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]<B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]<B.v[1]) ? A.v[1] : B.v[1];
    A.v[2] = (A.v[2]<B.v[2]) ? A.v[2] : B.v[2];
    A.v[3] = (A.v[3]<B.v[3]) ? A.v[3] : B.v[3];
    return A.m;
}

#define _mm256_cmplt_epi64_rpl(a,b) _mm256_cmpgt_epi64(b,a)

#define _mm256_srli_si256_rpl(a,imm) _mm256_or_si256(_mm256_slli_si256(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(3,0,0,1)), 16-imm), _mm256_srli_si256(a, imm))


static inline void arr_store_si256(
        int8_t *array,
        __m256i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[1LL*(i+0)*s2Len + (j-0)] = (int8_t)_mm256_extract_epi64_rpl(vWH, 3);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[1LL*(i+1)*s2Len + (j-1)] = (int8_t)_mm256_extract_epi64_rpl(vWH, 2);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[1LL*(i+2)*s2Len + (j-2)] = (int8_t)_mm256_extract_epi64_rpl(vWH, 1);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[1LL*(i+3)*s2Len + (j-3)] = (int8_t)_mm256_extract_epi64_rpl(vWH, 0);
    }
}

#define FNAME parasail_sw_trace_diag_avx2_256_64

parasail_result_t* FNAME(
        const char * const restrict _s1, const int _s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    /* declare local variables */
    int32_t N = 0;
    int32_t PAD = 0;
    int32_t PAD2 = 0;
    int32_t s1Len_PAD = 0;
    int32_t s2Len_PAD = 0;
    int64_t * restrict s1 = NULL;
    int64_t * restrict s2B = NULL;
    int64_t * restrict _H_pr = NULL;
    int64_t * restrict _F_pr = NULL;
    int64_t * restrict s2 = NULL;
    int64_t * restrict H_pr = NULL;
    int64_t * restrict F_pr = NULL;
    parasail_result_t *result = NULL;
    int32_t i = 0;
    int32_t j = 0;
    int32_t s1Len = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    int64_t NEG_LIMIT = 0;
    int64_t POS_LIMIT = 0;
    int64_t score = 0;
    __m256i vNegLimit;
    __m256i vPosLimit;
    __m256i vSaturationCheckMin;
    __m256i vSaturationCheckMax;
    __m256i vNegInf;
    __m256i vNegInf0;
    __m256i vOpen;
    __m256i vGap;
    __m256i vZero;
    __m256i vOne;
    __m256i vN;
    __m256i vNegOne;
    __m256i vI;
    __m256i vJreset;
    __m256i vMaxH;
    __m256i vEndI;
    __m256i vEndJ;
    __m256i vILimit;
    __m256i vJLimit;
    __m256i vTDiag;
    __m256i vTIns;
    __m256i vTDel;
    __m256i vTZero;
    __m256i vTDiagE;
    __m256i vTInsE;
    __m256i vTDiagF;
    __m256i vTDelF;

    /* validate inputs */
    PARASAIL_CHECK_NULL(_s2);
    PARASAIL_CHECK_GT0(s2Len);
    PARASAIL_CHECK_GE0(open);
    PARASAIL_CHECK_GE0(gap);
    PARASAIL_CHECK_NULL(matrix);
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        PARASAIL_CHECK_NULL(_s1);
        PARASAIL_CHECK_GT0(_s1Len);
    }
        
    /* initialize stack variables */
    N = 4; /* number of values in vector */
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
    vNegLimit = _mm256_set1_epi64x_rpl(NEG_LIMIT);
    vPosLimit = _mm256_set1_epi64x_rpl(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;
    vNegInf = _mm256_set1_epi64x_rpl(NEG_LIMIT);
    vNegInf0 = _mm256_srli_si256_rpl(vNegInf, 8); /* shift in a 0 */
    vOpen = _mm256_set1_epi64x_rpl(open);
    vGap  = _mm256_set1_epi64x_rpl(gap);
    vZero = _mm256_set1_epi64x_rpl(0);
    vOne = _mm256_set1_epi64x_rpl(1);
    vN = _mm256_set1_epi64x_rpl(N);
    vNegOne = _mm256_set1_epi64x_rpl(-1);
    vI = _mm256_set_epi64x_rpl(0,1,2,3);
    vJreset = _mm256_set_epi64x_rpl(0,-1,-2,-3);
    vMaxH = vNegInf;
    vEndI = vNegInf;
    vEndJ = vNegInf;
    vILimit = _mm256_set1_epi64x_rpl(s1Len);
    vJLimit = _mm256_set1_epi64x_rpl(s2Len);
    vTDiag = _mm256_set1_epi64x_rpl(PARASAIL_DIAG);
    vTIns = _mm256_set1_epi64x_rpl(PARASAIL_INS);
    vTDel = _mm256_set1_epi64x_rpl(PARASAIL_DEL);
    vTZero = _mm256_set1_epi64x_rpl(PARASAIL_ZERO);
    vTDiagE = _mm256_set1_epi64x_rpl(PARASAIL_DIAG_E);
    vTInsE = _mm256_set1_epi64x_rpl(PARASAIL_INS_E);
    vTDiagF = _mm256_set1_epi64x_rpl(PARASAIL_DIAG_F);
    vTDelF = _mm256_set1_epi64x_rpl(PARASAIL_DEL_F);

    /* initialize result */
    result = parasail_result_new_trace(s1Len, s2Len, 32, sizeof(int8_t));
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_SW | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_64 | PARASAIL_FLAG_LANES_4;

    /* initialize heap variables */
    s2B= parasail_memalign_int64_t(32, s2Len+PAD2);
    _H_pr = parasail_memalign_int64_t(32, s2Len+PAD2);
    _F_pr = parasail_memalign_int64_t(32, s2Len+PAD2);
    s2 = s2B+PAD; /* will allow later for negative indices */
    H_pr = _H_pr+PAD;
    F_pr = _F_pr+PAD;

    /* validate heap variables */
    if (!s2B) return NULL;
    if (!_H_pr) return NULL;
    if (!_F_pr) return NULL;

    /* convert _s1 from char to int in range 0-23 */
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        s1 = parasail_memalign_int64_t(32, s1Len+PAD);
        if (!s1) return NULL;
        for (i=0; i<s1Len; ++i) {
            s1[i] = matrix->mapper[(unsigned char)_s1[i]];
        }
        /* pad back of s1 with dummy values */
        for (i=s1Len; i<s1Len_PAD; ++i) {
            s1[i] = 0; /* point to first matrix row because we don't care */
        }
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
        F_pr[j] = NEG_LIMIT;
    }
    /* pad front of stored row values */
    for (j=-PAD; j<0; ++j) {
        H_pr[j] = NEG_LIMIT;
        F_pr[j] = NEG_LIMIT;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        H_pr[j] = NEG_LIMIT;
        F_pr[j] = NEG_LIMIT;
    }

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        __m256i vNH = vNegInf0;
        __m256i vWH = vNegInf0;
        __m256i vE = vNegInf;
        __m256i vE_opn = vNegInf;
        __m256i vE_ext = vNegInf;
        __m256i vF = vNegInf;
        __m256i vF_opn = vNegInf;
        __m256i vF_ext = vNegInf;
        __m256i vJ = vJreset;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+0] : ((i+0 >= s1Len) ? s1Len-1 : i+0))];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+1] : ((i+1 >= s1Len) ? s1Len-1 : i+1))];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+2] : ((i+2 >= s1Len) ? s1Len-1 : i+2))];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size * ((matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) ? s1[i+3] : ((i+3 >= s1Len) ? s1Len-1 : i+3))];
        __m256i vIltLimit = _mm256_cmplt_epi64_rpl(vI, vILimit);
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m256i vMat;
            __m256i vNWH = vNH;
            vNH = _mm256_srli_si256_rpl(vWH, 8);
            vNH = _mm256_insert_epi64_rpl(vNH, H_pr[j], 3);
            vF = _mm256_srli_si256_rpl(vF, 8);
            vF = _mm256_insert_epi64_rpl(vF, F_pr[j], 3);
            vF_opn = _mm256_sub_epi64(vNH, vOpen);
            vF_ext = _mm256_sub_epi64(vF, vGap);
            vF = _mm256_max_epi64_rpl(vF_opn, vF_ext);
            vE_opn = _mm256_sub_epi64(vWH, vOpen);
            vE_ext = _mm256_sub_epi64(vE, vGap);
            vE = _mm256_max_epi64_rpl(vE_opn, vE_ext);
            vMat = _mm256_set_epi64x_rpl(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]]
                    );
            vNWH = _mm256_add_epi64(vNWH, vMat);
            vNWH = _mm256_max_epi64_rpl(vNWH, vZero);
            vWH = _mm256_max_epi64_rpl(vNWH, vE);
            vWH = _mm256_max_epi64_rpl(vWH, vF);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m256i cond = _mm256_cmpeq_epi64(vJ,vNegOne);
                vWH = _mm256_andnot_si256(cond, vWH);
                vF = _mm256_blendv_epi8(vF, vNegInf, cond);
                vE = _mm256_blendv_epi8(vE, vNegInf, cond);
            }
            /* cannot start checking sat until after J clears boundary */
            if (j > PAD) {
                vSaturationCheckMin = _mm256_min_epi64_rpl(vSaturationCheckMin, vWH);
                vSaturationCheckMax = _mm256_max_epi64_rpl(vSaturationCheckMax, vWH);
            }
            /* trace table */
            {
                __m256i cond_zero = _mm256_cmpeq_epi64(vWH, vZero);
                __m256i case1 = _mm256_cmpeq_epi64(vWH, vNWH);
                __m256i case2 = _mm256_cmpeq_epi64(vWH, vF);
                __m256i vT = _mm256_blendv_epi8(
                        _mm256_blendv_epi8(vTIns, vTDel, case2),
                        _mm256_blendv_epi8(vTDiag, vTZero, cond_zero),
                        case1);
                __m256i condE = _mm256_cmpgt_epi64(vE_opn, vE_ext);
                __m256i condF = _mm256_cmpgt_epi64(vF_opn, vF_ext);
                __m256i vET = _mm256_blendv_epi8(vTInsE, vTDiagE, condE);
                __m256i vFT = _mm256_blendv_epi8(vTDelF, vTDiagF, condF);
                vT = _mm256_or_si256(vT, vET);
                vT = _mm256_or_si256(vT, vFT);
                arr_store_si256(result->trace->trace_table, vT, i, s1Len, j, s2Len);
            }
            H_pr[j-3] = (int64_t)_mm256_extract_epi64_rpl(vWH,0);
            F_pr[j-3] = (int64_t)_mm256_extract_epi64_rpl(vF,0);
            /* as minor diagonal vector passes across table, extract
             * max values within the i,j bounds */
            {
                __m256i cond_valid_J = _mm256_and_si256(
                        _mm256_cmpgt_epi64(vJ, vNegOne),
                        _mm256_cmplt_epi64_rpl(vJ, vJLimit));
                __m256i cond_valid_IJ = _mm256_and_si256(cond_valid_J, vIltLimit);
                __m256i cond_eq = _mm256_cmpeq_epi64(vWH, vMaxH);
                __m256i cond_max = _mm256_cmpgt_epi64(vWH, vMaxH);
                __m256i cond_all = _mm256_and_si256(cond_max, cond_valid_IJ);
                __m256i cond_Jlt = _mm256_cmplt_epi64_rpl(vJ, vEndJ);
                vMaxH = _mm256_blendv_epi8(vMaxH, vWH, cond_all);
                vEndI = _mm256_blendv_epi8(vEndI, vI, cond_all);
                vEndJ = _mm256_blendv_epi8(vEndJ, vJ, cond_all);
                cond_all = _mm256_and_si256(cond_Jlt, cond_eq);
                cond_all = _mm256_and_si256(cond_all, cond_valid_IJ);
                vEndI = _mm256_blendv_epi8(vEndI, vI, cond_all);
                vEndJ = _mm256_blendv_epi8(vEndJ, vJ, cond_all);
            }
            vJ = _mm256_add_epi64(vJ, vOne);
        }
        vI = _mm256_add_epi64(vI, vN);
    }

    /* alignment ending position */
    {
        int64_t *t = (int64_t*)&vMaxH;
        int64_t *i = (int64_t*)&vEndI;
        int64_t *j = (int64_t*)&vEndJ;
        int32_t k;
        for (k=0; k<N; ++k, ++t, ++i, ++j) {
            if (*t > score) {
                score = *t;
                end_query = *i;
                end_ref = *j;
            }
            else if (*t == score) {
                if (*j < end_ref) {
                    end_query = *i;
                    end_ref = *j;
                }
                else if (*j == end_ref && *i < end_query) {
                    end_query = *i;
                    end_ref = *j;
                }
            }
        }
    }

    if (_mm256_movemask_epi8(_mm256_or_si256(
            _mm256_cmplt_epi64_rpl(vSaturationCheckMin, vNegLimit),
            _mm256_cmpgt_epi64(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(_F_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        parasail_free(s1);
    }

    return result;
}


