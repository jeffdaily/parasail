/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <immintrin.h>

#include "parasail/memory.h"
#include "parasail/internal_knc.h"

__m512i * parasail_memalign___m512i(size_t alignment, size_t size)
{
    return (__m512i *) parasail_memalign(alignment, size*sizeof(__m512i));
}

void parasail_memset___m512i(__m512i *b, __m512i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm512_store_epi32(&b[i], c);
    }
}

parasail_profile_t * parasail_profile_create_knc_512_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 16; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m512i* const restrict vProfile = parasail_memalign___m512i(64, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            __m512i_32_t t;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t.v[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                j += segLen;
            }
            _mm512_store_epi32(&vProfile[index], t.m);
            ++index;
        }
    }

    profile->profile = vProfile;
    profile->free = &parasail_profile_free___m512i;
    return profile;
}

parasail_profile_t * parasail_profile_create_stats_knc_512_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 8; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m512i* const restrict vProfile = parasail_memalign___m512i(64, n * segLen);
    __m512i* const restrict vProfileM = parasail_memalign___m512i(64, n * segLen);
    __m512i* const restrict vProfileS = parasail_memalign___m512i(64, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            __m512i_32_t p;
            __m512i_32_t m;
            __m512i_32_t s;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                p.v[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                m.v[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                s.v[segNum] = p.v[segNum] > 0;
                j += segLen;
            }
            _mm512_store_epi32(&vProfile[index], p.m);
            _mm512_store_epi32(&vProfileM[index], m.m);
            _mm512_store_epi32(&vProfileS[index], s.m);
            ++index;
        }
    }

    profile->profile = vProfile;
    profile->profile_m = vProfileM;
    profile->profile_s = vProfileS;
    profile->free = &parasail_profile_free___m512i;
    return profile;
}

void parasail_profile_free___m512i(void *profile)
{
    free((__m512i*)profile);
}

