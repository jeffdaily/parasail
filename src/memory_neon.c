/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include "parasail/memory.h"
#include "parasail/internal_neon.h"

simde__m128i * parasail_memalign_simde__m128i(size_t alignment, size_t size)
{
    return (simde__m128i *) parasail_memalign(alignment, size*sizeof(simde__m128i));
}

void parasail_memset_simde__m128i(simde__m128i *b, simde__m128i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        simde_mm_store_si128(&b[i], c);
    }
}

parasail_profile_t * parasail_profile_create_neon_128_8(
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
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i_private t_;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t_.i8[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], simde__m128i_from_private(t_));
            ++index;
        }
    }

    profile->profile8.score = vProfile;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_neon_128_16(
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
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i_private t_;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t_.i16[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], simde__m128i_from_private(t_));
            ++index;
        }
    }

    profile->profile16.score = vProfile;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 4; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i_private t_;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t_.i32[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], simde__m128i_from_private(t_));
            ++index;
        }
    }

    profile->profile32.score = vProfile;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 2; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i_private t_;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                t_.i64[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], simde__m128i_from_private(t_));
            ++index;
        }
    }

    profile->profile64.score = vProfile;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t* parasail_profile_create_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile8 = parasail_profile_create_neon_128_8(s1, s1Len, matrix);
    parasail_profile_t *profile16 = parasail_profile_create_neon_128_16(s1, s1Len, matrix);
    parasail_profile_t *profile32 = parasail_profile_create_neon_128_32(s1, s1Len, matrix);
    profile8->profile16 = profile16->profile16;
    profile8->profile32 = profile32->profile32;
    free(profile16);
    free(profile32);

    return profile8;
}

parasail_profile_t * parasail_profile_create_stats_neon_128_8(
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
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileM = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileS = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i_private p_;
            simde__m128i_private m_;
            simde__m128i_private s_;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                p_.i8[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                m_.i8[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                s_.i8[segNum] = p_.i8[segNum] > 0;
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], simde__m128i_from_private(p_));
            simde_mm_store_si128(&vProfileM[index], simde__m128i_from_private(m_));
            simde_mm_store_si128(&vProfileS[index], simde__m128i_from_private(s_));
            ++index;
        }
    }

    profile->profile8.score = vProfile;
    profile->profile8.matches = vProfileM;
    profile->profile8.similar = vProfileS;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_stats_neon_128_16(
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
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileM = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileS = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i_private p_;
            simde__m128i_private m_;
            simde__m128i_private s_;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                p_.i16[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                m_.i16[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                s_.i16[segNum] = p_.i16[segNum] > 0;
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], simde__m128i_from_private(p_));
            simde_mm_store_si128(&vProfileM[index], simde__m128i_from_private(m_));
            simde_mm_store_si128(&vProfileS[index], simde__m128i_from_private(s_));
            ++index;
        }
    }

    profile->profile16.score = vProfile;
    profile->profile16.matches = vProfileM;
    profile->profile16.similar = vProfileS;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_stats_neon_128_32(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 4; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileM = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileS = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i_private p_;
            simde__m128i_private m_;
            simde__m128i_private s_;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                p_.i32[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                m_.i32[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                s_.i32[segNum] = p_.i32[segNum] > 0;
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], simde__m128i_from_private(p_));
            simde_mm_store_si128(&vProfileM[index], simde__m128i_from_private(m_));
            simde_mm_store_si128(&vProfileS[index], simde__m128i_from_private(s_));
            ++index;
        }
    }

    profile->profile32.score = vProfile;
    profile->profile32.matches = vProfileM;
    profile->profile32.similar = vProfileS;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t * parasail_profile_create_stats_neon_128_64(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    const int32_t n = matrix->size; /* number of amino acids in table */
    const int32_t segWidth = 2; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    simde__m128i* const restrict vProfile = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileM = parasail_memalign_simde__m128i(16, n * segLen);
    simde__m128i* const restrict vProfileS = parasail_memalign_simde__m128i(16, n * segLen);
    int32_t index = 0;

    parasail_profile_t *profile = parasail_profile_new(s1, s1Len, matrix);

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            simde__m128i_private p_;
            simde__m128i_private m_;
            simde__m128i_private s_;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                p_.i64[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                m_.i64[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                s_.i64[segNum] = p_.i64[segNum] > 0;
                j += segLen;
            }
            simde_mm_store_si128(&vProfile[index], simde__m128i_from_private(p_));
            simde_mm_store_si128(&vProfileM[index], simde__m128i_from_private(m_));
            simde_mm_store_si128(&vProfileS[index], simde__m128i_from_private(s_));
            ++index;
        }
    }

    profile->profile64.score = vProfile;
    profile->profile64.matches = vProfileM;
    profile->profile64.similar = vProfileS;
    profile->free = &parasail_free_simde__m128i;
    return profile;
}

parasail_profile_t* parasail_profile_create_stats_neon_128_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile8 = parasail_profile_create_stats_neon_128_8(s1, s1Len, matrix);
    parasail_profile_t *profile16 = parasail_profile_create_stats_neon_128_16(s1, s1Len, matrix);
    parasail_profile_t *profile32 = parasail_profile_create_stats_neon_128_32(s1, s1Len, matrix);
    profile8->profile16 = profile16->profile16;
    profile8->profile32 = profile32->profile32;
    free(profile16);
    free(profile32);

    return profile8;
}

void parasail_free_simde__m128i(void *ptr)
{
    parasail_free((simde__m128i*)ptr);
}

