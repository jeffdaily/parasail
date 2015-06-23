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

#include "parasail.h"
#include "parasail/cpuid.h"

typedef parasail_profile_t* parasail_proffunc_t(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix);

/* forward declare the dispatcher functions */
parasail_proffunc_t parasail_profile_create_64_dispatcher;
parasail_proffunc_t parasail_profile_create_32_dispatcher;
parasail_proffunc_t parasail_profile_create_16_dispatcher;
parasail_proffunc_t parasail_profile_create_8_dispatcher;
parasail_proffunc_t parasail_profile_create_stats_64_dispatcher;
parasail_proffunc_t parasail_profile_create_stats_32_dispatcher;
parasail_proffunc_t parasail_profile_create_stats_16_dispatcher;
parasail_proffunc_t parasail_profile_create_stats_8_dispatcher;

/* declare and initialize the pointer to the dispatcher function */
parasail_proffunc_t * parasail_profile_create_64_pointer = parasail_profile_create_64_dispatcher;
parasail_proffunc_t * parasail_profile_create_32_pointer = parasail_profile_create_32_dispatcher;
parasail_proffunc_t * parasail_profile_create_16_pointer = parasail_profile_create_16_dispatcher;
parasail_proffunc_t * parasail_profile_create_8_pointer = parasail_profile_create_8_dispatcher;
parasail_proffunc_t * parasail_profile_create_stats_64_pointer = parasail_profile_create_stats_64_dispatcher;
parasail_proffunc_t * parasail_profile_create_stats_32_pointer = parasail_profile_create_stats_32_dispatcher;
parasail_proffunc_t * parasail_profile_create_stats_16_pointer = parasail_profile_create_stats_16_dispatcher;
parasail_proffunc_t * parasail_profile_create_stats_8_pointer = parasail_profile_create_stats_8_dispatcher;

/* dispatcher function implementations */

parasail_profile_t* parasail_profile_create_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
#if HAVE_KNC
    if (1) {
        parasail_profile_create_64_pointer = parasail_profile_create_knc_512_32;
    }
    else
#endif
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_profile_create_64_pointer = parasail_profile_create_avx_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_profile_create_64_pointer = parasail_profile_create_sse_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_profile_create_64_pointer = parasail_profile_create_sse_128_64;
    }
    else
#endif
    {
        /* no fallback */
        parasail_profile_create_64_pointer = NULL;
    }
    return parasail_profile_create_64_pointer(s1, s1Len, matrix);
}

parasail_profile_t* parasail_profile_create_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
#if HAVE_KNC
    if (1) {
        parasail_profile_create_32_pointer = parasail_profile_create_knc_512_32;
    }
    else
#endif
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_profile_create_32_pointer = parasail_profile_create_avx_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_profile_create_32_pointer = parasail_profile_create_sse_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_profile_create_32_pointer = parasail_profile_create_sse_128_32;
    }
    else
#endif
    {
        /* no fallback */
        parasail_profile_create_32_pointer = NULL;
    }
    return parasail_profile_create_32_pointer(s1, s1Len, matrix);
}

parasail_profile_t* parasail_profile_create_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
#if HAVE_KNC
    if (1) {
        parasail_profile_create_16_pointer = parasail_profile_create_knc_512_32;
    }
    else
#endif
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_profile_create_16_pointer = parasail_profile_create_avx_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_profile_create_16_pointer = parasail_profile_create_sse_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_profile_create_16_pointer = parasail_profile_create_sse_128_16;
    }
    else
#endif
    {
        /* no fallback */
        parasail_profile_create_16_pointer = NULL;
    }
    return parasail_profile_create_16_pointer(s1, s1Len, matrix);
}

parasail_profile_t* parasail_profile_create_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
#if HAVE_KNC
    if (1) {
        parasail_profile_create_8_pointer = parasail_profile_create_knc_512_32;
    }
    else
#endif
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_profile_create_8_pointer = parasail_profile_create_avx_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_profile_create_8_pointer = parasail_profile_create_sse_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_profile_create_8_pointer = parasail_profile_create_sse_128_8;
    }
    else
#endif
    {
        /* no fallback */
        parasail_profile_create_8_pointer = NULL;
    }
    return parasail_profile_create_8_pointer(s1, s1Len, matrix);
}

parasail_profile_t* parasail_profile_create_stats_64_dispatcher(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
#if HAVE_KNC
    if (1) {
        parasail_profile_create_stats_64_pointer = parasail_profile_create_stats_knc_512_32;
    }
    else
#endif
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_profile_create_stats_64_pointer = parasail_profile_create_stats_avx_256_64;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_profile_create_stats_64_pointer = parasail_profile_create_stats_sse_128_64;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_profile_create_stats_64_pointer = parasail_profile_create_stats_sse_128_64;
    }
    else
#endif
    {
        /* no fallback */
        parasail_profile_create_stats_64_pointer = NULL;
    }
    return parasail_profile_create_stats_64_pointer(s1, s1Len, matrix);
}

parasail_profile_t* parasail_profile_create_stats_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
#if HAVE_KNC
    if (1) {
        parasail_profile_create_stats_32_pointer = parasail_profile_create_stats_knc_512_32;
    }
    else
#endif
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_profile_create_stats_32_pointer = parasail_profile_create_stats_avx_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_profile_create_stats_32_pointer = parasail_profile_create_stats_sse_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_profile_create_stats_32_pointer = parasail_profile_create_stats_sse_128_32;
    }
    else
#endif
    {
        /* no fallback */
        parasail_profile_create_stats_32_pointer = NULL;
    }
    return parasail_profile_create_stats_32_pointer(s1, s1Len, matrix);
}

parasail_profile_t* parasail_profile_create_stats_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
#if HAVE_KNC
    if (1) {
        parasail_profile_create_stats_16_pointer = parasail_profile_create_stats_knc_512_32;
    }
    else
#endif
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_profile_create_stats_16_pointer = parasail_profile_create_stats_avx_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_profile_create_stats_16_pointer = parasail_profile_create_stats_sse_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_profile_create_stats_16_pointer = parasail_profile_create_stats_sse_128_16;
    }
    else
#endif
    {
        /* no fallback */
        parasail_profile_create_stats_16_pointer = NULL;
    }
    return parasail_profile_create_stats_16_pointer(s1, s1Len, matrix);
}

parasail_profile_t* parasail_profile_create_stats_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
#if HAVE_KNC
    if (1) {
        parasail_profile_create_stats_8_pointer = parasail_profile_create_stats_knc_512_32;
    }
    else
#endif
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_profile_create_stats_8_pointer = parasail_profile_create_stats_avx_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_profile_create_stats_8_pointer = parasail_profile_create_stats_sse_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_profile_create_stats_8_pointer = parasail_profile_create_stats_sse_128_8;
    }
    else
#endif
    {
        /* no fallback */
        parasail_profile_create_stats_8_pointer = NULL;
    }
    return parasail_profile_create_stats_8_pointer(s1, s1Len, matrix);
}

