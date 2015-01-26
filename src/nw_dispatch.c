/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include "parasail.h"
#include "parasail_cpuid.h"
#include "parasail_internal.h"

/* typedef for the nw function */
typedef parasail_result_t* parasail_func(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24]);

/* forward declare the dispatcher function */
parasail_func parasail_nw_scan_dispatcher;
parasail_func parasail_nw_scan_32_dispatcher;
parasail_func parasail_nw_scan_16_dispatcher;
parasail_func parasail_nw_scan_8_dispatcher;
parasail_func parasail_nw_striped_dispatcher;
parasail_func parasail_nw_striped_32_dispatcher;
parasail_func parasail_nw_striped_16_dispatcher;
parasail_func parasail_nw_striped_8_dispatcher;
parasail_func parasail_nw_diag_dispatcher;
parasail_func parasail_nw_diag_32_dispatcher;
parasail_func parasail_nw_diag_16_dispatcher;
parasail_func parasail_nw_diag_8_dispatcher;

/* declare and initialize the nw pointer to the dispatcher function */
parasail_func * parasail_nw_scan_32_pointer = &parasail_nw_scan_32_dispatcher;
parasail_func * parasail_nw_scan_16_pointer = &parasail_nw_scan_16_dispatcher;
parasail_func * parasail_nw_scan_8_pointer  = &parasail_nw_scan_8_dispatcher;
parasail_func * parasail_nw_striped_32_pointer = &parasail_nw_striped_32_dispatcher;
parasail_func * parasail_nw_striped_16_pointer = &parasail_nw_striped_16_dispatcher;
parasail_func * parasail_nw_striped_8_pointer  = &parasail_nw_striped_8_dispatcher;
parasail_func * parasail_nw_diag_32_pointer = &parasail_nw_diag_32_dispatcher;
parasail_func * parasail_nw_diag_16_pointer = &parasail_nw_diag_16_dispatcher;
parasail_func * parasail_nw_diag_8_pointer  = &parasail_nw_diag_8_dispatcher;

/* dispatcher function implementation */
parasail_result_t* parasail_nw_scan_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#if HAVE_KNC
    parasail_nw_scan_32_pointer = nw_scan_knc_512_32;
#else
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_nw_scan_32_pointer = nw_scan_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_nw_scan_32_pointer = nw_scan_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_nw_scan_32_pointer = nw_scan_sse2_128_32;
    }
    else
#endif
#endif
    {
    }
    return parasail_nw_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* dispatcher function implementation */
parasail_result_t* parasail_nw_scan_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#if HAVE_KNC
    parasail_nw_scan_16_pointer = nw_scan_knc_512_32;
#else
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_nw_scan_16_pointer = nw_scan_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_nw_scan_16_pointer = nw_scan_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_nw_scan_16_pointer = nw_scan_sse2_128_16;
    }
    else
#endif
#endif
    {
    }
    return parasail_nw_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* dispatcher function implementation */
parasail_result_t* parasail_nw_scan_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#if HAVE_KNC
    parasail_nw_scan_8_pointer = nw_scan_knc_512_32;
#else
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_nw_scan_8_pointer = nw_scan_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_nw_scan_8_pointer = nw_scan_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_nw_scan_8_pointer = nw_scan_sse2_128_8;
    }
    else
#endif
#endif
    {
    }
    return parasail_nw_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* dispatcher function implementation */
parasail_result_t* parasail_nw_scan(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    parasail_result_t *result;

    result = parasail_nw_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

/* nw implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct nw impl */
parasail_result_t* parasail_nw_scan_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_nw_scan_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* nw implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct nw impl */
parasail_result_t* parasail_nw_scan_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_nw_scan_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* nw implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct nw impl */
parasail_result_t* parasail_nw_scan_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_nw_scan_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* dispatcher function implementation */
parasail_result_t* parasail_nw_striped_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#if HAVE_KNC
    parasail_nw_striped_32_pointer = nw_striped_knc_512_32;
#else
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_nw_striped_32_pointer = nw_striped_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_nw_striped_32_pointer = nw_striped_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_nw_striped_32_pointer = nw_striped_sse2_128_32;
    }
    else
#endif
#endif
    {
    }
    return parasail_nw_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* dispatcher function implementation */
parasail_result_t* parasail_nw_striped_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#if HAVE_KNC
    parasail_nw_striped_16_pointer = nw_striped_knc_512_32;
#else
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_nw_striped_16_pointer = nw_striped_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_nw_striped_16_pointer = nw_striped_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_nw_striped_16_pointer = nw_striped_sse2_128_16;
    }
    else
#endif
#endif
    {
    }
    return parasail_nw_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* dispatcher function implementation */
parasail_result_t* parasail_nw_striped_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#if HAVE_KNC
    parasail_nw_striped_8_pointer = nw_striped_knc_512_32;
#else
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_nw_striped_8_pointer = nw_striped_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_nw_striped_8_pointer = nw_striped_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_nw_striped_8_pointer = nw_striped_sse2_128_8;
    }
    else
#endif
#endif
    {
    }
    return parasail_nw_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* dispatcher function implementation */
parasail_result_t* parasail_nw_striped(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    parasail_result_t *result;

    result = parasail_nw_striped_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_striped_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_striped_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

/* nw implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct nw impl */
parasail_result_t* parasail_nw_striped_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_nw_striped_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* nw implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct nw impl */
parasail_result_t* parasail_nw_striped_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_nw_striped_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* nw implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct nw impl */
parasail_result_t* parasail_nw_striped_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_nw_striped_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* dispatcher function implementation */
parasail_result_t* parasail_nw_diag_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#if HAVE_KNC
    parasail_nw_diag_32_pointer = nw_diag_knc_512_32;
#else
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_nw_diag_32_pointer = nw_diag_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_nw_diag_32_pointer = nw_diag_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_nw_diag_32_pointer = nw_diag_sse2_128_32;
    }
    else
#endif
#endif
    {
    }
    return parasail_nw_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* dispatcher function implementation */
parasail_result_t* parasail_nw_diag_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#if HAVE_KNC
    parasail_nw_diag_16_pointer = nw_diag_knc_512_32;
#else
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_nw_diag_16_pointer = nw_diag_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_nw_diag_16_pointer = nw_diag_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_nw_diag_16_pointer = nw_diag_sse2_128_16;
    }
    else
#endif
#endif
    {
    }
    return parasail_nw_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* dispatcher function implementation */
parasail_result_t* parasail_nw_diag_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#if HAVE_KNC
    parasail_nw_diag_8_pointer = nw_diag_knc_512_32;
#else
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_nw_diag_8_pointer = nw_diag_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_nw_diag_8_pointer = nw_diag_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_nw_diag_8_pointer = nw_diag_sse2_128_8;
    }
    else
#endif
#endif
    {
    }
    return parasail_nw_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* dispatcher function implementation */
parasail_result_t* parasail_nw_diag(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    parasail_result_t *result;

    result = parasail_nw_diag_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_diag_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_nw_diag_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}

/* nw implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct nw impl */
parasail_result_t* parasail_nw_diag_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_nw_diag_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* nw implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct nw impl */
parasail_result_t* parasail_nw_diag_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_nw_diag_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* nw implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct nw impl */
parasail_result_t* parasail_nw_diag_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_nw_diag_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

