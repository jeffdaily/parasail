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

/* typedef for the sw function */
typedef parasail_result_t* parasail_func(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24]);

/* forward declare the dispatcher function */
static parasail_func parasail_sw_dispatcher;
static parasail_func parasail_sw_32_dispatcher;
static parasail_func parasail_sw_16_dispatcher;
static parasail_func parasail_sw_8_dispatcher;

/* forward declare the vector impls, if needed */
#if HAVE_AVX2
static parasail_func sw_scan_avx2;
#endif
#if HAVE_SSE41
static parasail_func sw_scan_sse41;
#endif
#if HAVE_SSE2
static parasail_func sw_scan_sse2;
#endif

/* declare and initialize the sw pointer to the dispatcher function */
static parasail_func * parasail_sw_pointer = &parasail_sw_dispatcher;
static parasail_func * parasail_sw_32_pointer = &parasail_sw_32_dispatcher;
static parasail_func * parasail_sw_16_pointer = &parasail_sw_16_dispatcher;
static parasail_func * parasail_sw_8_pointer = &parasail_sw_8_dispatcher;

/* dispatcher function implementation */
static parasail_result_t* parasail_sw_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#if HAVE_KNC
    parasail_sw_pointer = sw_scan_knc_512_32;
#else
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sw_pointer = sw_scan_avx2;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sw_pointer = sw_scan_sse41;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sw_pointer = sw_scan_sse2;
    }
    else
#endif
#endif
    {
    }
    return parasail_sw_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* dispatcher function implementation */
static parasail_result_t* parasail_sw_32_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#if HAVE_KNC
    parasail_sw_pointer = sw_scan_knc_512_32;
#else
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sw_pointer = sw_scan_avx2_256_32;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sw_pointer = sw_scan_sse41_128_32;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sw_pointer = sw_scan_sse2_128_32;
    }
    else
#endif
#endif
    {
    }
    return parasail_sw_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* dispatcher function implementation */
static parasail_result_t* parasail_sw_16_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#if HAVE_KNC
    parasail_sw_pointer = sw_scan_knc_512_32;
#else
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sw_pointer = sw_scan_avx2_256_16;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sw_pointer = sw_scan_sse41_128_16;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sw_pointer = sw_scan_sse2_128_16;
    }
    else
#endif
#endif
    {
    }
    return parasail_sw_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* dispatcher function implementation */
static parasail_result_t* parasail_sw_8_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#if HAVE_KNC
    parasail_sw_pointer = sw_scan_knc_512_32;
#else
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        parasail_sw_pointer = sw_scan_avx2_256_8;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        parasail_sw_pointer = sw_scan_sse41_128_8;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        parasail_sw_pointer = sw_scan_sse2_128_8;
    }
    else
#endif
#endif
    {
    }
    return parasail_sw_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

#if HAVE_AVX2
/* sw vec impl which checks for saturation */
static parasail_result_t* sw_scan_avx2(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    parasail_result_t *result;

    result = sw_scan_avx2_256_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = sw_scan_avx2_256_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = sw_scan_avx2_256_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    return result;
}
#endif

#if HAVE_SSE41
/* sw vec impl which checks for saturation */
static parasail_result_t* sw_scan_sse41(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    parasail_result_t *result;

    result = sw_scan_sse41_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = sw_scan_sse41_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = sw_scan_sse41_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    return result;
}
#endif

#if HAVE_SSE2
/* sw vec impl which checks for saturation */
static parasail_result_t* sw_scan_sse2(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    parasail_result_t *result;

    result = sw_scan_sse2_128_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = sw_scan_sse2_128_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = sw_scan_sse2_128_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    return result;
}
#endif

/* sw implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct sw impl */
parasail_result_t* parasail_sw(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_sw_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* sw implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct sw impl */
parasail_result_t* parasail_sw_32(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_sw_32_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* sw implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct sw impl */
parasail_result_t* parasail_sw_16(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_sw_16_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

/* sw implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct sw impl */
parasail_result_t* parasail_sw_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_sw_8_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}

