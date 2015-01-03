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

#include <assert.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if HAVE_SSE2 || HAVE_SSE41
#include <emmintrin.h>
#include "parasail_internal_sse.h"
#endif

#if HAVE_AVX2
#include <immintrin.h>
#include "parasail_internal_avx.h"
#endif

#include "parasail.h"
#include "parasail_internal.h"

void* parasail_memalign(size_t alignment, size_t size)
{
    void *ptr = NULL;
#ifdef HAVE_POSIX_MEMALIGN
    int retcode = posix_memalign(&ptr, alignment, size);
    assert(0 == retcode);
#elif defined(HAVE_ALIGNED_ALLOC)
    ptr = aligned_alloc(alignment, size);
#elif defined(HAVE_MEMALIGN)
    ptr = memalign(alignment, size);
#else
    ptr = malloc(size);
#endif
    assert(NULL != ptr);
    return ptr;
}

int * parasail_memalign_int(size_t alignment, size_t size)
{
    return (int *) parasail_memalign(alignment, size*sizeof(int));
}

int8_t * parasail_memalign_int8_t(size_t alignment, size_t size)
{
    return (int8_t *) parasail_memalign(alignment, size*sizeof(int8_t));
}

int16_t * parasail_memalign_int16_t(size_t alignment, size_t size)
{
    return (int16_t *) parasail_memalign(alignment, size*sizeof(int16_t));
}

int32_t * parasail_memalign_int32_t(size_t alignment, size_t size)
{
    return (int32_t *) parasail_memalign(alignment, size*sizeof(int32_t));
}

#if HAVE_SSE2 || HAVE_SSE41
__m128i * parasail_memalign_m128i(size_t alignment, size_t size)
{
    return (__m128i *) parasail_memalign(alignment, size*sizeof(__m128i));
}
#endif

#if HAVE_AVX2
__m256i * parasail_memalign_m256i(size_t alignment, size_t size)
{
    return (__m256i *) parasail_memalign(alignment, size*sizeof(__m256i));
}
#endif

void parasail_memset(void *b, int c, size_t len)
{
    (void)memset(b, c, len);
}

void parasail_memset_int(int *b, int c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

void parasail_memset_int8_t(int8_t *b, int8_t c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

void parasail_memset_int16_t(int16_t *b, int16_t c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

void parasail_memset_int32_t(int32_t *b, int32_t c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

#if HAVE_SSE2 || HAVE_SSE41
void parasail_memset_m128i(__m128i *b, __m128i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm_store_si128(&b[i], c);
    }
}
#endif

#if HAVE_AVX2
void parasail_memset_m256i(__m256i *b, __m256i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm256_store_si256(&b[i], c);
    }
}
#endif

parasail_result_t* parasail_result_new()
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    result = (parasail_result_t*)malloc(sizeof(parasail_result_t));
    assert(result);

    result->score = 0;
    result->matches = 0;
    result->length = 0;
    result->score_table = NULL;
    result->matches_table = NULL;
    result->length_table = NULL;

    return result;
}

parasail_result_t* parasail_result_new_table1(const int a, const int b)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    /* validate inputs */
    assert(a > 0);
    assert(b > 0);
    
    /* allocate struct to hold memory */
    result = parasail_result_new();

    /* allocate only score table */
    result->score_table = (int *)malloc(sizeof(int)*a*b);
    assert(result->score_table);

    return result;
}

parasail_result_t* parasail_result_new_table3(const int a, const int b)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    /* validate inputs */
    assert(a > 0);
    assert(b > 0);
    
    /* allocate struct to hold memory */
    result = parasail_result_new_table1(a, b);
    
    result->matches_table = (int *)malloc(sizeof(int)*a*b);
    assert(result->matches_table);
    result->length_table = (int *)malloc(sizeof(int)*a*b);
    assert(result->length_table);

    return result;
}

void parasail_result_free(parasail_result_t *result)
{
    /* validate inputs */
    assert(NULL != result);
    
    if (NULL != result->score_table) free(result->score_table);
    if (NULL != result->matches_table) free(result->matches_table);
    if (NULL != result->length_table) free(result->length_table);
    free(result);
}

