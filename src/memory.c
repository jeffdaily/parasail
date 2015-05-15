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
#include <ctype.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "parasail.h"
#include "parasail/memory.h"

void* parasail_memalign(size_t alignment, size_t size)
{
    void *ptr = NULL;
#if defined(HAVE__ALIGNED_MALLOC)
    ptr = _aligned_malloc(size, alignment);
#elif defined(HAVE_POSIX_MEMALIGN)
    int retcode = posix_memalign(&ptr, alignment, size);
    assert(0 == retcode);
#elif defined(HAVE_ALIGNED_ALLOC)
    ptr = aligned_alloc(alignment, size);
#elif defined(HAVE_MEMALIGN)
    ptr = memalign(alignment, size);
#else
#error "No suitable memory alignment routine found."
#endif
    assert(NULL != ptr);
    return ptr;
}

void parasail_free(void *ptr)
{
#if defined(HAVE__ALIGNED_MALLOC)
     _aligned_free(ptr);
#else
    free(ptr);
#endif
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

int64_t * parasail_memalign_int64_t(size_t alignment, size_t size)
{
    return (int64_t *) parasail_memalign(alignment, size*sizeof(int64_t));
}

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

void parasail_memset_int64_t(int64_t *b, int64_t c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

parasail_result_t* parasail_result_new()
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    result = (parasail_result_t*)malloc(sizeof(parasail_result_t));
    assert(result);

    result->saturated = 0;
    result->score = 0;
    result->matches = 0;
    result->similar = 0;
    result->length = 0;
    result->score_table = NULL;
    result->matches_table = NULL;
    result->similar_table = NULL;
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
    result->similar_table = (int *)malloc(sizeof(int)*a*b);
    assert(result->similar_table);
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
    if (NULL != result->similar_table) free(result->similar_table);
    if (NULL != result->length_table) free(result->length_table);
    free(result);
}

void parasail_version(int *major, int *minor, int *patch)
{
    *major = PARASAIL_VERSION_MAJOR;
    *minor = PARASAIL_VERSION_MINOR;
    *patch = PARASAIL_VERSION_PATCH;
}

#include <stdio.h>
parasail_matrix_t* parasail_matrix_create(
        const char *alphabet, const int match, const int mismatch)
{
    parasail_matrix_t *retval = NULL;
    int *matrix = NULL;
    int *mapper = NULL;
    size_t size = 0;
    size_t size1 = 0;
    size_t i = 0;
    size_t j = 0;
    size_t c = 0;

    size = strlen(alphabet);
    size1 = size + 1;

    matrix = (int*)malloc(sizeof(int)*size1*size1);
    assert(matrix);
    for (i=0; i<size; ++i) {
        for (j=0; j<size; ++j) {
            if (i == j) {
                matrix[c++] = match;
            }
            else {
                matrix[c++] = mismatch;
            }
        }
        matrix[c++] = mismatch;
    }
    for (j=0; j<size1; ++j) {
        matrix[c++] = mismatch;
    }

    mapper = (int*)malloc(sizeof(int)*256);
    assert(mapper);
    parasail_memset_int(mapper, size, 256);
    for (i=0; i<size; ++i) {
        mapper[toupper((unsigned char)alphabet[i])] = (int)i;
        mapper[tolower((unsigned char)alphabet[i])] = (int)i;
    }

    retval = (parasail_matrix_t*)malloc(sizeof(parasail_matrix_t));
    assert(retval);
    retval->name = "";
    retval->matrix = matrix;
    retval->mapper = mapper;
    retval->size = (int)size;
    retval->need_free = 1;
    return retval;
}

void parasail_matrix_free(parasail_matrix_t *matrix)
{
    /* validate inputs */
    assert(NULL != matrix);
    if (matrix->need_free) {
        free((void*)matrix->matrix);
        free((void*)matrix->mapper);
        free(matrix);
    }
    else {
        fprintf(stderr, "attempted to free built-in matrix '%s'\n",
                matrix->name);
    }
}

