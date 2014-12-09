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
#include <malloc.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include "parasail.h"
#include "parasail_internal.h"

inline void* parasail_memalign(size_t alignment, size_t size)
{
    void *ptr = NULL;
#if HAVE_POSIX_MEMALIGN
    int retcode = posix_memalign(&ptr, alignment, size);
    assert(0 == retcode);
#elif HAVE_ALIGNED_ALLOC
    ptr = aligned_alloc(alignment, size);
#elif HAVE_MEMALIGN
    ptr = memalign(alignment, size);
#else
    ptr = malloc(size);
#endif
    assert(NULL != ptr);
    return ptr;
}

inline int * restrict parasail_memalign_int(size_t alignment, size_t size)
{
    return (int * restrict) parasail_memalign(alignment, size*sizeof(int));
}

parasail_result_t* parasail_result_new()
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    result = (parasail_result_t*)malloc(sizeof(parasail_result_t));
    assert(result);

    result->score = 0;
    result->matches = 0;
    result->similarities = 0;
    result->length = 0;
    result->score_table = NULL;
    result->matches_table = NULL;
    result->similarities_table = NULL;
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
    result->score_table = (int * restrict)malloc(sizeof(int)*a*b);
    assert(result->score_table);

    return result;
}

parasail_result_t* parasail_result_new_table4(const int a, const int b)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    /* validate inputs */
    assert(a > 0);
    assert(b > 0);
    
    /* allocate struct to hold memory */
    result = parasail_result_new_table1(a, b);
    
    result->matches_table = (int * restrict)malloc(sizeof(int)*a*b);
    assert(result->matches_table);
    result->similarities_table = (int * restrict)malloc(sizeof(int)*a*b);
    assert(result->similarities_table);
    result->length_table = (int * restrict)malloc(sizeof(int)*a*b);
    assert(result->length_table);

    return result;
}

void parasail_result_free(parasail_result_t *result)
{
    /* validate inputs */
    assert(NULL != result);
    
    if (NULL != result->score_table) free(result->score_table);
    if (NULL != result->matches_table) free(result->matches_table);
    if (NULL != result->similarities_table) free(result->similarities_table);
    if (NULL != result->length_table) free(result->length_table);
    free(result);
}

