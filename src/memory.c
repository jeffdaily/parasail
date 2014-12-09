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
#include <stddef.h>
#include <stdint.h>

#include "parasail.h"
#include "parasail_internal.h"

static inline void* my_memalign(size_t alignment, size_t size)
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
    assert(ptr);
    return ptr;
}

inline void parasail_workspace_allocate_s(
        parasail_workspace_t *workspace,
        const int length)
{
    /* declare all variables */
    int32_t padding = 0;

    /* validate inputs */
    assert(length > 0);
    
    /* The wozniak anti-diagonal algorithm needs to pad the input
     * sequences with the number of elements in the vector. We choose
     * the largest possible value here. 512-bit vector unit with 32-bit
     * elements is 16. We need twice that (pad both ends of sequence). */
    padding = 16*2;

    workspace->s1 = (int * restrict)malloc(sizeof(int)*(length+padding));
    assert(workspace->s1);
    workspace->s2 = (int * restrict)malloc(sizeof(int)*(length+padding));
    assert(workspace->s2);
}

inline void parasail_workspace_allocate_serial2(
        parasail_workspace_t *workspace,
        const int length)
{
    workspace->a1 = (int * restrict)malloc(sizeof(int)*length);
    assert(workspace->a1);
    workspace->a2 = (int * restrict)malloc(sizeof(int)*length);
    assert(workspace->a2);
}

inline void parasail_workspace_allocate_serial4(
        parasail_workspace_t *workspace,
        const int length)
{
    parasail_workspace_allocate_serial2(workspace, length);
    workspace->a3 = (int * restrict)malloc(sizeof(int)*length);
    assert(workspace->a1);
    workspace->a4 = (int * restrict)malloc(sizeof(int)*length);
    assert(workspace->a2);
}

inline void parasail_workspace_allocate_serial9(
        parasail_workspace_t *workspace,
        const int length)
{
    parasail_workspace_allocate_serial4(workspace, length);
    workspace->a5 = (int * restrict)malloc(sizeof(int)*length);
    assert(workspace->a5);
    workspace->a6 = (int * restrict)malloc(sizeof(int)*length);
    assert(workspace->a6);
    workspace->a7 = (int * restrict)malloc(sizeof(int)*length);
    assert(workspace->a7);
    workspace->a8 = (int * restrict)malloc(sizeof(int)*length);
    assert(workspace->a8);
    workspace->a9 = (int * restrict)malloc(sizeof(int)*length);
    assert(workspace->a9);
}

inline parasail_workspace_t* parasail_workspace_new() {
    return (parasail_workspace_t*)malloc(sizeof(parasail_workspace_t));
}

parasail_workspace_t* parasail_workspace_allocate(const int length)
{
    /* declare all variables */
    parasail_workspace_t *workspace = NULL;
    int32_t vector_count = 0;
    int32_t vector_size = 0;

    /* validate inputs */
    assert(length > 0);
    
    /* allocate struct to hold memory */
    workspace = parasail_workspace_new();
    assert(workspace);

    parasail_workspace_allocate_s(workspace, length);
    parasail_workspace_allocate_serial9(workspace, length);

    workspace->boundary = (int * restrict)malloc(sizeof(int)*length);
    assert(workspace->boundary);

#if HAVE_SSE2 || HAVE_SSE41
    vector_count = (length + 7) / 8;
    vector_size = vector_count * sizeof(__m128i);
    workspace->sse_1 = (__m128i*)my_memalign(16, vector_size);
    workspace->sse_2 = (__m128i*)my_memalign(16, vector_size);
    workspace->sse_3 = (__m128i*)my_memalign(16, vector_size);
    workspace->sse_4 = (__m128i*)my_memalign(16, vector_size);
    workspace->sse_5 = (__m128i*)my_memalign(16, vector_size);
    workspace->sse_6 = (__m128i*)my_memalign(16, vector_size);
    workspace->sse_7 = (__m128i*)my_memalign(16, vector_size);
    workspace->sse_8 = (__m128i*)my_memalign(16, vector_size);
    workspace->sse_9 = (__m128i*)my_memalign(16, vector_size);
    workspace->sse_10 = (__m128i*)my_memalign(16, vector_size);
    workspace->sse_11 = (__m128i*)my_memalign(16, vector_size);
    workspace->sse_12 = (__m128i*)my_memalign(16, vector_size);
#endif

#if HAVE_AVX2
    vector_count = (length + 7) / 8;
    vector_size = vector_count * sizeof(__m256i);
    workspace->avx_1 = (__m256i*)my_memalign(32, vector_size);
    workspace->avx_2 = (__m256i*)my_memalign(32, vector_size);
    workspace->avx_3 = (__m256i*)my_memalign(32, vector_size);
    workspace->avx_4 = (__m256i*)my_memalign(32, vector_size);
    workspace->avx_5 = (__m256i*)my_memalign(32, vector_size);
    workspace->avx_6 = (__m256i*)my_memalign(32, vector_size);
    workspace->avx_7 = (__m256i*)my_memalign(32, vector_size);
    workspace->avx_8 = (__m256i*)my_memalign(32, vector_size);
    workspace->avx_9 = (__m256i*)my_memalign(32, vector_size);
    workspace->avx_10 = (__m256i*)my_memalign(32, vector_size);
    workspace->avx_11 = (__m256i*)my_memalign(32, vector_size);
    workspace->avx_12 = (__m256i*)my_memalign(32, vector_size);
#endif

#if HAVE_KNC
    vector_count = (length + 15) / 16;
    vector_size = vector_count * sizeof(__m512i);
    workspace->knc_1 = (__m512i*)my_memalign(64, vector_size);
    workspace->knc_2 = (__m512i*)my_memalign(64, vector_size);
    workspace->knc_3 = (__m512i*)my_memalign(64, vector_size);
    workspace->knc_4 = (__m512i*)my_memalign(64, vector_size);
    workspace->knc_5 = (__m512i*)my_memalign(64, vector_size);
    workspace->knc_6 = (__m512i*)my_memalign(64, vector_size);
    workspace->knc_7 = (__m512i*)my_memalign(64, vector_size);
    workspace->knc_8 = (__m512i*)my_memalign(64, vector_size);
    workspace->knc_9 = (__m512i*)my_memalign(64, vector_size);
    workspace->knc_10 = (__m512i*)my_memalign(64, vector_size);
    workspace->knc_11 = (__m512i*)my_memalign(64, vector_size);
    workspace->knc_12 = (__m512i*)my_memalign(64, vector_size);
#endif

    return workspace;
}

void parasail_workspace_free(parasail_workspace_t *workspace)
{
    /* validate inputs */
    assert(workspace);
    
    free(workspace->s1);
    free(workspace->s2);

    free(workspace->a1);
    free(workspace->a2);
    free(workspace->a3);
    free(workspace->a4);
    free(workspace->a5);
    free(workspace->a6);
    free(workspace->a7);
    free(workspace->a8);
    free(workspace->a9);

    free(workspace->boundary);

#if HAVE_SSE2 || HAVE_SSE41
    free(workspace->sse_1);
    free(workspace->sse_2);
    free(workspace->sse_3);
    free(workspace->sse_4);
    free(workspace->sse_5);
    free(workspace->sse_6);
    free(workspace->sse_7);
    free(workspace->sse_8);
    free(workspace->sse_9);
    free(workspace->sse_10);
    free(workspace->sse_11);
    free(workspace->sse_12);
#endif

#if HAVE_AVX2
    free(workspace->avx_1);
    free(workspace->avx_2);
    free(workspace->avx_3);
    free(workspace->avx_4);
    free(workspace->avx_5);
    free(workspace->avx_6);
    free(workspace->avx_7);
    free(workspace->avx_8);
    free(workspace->avx_9);
    free(workspace->avx_10);
    free(workspace->avx_11);
    free(workspace->avx_12);
#endif

#if HAVE_KNC
    free(workspace->knc_1);
    free(workspace->knc_2);
    free(workspace->knc_3);
    free(workspace->knc_4);
    free(workspace->knc_5);
    free(workspace->knc_6);
    free(workspace->knc_7);
    free(workspace->knc_8);
    free(workspace->knc_9);
    free(workspace->knc_10);
    free(workspace->knc_11);
    free(workspace->knc_12);
#endif

    free(workspace);
}

parasail_result_t* parasail_result_allocate(const int length)
{
    /* declare all variables */
    parasail_result_t *retval = NULL;

    /* validate inputs */
    assert(length > 0);
    
    /* allocate struct to hold memory */
    retval = (parasail_result_t*)malloc(sizeof(parasail_result_t));
    assert(retval);

    retval->score = 0;
    retval->matches = 0;
    retval->similarities = 0;
    retval->length = 0;

    retval->score_table = (int * restrict)malloc(sizeof(int)*length);
    assert(retval->score_table);
    retval->matches_table = (int * restrict)malloc(sizeof(int)*length);
    assert(retval->matches_table);
    retval->similarities_table = (int * restrict)malloc(sizeof(int)*length);
    assert(retval->similarities_table);
    retval->length_table = (int * restrict)malloc(sizeof(int)*length);
    assert(retval->length_table);

    return retval;
}

void parasail_result_free(parasail_result_t *result)
{
    /* validate inputs */
    assert(NULL != result);
    
    free(result->score_table);
    free(result->matches_table);
    free(result->similarities_table);
    free(result->length_table);
    free(result);
}

