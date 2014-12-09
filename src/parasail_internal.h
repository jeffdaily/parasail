/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#ifndef _PARASAIL_INTERNAL_H_
#define _PARASAIL_INTERNAL_H_

#include "config.h"

#if HAVE_SSE2
#include <emmintrin.h>
#endif

#if HAVE_SSE41
#include <smmintrin.h>
#endif

#if HAVE_AVX2
#include <immintrin.h>
#endif

#if HAVE_KNC
#include <zmmintrin.h>
#endif

#include "parasail.h"

void parasail_workspace_allocate_s(parasail_workspace_t *workspace, const int length);
void parasail_workspace_allocate_serial2(parasail_workspace_t *workspace, const int length);
void parasail_workspace_allocate_serial4(parasail_workspace_t *workspace, const int length);
void parasail_workspace_allocate_serial9(parasail_workspace_t *workspace, const int length);

void parasail_workspace_free_s(parasail_workspace_t *workspace, const int length);
void parasail_workspace_free_serial2(parasail_workspace_t *workspace, const int length);
void parasail_workspace_free_serial4(parasail_workspace_t *workspace, const int length);
void parasail_workspace_free_serial9(parasail_workspace_t *workspace, const int length);

struct parasail_workspace {
    /* Used by some methods to convert character arrays into integer
     * arrays to speed up later processing. */
    int * restrict s1;
    int * restrict s2;

    /* Serial methods use at most nine arrays of integers. */
    int * restrict a1;
    int * restrict a2;
    int * restrict a3;
    int * restrict a4;
    int * restrict a5;
    int * restrict a6;
    int * restrict a7;
    int * restrict a8;
    int * restrict a9;

    /* Used by some nw methods to set table boundary conditions. */
    int * restrict boundary;

    /* Vector methods use at most 12 arrays of integer vectors. */
#if HAVE_SSE2 || HAVE_SSE41
    __m128i * restrict sse_1;
    __m128i * restrict sse_2;
    __m128i * restrict sse_3;
    __m128i * restrict sse_4;
    __m128i * restrict sse_5;
    __m128i * restrict sse_6;
    __m128i * restrict sse_7;
    __m128i * restrict sse_8;
    __m128i * restrict sse_9;
    __m128i * restrict sse_10;
    __m128i * restrict sse_11;
    __m128i * restrict sse_12;
#endif
#if HAVE_AVX2
    __m256i * restrict avx_1;
    __m256i * restrict avx_2;
    __m256i * restrict avx_3;
    __m256i * restrict avx_4;
    __m256i * restrict avx_5;
    __m256i * restrict avx_6;
    __m256i * restrict avx_7;
    __m256i * restrict avx_8;
    __m256i * restrict avx_9;
    __m256i * restrict avx_10;
    __m256i * restrict avx_11;
    __m256i * restrict avx_12;
#endif
#if HAVE_KNC
    __m512i * restrict knc_1;
    __m512i * restrict knc_2;
    __m512i * restrict knc_3;
    __m512i * restrict knc_4;
    __m512i * restrict knc_5;
    __m512i * restrict knc_6;
    __m512i * restrict knc_7;
    __m512i * restrict knc_8;
    __m512i * restrict knc_9;
    __m512i * restrict knc_10;
    __m512i * restrict knc_11;
    __m512i * restrict knc_12;
#endif
};

#endif /* _PARASAIL_INTERNAL_H_ */
