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

#include <stdint.h>

#include "parasail.h"

void* parasail_memalign(size_t alignment, size_t size);
int * restrict parasail_memalign_int(size_t alignment, size_t size);
int8_t * restrict parasail_memalign_int8_t(size_t alignment, size_t size);
int16_t * restrict parasail_memalign_int16_t(size_t alignment, size_t size);
int32_t * restrict parasail_memalign_int32_t(size_t alignment, size_t size);
parasail_result_t* parasail_result_new();
parasail_result_t* parasail_result_new_table1(const int a, const int b);
parasail_result_t* parasail_result_new_table4(const int a, const int b);

#if HAVE_SSE2 || HAVE_SSE41

#include <emmintrin.h>

typedef union __m128i_8 {
    __m128i m;
    int8_t v[16];
} __m128i_8_t;

typedef union __m128i_16 {
    __m128i m;
    int16_t v[8];
} __m128i_16_t;

typedef union __m128i_32 {
    __m128i m;
    int32_t v[4];
} __m128i_32_t;

__m128i * restrict parasail_memalign_m128i(size_t alignment, size_t size);

#endif /* HAVE_SSE2 || HAVE_SSE41 */

#endif /* _PARASAIL_INTERNAL_H_ */
