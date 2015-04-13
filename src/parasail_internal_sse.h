/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#ifndef _PARASAIL_INTERNAL_SSE_H_
#define _PARASAIL_INTERNAL_SSE_H_

#include <stdint.h>

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

typedef union __m128i_64 {
    __m128i m;
    int64_t v[2];
} __m128i_64_t;

__m128i * parasail_memalign___m128i(size_t alignment, size_t size);
void parasail_memset___m128i(__m128i *b, __m128i c, size_t len);

#endif /* _PARASAIL_INTERNAL_SSE_H_ */
