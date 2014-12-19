/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#ifndef _PARASAIL_INTERNAL_AVX_H_
#define _PARASAIL_INTERNAL_AVX_H_

#include <stdint.h>

#include <immintrin.h>

typedef union __m256i_8 {
    __m256i m;
    int8_t v[32];
} __m256i_8_t;

typedef union __m256i_16 {
    __m256i m;
    int16_t v[16];
} __m256i_16_t;

typedef union __m256i_32 {
    __m256i m;
    int32_t v[8];
} __m256i_32_t;

__m256i * parasail_memalign_m256i(size_t alignment, size_t size);
void parasail_memset_m256i(__m256i *b, __m256i c, size_t len);

#endif /* _PARASAIL_INTERNAL_AVX_H_ */
