/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#ifndef _PARASAIL_INTERNAL_KNC_H_
#define _PARASAIL_INTERNAL_KNC_H_

#include <stdint.h>

#include <immintrin.h>

typedef union __m512i_32 {
    __m512i m;
    int32_t v[16] __attribute__((aligned(64)));
} __m512i_32_t;

__m512i * parasail_memalign_m512i(size_t alignment, size_t size);
void parasail_memset_m512i(__m512i *b, __m512i c, size_t len);

#endif /* _PARASAIL_INTERNAL_KNC_H_ */
