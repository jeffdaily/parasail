/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <immintrin.h>

#include "parasail/memory.h"
#include "parasail/internal_avx.h"

__m256i * parasail_memalign___m256i(size_t alignment, size_t size)
{
    return (__m256i *) parasail_memalign(alignment, size*sizeof(__m256i));
}

void parasail_memset___m256i(__m256i *b, __m256i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm256_store_si256(&b[i], c);
    }
}

