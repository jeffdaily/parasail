/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <emmintrin.h>

#include "parasail/memory.h"
#include "parasail/internal_sse.h"

__m128i * parasail_memalign___m128i(size_t alignment, size_t size)
{
    return (__m128i *) parasail_memalign(alignment, size*sizeof(__m128i));
}

void parasail_memset___m128i(__m128i *b, __m128i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm_store_si128(&b[i], c);
    }
}

