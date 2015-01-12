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

#include <emmintrin.h>

#include "parasail_internal_sse.h"

__m128i * parasail_memalign_m128i(size_t alignment, size_t size)
{
    return (__m128i *) parasail_memalign(alignment, size*sizeof(__m128i));
}

void parasail_memset_m128i(__m128i *b, __m128i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm_store_si128(&b[i], c);
    }
}

