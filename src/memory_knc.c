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

#include <immintrin.h>

#include "parasail_internal.h"
#include "parasail_internal_knc.h"

__m512i * parasail_memalign_m512i(size_t alignment, size_t size)
{
    return (__m512i *) parasail_memalign(alignment, size*sizeof(__m512i));
}

void parasail_memset_m512i(__m512i *b, __m512i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm512_store_epi32(&b[i], c);
    }
}

