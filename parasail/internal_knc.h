/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_INTERNAL_KNC_H_
#define _PARASAIL_INTERNAL_KNC_H_

#include <stdint.h>

#include <immintrin.h>

#include "parasail.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef union __m512i_32 {
    __m512i m;
    int32_t v[16] __attribute__((aligned(64)));
} __m512i_32_t;

extern PARASAIL_LOCAL
__m512i * parasail_memalign___m512i(size_t alignment, size_t size);

extern PARASAIL_LOCAL
void parasail_memset___m512i(__m512i *b, __m512i c, size_t len);

extern PARASAIL_LOCAL
void parasail_profile_free___m512i(void *profile);

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_INTERNAL_KNC_H_ */
