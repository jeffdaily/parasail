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

void * parasail_memalign(size_t alignment, size_t size);
int * parasail_memalign_int(size_t alignment, size_t size);
int8_t * parasail_memalign_int8_t(size_t alignment, size_t size);
int16_t * parasail_memalign_int16_t(size_t alignment, size_t size);
int32_t * parasail_memalign_int32_t(size_t alignment, size_t size);

void parasail_memset(void *b, int c, size_t len);
void parasail_memset_int(int *b, int c, size_t len);
void parasail_memset_int8_t(int8_t *b, int8_t c, size_t len);
void parasail_memset_int16_t(int16_t *b, int16_t c, size_t len);
void parasail_memset_int32_t(int32_t *b, int32_t c, size_t len);

parasail_result_t* parasail_result_new();
parasail_result_t* parasail_result_new_table1(const int a, const int b);
parasail_result_t* parasail_result_new_table4(const int a, const int b);

#endif /* _PARASAIL_INTERNAL_H_ */
