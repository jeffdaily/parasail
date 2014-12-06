/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#ifndef _PARASAIL_ALIGN_STRIPED_128_8_H_
#define _PARASAIL_ALIGN_STRIPED_128_8_H_

#include <limits.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX(a,b) ((a)>(b)?(a):(b))

#define NEG_INF_8 INT8_MIN

extern int nw_striped_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix);

extern int sg_striped_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix);

extern int sw_striped_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix);

extern int nw_stats_striped_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict matches, int * const restrict length);

extern int sg_stats_striped_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict matches, int * const restrict length);

extern int sw_stats_striped_128_8(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict matches, int * const restrict length);

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_ALIGN_STRIPED_128_16_H_ */
