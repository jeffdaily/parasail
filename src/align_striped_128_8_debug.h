/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#ifndef _PARASAIL_ALIGN_STRIPED_128_8_DEBUG_H_
#define _PARASAIL_ALIGN_STRIPED_128_8_DEBUG_H_

#include <limits.h>
#include <stdint.h>

#define MAX(a,b) ((a)>(b)?(a):(b))

#define NEG_INF_8 INT8_MIN

extern int nw_striped_128_8_debug(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict score_table);

extern int sg_striped_128_8_debug(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict score_table);

extern int sw_striped_128_8_debug(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict score_table);

extern int nw_stats_striped_128_8_debug(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict matches, int * const restrict length,
        int * const restrict score_table,
        int * const restrict match_table,
        int * const restrict length_table);

extern int sg_stats_striped_128_8_debug(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict matches, int * const restrict length,
        int * const restrict score_table,
        int * const restrict match_table,
        int * const restrict length_table);

extern int sw_stats_striped_128_8_debug(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict matches, int * const restrict length,
        int * const restrict score_table,
        int * const restrict match_table,
        int * const restrict length_table);

#endif /* _PARASAIL_ALIGN_STRIPED_128_8_DEBUG_H_ */
