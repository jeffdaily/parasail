#ifndef _PGRAPH_ALIGN_SCAN_512_32_H_
#define _PGRAPH_ALIGN_SCAN_512_32_H_

#include <limits.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX(a,b) ((a)>(b)?(a):(b))

#define NEG_INF_32 (INT32_MIN/(int32_t)2)

extern int32_t nw_scan_512_32(
        const char * const restrict s1, const int32_t s1Len,
        const char * const restrict s2, const int32_t s2Len,
        const int32_t open, const int32_t gap,
        const int8_t * const restrict matrix);

extern int32_t sg_scan_512_32(
        const char * const restrict s1, const int32_t s1Len,
        const char * const restrict s2, const int32_t s2Len,
        const int32_t open, const int32_t gap,
        const int8_t * const restrict matrix);

extern int32_t sw_scan_512_32(
        const char * const restrict s1, const int32_t s1Len,
        const char * const restrict s2, const int32_t s2Len,
        const int32_t open, const int32_t gap,
        const int8_t * const restrict matrix);

extern int32_t nw_stats_scan_512_32(
        const char * const restrict s1, const int32_t s1Len,
        const char * const restrict s2, const int32_t s2Len,
        const int32_t open, const int32_t gap,
        const int8_t * const restrict matrix,
        int32_t * const restrict matches, int32_t * const restrict length);

extern int32_t sg_stats_scan_512_32(
        const char * const restrict s1, const int32_t s1Len,
        const char * const restrict s2, const int32_t s2Len,
        const int32_t open, const int32_t gap,
        const int8_t * const restrict matrix,
        int32_t * const restrict matches, int32_t * const restrict length);

extern int32_t sw_stats_scan_512_32(
        const char * const restrict s1, const int32_t s1Len,
        const char * const restrict s2, const int32_t s2Len,
        const int32_t open, const int32_t gap,
        const int8_t * const restrict matrix,
        int32_t * const restrict matches, int32_t * const restrict length);

#ifdef __cplusplus
}
#endif

#endif /* _PGRAPH_ALIGN_SCAN_512_32_H_ */
