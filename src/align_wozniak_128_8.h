#ifndef _PARASAIL_ALIGN_WOZNIAK_128_8_H_
#define _PARASAIL_ALIGN_WOZNIAK_128_8_H_

#include <limits.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX(a,b) ((a)>(b)?(a):(b))

/*#define NEG_INF_8 (INT8_MIN/(int8_t)2)*/
#define NEG_INF_8 INT8_MIN

extern int nw_wozniak_128_8(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr);

extern int sg_wozniak_128_8(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr);

extern int sw_wozniak_128_8(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr);

extern int nw_stats_wozniak_128_8(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict matches, int * const restrict length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr);

extern int sg_stats_wozniak_128_8(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict matches, int * const restrict length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr);

extern int sw_stats_wozniak_128_8(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict matches, int * const restrict length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr);

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_ALIGN_WOZNIAK_128_8_H_ */
