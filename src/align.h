#ifndef _PARASAIL_ALIGN_H_
#define _PARASAIL_ALIGN_H_

#include <limits.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX(a,b) ((a)>(b)?(a):(b))

#define NEG_INF_32 (INT32_MIN/(int32_t)2)

extern int nw(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr);

extern int nw_scan(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr);

extern int nw_scan_row(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr);

extern int sg(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr);

extern int sg_scan(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr);

extern int sw(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr);

extern int sw_scan(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr);

extern int nw_stats(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict matches, int * const restrict length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr);

extern int nw_stats_scan(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict matches, int * const restrict length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr);

extern int sg_stats(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict matches, int * const restrict length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr);

extern int sg_stats_scan(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict matches, int * const restrict length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr);

extern int sw_stats(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict matches, int * const restrict length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr);

extern int sw_stats_scan(
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

#endif /* _PARASAIL_ALIGN_H_ */
