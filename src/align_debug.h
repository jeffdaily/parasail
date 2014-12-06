#ifndef _PARASAIL_ALIGN_DEBUG_H_
#define _PARASAIL_ALIGN_DEBUG_H_

#include <limits.h>
#include <stdint.h>

#define MAX(a,b) ((a)>(b)?(a):(b))

#define NEG_INF_32 (INT32_MIN/(int32_t)2)

extern int nw_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict score_table);

extern int nw_scan_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict score_table);

extern int nw_scan_row_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict score_table);

extern int sg_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict score_table);

extern int sg_scan_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict score_table);

extern int sw_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict score_table);

extern int sw_scan_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict score_table);

extern int nw_stats_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict matches, int * const restrict length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr,
        int * const restrict score_table,
        int * const restrict match_table,
        int * const restrict length_table);

extern int nw_stats_scan_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict matches, int * const restrict length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr,
        int * const restrict score_table,
        int * const restrict match_table,
        int * const restrict length_table);

extern int sg_stats_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict matches, int * const restrict length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr,
        int * const restrict score_table,
        int * const restrict match_table,
        int * const restrict length_table);

extern int sg_stats_scan_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict matches, int * const restrict length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr,
        int * const restrict score_table,
        int * const restrict match_table,
        int * const restrict length_table);

extern int sw_stats_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict matches, int * const restrict length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr,
        int * const restrict score_table,
        int * const restrict match_table,
        int * const restrict length_table);

extern int sw_stats_scan_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict matches, int * const restrict length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr,
        int * const restrict score_table,
        int * const restrict match_table,
        int * const restrict length_table);

#endif /* _PARASAIL_ALIGN_DEBUG_H_ */
