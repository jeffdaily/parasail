#ifndef _PGRAPH_ALIGN_H_
#define _PGRAPH_ALIGN_H_

#include <limits.h>
#include <stdint.h>

#define MAX(a,b) ((a)>(b)?(a):(b))

static const int16_t NEG_INF = SHRT_MIN / 2;

extern int nw_debug(
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

extern int sw_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict score_table);

extern int nw_woz_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict score_table);

extern int sg_woz_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict score_table);

extern int sw_woz_debug(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict score_table);

extern int nw_striped_debug(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict score_table);

extern int sg_striped_debug(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict score_table);

extern int sw_striped_debug(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
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

extern int nw_stats_woz_debug(
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

extern int sg_stats_woz_debug(
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

extern int sw_stats_woz_debug(
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

extern int nw_stats_striped_debug(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict matches, int * const restrict length,
        int * const restrict score_table,
        int * const restrict match_table,
        int * const restrict length_table);

extern int sg_stats_striped_debug(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict matches, int * const restrict length,
        int * const restrict score_table,
        int * const restrict match_table,
        int * const restrict length_table);

extern int sw_stats_striped_debug(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int8_t * const restrict matrix,
        int * const restrict matches, int * const restrict length,
        int * const restrict score_table,
        int * const restrict match_table,
        int * const restrict length_table);

#endif /* _PGRAPH_ALIGN_H_ */
