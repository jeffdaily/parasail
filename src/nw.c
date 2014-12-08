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

#include <stdint.h>
#include <stdlib.h>

#include "parasail.h"
#include "parasail_internal.h"
#include "blosum/blosum_map.h"

#define NEG_INF_32 (INT32_MIN/2)
#define MAX(a,b) ((a)>(b)?(a):(b))

#ifdef PARASAIL_TABLE
#define ENAME nw_table_ext
#else
#define ENAME nw_ext
#endif

void ENAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24],
        parasail_result_t *result, parasail_workspace_t *workspace)
{
    int * const restrict s1 = workspace->s1;
    int * const restrict s2 = workspace->s2;
    int * const restrict tbl_pr = workspace->a1;
    int * const restrict del_pr = workspace->a2;
    int i = 0;
    int j = 0;

    for (i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[(unsigned char)_s1[i]];
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[(unsigned char)_s2[j]];
    }

    /* upper left corner */
    tbl_pr[0] = 0;
    del_pr[0] = NEG_INF_32;
    
    /* first row */
    for (j=1; j<=s2Len; ++j) {
        tbl_pr[j] = -open -(j-1)*gap;
        del_pr[j] = NEG_INF_32;
    }

    /* iter over first sequence */
    for (i=1; i<=s1Len; ++i) {
        const int * const restrict matrow = matrix[s1[i-1]];
        /* init first column */
        int Nscore = tbl_pr[0];
        int Wscore = -open - (i-1)*gap;
        int ins_cr = NEG_INF_32;
        tbl_pr[0] = Wscore;
        for (j=1; j<=s2Len; ++j) {
            int NWscore = Nscore;
            Nscore = tbl_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            tbl_pr[j] = NWscore + matrow[s2[j-1]];
            Wscore = tbl_pr[j] = MAX(tbl_pr[j],MAX(ins_cr,del_pr[j]));
#ifdef PARASAIL_TABLE
            workspace->score_table[(i-1)*s2Len + (j-1)] = Wscore;
#endif
        }
    }

    result->score = tbl_pr[s2Len];
}

#ifdef PARASAIL_TABLE
#define NAME nw_table
#else
#define NAME nw
#endif

void NAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24],
        parasail_result_t *result)
{
    const int longest = MAX(s1Len, s2Len);
    parasail_workspace_t *workspace = parasail_workspace_allocate(longest);
    ENAME(_s1, s1Len, _s2, s2Len, open, gap, matrix, result, workspace);
    parasail_workspace_free(workspace);
}

