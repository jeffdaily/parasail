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
#define ENAME sg_stats_table
#else
#define ENAME sg_stats
#endif

parasail_result_t* ENAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#if PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table4(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
    int * const restrict s1 = parasail_memalign_int(16, s1Len);
    int * const restrict s2 = parasail_memalign_int(16, s2Len);
    int * const restrict tbl_pr = parasail_memalign_int(16, s2Len+1);
    int * const restrict del_pr = parasail_memalign_int(16, s2Len+1);
    int * const restrict mch_pr = parasail_memalign_int(16, s2Len+1);
    int * const restrict len_pr = parasail_memalign_int(16, s2Len+1);
    int i = 0;
    int j = 0;
    int score = NEG_INF_32;
    int matches = NEG_INF_32;
    int length = NEG_INF_32;

    for (i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[(unsigned char)_s1[i]];
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[(unsigned char)_s2[j]];
    }

    /* upper left corner */
    tbl_pr[0] = 0;
    del_pr[0] = NEG_INF_32;
    mch_pr[0] = 0;
    len_pr[0] = 0;
    
    /* first row */
    for (j=1; j<=s2Len; ++j) {
        tbl_pr[j] = 0;
        del_pr[j] = NEG_INF_32;
        mch_pr[j] = 0;
        len_pr[j] = 0;
    }

    /* iter over first sequence */
    for (i=1; i<s1Len; ++i) {
        const int * const restrict matrow = matrix[s1[i-1]];
        /* init first column */
        int Nscore = tbl_pr[0];
        int Nmatches = mch_pr[0];
        int Nlength = len_pr[0];
        int Wscore = 0;
        int Wmatches = 0;
        int Wlength = 0;
        int ins_cr = NEG_INF_32;
        tbl_pr[0] = Wscore;
        mch_pr[0] = Wmatches;
        len_pr[0] = Wlength;
        for (j=1; j<=s2Len; ++j) {
            int NWscore = Nscore;
            int NWmatches = Nmatches;
            int NWlength = Nlength;
            Nscore = tbl_pr[j];
            Nmatches = mch_pr[j];
            Nlength = len_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            tbl_pr[j] = NWscore + matrow[s2[j-1]];
            if ((tbl_pr[j] >= del_pr[j]) && (tbl_pr[j] >= ins_cr)) {
                Wscore = tbl_pr[j];
                Wmatches  = NWmatches + (s1[i-1] == s2[j-1]);
                Wlength  = NWlength + 1;
            } else if (del_pr[j] >= ins_cr) {
                Wscore = del_pr[j];
                Wmatches  = Nmatches;
                Wlength  = Nlength + 1;
            } else {
                Wscore = ins_cr;
                Wmatches  = Wmatches;
                Wlength  = Wlength + 1;
            }
            tbl_pr[j] = Wscore;
            mch_pr[j] = Wmatches;
            len_pr[j] = Wlength;
#ifdef PARASAIL_TABLE
            result->score_table[(i-1)*s2Len + (j-1)] = Wscore;
            result->matches_table[(i-1)*s2Len + (j-1)] = Wmatches;
            result->length_table[(i-1)*s2Len + (j-1)] = Wlength;
#endif
        }
        if (Wscore > score) {
            score = Wscore;
            matches = Wmatches;
            length = Wlength;
        }
    }
    {
        /* i == s1Len */
        const int * const restrict matrow = matrix[s1[i-1]];
        /* init first column */
        int Nscore = tbl_pr[0];
        int Nmatches = mch_pr[0];
        int Nlength = len_pr[0];
        int Wscore = 0;
        int Wmatches = 0;
        int Wlength = 0;
        int ins_cr = NEG_INF_32;
        tbl_pr[0] = Wscore;
        mch_pr[0] = Wmatches;
        len_pr[0] = Wlength;
        for (j=1; j<=s2Len; ++j) {
            int NWscore = Nscore;
            int NWmatches = Nmatches;
            int NWlength = Nlength;
            Nscore = tbl_pr[j];
            Nmatches = mch_pr[j];
            Nlength = len_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            tbl_pr[j] = NWscore + matrow[s2[j-1]];
            if ((tbl_pr[j] >= del_pr[j]) && (tbl_pr[j] >= ins_cr)) {
                Wscore = tbl_pr[j];
                Wmatches  = NWmatches + (s1[i-1] == s2[j-1]);
                Wlength  = NWlength + 1;
            } else if (del_pr[j] >= ins_cr) {
                Wscore = del_pr[j];
                Wmatches  = Nmatches;
                Wlength  = Nlength + 1;
            } else {
                Wscore = ins_cr;
                Wmatches  = Wmatches;
                Wlength  = Wlength + 1;
            }
            tbl_pr[j] = Wscore;
            mch_pr[j] = Wmatches;
            len_pr[j] = Wlength;
#ifdef PARASAIL_TABLE
            result->score_table[(i-1)*s2Len + (j-1)] = Wscore;
            result->matches_table[(i-1)*s2Len + (j-1)] = Wmatches;
            result->length_table[(i-1)*s2Len + (j-1)] = Wlength;
#endif
            if (Wscore > score) {
                score = Wscore;
                matches = Wmatches;
                length = Wlength;
            }
        }
    }

    result->score = score;
    result->matches = matches;
    result->length = length;

    free(len_pr);
    free(mch_pr);
    free(del_pr);
    free(tbl_pr);
    free(s2);
    free(s1);

    return result;
}

