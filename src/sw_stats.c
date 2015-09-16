/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include "parasail.h"
#include "parasail/memory.h"

#define NEG_INF_32 (INT32_MIN/2)
#define MAX(a,b) ((a)>(b)?(a):(b))

#ifdef PARASAIL_TABLE
#define ENAME parasail_sw_stats_table
#else
#ifdef PARASAIL_ROWCOL
#define ENAME parasail_sw_stats_rowcol
#else
#define ENAME parasail_sw_stats
#endif
#endif

parasail_result_t* ENAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif
    int * const restrict s1 = parasail_memalign_int(16, s1Len);
    int * const restrict s2 = parasail_memalign_int(16, s2Len);
    int * const restrict tbl_pr = parasail_memalign_int(16, s2Len+1);
    int * const restrict del_pr = parasail_memalign_int(16, s2Len+1);
    int * const restrict mch_pr = parasail_memalign_int(16, s2Len+1);
    int * const restrict sim_pr = parasail_memalign_int(16, s2Len+1);
    int * const restrict len_pr = parasail_memalign_int(16, s2Len+1);
    int i = 0;
    int j = 0;
    int score = NEG_INF_32;
    int matches = NEG_INF_32;
    int similar = NEG_INF_32;
    int length = NEG_INF_32;
    int end_query = s1Len;
    int end_ref = s2Len;

    for (i=0; i<s1Len; ++i) {
        s1[i] = matrix->mapper[(unsigned char)_s1[i]];
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = matrix->mapper[(unsigned char)_s2[j]];
    }

    /* upper left corner */
    tbl_pr[0] = 0;
    del_pr[0] = NEG_INF_32;
    mch_pr[0] = 0;
    sim_pr[0] = 0;
    len_pr[0] = 0;
    
    /* first row */
    for (j=1; j<=s2Len; ++j) {
        tbl_pr[j] = 0;
        del_pr[j] = NEG_INF_32;
        mch_pr[j] = 0;
        sim_pr[j] = 0;
        len_pr[j] = 0;
    }

    /* iter over first sequence */
    for (i=1; i<=s1Len; ++i) {
        const int * const restrict matrow = &matrix->matrix[matrix->size*s1[i-1]];
        /* init first column */
        int Nscore = tbl_pr[0];
        int Nmatches = mch_pr[0];
        int Nsimilar = sim_pr[0];
        int Nlength = len_pr[0];
        int Wscore = 0;
        int Wmatches = 0;
        int Wsimilar = 0;
        int Wlength = 0;
        int ins_cr = NEG_INF_32;
        tbl_pr[0] = Wscore;
        mch_pr[0] = Wmatches;
        sim_pr[0] = Wsimilar;
        len_pr[0] = Wlength;
        for (j=1; j<=s2Len; ++j) {
            int NWscore = Nscore;
            int NWmatches = Nmatches;
            int NWsimilar = Nsimilar;
            int NWlength = Nlength;
            Nscore = tbl_pr[j];
            Nmatches = mch_pr[j];
            Nsimilar = sim_pr[j];
            Nlength = len_pr[j];
            del_pr[j] = MAX(Nscore - open, del_pr[j] - gap);
            ins_cr    = MAX(Wscore - open, ins_cr    - gap);
            tbl_pr[j] = MAX(NWscore + matrow[s2[j-1]], 0);
            if ((tbl_pr[j] >= del_pr[j]) && (tbl_pr[j] >= ins_cr)) {
                Wscore = tbl_pr[j];
                Wmatches  = NWmatches + (s1[i-1] == s2[j-1]);
                Wsimilar  = NWsimilar + (matrow[s2[j-1]] > 0);
                Wlength  = NWlength + 1;
            } else if (del_pr[j] >= ins_cr) {
                Wscore = del_pr[j];
                Wmatches  = Nmatches;
                Wsimilar  = Nsimilar;
                Wlength  = Nlength + 1;
            } else {
                Wscore = ins_cr;
                /*Wmatches  = Wmatches;*/
                /*Wsimilar  = Wsimilar;*/
                Wlength  = Wlength + 1;
            }
            if (Wscore <= 0) {
                Wscore = 0;
                Wmatches = 0;
                Wsimilar = 0;
                Wlength = 0;
            }
            tbl_pr[j] = Wscore;
            mch_pr[j] = Wmatches;
            sim_pr[j] = Wsimilar;
            len_pr[j] = Wlength;
#ifdef PARASAIL_TABLE
            result->score_table[(i-1)*s2Len + (j-1)] = Wscore;
            result->matches_table[(i-1)*s2Len + (j-1)] = Wmatches;
            result->similar_table[(i-1)*s2Len + (j-1)] = Wsimilar;
            result->length_table[(i-1)*s2Len + (j-1)] = Wlength;
#endif
            if (Wscore > score) {
                score = Wscore;
                matches = Wmatches;
                similar = Wsimilar;
                length = Wlength;
                end_query = i-1;
                end_ref = j-1;
            }
            else if (score == Wscore && j-1 < end_ref) {
                matches = Wmatches;
                similar = Wsimilar;
                length = Wlength;
                end_query = i-1;
                end_ref = j-1;
            }
        }
#ifdef PARASAIL_ROWCOL
        result->score_col[i-1] = Wscore;
        result->matches_col[i-1] = Wmatches;
        result->similar_col[i-1] = Wsimilar;
        result->length_col[i-1] = Wlength;
#endif
    }
#ifdef PARASAIL_ROWCOL
    for (j=1; j<=s2Len; ++j) {
        result->score_row[j-1] = tbl_pr[j];
        result->matches_row[j-1] = mch_pr[j];
        result->similar_row[j-1] = sim_pr[j];
        result->length_row[j-1] = len_pr[j];
    }
#endif

    result->score = score;
    result->matches = matches;
    result->similar = similar;
    result->length = length;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(len_pr);
    parasail_free(sim_pr);
    parasail_free(mch_pr);
    parasail_free(del_pr);
    parasail_free(tbl_pr);
    parasail_free(s2);
    parasail_free(s1);

    return result;
}

