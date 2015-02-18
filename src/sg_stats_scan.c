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
#define ENAME sg_stats_table_scan
#else
#define ENAME sg_stats_scan
#endif

parasail_result_t* ENAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new();
#endif
    int * const restrict s1 = parasail_memalign_int(16, s1Len);
    int * const restrict s2 = parasail_memalign_int(16, s2Len);
    int * const restrict HB = parasail_memalign_int(16, s1Len+1);
    int * const restrict H  = HB+1;
    int * const restrict E  = parasail_memalign_int(16, s1Len);
    int * const restrict HtB= parasail_memalign_int(16, s1Len+1);
    int * const restrict Ht = HtB+1;
    int * const restrict FtB= parasail_memalign_int(16, s1Len+1);
    int * const restrict Ft = FtB+1;
    int * const restrict MB = parasail_memalign_int(16, s1Len+1);
    int * const restrict M  = MB+1;
    int * const restrict SB = parasail_memalign_int(16, s1Len+1);
    int * const restrict S  = SB+1;
    int * const restrict LB = parasail_memalign_int(16, s1Len+1);
    int * const restrict L  = LB+1;
    int * const restrict Mt = parasail_memalign_int(16, s1Len);
    int * const restrict St = parasail_memalign_int(16, s1Len);
    int * const restrict Lt = parasail_memalign_int(16, s1Len);
    int * const restrict Ex = parasail_memalign_int(16, s1Len);
    int i = 0;
    int j = 0;
    int score = NEG_INF_32;
    int matches = 0;
    int similar = 0;
    int length = 0;

    for (i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[(unsigned char)_s1[i]];
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[(unsigned char)_s2[j]];
    }

    /* initialize H */
    H[-1] = 0;
    Ht[-1] = 0;
    for (i=0; i<s1Len; ++i) {
        H[i] = 0;
    }

    /* initialize M */
    M[-1] = 0;
    for (i=0; i<s1Len; ++i) {
        M[i] = 0;
    }

    /* initialize S */
    S[-1] = 0;
    for (i=0; i<s1Len; ++i) {
        S[i] = 0;
    }

    /* initialize L */
    L[-1] = 0;
    for (i=0; i<s1Len; ++i) {
        L[i] = 0;
    }

    /* initialize E */
    for (i=0; i<s1Len; ++i) {
        E[i] = NEG_INF_32;
    }

    /* iterate over database */
    for (j=0; j<s2Len; ++j) {
        const int * const restrict matcol = matrix[s2[j]];
        int FM = 0;
        int FS = 0;
        int FL = 0;
        /* calculate E */
        for (i=0; i<s1Len; ++i) {
            E[i] = MAX(E[i]-gap, H[i]-open);
        }
        /* calculate Ht */
        for (i=0; i<s1Len; ++i) {
            int tmp = H[i-1]+matcol[s1[i]];
            Ht[i] = MAX(tmp, E[i]);
            Mt[i] = tmp >= E[i] ? M[i-1] + (s1[i]==s2[j]) : M[i];
            St[i] = tmp >= E[i] ? S[i-1] + (matcol[s1[i]] > 0) : S[i];
            Lt[i] = tmp >= E[i] ? L[i-1] + 1 : L[i] + 1;
            Ex[i] = (E[i] > tmp);
        }
        Ht[-1] = 0;
        Ft[-1] = NEG_INF_32;
        /* calculate Ft */
        for (i=0; i<s1Len; ++i) {
            Ft[i] = MAX(Ft[i-1]-gap, Ht[i-1]);
        }
        /* calculate H,M,S,L */
        for (i=0; i<s1Len; ++i) {
            int tmp = Ft[i]-open;
            H[i] = MAX(Ht[i], tmp);
            if ((Ht[i] == tmp && Ex[i]) || (Ht[i] < tmp)) {
                /* we favor F/up/del when F and E scores tie */
                M[i] = FM;
                S[i] = FS;
                L[i] = FL + 1;
            }
            else {
                M[i] = Mt[i];
                S[i] = St[i];
                L[i] = Lt[i];
            }
#ifdef PARASAIL_TABLE
            result->score_table[i*s2Len + j] = H[i];
            result->matches_table[i*s2Len + j] = M[i];
            result->similar_table[i*s2Len + j] = S[i];
            result->length_table[i*s2Len + j] = L[i];
#endif
            FM = M[i];
            FS = S[i];
            FL = L[i];
        }
        H[-1] = 0;
        /* last value from column */
        if (H[s1Len-1] > score) {
            score = H[s1Len-1];
            matches = M[s1Len-1];
            similar = S[s1Len-1];
            length = L[s1Len-1];
        }
    }
    /* max of last column */
    for (i=0; i<s1Len; ++i) {
        if (H[i] > score) {
            score = H[i];
            matches = M[i];
            similar = S[i];
            length = L[i];
        }
    }

    result->score = score;
    result->matches = matches;
    result->similar = similar;
    result->length = length;

    parasail_free(Ex);
    parasail_free(Lt);
    parasail_free(St);
    parasail_free(Mt);
    parasail_free(LB);
    parasail_free(SB);
    parasail_free(MB);
    parasail_free(FtB);
    parasail_free(HtB);
    parasail_free(E);
    parasail_free(HB);
    parasail_free(s2);
    parasail_free(s1);

    return result;
}

