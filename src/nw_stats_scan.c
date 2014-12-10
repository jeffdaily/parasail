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

#include <stdlib.h>

#ifdef PARASAIL_TABLE
#include "align_table.h"
#else
#include "align.h"
#endif
#include "blosum/blosum_map.h"

#ifdef PARASAIL_TABLE
#define FNAME nw_stats_scan_table
#else
#define FNAME nw_stats_scan
#endif

int FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * matches, int * length,
        int * const restrict _H, int * const restrict E,
        int * const restrict _M, int * const restrict _L
#ifdef PARASAIL_TABLE
        , int * const restrict result->score_table
        , int * const restrict result->matches_table
        , int * const restrict result->length_table
#endif
        )
{
    int * const restrict s1 = (int * const restrict)malloc(sizeof(int)*s1Len);
    int * const restrict s2 = (int * const restrict)malloc(sizeof(int)*s2Len);
    int * const restrict HtB= (int * const restrict)malloc(sizeof(int)*(s1Len+1));
    int * const restrict Ht = HtB+1;
    int * const restrict FtB= (int * const restrict)malloc(sizeof(int)*(s1Len+1));
    int * const restrict Ft = FtB+1;
    int i = 0;
    int j = 0;
    int * const restrict H  = _H+1;
    int * const restrict M  = _M+1;
    int * const restrict L  = _L+1;
    int * const restrict Mt = (int * const restrict)malloc(sizeof(int)*s1Len);
    int * const restrict Lt = (int * const restrict)malloc(sizeof(int)*s1Len);
    int * const restrict Ex = (int * const restrict)malloc(sizeof(int)*s1Len);

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
        H[i] = -open - i*gap;
    }

    /* initialize M */
    M[-1] = 0;
    for (i=0; i<s1Len; ++i) {
        M[i] = 0;
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
            Lt[i] = tmp >= E[i] ? L[i-1] + 1 : L[i] + 1;
            Ex[i] = (E[i] > tmp);
        }
        Ht[-1] = -open -j*gap;
        Ft[-1] = NEG_INF_32;
        /* calculate Ft */
        for (i=0; i<s1Len; ++i) {
            Ft[i] = MAX(Ft[i-1]-gap, Ht[i-1]);
        }
        /* calculate H,M,L */
        for (i=0; i<s1Len; ++i) {
            int tmp = Ft[i]-open;
            H[i] = MAX(Ht[i], tmp);
            if ((Ht[i] == tmp && Ex[i]) || (Ht[i] < tmp)) {
                /* we favor F/up/del when F and E scores tie */
                M[i] = FM;
                L[i] = FL + 1;
            }
            else {
                M[i] = Mt[i];
                L[i] = Lt[i];
            }
#ifdef PARASAIL_TABLE
            result->score_table[i*s2Len + j] = H[i];
            result->matches_table[i*s2Len + j] = M[i];
            result->length_table[i*s2Len + j] = L[i];
#endif
            FM = M[i];
            FL = L[i];
        }
        H[-1] = -open - j*gap;
    }

    free(s1);
    free(s2);
    free(HtB);
    free(FtB);
    free(Mt);
    free(Lt);
    free(Ex);

    *matches = M[s1Len-1];
    *length = L[s1Len-1];
    return H[s1Len-1];
}
