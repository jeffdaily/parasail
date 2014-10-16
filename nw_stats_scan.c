#include "config.h"

#include <stdlib.h>

#ifdef ALIGN_EXTRA
#include "align/align_debug.h"
#else
#include "align/align.h"
#endif
#include "blosum/blosum_map.h"

#ifdef ALIGN_EXTRA
#define FNAME nw_stats_scan_debug
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
        int * const restrict _M, int * const restrict L
#ifdef ALIGN_EXTRA
        , int * const restrict score_table
        , int * const restrict match_table
        , int * const restrict length_table
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
    int * const restrict H = _H+1;
    int * const restrict M = _M+1;

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
    for (i=0; i<s1Len; ++i) {
        L[i] = 0;
    }

    /* initialize E */
    for (i=0; i<s1Len; ++i) {
        E[i] = NEG_INF_32;
    }

    /* iterate over database */
    for (j=0; j<s2Len; ++j) {
        int Wmatches = M[-1];
        int Nmatches = 0;
        const int * const restrict matcol = matrix[s2[j]];
        /* calculate E */
        for (i=0; i<s1Len; ++i) {
            E[i] = MAX(E[i]-gap, H[i]-open);
        }
        /* calculate Ht */
        for (i=0; i<s1Len; ++i) {
            Ht[i] = MAX(H[i-1]+matcol[s1[i]], E[i]);
        }
        Ht[-1] = -open -j*gap;
        Ft[-1] = NEG_INF_32;
        /* calculate Ft */
        for (i=0; i<s1Len; ++i) {
            Ft[i] = MAX(Ft[i-1]-gap, Ht[i-1]);
        }
        /* calculate H,M,L */
        for (i=0; i<s1Len; ++i) {
            int NWmatches = Wmatches;
            Wmatches = M[i];
            H[i] = MAX(Ht[i], Ft[i]-open);
#ifdef ALIGN_EXTRA
            score_table[i*s2Len + j] = H[i];
#endif
            //if (H[i] >= Ft[i]-open && H[i] >= E[i]) {
            //    Nmatches = NWmatches + (s1[i] == s2[j]);
            //}
            //else if (Ft[i]-open >= E[i]) {
            //    Nmatches = Nmatches;
            //}
            //else {
            //    Nmatches = Wmatches;
            //}
            M[i] = Nmatches;
#ifdef ALIGN_EXTRA
            match_table[i*s2Len + j] = Nmatches;
#endif
        }
        H[-1] = -open - j*gap;
        /* calculate M, L */
        for (i=0; i<s1Len; ++i) {
        }
    }

    free(s1);
    free(s2);
    free(HtB);
    free(FtB);

    *matches = M[s2Len];
    *length = L[s2Len];
    return H[s1Len-1];
}
