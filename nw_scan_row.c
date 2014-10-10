#include "config.h"

#include <stdlib.h>

#ifdef ALIGN_EXTRA
#include "align/align_debug.h"
#else
#include "align/align.h"
#endif
#include "blosum/blosum_map.h"

#ifdef ALIGN_EXTRA
#define FNAME nw_scan_row_debug
#else
#define FNAME nw_scan_row
#endif

int FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict _H, int * const restrict F
#ifdef ALIGN_EXTRA
        , int * const restrict score_table
#endif
        )
{
    int * const restrict s1 = (int * const restrict)malloc(sizeof(int)*s1Len);
    int * const restrict s2 = (int * const restrict)malloc(sizeof(int)*s2Len);
    int * const restrict HtB= (int * const restrict)malloc(sizeof(int)*(s2Len+1));
    int * const restrict Ht = HtB+1;
    int * const restrict EtB= (int * const restrict)malloc(sizeof(int)*(s2Len+1));
    int * const restrict Et = EtB+1;
    int i = 0;
    int j = 0;
    int * const restrict H = _H+1;

    for (i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[(unsigned char)_s1[i]];
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[(unsigned char)_s2[j]];
    }

    /* initialize H */
    H[-1] = 0;
    for (j=0; j<s2Len; ++j) {
        H[j] = -open - j*gap;
    }

    /* initialize F */
    for (j=0; j<s2Len; ++j) {
        F[j] = NEG_INF_32;
    }

    /* iterate over query */
    for (i=0; i<s1Len; ++i) {
        const int * const restrict matrow = matrix[s1[i]];
        /* calculate F */
        for (j=0; j<s2Len; ++j) {
            F[j] = MAX(F[j]-gap, H[j]-open);
        }
        /* calculate Ht */
        for (j=0; j<s2Len; ++j) {
            Ht[j] = MAX(H[j-1]+matrow[s2[j]], F[j]);
        }
        Ht[-1] = -open -i*gap;
        Et[-1] = NEG_INF_32;
        /* calculate Et */
        for (j=0; j<s2Len; ++j) {
            Et[j] = MAX(Et[j-1]-gap, Ht[j-1]);
        }
        /* calculate H */
        for (j=0; j<s2Len; ++j) {
            H[j] = MAX(Ht[j], Et[j]-open);
#ifdef ALIGN_EXTRA
            score_table[i*s2Len + j] = H[j];
#endif
        }
        H[-1] = -open - i*gap;
    }

    free(s1);
    free(s2);
    free(HtB);
    free(EtB);

    return H[s2Len-1];
}
