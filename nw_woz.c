#include "config.h"

#include <stdlib.h>

#include <emmintrin.h>

#ifdef ALIGN_EXTRA
#include "align/align_debug.h"
#else
#include "align/align.h"
#endif
#include "blosum/blosum_map.h"


#define BLOSUM0_(ch1, ch2) (matrow0[(ch2)])
#define BLOSUM1_(ch1, ch2) (matrow1[(ch2)])
#define BLOSUM2_(ch1, ch2) (matrow2[(ch2)])
#define BLOSUM3_(ch1, ch2) (matrow3[(ch2)])
#define BLOSUM4_(ch1, ch2) (matrow4[(ch2)])
#define BLOSUM5_(ch1, ch2) (matrow5[(ch2)])
#define BLOSUM6_(ch1, ch2) (matrow6[(ch2)])
#define BLOSUM7_(ch1, ch2) (matrow7[(ch2)])

/* on OSX, _mm_extract_epi16 got wrong answer, but taking the union of
 * the vector and extracting that way seemed to work... */
#define EXTRACT extract
//#define EXTRACT _mm_extract_epi16

/* shift given vector v, insert val, return shifted val */
static inline __m128i vshift16(const __m128i v, const int val)
{
    __m128i ret = _mm_srli_si128(v, 2);
    ret = _mm_insert_epi16(ret, val, 7);
    return ret;
}


static inline int extract(const __m128i m, const int pos)
{
    union {
        __m128i m;
        int16_t v[8];
    } tmp;
    tmp.m = m;
    return tmp.v[pos];
}


#ifdef ALIGN_EXTRA
static inline void arr_store_si128(
        int *array,
        __m128i vWscore,
        int i,
        int s1Len,
        int j,
        int s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[(i+0)*s2Len + (j-0)] = _mm_extract_epi16(vWscore, 7);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[(i+1)*s2Len + (j-1)] = _mm_extract_epi16(vWscore, 6);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[(i+2)*s2Len + (j-2)] = _mm_extract_epi16(vWscore, 5);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[(i+3)*s2Len + (j-3)] = _mm_extract_epi16(vWscore, 4);
    }
    if (0 <= i+4 && i+4 < s1Len && 0 <= j-4 && j-4 < s2Len) {
        array[(i+4)*s2Len + (j-4)] = _mm_extract_epi16(vWscore, 3);
    }
    if (0 <= i+5 && i+5 < s1Len && 0 <= j-5 && j-5 < s2Len) {
        array[(i+5)*s2Len + (j-5)] = _mm_extract_epi16(vWscore, 2);
    }
    if (0 <= i+6 && i+6 < s1Len && 0 <= j-6 && j-6 < s2Len) {
        array[(i+6)*s2Len + (j-6)] = _mm_extract_epi16(vWscore, 1);
    }
    if (0 <= i+7 && i+7 < s1Len && 0 <= j-7 && j-7 < s2Len) {
        array[(i+7)*s2Len + (j-7)] = _mm_extract_epi16(vWscore, 0);
    }
}
#endif


#ifdef ALIGN_EXTRA
#define FNAME nw_woz_debug
#else
#define FNAME nw_woz
#endif

#if 1

int FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict _tbl_pr, int * const restrict _del_pr
#ifdef ALIGN_EXTRA
        , int * const restrict score_table
#endif
        )
{
    const int N = 8; /* number of values in vector */
    const int PAD2 = N-1; /* N 16-byte values in vector, so N - 1 */
    const int PAD = PAD2*2;
    int * const restrict s1 = (int * const restrict)malloc(sizeof(int)*(s1Len+PAD2));
    int * const restrict s2B= (int * const restrict)malloc(sizeof(int)*(s2Len+PAD));
    int * const restrict s2 = s2B+PAD2; /* will allow later for negative indices */
    int i = 0;
    int j = 0;
    int score = NEG_INF;
    int * const restrict tbl_pr = _tbl_pr+PAD2;
    int * const restrict del_pr = _del_pr+PAD2;
    __m128i vNegInf = _mm_set1_epi16(NEG_INF);
    __m128i vNegInf0 = _mm_srli_si128(vNegInf, 2); /* shift in a 0 */
    __m128i vOpen = _mm_set1_epi16(open);
    __m128i vGap  = _mm_set1_epi16(gap);
    __m128i vOne = _mm_set1_epi16(1);
    __m128i vN = _mm_set1_epi16(N);
    __m128i vNegOne = _mm_set1_epi16(-1);
    __m128i vI = _mm_set_epi16(0,1,2,3,4,5,6,7);
    __m128i vJreset = _mm_set_epi16(0,-1,-2,-3,-4,-5,-6,-7);
    __m128i vMax = vNegInf;
    __m128i vILimit = _mm_set1_epi16(s1Len);
    __m128i vILimit1 = _mm_sub_epi16(vILimit, vOne);
    __m128i vJLimit = _mm_set1_epi16(s2Len);
    __m128i vJLimit1 = _mm_sub_epi16(vJLimit, vOne);
    __m128i vIBoundary = _mm_set_epi16(
            -open-0*gap,
            -open-1*gap,
            -open-2*gap,
            -open-3*gap,
            -open-4*gap,
            -open-5*gap,
            -open-6*gap,
            -open-7*gap);

    /* convert _s1 from char to int in range 0-23 */
    for (i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[(unsigned char)_s1[i]];
    }
    /* pad back of s1 with dummy values */
    for (i=s1Len; i<s1Len+PAD2; ++i) {
        s1[i] = 0; /* point to first matrix row because we don't care */
    }

    /* convert _s2 from char to int in range 0-23 */
    for (j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[(unsigned char)_s2[j]];
    }
    /* pad front of s2 with dummy values */
    for (j=-PAD2; j<0; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }
    /* pad back of s2 with dummy values */
    for (j=s2Len; j<s2Len+PAD2; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }

    /* set initial values for stored row */
    for (j=0; j<s2Len; ++j) {
        tbl_pr[j] = -open - j*gap;
        del_pr[j] = NEG_INF;
    }
    /* pad front of stored row values */
    for (j=-PAD2; j<0; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD2; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
    }
    tbl_pr[-1] = 0; /* upper left corner */

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        __m128i vNscore = vNegInf;
        __m128i vWscore = vNegInf;
        __m128i vIns = vNegInf;
        __m128i vDel = vNegInf;
        __m128i vJ = vJreset;
        const int * const restrict matrow0 = matrix[s1[i+0]];
        const int * const restrict matrow1 = matrix[s1[i+1]];
        const int * const restrict matrow2 = matrix[s1[i+2]];
        const int * const restrict matrow3 = matrix[s1[i+3]];
        const int * const restrict matrow4 = matrix[s1[i+4]];
        const int * const restrict matrow5 = matrix[s1[i+5]];
        const int * const restrict matrow6 = matrix[s1[i+6]];
        const int * const restrict matrow7 = matrix[s1[i+7]];
        vNscore = vshift16(vNscore, tbl_pr[-1]);
        vWscore = vshift16(vWscore, -open - i*gap);
        tbl_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD2; ++j) {
            __m128i vMat;
            __m128i vNWscore = vNscore;
            vNscore = vshift16(vWscore, tbl_pr[j]);
            vDel = vshift16(vDel, del_pr[j]);
            vDel = _mm_max_epi16(
                    _mm_sub_epi16(vNscore, vOpen),
                    _mm_sub_epi16(vDel, vGap));
            vIns = _mm_max_epi16(
                    _mm_sub_epi16(vWscore, vOpen),
                    _mm_sub_epi16(vIns, vGap));
            vMat = _mm_set_epi16(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]],
                    matrow4[s2[j-4]],
                    matrow5[s2[j-5]],
                    matrow6[s2[j-6]],
                    matrow7[s2[j-7]]
                    );
            vWscore = _mm_add_epi16(vNWscore, vMat);
            vWscore = _mm_max_epi16(vWscore, vIns);
            vWscore = _mm_max_epi16(vWscore, vDel);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m128i cond = _mm_cmpeq_epi16(vJ,vNegOne);
                vWscore = _mm_andnot_si128(cond, vWscore); /* all but j=-1 */
                vWscore = _mm_or_si128(vWscore,
                        _mm_and_si128(cond, vIBoundary));
                vDel = _mm_andnot_si128(cond, vDel);
                vDel = _mm_or_si128(vDel, _mm_and_si128(cond, vNegInf));
                vIns = _mm_andnot_si128(cond, vIns);
                vIns = _mm_or_si128(vIns, _mm_and_si128(cond, vNegInf));
            }
#ifdef ALIGN_EXTRA
            arr_store_si128(score_table, vWscore, i, s1Len, j, s2Len);
#endif
            tbl_pr[j-7] = _mm_extract_epi16(vWscore,0);
            del_pr[j-7] = _mm_extract_epi16(vDel,0);
            /* as minor diagonal vector passes across table, extract
               last table value at the i,j bound */
            {
                __m128i cond_valid_I = _mm_cmpeq_epi16(vI, vILimit1);
                __m128i cond_valid_J = _mm_cmpeq_epi16(vJ, vJLimit1);
                __m128i cond_max = _mm_cmpgt_epi16(vWscore, vMax);
                __m128i cond_all = _mm_and_si128(cond_max,
                        _mm_and_si128(cond_valid_I, cond_valid_J));
                vMax = _mm_andnot_si128(cond_all, vMax); /* keep old */
                vMax = _mm_or_si128(vMax,
                        _mm_and_si128(cond_all, vWscore));
            }
            vJ = _mm_add_epi16(vJ, vOne);
        }
        vI = _mm_add_epi16(vI, vN);
        vIBoundary = _mm_sub_epi16(vIBoundary, _mm_mullo_epi16(vN, vGap));
    }

    /* max in vMax */
    {
        int16_t value;
        value = (int16_t) _mm_extract_epi16(vMax, 0);
        if (value > score) {
            score = value;
        }
        value = (int16_t) _mm_extract_epi16(vMax, 1);
        if (value > score) {
            score = value;
        }
        value = (int16_t) _mm_extract_epi16(vMax, 2);
        if (value > score) {
            score = value;
        }
        value = (int16_t) _mm_extract_epi16(vMax, 3);
        if (value > score) {
            score = value;
        }
        value = (int16_t) _mm_extract_epi16(vMax, 4);
        if (value > score) {
            score = value;
        }
        value = (int16_t) _mm_extract_epi16(vMax, 5);
        if (value > score) {
            score = value;
        }
        value = (int16_t) _mm_extract_epi16(vMax, 6);
        if (value > score) {
            score = value;
        }
        value = (int16_t) _mm_extract_epi16(vMax, 7);
        if (value > score) {
            score = value;
        }
    }

    free(s1);
    free(s2B);

    return score;
}

#else

int FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict tbl_pr, int * const restrict del_pr
#ifdef ALIGN_EXTRA
        , int * const restrict score_table
#endif
        )
{
    int * const restrict s1 = (int * const restrict)malloc(sizeof(int)*(s1Len+8));
    int * const restrict s2 = (int * const restrict)malloc(sizeof(int)*s2Len);
    int i = 0;
    int j = 0;
    int score = 0;

    for (i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[(unsigned char)_s1[i]];
    }
    for (i=s1Len; i<s1Len+8; ++i) {
        s1[i] = 23;
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[(unsigned char)_s2[j]];
    }

    /* dummy padding */
    for (j=0; j<7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
    }

    /* upper left corner */
    tbl_pr[7] = 0;
    del_pr[7] = NEG_INF;

    /* first row */
    for (j=8; j<s2Len+8; ++j) {
        tbl_pr[j] = -open -(j-8)*gap;
        del_pr[j] = NEG_INF;
    }

    /* dummy padding */
    for (j=s2Len+8; j<s2Len+8+7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
    }

    __m128i vOpen = _mm_set1_epi16(open);
    __m128i vGap  = _mm_set1_epi16(gap);

    /* iter over first sequence */
    for (i=1; i<=s1Len; i+=8) {
        int j;
        __m128i NWscore = _mm_set1_epi16(NEG_INF);
        __m128i Nscore  = _mm_set1_epi16(NEG_INF);
        __m128i Wscore  = _mm_set1_epi16(NEG_INF);
        __m128i vTbl    = _mm_set1_epi16(NEG_INF);
        __m128i vDel    = _mm_set1_epi16(NEG_INF);
        __m128i vIns    = _mm_set1_epi16(NEG_INF);
        __m128i vMat;

        const int * const restrict matrow0 = matrix[s1[i-1+0]];
        const int * const restrict matrow1 = matrix[s1[i-1+1]];
        const int * const restrict matrow2 = matrix[s1[i-1+2]];
        const int * const restrict matrow3 = matrix[s1[i-1+3]];
        const int * const restrict matrow4 = matrix[s1[i-1+4]];
        const int * const restrict matrow5 = matrix[s1[i-1+5]];
        const int * const restrict matrow6 = matrix[s1[i-1+6]];
        const int * const restrict matrow7 = matrix[s1[i-1+7]];

        /* j = 0 */
        j = 0;
        Nscore = _mm_insert_epi16(Nscore, tbl_pr[j+7], 7);
        Wscore = _mm_insert_epi16(Wscore, -open-(i-1)*gap, 7);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 7);
#ifdef ALIGN_EXTRA
        score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wscore, 7);
#endif

        /* j = 1 */
        j = 1;
        NWscore = Nscore;
        Nscore  = vshift16(Wscore, tbl_pr[j+7]);
        vDel    = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1],s2[j-1]),
            0,
            0,
            0,
            0,
            0,
            0,
            0
        );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+0)*gap, 6);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 6);
#ifdef ALIGN_EXTRA
        score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wscore, 6);
#endif

        /* j = 2 */
        j = 2;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            0,
            0,
            0,
            0,
            0,
            0
        );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+1)*gap, 5);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 5);
#ifdef ALIGN_EXTRA
        score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wscore, 5);
#endif

        /* j = 3 */
        j = 3;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            0,
            0,
            0,
            0,
            0
        );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+2)*gap, 4);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 4);
#ifdef ALIGN_EXTRA
        score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wscore, 4);
#endif

        /* j = 4 */
        j = 4;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            0,
            0,
            0,
            0
        );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+3)*gap, 3);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 3);
#ifdef ALIGN_EXTRA
        score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wscore, 3);
#endif

        /* j = 5 */
        j = 5;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            BLOSUM4_(s1[i-1+4],s2[j-1-4]),
            0,
            0,
            0
        );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+4)*gap, 2);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 2);
#ifdef ALIGN_EXTRA
        score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wscore, 2);
#endif

        /* j = 6 */
        j = 6;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            BLOSUM4_(s1[i-1+4],s2[j-1-4]),
            BLOSUM5_(s1[i-1+5],s2[j-1-5]),
            0,
            0
        );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+5)*gap, 1);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 1);
#ifdef ALIGN_EXTRA
        score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wscore, 1);
#endif

        /* j = 7 */
        j = 7;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
            BLOSUM0_(s1[i-1+0],s2[j-1-0]),
            BLOSUM1_(s1[i-1+1],s2[j-1-1]),
            BLOSUM2_(s1[i-1+2],s2[j-1-2]),
            BLOSUM3_(s1[i-1+3],s2[j-1-3]),
            BLOSUM4_(s1[i-1+4],s2[j-1-4]),
            BLOSUM5_(s1[i-1+5],s2[j-1-5]),
            BLOSUM6_(s1[i-1+6],s2[j-1-6]),
            0
        );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, -open-(i+6)*gap, 0);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 0);
        tbl_pr[7] = -open-(i+6)*gap;
#ifdef ALIGN_EXTRA
        score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wscore, 0);
#endif

        for (j=8; j<=s2Len; ++j) {
            NWscore = Nscore;
            Nscore = vshift16(Wscore, tbl_pr[j+7]);
            vDel = vshift16(vDel, del_pr[j+7]);
            vDel = _mm_max_epi16(
                    _mm_sub_epi16(Nscore, vOpen),
                    _mm_sub_epi16(vDel, vGap));
            vIns = _mm_max_epi16(
                    _mm_sub_epi16(Wscore,vOpen),
                    _mm_sub_epi16(vIns,vGap));
            vMat = _mm_set_epi16(
                    BLOSUM0_(s1[i-1+0],s2[j-1-0]),
                    BLOSUM1_(s1[i-1+1],s2[j-1-1]),
                    BLOSUM2_(s1[i-1+2],s2[j-1-2]),
                    BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                    BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                    BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                    BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                    BLOSUM7_(s1[i-1+7],s2[j-1-7])
                    );
            vTbl = _mm_add_epi16(NWscore, vMat);
            vTbl = _mm_max_epi16(vTbl, vDel);
            vTbl = _mm_max_epi16(vTbl, vIns);
            Wscore = vTbl;
            tbl_pr[j] = EXTRACT(vTbl, 0);
            del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
            score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
            score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
            score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
            score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
            score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
            score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
            score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
            score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        }
        if (i+0 == s1Len) score = EXTRACT(vTbl, 7);

        /* j = s2Len + 1 */
        j = s2Len + 1;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
                0,
                BLOSUM1_(s1[i-1+1],s2[j-1-1]),
                BLOSUM2_(s1[i-1+2],s2[j-1-2]),
                BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+1 == s1Len) score = EXTRACT(vTbl, 6);

        /* j = s2Len + 2 */
        j = s2Len + 2;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
                0,
                0,
                BLOSUM2_(s1[i-1+2],s2[j-1-2]),
                BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+2 == s1Len) score = EXTRACT(vTbl, 5);

        /* j = s2Len + 3 */
        j = s2Len + 3;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                BLOSUM3_(s1[i-1+3],s2[j-1-3]),
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+3 == s1Len) score = EXTRACT(vTbl, 4);

        /* j = s2Len + 4 */
        j = s2Len + 4;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                BLOSUM4_(s1[i-1+4],s2[j-1-4]),
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+4 == s1Len) score = EXTRACT(vTbl, 3);

        /* j = s2Len + 5 */
        j = s2Len + 5;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                BLOSUM5_(s1[i-1+5],s2[j-1-5]),
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+5 == s1Len) score = EXTRACT(vTbl, 2);

        /* j = s2Len + 6 */
        j = s2Len + 6;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                0,
                BLOSUM6_(s1[i-1+6],s2[j-1-6]),
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+6 == s1Len) score = EXTRACT(vTbl, 1);

        /* j = s2Len + 7 */
        j = s2Len + 7;
        NWscore = Nscore;
        Nscore = vshift16(Wscore, tbl_pr[j+7]);
        vDel = vshift16(vDel, del_pr[j+7]);
        vDel = _mm_max_epi16(
                _mm_sub_epi16(Nscore, vOpen),
                _mm_sub_epi16(vDel, vGap));
        vIns = _mm_max_epi16(
                _mm_sub_epi16(Wscore,vOpen),
                _mm_sub_epi16(vIns,vGap));
        vMat = _mm_set_epi16(
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                BLOSUM7_(s1[i-1+7],s2[j-1-7])
                );
        vTbl = _mm_add_epi16(NWscore, vMat);
        vTbl = _mm_max_epi16(vTbl, vDel);
        vTbl = _mm_max_epi16(vTbl, vIns);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+7 == s1Len) score = EXTRACT(vTbl, 0);
    }

    free(s1);
    free(s2);
    return score;
}
#endif

