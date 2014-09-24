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

/* shift given vector v, return shifted val */
static inline __m128i vshift16_(const __m128i v)
{
    return _mm_srli_si128(v, 2);
}


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
#define FNAME nw_stats_woz_debug
#else
#define FNAME nw_stats_woz
#endif

int FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24],
        int * const restrict _matches, int * const restrict _length,
        int * const restrict tbl_pr, int * const restrict del_pr,
        int * const restrict mch_pr, int * const restrict len_pr
#ifdef ALIGN_EXTRA
        , int * const restrict score_table
        , int * const restrict match_table
        , int * const restrict length_table
#endif
        )
{
    int * const restrict s1 = (int * const restrict)malloc(sizeof(int)*(s1Len+8));
    int * const restrict s2 = (int * const restrict)malloc(sizeof(int)*s2Len);
    int i = 0;
    int j = 0;
    int score = 0;
    int match = 0;
    int length = 0;

    for (i=0; i<s1Len; ++i) {
        s1[i] = MAP_BLOSUM_[_s1[i]];
    }
    for (i=s1Len; i<s1Len+8; ++i) {
        s1[i] = 23;
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = MAP_BLOSUM_[_s2[j]];
    }

    /* dummy padding */
    for (j=0; j<7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
    }

    /* upper left corner */
    tbl_pr[7] = 0;
    del_pr[7] = NEG_INF;
    mch_pr[7] = 0;
    len_pr[7] = 0;

    /* first row */
    for (j=8; j<s2Len+8; ++j) {
        tbl_pr[j] = -open -(j-8)*gap;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
    }

    /* dummy padding */
    for (j=s2Len+8; j<s2Len+8+7; ++j) {
        tbl_pr[j] = NEG_INF;
        del_pr[j] = NEG_INF;
        mch_pr[j] = 0;
        len_pr[j] = 0;
    }

    __m128i vOpen = _mm_set1_epi16(open);
    __m128i vGap  = _mm_set1_epi16(gap);

    /* iter over first sequence */
    for (i=1; i<=s1Len; i+=8) {
        int j;
        __m128i NWscore  = _mm_set1_epi16(NEG_INF);
        __m128i NWmatch  = _mm_set1_epi16(NEG_INF);
        __m128i NWlength = _mm_set1_epi16(NEG_INF);
        __m128i Nscore   = _mm_set1_epi16(NEG_INF);
        __m128i Nmatch   = _mm_set1_epi16(NEG_INF);
        __m128i Nlength  = _mm_set1_epi16(NEG_INF);
        __m128i Wscore   = _mm_set1_epi16(NEG_INF);
        __m128i Wmatch   = _mm_set1_epi16(NEG_INF);
        __m128i Wlength  = _mm_set1_epi16(NEG_INF);
        __m128i vTbl     = _mm_set1_epi16(NEG_INF);
        __m128i vDel     = _mm_set1_epi16(NEG_INF);
        __m128i vIns     = _mm_set1_epi16(NEG_INF);
        __m128i vMat;
        __m128i vs1;
        __m128i vs2;
        __m128i vOne     = _mm_set1_epi16(1);
        __m128i case1not;
        __m128i case2not;
        __m128i case2;
        __m128i case3;
        __m128i Cscore;
        __m128i Cmatch;
        __m128i Clength;

        const int * const restrict matrow0 = matrix[s1[i-1+0]];
        const int * const restrict matrow1 = matrix[s1[i-1+1]];
        const int * const restrict matrow2 = matrix[s1[i-1+2]];
        const int * const restrict matrow3 = matrix[s1[i-1+3]];
        const int * const restrict matrow4 = matrix[s1[i-1+4]];
        const int * const restrict matrow5 = matrix[s1[i-1+5]];
        const int * const restrict matrow6 = matrix[s1[i-1+6]];
        const int * const restrict matrow7 = matrix[s1[i-1+7]];

        vs1 = _mm_set_epi16(
                s1[i-1+0],
                s1[i-1+1],
                s1[i-1+2],
                s1[i-1+3],
                s1[i-1+4],
                s1[i-1+5],
                s1[i-1+6],
                s1[i-1+7]
                );

        /* j = 0 */
        j = 0;
        Nscore = _mm_insert_epi16(Nscore, tbl_pr[j+7], 7);
        Nmatch = _mm_insert_epi16(Nscore, mch_pr[j+7], 7);
        Nlength= _mm_insert_epi16(Nscore, len_pr[j+7], 7);
        Wscore = _mm_insert_epi16(Wscore, -open-(i-1)*gap, 7);
        Wmatch = _mm_insert_epi16(Wscore, 0, 7);
        Wlength= _mm_insert_epi16(Wscore, 0, 7);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 7);
#ifdef ALIGN_EXTRA
        match_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        length_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength,7);
        score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wscore, 7);
#endif

        /* j = 1 */
        j = 1;
        /* this block never changes, so let's not copy-and-paste */
#ifdef SETUP_BLOCK
#undef SETUP_BLOCK
#endif
#define SETUP_BLOCK                                 \
        NWscore = Nscore;                           \
        NWmatch = Nmatch;                           \
        NWlength= Nlength;                          \
        Nscore  = vshift16(Wscore, tbl_pr[j+7]);    \
        Nmatch  = vshift16(Wmatch, mch_pr[j+7]);    \
        Nlength = vshift16(Wlength,len_pr[j+7]);    \
        vDel    = vshift16(vDel,   del_pr[j+7]);    \
        vDel = _mm_max_epi16(                       \
                _mm_sub_epi16(Nscore, vOpen),       \
                _mm_sub_epi16(vDel, vGap));         \
        vIns = _mm_max_epi16(                       \
                _mm_sub_epi16(Wscore,vOpen),        \
                _mm_sub_epi16(vIns,vGap));
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
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
        /* we reuse this logic throughout this function,
         * best to not copy-and-paste everywhere */
#ifdef CONDITIONAL_BLOCK
#undef CONDITIONAL_BLOCK
#endif
#define CONDITIONAL_BLOCK                                                   \
        vTbl = _mm_add_epi16(NWscore, vMat);                                \
        case1not = _mm_or_si128(                                            \
                _mm_cmplt_epi16(vTbl,vDel),_mm_cmplt_epi16(vTbl,vIns));     \
        case2not = _mm_cmplt_epi16(vDel,vIns);                              \
        case2 = _mm_andnot_si128(case2not,case1not);                        \
        case3 = _mm_and_si128(case1not,case2not);                           \
        Cscore = _mm_andnot_si128(case1not, vTbl);                          \
        Cmatch = _mm_andnot_si128(case1not,                                 \
                    _mm_add_epi16(NWmatch, _mm_and_si128(                   \
                            _mm_cmpeq_epi16(vs1,vs2),vOne)));               \
        Clength= _mm_andnot_si128(case1not, _mm_add_epi16(NWlength, vOne)); \
        Cscore = _mm_or_si128(Cscore, _mm_and_si128(case2, vDel));          \
        Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(case2, Nmatch));        \
        Clength= _mm_or_si128(Clength,_mm_and_si128(case2,                  \
                    _mm_add_epi16(Nlength, vOne)));                         \
        Cscore = _mm_or_si128(Cscore, _mm_and_si128(case3, vIns));          \
        Cmatch = _mm_or_si128(Cmatch, _mm_and_si128(case3, Wmatch));        \
        Clength= _mm_or_si128(Clength,_mm_and_si128(case3,                  \
                    _mm_add_epi16(Wlength, vOne)));                         \
        Wscore = vTbl = Cscore;                                             \
        Wmatch = Cmatch;                                                    \
        Wlength= Clength;
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, -open-(i+0)*gap, 6);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 6);
        Wlength= _mm_insert_epi16(Wlength,0, 6);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 6);
#ifdef ALIGN_EXTRA
        match_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        match_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        length_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        length_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wscore, 6);
#endif

        /* j = 2 */
        j = 2;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
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
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, -open-(i+1)*gap, 5);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 5);
        Wlength= _mm_insert_epi16(Wlength,0, 5);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 5);
#ifdef ALIGN_EXTRA
        match_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        match_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        match_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        length_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        length_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        length_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wscore, 5);
#endif

        /* j = 3 */
        j = 3;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
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
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, -open-(i+2)*gap, 4);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 4);
        Wlength= _mm_insert_epi16(Wlength,0, 4);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 4);
#ifdef ALIGN_EXTRA
        match_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        match_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        match_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        match_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        length_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        length_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        length_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        length_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wscore, 4);
#endif

        /* j = 4 */
        j = 4;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
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
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, -open-(i+3)*gap, 3);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 3);
        Wlength= _mm_insert_epi16(Wlength,0, 3);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 3);
#ifdef ALIGN_EXTRA
        match_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        match_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        match_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        match_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        match_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        length_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        length_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        length_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        length_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        length_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wscore, 3);
#endif

        /* j = 5 */
        j = 5;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
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
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, -open-(i+4)*gap, 2);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 2);
        Wlength= _mm_insert_epi16(Wlength,0, 2);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 2);
#ifdef ALIGN_EXTRA
        match_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        match_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        match_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        match_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        match_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        match_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        length_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        length_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        length_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        length_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        length_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        length_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        score_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(vTbl, 7);
        score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wscore, 2);
#endif

        /* j = 6 */
        j = 6;
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
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
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, -open-(i+5)*gap, 1);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 1);
        Wlength= _mm_insert_epi16(Wlength,0, 1);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 1);
#ifdef ALIGN_EXTRA
        match_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        match_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        match_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        match_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        match_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        match_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        match_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        length_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        length_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        length_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        length_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        length_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        length_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        length_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
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
        SETUP_BLOCK
        vs2 = vshift16(vs2, s2[j-1]);
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
        CONDITIONAL_BLOCK
        Wscore = _mm_insert_epi16(Wscore, -open-(i+6)*gap, 0);
        Wmatch = _mm_insert_epi16(Wmatch, 0, 0);
        Wlength= _mm_insert_epi16(Wlength,0, 0);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 0);
        tbl_pr[7] = -open-(i+6)*gap;
        mch_pr[7] = 0;
        len_pr[7] = 0;
#ifdef ALIGN_EXTRA
        match_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
        match_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        match_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        match_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        match_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        match_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        match_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        match_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
        length_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
        length_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        length_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        length_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        length_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        length_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        length_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        length_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
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
            SETUP_BLOCK
            vs2 = vshift16(vs2, s2[j-1]);
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
            CONDITIONAL_BLOCK
            tbl_pr[j] = EXTRACT(vTbl, 0);
            mch_pr[j] = EXTRACT(Wmatch, 0);
            len_pr[j] = EXTRACT(Wlength, 0);
            del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
            match_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wmatch, 7);
            match_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
            match_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
            match_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
            match_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
            match_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
            match_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
            match_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
            length_table[(i+0)*(s2Len+14+1) + j + 7] = EXTRACT(Wlength, 7);
            length_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
            length_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
            length_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
            length_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
            length_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
            length_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
            length_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
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
        if (i+0 == s1Len) {
            score = EXTRACT(vTbl, 7);
            match = EXTRACT(Wmatch, 7);
            length= EXTRACT(Wlength, 7);
        }

        /* j = s2Len + 1 */
        j = s2Len + 1;
        SETUP_BLOCK
        vs2 = vshift16_(vs2);
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
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        match_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wmatch, 6);
        match_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        match_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        match_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        match_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        match_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        match_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
        length_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(Wlength, 6);
        length_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        length_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        length_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        length_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        length_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        length_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
        score_table[(i+1)*(s2Len+14+1) + j + 6] = EXTRACT(vTbl, 6);
        score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+1 == s1Len) {
            score = EXTRACT(vTbl, 6);
            match = EXTRACT(Wmatch, 6);
            length= EXTRACT(Wlength, 6);
        }

        /* j = s2Len + 2 */
        j = s2Len + 2;
        SETUP_BLOCK
        vs2 = vshift16_(vs2);
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
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        match_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wmatch, 5);
        match_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        match_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        match_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        match_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        match_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
        length_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(Wlength, 5);
        length_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        length_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        length_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        length_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        length_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
        score_table[(i+2)*(s2Len+14+1) + j + 5] = EXTRACT(vTbl, 5);
        score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+2 == s1Len) {
            score = EXTRACT(vTbl, 5);
            match = EXTRACT(Wmatch, 5);
            length= EXTRACT(Wlength, 5);
        }

        /* j = s2Len + 3 */
        j = s2Len + 3;
        SETUP_BLOCK
        vs2 = vshift16_(vs2);
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
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        match_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wmatch, 4);
        match_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        match_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        match_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        match_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
        length_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(Wlength, 4);
        length_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        length_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        length_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        length_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
        score_table[(i+3)*(s2Len+14+1) + j + 4] = EXTRACT(vTbl, 4);
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+3 == s1Len) {
            score = EXTRACT(vTbl, 4);
            match = EXTRACT(Wmatch, 4);
            length= EXTRACT(Wlength, 4);
        }

        /* j = s2Len + 4 */
        j = s2Len + 4;
        SETUP_BLOCK
        vs2 = vshift16_(vs2);
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
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        match_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wmatch, 3);
        match_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        match_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        match_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
        length_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(Wlength, 3);
        length_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        length_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        length_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+4 == s1Len) {
            score = EXTRACT(vTbl, 3);
            match = EXTRACT(Wmatch, 3);
            length= EXTRACT(Wlength, 3);
        }

        /* j = s2Len + 5 */
        j = s2Len + 5;
        SETUP_BLOCK
        vs2 = vshift16_(vs2);
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
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        match_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wmatch, 2);
        match_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        match_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
        length_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(Wlength, 2);
        length_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        length_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+5 == s1Len) {
            score = EXTRACT(vTbl, 2);
            match = EXTRACT(Wmatch, 2);
            length= EXTRACT(Wlength, 2);
        }

        /* j = s2Len + 6 */
        j = s2Len + 6;
        SETUP_BLOCK
        vs2 = vshift16_(vs2);
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
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        match_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wmatch, 1);
        match_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
        length_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(Wlength, 1);
        length_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+6 == s1Len) {
            score = EXTRACT(vTbl, 1);
            match = EXTRACT(Wmatch, 1);
            length= EXTRACT(Wlength, 1);
        }

        /* j = s2Len + 7 */
        j = s2Len + 7;
        SETUP_BLOCK
        vs2 = vshift16_(vs2);
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
        CONDITIONAL_BLOCK
        tbl_pr[j] = EXTRACT(vTbl, 0);
        mch_pr[j] = EXTRACT(Wmatch, 0);
        len_pr[j] = EXTRACT(Wlength, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        match_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wmatch, 0);
        length_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(Wlength, 0);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
        if (i+7 == s1Len) {
            score = EXTRACT(vTbl, 0);
            match = EXTRACT(Wmatch, 0);
            length= EXTRACT(Wlength, 0);
        }
    }

    free(s1);
    free(s2);
    *_matches = match;
    *_length = length;
    return score;
}

