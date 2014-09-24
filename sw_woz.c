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

#define max8(m, vm) (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 8)); \
                    (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 4)); \
                    (vm) = _mm_max_epi16((vm), _mm_srli_si128((vm), 2)); \
                    (m) = EXTRACT((vm), 0)

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
#define FNAME sw_woz_debug
#else
#define FNAME sw_woz
#endif

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
    int * const restrict s2 = (int * const restrict)malloc(sizeof(int)*(s2Len));
    int i = 0;
    int j = 0;
    int score = NEG_INF;
    __m128i vScore = _mm_setzero_si128();

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
    }

    /* upper left corner */
    tbl_pr[7] = 0;
    del_pr[7] = NEG_INF;

    /* first row */
    for (j=8; j<s2Len+8; ++j) {
        tbl_pr[j] = 0;
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
        __m128i vZero   = _mm_set1_epi16(0);
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
        Wscore = _mm_insert_epi16(Wscore, 0, 7);
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
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 6);
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
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 5);
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
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 4);
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
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 3);
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
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 2);
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
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 1);
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
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        Wscore = _mm_insert_epi16(Wscore, 0, 0);
        vIns   = _mm_insert_epi16(vIns, NEG_INF, 0);
        tbl_pr[7] = 0;
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
            vTbl = _mm_max_epi16(vTbl, vZero);
            vScore = _mm_max_epi16(vScore, vTbl);
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
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
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
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
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
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
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
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        score_table[(i+4)*(s2Len+14+1) + j + 3] = EXTRACT(vTbl, 3);
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif

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
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        score_table[(i+5)*(s2Len+14+1) + j + 2] = EXTRACT(vTbl, 2);
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif

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
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        score_table[(i+6)*(s2Len+14+1) + j + 1] = EXTRACT(vTbl, 1);
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif

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
        vTbl = _mm_max_epi16(vTbl, vZero);
        vScore = _mm_max_epi16(vScore, vTbl);
        Wscore = vTbl;
        tbl_pr[j] = EXTRACT(vTbl, 0);
        del_pr[j] = EXTRACT(vDel, 0);
#ifdef ALIGN_EXTRA
        score_table[(i+7)*(s2Len+14+1) + j + 0] = EXTRACT(vTbl, 0);
#endif
    }

    max8(score, vScore);
    free(s1);
    free(s2);
    return score;
}

