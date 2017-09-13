/* this file is only included by traceback.c, multiple times */

#define D CONCAT3(int, T, _t)

#if defined(STRIPED)
#define NAME parasail_traceback_striped_
#define LOC LOC_STRIPED
#else
#define NAME parasail_traceback_
#define LOC LOC_NOVEC
#endif

static inline void CONCAT(NAME, T) (
        const char *seqA,
        int lena,
        const char *seqB,
        int lenb,
        const char *nameA,
        const char *nameB,
        const parasail_matrix_t *matrix,
        parasail_result_t *result,
        char match, char pos, char neg,
        int width,
        int name_width)
{
    char *q = malloc(sizeof(char)*(lena+lenb));
    char *d = malloc(sizeof(char)*(lena+lenb));
    char *a = malloc(sizeof(char)*(lena+lenb));
    char *qc = q;
    char *dc = d;
    char *ac = a;
    int i = result->end_query;
    int j = result->end_ref;
    int where = PARASAIL_DIAG;
    int c_ins = 0;
    int c_del = 0;
    D *HT = (D*)result->trace->trace_table;
    D *ET = (D*)result->trace->trace_ins_table;
    D *FT = (D*)result->trace->trace_del_table;
    int namelenA = (NULL == nameA) ? 0 : (int)strlen(nameA);
    int namelenB = (NULL == nameB) ? 0 : (int)strlen(nameB);
#if defined(STRIPED)
    int32_t segWidth = 0;
    int32_t segLen = 0;
    if (result->flag & PARASAIL_FLAG_LANES_1) {
        segWidth = 1;
    }
    if (result->flag & PARASAIL_FLAG_LANES_2) {
        segWidth = 2;
    }
    if (result->flag & PARASAIL_FLAG_LANES_4) {
        segWidth = 4;
    }
    if (result->flag & PARASAIL_FLAG_LANES_8) {
        segWidth = 8;
    }
    if (result->flag & PARASAIL_FLAG_LANES_16) {
        segWidth = 16;
    }
    if (result->flag & PARASAIL_FLAG_LANES_32) {
        segWidth = 32;
    }
    if (result->flag & PARASAIL_FLAG_LANES_64) {
        segWidth = 64;
    }
    segLen = (lena + segWidth - 1) / segWidth;
#endif
    /* semi-global alignment includes the end gaps */
    if (result->flag & PARASAIL_FLAG_SG) {
        int k;
        if (result->end_query+1 == lena) {
            k = lenb-1;
            while (k > j) {
                ++c_ins;
                *(qc++) = '-';
                *(dc++) = seqB[k];
                *(ac++) = ' ';
                --k;
            }
        }
        else if (result->end_ref+1 == lenb) {
            k = lena-1;
            while (k > i) {
                ++c_del;
                *(qc++) = seqA[k];
                *(dc++) = '-';
                *(ac++) = ' ';
                --k;
            }
        }
        else {
            assert(0);
        }
    }
    //while (i >= 0 && j >= 0) {
    while (i >= 0 || j >= 0) {
        LOC
        //assert(i >= 0 && j >= 0);
        if (i < 0) {
            while (j >= 0) {
                ++c_ins;
                *(qc++) = '-';
                *(dc++) = seqB[j];
                *(ac++) = ' ';
                --j;
            }
            break;
        }
        if (j < 0) {
            while (i >= 0) {
                ++c_del;
                *(qc++) = seqA[i];
                *(dc++) = '-';
                *(ac++) = ' ';
                --i;
            }
            break;
        }
        if (PARASAIL_DIAG == where) {
            if (HT[loc] == PARASAIL_DIAG) {
                *(qc++) = seqA[i];
                *(dc++) = seqB[j];
                *(ac++) = match_char(seqA[i], seqB[j], matrix, match, pos, neg);
                --i;
                --j;
            }
            else if (HT[loc] == PARASAIL_INS) {
                where = PARASAIL_INS;
            }
            else if (HT[loc] == PARASAIL_DEL) {
                where = PARASAIL_DEL;
            }
            else if (HT[loc] == PARASAIL_ZERO) {
                break;
            }
            else {
                assert(0);
            }
        }
        else if (PARASAIL_INS == where) {
            ++c_ins;
            *(qc++) = '-';
            *(dc++) = seqB[j];
            *(ac++) = ' ';
            --j;
            if (ET[loc] == PARASAIL_DIAG) {
                where = PARASAIL_DIAG;
            }
            else if (ET[loc] == PARASAIL_INS) {
                where = PARASAIL_INS;
            }
            else {
                assert(0);
            }
        }
        else if (PARASAIL_DEL == where) {
            ++c_del;
            *(qc++) = seqA[i];
            *(dc++) = '-';
            *(ac++) = ' ';
            --i;
            if (FT[loc] == PARASAIL_DIAG) {
                where = PARASAIL_DIAG;
            }
            else if (FT[loc] == PARASAIL_DEL) {
                where = PARASAIL_DEL;
            }
            else {
                assert(0);
            }
        }
        else if (PARASAIL_ZERO == where) {
            break;
        }
        else {
            assert(0);
        }
    }
    *(qc++) = '\0';
    *(dc++) = '\0';
    *(ac++) = '\0';

    if (1) {
        char *qr = NULL;
        char *ar = NULL;
        char *dr = NULL;
        int mch = 0;
        int sim = 0;
        int gap = 0;
        int len = strlen(a);
        int q_pindex = 0;
        int d_pindex = 0;
        int qi = 0;
        int ai = 0;
        int di = 0;

        if (result->flag & PARASAIL_FLAG_SW) {
            q_pindex = result->end_query + 1 - len + c_ins;
            d_pindex = result->end_ref + 1 - len + c_del;
        }
        else {
            q_pindex = 0;
            d_pindex = 0;
        }
        qr = parasail_reverse(q, strlen(q));
        ar = parasail_reverse(a, strlen(a));
        dr = parasail_reverse(d, strlen(d));
        for (i=0; i<len; i+=width) {
            printf("\n");
            for (j=0; j<name_width; ++j) {
                printf("%c", nameA[j]);
                if (j>=namelenA) break;
            }
            for (; j<name_width; ++j) {
                printf(" ");
            }
            printf(" %7d ", q_pindex+1);
            for (j=0; j<len&&j<width&&qi<len; ++j) {
                if (qr[qi] != '-') ++q_pindex;
                printf("%c", qr[qi]);
                ++qi;
            }
            printf(" %7d\n", q_pindex);
            for (j=0; j<name_width+1; ++j) {
                printf(" ");
            }
            printf("        ");
            for (j=0; j<len&&j<width&&ai<len; ++j) {
                if (ar[ai] == match) { ++mch; ++sim; }
                else if (ar[ai] == pos) ++sim;
                else if (ar[ai] == neg) ;
                else if (ar[ai] == ' ') ++gap;
                else {
                    printf("bad char in traceback '%c'\n", ar[ai]);
                    assert(0);
                }
                printf("%c", ar[ai]);
                ++ai;
            }
            printf("\n");
            for (j=0; j<name_width; ++j) {
                printf("%c", nameB[j]);
                if (j>=namelenB) break;
            }
            for (; j<name_width; ++j) {
                printf(" ");
            }
            printf(" %7d ", d_pindex+1);
            for (j=0; j<len&&j<width&&di<len; ++j) {
                if (dr[di] != '-') ++d_pindex;
                printf("%c", dr[di]);
                ++di;
            }
            printf(" %7d\n", d_pindex);
        }
        printf("\n");
        printf("Length: %d\n", len);
        printf("Identity:   %7d/%d (%4.1f%%)\n", mch, len, 100.0*mch/len);
        printf("Similarity: %7d/%d (%4.1f%%)\n", sim, len, 100.0*sim/len);
        printf("Gaps:       %7d/%d (%4.1f%%)\n", gap, len, 100.0*gap/len);
        printf("Score: %d\n", result->score);
        free(qr);
        free(ar);
        free(dr);
    }
    else {
        printf("%s\n", q);
        printf("%s\n", a);
        printf("%s\n", d);
    }

    free(q);
    free(d);
    free(a);
}

#undef D
#undef NAME
#undef LOC

