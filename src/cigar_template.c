/* this file is only included by cigar.c, multiple times */

#define D CONCAT3(int, T, _t)

#define UNUSED(expr) do { (void)(expr); } while (0)

#if defined(STRIPED)
#define NAME parasail_cigar_striped_
#define LOC LOC_STRIPED
#else
#define NAME parasail_cigar_
#define LOC LOC_NOVEC
#endif

#define INC                                      \
do {                                             \
    cigar->len += 1;                             \
    if ((size_t)cigar->len >= size) {            \
        size = size * 2;                         \
        cigar->seq = realloc(cigar->seq, size);  \
    }                                            \
} while (0);

#define RESET  \
do {           \
    c_mat = 0; \
    c_mis = 0; \
    c_del = 0; \
    c_ins = 0; \
} while (0)

#define WRITE(VAL,CHAR)                                         \
do {                                                            \
    INC;                                                        \
    cigar->seq[cigar->len-1] = parasail_cigar_encode(VAL,CHAR); \
} while (0)

#define WRITE_ANY         \
do {                      \
    if (c_mat) {          \
        WRITE(c_mat,'='); \
    }                     \
    else if (c_mis) {     \
        WRITE(c_mis,'X'); \
    }                     \
    else if (c_del) {     \
        WRITE(c_del,'D'); \
    }                     \
    else if (c_ins) {     \
        WRITE(c_ins,'I'); \
    }                     \
    RESET;                \
} while (0)

static inline parasail_cigar_t* CONCAT(NAME, T) (
        const char *seqA,
        int lena,
        const char *seqB,
        int lenb,
        const parasail_matrix_t *matrix,
        parasail_result_t *result)
{
    size_t size = sizeof(char)*(lena+lenb);
    parasail_cigar_t *cigar = malloc(sizeof(parasail_cigar_t));
    uint32_t *cigar_reverse = NULL;
    uint32_t c_mat = 0;
    uint32_t c_mis = 0;
    uint32_t c_del = 0;
    uint32_t c_ins = 0;
    int i = result->end_query;
    int j = result->end_ref;
    int where = PARASAIL_DIAG;
    D *HT = (D*)result->trace->trace_table;
    D *ET = (D*)result->trace->trace_ins_table;
    D *FT = (D*)result->trace->trace_del_table;
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
    cigar->seq = malloc(size);
    cigar->len = 0;
    cigar->beg_query = 0;
    cigar->beg_ref = 0;
    UNUSED(matrix);
    /* semi-global alignment includes the end gaps */
    if (result->flag & PARASAIL_FLAG_SG) {
        int k;
        if (result->end_query+1 == lena) {
            k = lenb-1;
            while (k > j) {
                ++c_ins;
                --k;
            }
        }
        else if (result->end_ref+1 == lenb) {
            k = lena-1;
            while (k > i) {
                ++c_del;
                --k;
            }
        }
        else {
            assert(0);
        }
    }
    while (i >= 0 || j >= 0) {
        LOC
        /*assert(i >= 0 && j >= 0);*/
        if (i < 0) {
            if (0 == c_ins) {
                WRITE_ANY;
            }
            while (j >= 0) {
                ++c_ins;
                --j;
            }
            break;
        }
        if (j < 0) {
            if (0 == c_del) {
                WRITE_ANY;
            }
            while (i >= 0) {
                ++c_del;
                --i;
            }
            break;
        }
        if (PARASAIL_DIAG == where) {
            if (HT[loc] == PARASAIL_DIAG) {
                if (seqA[i] == seqB[j]) {
                    if (0 == c_mat) {
                        WRITE_ANY;
                    }
                    c_mat += 1;
                }
                else {
                    if (0 == c_mis) {
                        WRITE_ANY;
                    }
                    c_mis += 1;
                }
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
            if (0 == c_ins) {
                WRITE_ANY;
            }
            c_ins += 1;
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
            if (0 == c_del) {
                WRITE_ANY;
            }
            c_del += 1;
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

    /* in case we missed the last write */
    WRITE_ANY;

    cigar_reverse = parasail_reverse_uint32_t(cigar->seq, cigar->len);
    free(cigar->seq);
    cigar->seq = cigar_reverse;
    cigar->beg_query = i+1;
    cigar->beg_ref = j+1;

    return cigar;
}

#undef D
#undef NAME
#undef LOC

