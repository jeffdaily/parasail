#if defined(T)

#define D CONCAT3(int, T, _t)

#if defined(STRIPED)
#define NAME parasail_cigar_striped_
#define LOC LOC_STRIPED
#else
#define NAME parasail_cigar_
#define LOC LOC_NOVEC
#endif

#define INC                            \
do {                                   \
    c_len += 1;                        \
    if (c_len >= size) {               \
        size = size * 2;               \
        cigar = realloc(cigar, size);  \
    }                                  \
} while (0);

#define RESET \
do {          \
c_mat = 0;    \
c_mis = 0;    \
c_del = 0;    \
c_ins = 0;    \
} while (0)

#define WRITE(VAL,CHAR)                                    \
do {                                                       \
    int c = sprintf(tmp, "%llu", (unsigned long long)VAL); \
    INC;                                                   \
    cigar[c_len-1] = CHAR;                                 \
    for (; c > 0; --c) {                                   \
        INC;                                               \
        cigar[c_len-1] = tmp[c-1];                         \
    }                                                      \
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

static inline char* CONCAT(NAME, T) (
        const char *seqA,
        int lena,
        const char *seqB,
        int lenb,
        const parasail_matrix_t *matrix,
        parasail_result_t *result)
{
    size_t size = sizeof(char)*(lena+lenb);
    char *cigar = malloc(size);
    char *cigar_reverse = NULL;
    size_t c_len = 0;
    size_t c_mat = 0;
    size_t c_mis = 0;
    size_t c_del = 0;
    size_t c_ins = 0;
    int i = result->end_query;
    int j = result->end_ref;
    int where = PARASAIL_DIAG;
    D *HT = (D*)result->trace_table;
    D *ET = (D*)result->trace_ins_table;
    D *FT = (D*)result->trace_del_table;
    char tmp[20];
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
    while (i >= 0 && j >= 0) {
        LOC
        assert(i >= 0 && j >= 0);
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

    cigar_reverse = parasail_reverse(cigar, c_len);
    free(cigar);

    return cigar_reverse;
}

#undef D
#undef NAME
#undef LOC

#else

#include "config.h"

#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parasail.h"
#include "parasail/memory.h"

static inline int weight(
        const char a,
        const char b,
        const parasail_matrix_t *matrix)
{
    return matrix->matrix[matrix->mapper[(unsigned char)a]*matrix->size + matrix->mapper[(unsigned char)b]];
}

static inline char match_char(
        const char a,
        const char b,
        const parasail_matrix_t *matrix)
{
    if (a == b) {
        return '|';
    }
    else {
        int sub = weight(a, b, matrix);
        if (sub > 0) {
            return ':';
        }
        else {
            return '.';
        }
    }

    return 'X'; /* shouldn't happen */
}

#define CONCAT_(X, Y) X##Y
#define CONCAT(X, Y) CONCAT_(X, Y)
#define CONCAT3_(X, Y, Z) X##Y##Z
#define CONCAT3(X, Y, Z) CONCAT3_(X, Y, Z)
#define LOC_NOVEC int loc = i*lenb + j;
#define LOC_STRIPED int loc = j*segLen*segWidth + (i%segLen)*segWidth + (i/segLen);

#define T 8
#include "cigar.c"
#undef T

#define T 8
#define STRIPED
#include "cigar.c"
#undef T
#undef STRIPED

#define T 16
#include "cigar.c"
#undef T

#define T 16
#define STRIPED
#include "cigar.c"
#undef T
#undef STRIPED

#define T 32
#include "cigar.c"
#undef T

#define T 32
#define STRIPED
#include "cigar.c"
#undef T
#undef STRIPED

#define T 64
#include "cigar.c"
#undef T

#define T 64
#define STRIPED
#include "cigar.c"
#undef T
#undef STRIPED

char* parasail_cigar(
        const char *seqA,
        int lena,
        const char *seqB,
        int lenb,
        const parasail_matrix_t *matrix,
        parasail_result_t *result)
{
    if (result->flag & PARASAIL_FLAG_STRIPED) {
        if (result->flag & PARASAIL_FLAG_BITS_8) {
            return parasail_cigar_striped_8(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_16) {
            return parasail_cigar_striped_16(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_32) {
            return parasail_cigar_striped_32(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_64) {
            return parasail_cigar_striped_64(seqA, lena, seqB, lenb, matrix, result);
        }
    }
    else {
        if (result->flag & PARASAIL_FLAG_BITS_8) {
            return parasail_cigar_8(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_16) {
            return parasail_cigar_16(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_32) {
            return parasail_cigar_32(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_64) {
            return parasail_cigar_64(seqA, lena, seqB, lenb, matrix, result);
        }
    }
    
    /* should not get here, but to silence warnings */
    return NULL;
}

#endif

