#if defined(T)

#include <stdint.h>

#define D CONCAT3(int, T, _t)

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
    if (cigar->len >= size) {                    \
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
    D *HT = (D*)result->trace_table;
    D *ET = (D*)result->trace_ins_table;
    D *FT = (D*)result->trace_del_table;
    cigar->seq = malloc(size);
    cigar->len = 0;
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

    cigar_reverse = parasail_reverse_uint32_t(cigar->seq, cigar->len);
    free(cigar->seq);
    cigar->seq = cigar_reverse;

    return cigar;
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

#if defined(_MSC_VER)
#define snprintf _snprintf
#endif

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

/* array index is an ASCII character value from a CIGAR,
   element value is the corresponding integer opcode between 0 and 9 */
const uint8_t parasail_cigar_encoded_ops[] = {
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0,         0,         0,         0,
    0 /*   */, 0 /* ! */, 0 /* " */, 0 /* # */,
    0 /* $ */, 0 /* % */, 0 /* & */, 0 /* ' */,
    0 /* ( */, 0 /* ) */, 0 /* * */, 0 /* + */,
    0 /* , */, 0 /* - */, 0 /* . */, 0 /* / */,
    0 /* 0 */, 0 /* 1 */, 0 /* 2 */, 0 /* 3 */,
    0 /* 4 */, 0 /* 5 */, 0 /* 6 */, 0 /* 7 */,
    0 /* 8 */, 0 /* 9 */, 0 /* : */, 0 /* ; */,
    0 /* < */, 7 /* = */, 0 /* > */, 0 /* ? */,
    0 /* @ */, 0 /* A */, 0 /* B */, 0 /* C */,
    2 /* D */, 0 /* E */, 0 /* F */, 0 /* G */,
    5 /* H */, 1 /* I */, 0 /* J */, 0 /* K */,
    0 /* L */, 0 /* M */, 3 /* N */, 0 /* O */,
    6 /* P */, 0 /* Q */, 0 /* R */, 4 /* S */,
    0 /* T */, 0 /* U */, 0 /* V */, 0 /* W */,
    8 /* X */, 0 /* Y */, 0 /* Z */, 0 /* [ */,
    0 /* \ */, 0 /* ] */, 0 /* ^ */, 0 /* _ */,
    0 /* ` */, 0 /* a */, 0 /* b */, 0 /* c */,
    0 /* d */, 0 /* e */, 0 /* f */, 0 /* g */,
    0 /* h */, 0 /* i */, 0 /* j */, 0 /* k */,
    0 /* l */, 0 /* m */, 0 /* n */, 0 /* o */,
    0 /* p */, 0 /* q */, 0 /* r */, 0 /* s */,
    0 /* t */, 0 /* u */, 0 /* v */, 0 /* w */,
    0 /* x */, 0 /* y */, 0 /* z */, 0 /* { */,
    0 /* | */, 0 /* } */, 0 /* ~ */, 0 /*  */
};

uint32_t parasail_cigar_encode(uint32_t length, char op_letter)
{
    return (length << BAM_CIGAR_SHIFT) | (parasail_cigar_encoded_ops[(int)op_letter]);
}

parasail_cigar_t* parasail_cigar_encode_string(const char *cigar)
{
    int sscanf_retcode = 0;
    size_t offset = 0;
    int chars_read = 0;
    unsigned int len = 0;
    char op = 'M';
    int done = 0;
    size_t string_length = 0;
    size_t size = 0;
    parasail_cigar_t *ret = NULL;

    string_length = strlen(cigar);
    size = sizeof(uint32_t)*string_length;
    ret = malloc(sizeof(parasail_cigar_t));
    ret->seq = malloc(size);
    ret->len = 0;

    while (!done) {
        sscanf_retcode = sscanf(
                &cigar[offset], "%u%c%n", &len, &op, &chars_read);
        if (2 != sscanf_retcode) {
            fprintf(stderr, "invalid CIGAR string\n");
            parasail_cigar_free(ret);
            return NULL;
        }
        offset += chars_read;
        ret->len += 1;
        if (ret->len >= size) {
            size *= 2;
            ret->seq = realloc(ret->seq, size);
        }
        ret->seq[ret->len-1] = parasail_cigar_encode(len, op);
        if (offset >= string_length) {
            done = 1;
        }
    }

    return ret;
}

char parasail_cigar_decode_op(uint32_t cigar_int) {
    return (cigar_int & 0xfU) > 9 ? 'M': BAM_CIGAR_STR[cigar_int & 0xfU];
}

uint32_t parasail_cigar_decode_len(uint32_t cigar_int) {
    return cigar_int >> BAM_CIGAR_SHIFT;
}

char* parasail_cigar_decode(parasail_cigar_t *cigar)
{
#define SIZE 40
    char *ret = NULL;
    size_t retlen = 0;
    size_t size = 0;
    int i = 0;

    /* initial allocation for 1 op and 3 number characters per cigar int */
    size = sizeof(char)*cigar->len*4;
    ret = malloc(size+1);
    ret[0] = '\0';

    for (i=0; i<cigar->len; ++i) {
        char tmp[SIZE];
        char op = parasail_cigar_decode_op(cigar->seq[i]);
        uint32_t len = parasail_cigar_decode_len(cigar->seq[i]);
        int snprintf_retcode = snprintf(tmp, SIZE, "%u%c", len, op);
        if (snprintf_retcode >= SIZE) {
            fprintf(stderr, "invalid CIGAR\n");
            free(ret);
            return NULL;
        }
        retlen += snprintf_retcode;
        if (retlen >= size) {
            size *= 2;
            ret = realloc(ret, size+1);
        }
        strcat(ret, tmp);
    }

    return ret;
}

void parasail_cigar_free(parasail_cigar_t *cigar)
{
	free(cigar->seq);
	free(cigar);
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

parasail_cigar_t* parasail_cigar(
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

