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
        const parasail_matrix_t *matrix,
        char match,
        char pos,
        char neg)
{
    if (a == b) {
        return match;
    }
    else {
        int sub = weight(a, b, matrix);
        if (sub > 0) {
            return pos;
        }
        else {
            return neg;
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
#include "traceback_template.c"
#undef T

#define T 8
#define STRIPED
#include "traceback_template.c"
#undef T
#undef STRIPED

#define T 16
#include "traceback_template.c"
#undef T

#define T 16
#define STRIPED
#include "traceback_template.c"
#undef T
#undef STRIPED

#define T 32
#include "traceback_template.c"
#undef T

#define T 32
#define STRIPED
#include "traceback_template.c"
#undef T
#undef STRIPED

#define T 64
#include "traceback_template.c"
#undef T

#define T 64
#define STRIPED
#include "traceback_template.c"
#undef T
#undef STRIPED

void parasail_traceback_generic(
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
        int name_width,
        int use_stats)
{
    assert(parasail_result_is_trace(result));

    if (result->flag & PARASAIL_FLAG_STRIPED || result->flag & PARASAIL_FLAG_SCAN) {
        if (result->flag & PARASAIL_FLAG_BITS_8) {
            parasail_traceback_striped_8(seqA, lena, seqB, lenb, nameA, nameB, matrix, result, match, pos, neg, width, name_width, use_stats);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_16) {
            parasail_traceback_striped_16(seqA, lena, seqB, lenb, nameA, nameB, matrix, result, match, pos, neg, width, name_width, use_stats);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_32) {
            parasail_traceback_striped_32(seqA, lena, seqB, lenb, nameA, nameB, matrix, result, match, pos, neg, width, name_width, use_stats);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_64) {
            parasail_traceback_striped_64(seqA, lena, seqB, lenb, nameA, nameB, matrix, result, match, pos, neg, width, name_width, use_stats);
        }
    }
    else {
#if SIZEOF_INT == 2
        parasail_traceback_16(seqA, lena, seqB, lenb, nameA, nameB, matrix, result, match, pos, neg, width, name_width, use_stats);
#elif SIZEOF_INT == 4
        parasail_traceback_32(seqA, lena, seqB, lenb, nameA, nameB, matrix, result, match, pos, neg, width, name_width, use_stats);
#elif SIZEOF_INT == 8
        parasail_traceback_64(seqA, lena, seqB, lenb, nameA, nameB, matrix, result, match, pos, neg, width, name_width, use_stats);
#endif
    }
}

