/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 *
 * This file was converted to C code from the raw file found at
 * ftp://ftp.cbi.pku.edu.cn/pub/software/blast/matrices/PAM120, the
 * Center for Bioinformatics, Peking University, China.
 */
#ifndef _PARASAIL_PAM120_H_
#define _PARASAIL_PAM120_H_

#include "parasail.h"
#include "pam_map.h"

#ifdef __cplusplus
extern "C" {
#endif

/* # */
/* # This matrix was produced by "pam" Version 1.0.6 [28-Jul-93] */
/* # */
/* # PAM 120 substitution matrix, scale = ln(2)/2 = 0.346574 */
/* # */
/* # Expected score = -1.64, Entropy = 0.979 bits */
/* # */
/* # Lowest score = -8, Highest score = 12 */
/* # */

static const int parasail_pam120_[] = {
/*        A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   * */
/* A */   3, -3, -1,  0, -3, -1,  0,  1, -3, -1, -3, -2, -2, -4,  1,  1,  1, -7, -4,  0,  0, -1, -1, -8,
/* R */  -3,  6, -1, -3, -4,  1, -3, -4,  1, -2, -4,  2, -1, -5, -1, -1, -2,  1, -5, -3, -2, -1, -2, -8,
/* N */  -1, -1,  4,  2, -5,  0,  1,  0,  2, -2, -4,  1, -3, -4, -2,  1,  0, -4, -2, -3,  3,  0, -1, -8,
/* D */   0, -3,  2,  5, -7,  1,  3,  0,  0, -3, -5, -1, -4, -7, -3,  0, -1, -8, -5, -3,  4,  3, -2, -8,
/* C */  -3, -4, -5, -7,  9, -7, -7, -4, -4, -3, -7, -7, -6, -6, -4,  0, -3, -8, -1, -3, -6, -7, -4, -8,
/* Q */  -1,  1,  0,  1, -7,  6,  2, -3,  3, -3, -2,  0, -1, -6,  0, -2, -2, -6, -5, -3,  0,  4, -1, -8,
/* E */   0, -3,  1,  3, -7,  2,  5, -1, -1, -3, -4, -1, -3, -7, -2, -1, -2, -8, -5, -3,  3,  4, -1, -8,
/* G */   1, -4,  0,  0, -4, -3, -1,  5, -4, -4, -5, -3, -4, -5, -2,  1, -1, -8, -6, -2,  0, -2, -2, -8,
/* H */  -3,  1,  2,  0, -4,  3, -1, -4,  7, -4, -3, -2, -4, -3, -1, -2, -3, -3, -1, -3,  1,  1, -2, -8,
/* I */  -1, -2, -2, -3, -3, -3, -3, -4, -4,  6,  1, -3,  1,  0, -3, -2,  0, -6, -2,  3, -3, -3, -1, -8,
/* L */  -3, -4, -4, -5, -7, -2, -4, -5, -3,  1,  5, -4,  3,  0, -3, -4, -3, -3, -2,  1, -4, -3, -2, -8,
/* K */  -2,  2,  1, -1, -7,  0, -1, -3, -2, -3, -4,  5,  0, -7, -2, -1, -1, -5, -5, -4,  0, -1, -2, -8,
/* M */  -2, -1, -3, -4, -6, -1, -3, -4, -4,  1,  3,  0,  8, -1, -3, -2, -1, -6, -4,  1, -4, -2, -2, -8,
/* F */  -4, -5, -4, -7, -6, -6, -7, -5, -3,  0,  0, -7, -1,  8, -5, -3, -4, -1,  4, -3, -5, -6, -3, -8,
/* P */   1, -1, -2, -3, -4,  0, -2, -2, -1, -3, -3, -2, -3, -5,  6,  1, -1, -7, -6, -2, -2, -1, -2, -8,
/* S */   1, -1,  1,  0,  0, -2, -1,  1, -2, -2, -4, -1, -2, -3,  1,  3,  2, -2, -3, -2,  0, -1, -1, -8,
/* T */   1, -2,  0, -1, -3, -2, -2, -1, -3,  0, -3, -1, -1, -4, -1,  2,  4, -6, -3,  0,  0, -2, -1, -8,
/* W */  -7,  1, -4, -8, -8, -6, -8, -8, -3, -6, -3, -5, -6, -1, -7, -2, -6, 12, -2, -8, -6, -7, -5, -8,
/* Y */  -4, -5, -2, -5, -1, -5, -5, -6, -1, -2, -2, -5, -4,  4, -6, -3, -3, -2,  8, -3, -3, -5, -3, -8,
/* V */   0, -3, -3, -3, -3, -3, -3, -2, -3,  3,  1, -4,  1, -3, -2, -2,  0, -8, -3,  5, -3, -3, -1, -8,
/* B */   0, -2,  3,  4, -6,  0,  3,  0,  1, -3, -4,  0, -4, -5, -2,  0,  0, -6, -3, -3,  4,  2, -1, -8,
/* Z */  -1, -1,  0,  3, -7,  4,  4, -2,  1, -3, -3, -1, -2, -6, -1, -1, -2, -7, -5, -3,  2,  4, -1, -8,
/* X */  -1, -2, -1, -2, -4, -1, -1, -2, -2, -1, -2, -2, -2, -3, -2, -1, -1, -5, -3, -1, -1, -1, -2, -8,
/* * */  -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1
};

static const parasail_matrix_t parasail_pam120 = {
    "pam120",
    parasail_pam120_,
    parasail_pam_map,
    24,
    12,
    -8,
    NULL,
    PARASAIL_MATRIX_TYPE_SQUARE,
    24
};

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_PAM120_H_ */

