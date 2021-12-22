/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 *
 * This file was converted to C code from the raw file found at
 * ftp://ftp.cbi.pku.edu.cn/pub/software/blast/matrices/PAM70, the
 * Center for Bioinformatics, Peking University, China.
 */
#ifndef _PARASAIL_PAM70_H_
#define _PARASAIL_PAM70_H_

#include "parasail.h"
#include "pam_map.h"

#ifdef __cplusplus
extern "C" {
#endif

/* # */
/* # This matrix was produced by "pam" Version 1.0.6 [28-Jul-93] */
/* # */
/* # PAM 70 substitution matrix, scale = ln(2)/2 = 0.346574 */
/* # */
/* # Expected score = -2.77, Entropy = 1.60 bits */
/* # */
/* # Lowest score = -11, Highest score = 13 */
/* # */

static const int parasail_pam70_[] = {
/*        A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   * */
/* A */   5, -4, -2, -1, -4, -2, -1,  0, -4, -2, -4, -4, -3, -6,  0,  1,  1, -9, -5, -1, -1, -1, -2,-11,
/* R */  -4,  8, -3, -6, -5,  0, -5, -6,  0, -3, -6,  2, -2, -7, -2, -1, -4,  0, -7, -5, -4, -2, -3,-11,
/* N */  -2, -3,  6,  3, -7, -1,  0, -1,  1, -3, -5,  0, -5, -6, -3,  1,  0, -6, -3, -5,  5, -1, -2,-11,
/* D */  -1, -6,  3,  6, -9,  0,  3, -1, -1, -5, -8, -2, -7,-10, -4, -1, -2,-10, -7, -5,  5,  2, -3,-11,
/* C */  -4, -5, -7, -9,  9, -9, -9, -6, -5, -4,-10, -9, -9, -8, -5, -1, -5,-11, -2, -4, -8, -9, -6,-11,
/* Q */  -2,  0, -1,  0, -9,  7,  2, -4,  2, -5, -3, -1, -2, -9, -1, -3, -3, -8, -8, -4, -1,  5, -2,-11,
/* E */  -1, -5,  0,  3, -9,  2,  6, -2, -2, -4, -6, -2, -4, -9, -3, -2, -3,-11, -6, -4,  2,  5, -3,-11,
/* G */   0, -6, -1, -1, -6, -4, -2,  6, -6, -6, -7, -5, -6, -7, -3,  0, -3,-10, -9, -3, -1, -3, -3,-11,
/* H */  -4,  0,  1, -1, -5,  2, -2, -6,  8, -6, -4, -3, -6, -4, -2, -3, -4, -5, -1, -4,  0,  1, -3,-11,
/* I */  -2, -3, -3, -5, -4, -5, -4, -6, -6,  7,  1, -4,  1,  0, -5, -4, -1, -9, -4,  3, -4, -4, -3,-11,
/* L */  -4, -6, -5, -8,-10, -3, -6, -7, -4,  1,  6, -5,  2, -1, -5, -6, -4, -4, -4,  0, -6, -4, -4,-11,
/* K */  -4,  2,  0, -2, -9, -1, -2, -5, -3, -4, -5,  6,  0, -9, -4, -2, -1, -7, -7, -6, -1, -2, -3,-11,
/* M */  -3, -2, -5, -7, -9, -2, -4, -6, -6,  1,  2,  0, 10, -2, -5, -3, -2, -8, -7,  0, -6, -3, -3,-11,
/* F */  -6, -7, -6,-10, -8, -9, -9, -7, -4,  0, -1, -9, -2,  8, -7, -4, -6, -2,  4, -5, -7, -9, -5,-11,
/* P */   0, -2, -3, -4, -5, -1, -3, -3, -2, -5, -5, -4, -5, -7,  7,  0, -2, -9, -9, -3, -4, -2, -3,-11,
/* S */   1, -1,  1, -1, -1, -3, -2,  0, -3, -4, -6, -2, -3, -4,  0,  5,  2, -3, -5, -3,  0, -2, -1,-11,
/* T */   1, -4,  0, -2, -5, -3, -3, -3, -4, -1, -4, -1, -2, -6, -2,  2,  6, -8, -4, -1, -1, -3, -2,-11,
/* W */  -9,  0, -6,-10,-11, -8,-11,-10, -5, -9, -4, -7, -8, -2, -9, -3, -8, 13, -3,-10, -7,-10, -7,-11,
/* Y */  -5, -7, -3, -7, -2, -8, -6, -9, -1, -4, -4, -7, -7,  4, -9, -5, -4, -3,  9, -5, -4, -7, -5,-11,
/* V */  -1, -5, -5, -5, -4, -4, -4, -3, -4,  3,  0, -6,  0, -5, -3, -3, -1,-10, -5,  6, -5, -4, -2,-11,
/* B */  -1, -4,  5,  5, -8, -1,  2, -1,  0, -4, -6, -1, -6, -7, -4,  0, -1, -7, -4, -5,  5,  1, -2,-11,
/* Z */  -1, -2, -1,  2, -9,  5,  5, -3,  1, -4, -4, -2, -3, -9, -2, -2, -3,-10, -7, -4,  1,  5, -3,-11,
/* X */  -2, -3, -2, -3, -6, -2, -3, -3, -3, -3, -4, -3, -3, -5, -3, -1, -2, -7, -5, -2, -2, -3, -3,-11,
/* * */ -11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,  1
};

static const parasail_matrix_t parasail_pam70 = {
    "pam70",
    parasail_pam70_,
    parasail_pam_map,
    24,
    13,
    -11,
    NULL,
    PARASAIL_MATRIX_TYPE_SQUARE,
    24,
    "ARNDCQEGHILKMFPSTWYVBZX*",
    NULL
};

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_PAM70_H_ */

