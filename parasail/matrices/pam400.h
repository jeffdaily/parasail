/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 *
 * This file was converted to C code from the raw file found at
 * ftp://ftp.cbi.pku.edu.cn/pub/software/blast/matrices/PAM400, the
 * Center for Bioinformatics, Peking University, China.
 */
#ifndef _PARASAIL_PAM400_H_
#define _PARASAIL_PAM400_H_

#include "parasail.h"
#include "pam_map.h"

#ifdef __cplusplus
extern "C" {
#endif

/* # */
/* # This matrix was produced by "pam" Version 1.0.6 [28-Jul-93] */
/* # */
/* # PAM 400 substitution matrix, scale = ln(2)/5 = 0.138629 */
/* # */
/* # Expected score = -0.521, Entropy = 0.139 bits */
/* # */
/* # Lowest score = -8, Highest score = 26 */
/* # */

static const int parasail_pam400_[] = {
/*        A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   * */
/* A */   1, -1,  0,  1, -2,  0,  1,  1, -1,  0, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0,  1,  0,  0, -8,
/* R */  -1,  5,  1,  0, -4,  1,  0, -2,  2, -2, -3,  4,  0, -4,  0,  0,  0,  3, -4, -2,  0,  1,  0, -8,
/* N */   0,  1,  1,  2, -3,  1,  2,  1,  1, -1, -3,  1, -1, -3,  0,  1,  0, -4, -3, -1,  2,  1,  0, -8,
/* D */   1,  0,  2,  3, -5,  2,  3,  1,  1, -2, -3,  1, -2, -5,  0,  1,  0, -7, -4, -2,  2,  2,  0, -8,
/* C */  -2, -4, -3, -5, 17, -5, -5, -3, -3, -2, -6, -5, -5, -4, -2,  0, -2, -8,  1, -2, -4, -5, -3, -8,
/* Q */   0,  1,  1,  2, -5,  3,  2,  0,  3, -2, -2,  1, -1, -4,  0,  0,  0, -5, -4, -1,  2,  2,  0, -8,
/* E */   1,  0,  2,  3, -5,  2,  3,  1,  1, -2, -3,  1, -2, -5,  0,  0,  0, -7, -4, -1,  2,  3,  0, -8,
/* G */   1, -2,  1,  1, -3,  0,  1,  4, -1, -2, -3, -1, -2, -5,  0,  1,  1, -7, -5, -1,  1,  0,  0, -8,
/* H */  -1,  2,  1,  1, -3,  3,  1, -1,  5, -2, -2,  1, -1, -2,  0,  0, -1, -2,  0, -2,  1,  2,  0, -8,
/* I */   0, -2, -1, -2, -2, -2, -2, -2, -2,  4,  3, -2,  2,  2, -1, -1,  0, -5,  0,  3, -2, -2,  0, -8,
/* L */  -2, -3, -3, -3, -6, -2, -3, -3, -2,  3,  7, -2,  4,  3, -2, -2, -1, -2,  0,  3, -3, -2, -1, -8,
/* K */  -1,  4,  1,  1, -5,  1,  1, -1,  1, -2, -2,  4,  0, -5,  0,  0,  0, -3, -4, -2,  1,  1,  0, -8,
/* M */  -1,  0, -1, -2, -5, -1, -2, -2, -1,  2,  4,  0,  5,  1, -1, -1,  0, -4, -1,  2, -2, -1,  0, -8,
/* F */  -3, -4, -3, -5, -4, -4, -5, -5, -2,  2,  3, -5,  1, 11, -4, -3, -3,  2, 10,  0, -4, -5, -2, -8,
/* P */   1,  0,  0,  0, -2,  0,  0,  0,  0, -1, -2,  0, -1, -4,  5,  1,  1, -6, -5, -1,  0,  0,  0, -8,
/* S */   1,  0,  1,  1,  0,  0,  0,  1,  0, -1, -2,  0, -1, -3,  1,  1,  1, -3, -3, -1,  1,  0,  0, -8,
/* T */   1,  0,  0,  0, -2,  0,  0,  1, -1,  0, -1,  0,  0, -3,  1,  1,  1, -5, -3,  0,  0,  0,  0, -8,
/* W */  -6,  3, -4, -7, -8, -5, -7, -7, -2, -5, -2, -3, -4,  2, -6, -3, -5, 26,  1, -6, -5, -6, -4, -8,
/* Y */  -3, -4, -3, -4,  1, -4, -4, -5,  0,  0,  0, -4, -1, 10, -5, -3, -3,  1, 13, -2, -3, -4, -2, -8,
/* V */   0, -2, -1, -2, -2, -1, -1, -1, -2,  3,  3, -2,  2,  0, -1, -1,  0, -6, -2,  4, -1, -1,  0, -8,
/* B */   1,  0,  2,  2, -4,  2,  2,  1,  1, -2, -3,  1, -2, -4,  0,  1,  0, -5, -3, -1,  2,  2,  0, -8,
/* Z */   0,  1,  1,  2, -5,  2,  3,  0,  2, -2, -2,  1, -1, -5,  0,  0,  0, -6, -4, -1,  2,  3,  0, -8,
/* X */   0,  0,  0,  0, -3,  0,  0,  0,  0,  0, -1,  0,  0, -2,  0,  0,  0, -4, -2,  0,  0,  0, -1, -8,
/* * */  -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1
};

static const parasail_matrix_t parasail_pam400 = {
    "pam400",
    parasail_pam400_,
    parasail_pam_map,
    24,
    26,
    -8,
    NULL,
    PARASAIL_MATRIX_TYPE_SQUARE,
    24,
    "ARNDCQEGHILKMFPSTWYVBZX*"
};

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_PAM400_H_ */

