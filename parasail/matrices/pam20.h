/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 *
 * This file was converted to C code from the raw file found at
 * ftp://ftp.cbi.pku.edu.cn/pub/software/blast/matrices/PAM20, the
 * Center for Bioinformatics, Peking University, China.
 */
#ifndef _PARASAIL_PAM20_H_
#define _PARASAIL_PAM20_H_

#include "parasail.h"
#include "pam_map.h"

#ifdef __cplusplus
extern "C" {
#endif

/* # */
/* # This matrix was produced by "pam" Version 1.0.6 [28-Jul-93] */
/* # */
/* # PAM 20 substitution matrix, scale = ln(2)/2 = 0.346574 */
/* # */
/* # Expected score = -6.18, Entropy = 2.95 bits */
/* # */
/* # Lowest score = -19, Highest score = 13 */
/* # */

static const int parasail_pam20_[] = {
/*        A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   * */
/* A */   6, -8, -5, -4, -8, -5, -3, -3, -8, -6, -7, -8, -6, -9, -2, -1, -1,-16, -9, -3, -5, -4, -4,-19,
/* R */  -8,  9, -7,-12, -9, -2,-11,-11, -3, -6,-10, -1, -5,-10, -5, -4, -8, -3,-11, -9, -9, -5, -7,-19,
/* N */  -5, -7,  8,  1,-13, -5, -3, -4, -1, -6, -8, -2,-11,-10, -7, -1, -3, -9, -5, -9,  6, -4, -4,-19,
/* D */  -4,-12,  1,  8,-16, -4,  2, -4, -5, -9,-15, -6,-13,-17, -9, -5, -6,-17,-13, -9,  6,  0, -7,-19,
/* C */  -8, -9,-13,-16, 10,-16,-16,-11, -8, -7,-17,-16,-16,-15, -9, -4, -9,-18, -5, -7,-14,-16,-11,-19,
/* Q */  -5, -2, -5, -4,-16,  9,  0, -8,  0, -9, -6, -4, -5,-15, -4, -6, -7,-15,-14, -8, -4,  7, -6,-19,
/* E */  -3,-11, -3,  2,-16,  0,  8, -5, -6, -6,-10, -5, -8,-16, -7, -5, -7,-19, -9, -8,  0,  6, -6,-19,
/* G */  -3,-11, -4, -4,-11, -8, -5,  7,-10,-13,-12, -8,-10,-10, -7, -3, -7,-17,-16, -7, -4, -6, -6,-19,
/* H */  -8, -3, -1, -5, -8,  0, -6,-10,  9,-11, -7, -8,-13, -7, -5, -7, -8, -8, -4, -7, -2, -2, -6,-19,
/* I */  -6, -6, -6, -9, -7, -9, -6,-13,-11,  9, -2, -7, -2, -3,-10, -8, -3,-16, -7,  1, -7, -7, -6,-19,
/* L */  -7,-10, -8,-15,-17, -6,-10,-12, -7, -2,  7, -9,  0, -4, -8, -9, -8, -7, -8, -3,-10, -8, -7,-19,
/* K */  -8, -1, -2, -6,-16, -4, -5, -8, -8, -7, -9,  7, -3,-16, -8, -5, -4,-14,-10,-10, -3, -5, -6,-19,
/* M */  -6, -5,-11,-13,-16, -5, -8,-10,-13, -2,  0, -3, 11, -5, -9, -6, -5,-15,-13, -2,-12, -6, -6,-19,
/* F */  -9,-10,-10,-17,-15,-15,-16,-10, -7, -3, -4,-16, -5,  9,-11, -7,-10, -6,  1, -9,-12,-16, -9,-19,
/* P */  -2, -5, -7, -9, -9, -4, -7, -7, -5,-10, -8, -8, -9,-11,  8, -3, -5,-16,-16, -7, -8, -5, -6,-19,
/* S */  -1, -4, -1, -5, -4, -6, -5, -3, -7, -8, -9, -5, -6, -7, -3,  7,  0, -6, -8, -8, -2, -6, -4,-19,
/* T */  -1, -8, -3, -6, -9, -7, -7, -7, -8, -3, -8, -4, -5,-10, -5,  0,  7,-15, -7, -4, -4, -7, -5,-19,
/* W */ -16, -3, -9,-17,-18,-15,-19,-17, -8,-16, -7,-14,-15, -6,-16, -6,-15, 13, -6,-18,-11,-17,-13,-19,
/* Y */  -9,-11, -5,-13, -5,-14, -9,-16, -4, -7, -8,-10,-13,  1,-16, -8, -7, -6, 10, -8, -7,-11, -9,-19,
/* V */  -3, -9, -9, -9, -7, -8, -8, -7, -7,  1, -3,-10, -2, -9, -7, -8, -4,-18, -8,  7, -9, -8, -6,-19,
/* B */  -5, -9,  6,  6,-14, -4,  0, -4, -2, -7,-10, -3,-12,-12, -8, -2, -4,-11, -7, -9,  6, -1, -6,-19,
/* Z */  -4, -5, -4,  0,-16,  7,  6, -6, -2, -7, -8, -5, -6,-16, -5, -6, -7,-17,-11, -8, -1,  6, -6,-19,
/* X */  -4, -7, -4, -7,-11, -6, -6, -6, -6, -6, -7, -6, -6, -9, -6, -4, -5,-13, -9, -6, -6, -6, -6,-19,
/* * */ -19,-19,-19,-19,-19,-19,-19,-19,-19,-19,-19,-19,-19,-19,-19,-19,-19,-19,-19,-19,-19,-19,-19,  1
};

static const parasail_matrix_t parasail_pam20 = {
    "pam20",
    parasail_pam20_,
    parasail_pam_map,
    24,
    13,
    -19,
    NULL,
    PARASAIL_MATRIX_TYPE_SQUARE,
    24
};

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_PAM20_H_ */

