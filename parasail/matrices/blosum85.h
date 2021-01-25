/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 *
 * This file was converted to C code from the raw file found at
 * ftp://ftp.cbi.pku.edu.cn/pub/software/blast/matrices/BLOSUM85, the
 * Center for Bioinformatics, Peking University, China.
 */
#ifndef _PARASAIL_BLOSUM85_H_
#define _PARASAIL_BLOSUM85_H_

#include "parasail.h"
#include "blosum_map.h"

#ifdef __cplusplus
extern "C" {
#endif

/* #  Matrix made by matblas from blosum85.iij */
/* #  * column uses minimum score */
/* #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units */
/* #  Blocks Database = /data/blocks_5.0/blocks.dat */
/* #  Cluster Percentage: >= 85 */
/* #  Entropy =   1.0805, Expected =  -0.8153 */

static const int parasail_blosum85_[] = {
/*        A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   * */
/* A */   5, -2, -2, -2, -1, -1, -1,  0, -2, -2, -2, -1, -2, -3, -1,  1,  0, -3, -3, -1, -2, -1, -1, -6,
/* R */  -2,  6, -1, -2, -4,  1, -1, -3,  0, -4, -3,  2, -2, -4, -2, -1, -2, -4, -3, -3, -2,  0, -2, -6,
/* N */  -2, -1,  7,  1, -4,  0, -1, -1,  0, -4, -4,  0, -3, -4, -3,  0,  0, -5, -3, -4,  4, -1, -2, -6,
/* D */  -2, -2,  1,  7, -5, -1,  1, -2, -2, -5, -5, -1, -4, -4, -2, -1, -2, -6, -4, -4,  4,  1, -2, -6,
/* C */  -1, -4, -4, -5,  9, -4, -5, -4, -5, -2, -2, -4, -2, -3, -4, -2, -2, -4, -3, -1, -4, -5, -3, -6,
/* Q */  -1,  1,  0, -1, -4,  6,  2, -3,  1, -4, -3,  1,  0, -4, -2, -1, -1, -3, -2, -3, -1,  4, -1, -6,
/* E */  -1, -1, -1,  1, -5,  2,  6, -3, -1, -4, -4,  0, -3, -4, -2, -1, -1, -4, -4, -3,  0,  4, -1, -6,
/* G */   0, -3, -1, -2, -4, -3, -3,  6, -3, -5, -5, -2, -4, -4, -3, -1, -2, -4, -5, -4, -1, -3, -2, -6,
/* H */  -2,  0,  0, -2, -5,  1, -1, -3,  8, -4, -3, -1, -3, -2, -3, -1, -2, -3,  2, -4, -1,  0, -2, -6,
/* I */  -2, -4, -4, -5, -2, -4, -4, -5, -4,  5,  1, -3,  1, -1, -4, -3, -1, -3, -2,  3, -5, -4, -2, -6,
/* L */  -2, -3, -4, -5, -2, -3, -4, -5, -3,  1,  4, -3,  2,  0, -4, -3, -2, -3, -2,  0, -5, -4, -2, -6,
/* K */  -1,  2,  0, -1, -4,  1,  0, -2, -1, -3, -3,  6, -2, -4, -2, -1, -1, -5, -3, -3, -1,  1, -1, -6,
/* M */  -2, -2, -3, -4, -2,  0, -3, -4, -3,  1,  2, -2,  7, -1, -3, -2, -1, -2, -2,  0, -4, -2, -1, -6,
/* F */  -3, -4, -4, -4, -3, -4, -4, -4, -2, -1,  0, -4, -1,  7, -4, -3, -3,  0,  3, -1, -4, -4, -2, -6,
/* P */  -1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -4, -2, -3, -4,  8, -1, -2, -5, -4, -3, -3, -2, -2, -6,
/* S */   1, -1,  0, -1, -2, -1, -1, -1, -1, -3, -3, -1, -2, -3, -1,  5,  1, -4, -2, -2,  0, -1, -1, -6,
/* T */   0, -2,  0, -2, -2, -1, -1, -2, -2, -1, -2, -1, -1, -3, -2,  1,  5, -4, -2,  0, -1, -1, -1, -6,
/* W */  -3, -4, -5, -6, -4, -3, -4, -4, -3, -3, -3, -5, -2,  0, -5, -4, -4, 11,  2, -3, -5, -4, -3, -6,
/* Y */  -3, -3, -3, -4, -3, -2, -4, -5,  2, -2, -2, -3, -2,  3, -4, -2, -2,  2,  7, -2, -4, -3, -2, -6,
/* V */  -1, -3, -4, -4, -1, -3, -3, -4, -4,  3,  0, -3,  0, -1, -3, -2,  0, -3, -2,  5, -4, -3, -1, -6,
/* B */  -2, -2,  4,  4, -4, -1,  0, -1, -1, -5, -5, -1, -4, -4, -3,  0, -1, -5, -4, -4,  4,  0, -2, -6,
/* Z */  -1,  0, -1,  1, -5,  4,  4, -3,  0, -4, -4,  1, -2, -4, -2, -1, -1, -4, -3, -3,  0,  4, -1, -6,
/* X */  -1, -2, -2, -2, -3, -1, -1, -2, -2, -2, -2, -1, -1, -2, -2, -1, -1, -3, -2, -1, -2, -1, -2, -6,
/* * */  -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,  1
};

static const parasail_matrix_t parasail_blosum85 = {
    "blosum85",
    parasail_blosum85_,
    parasail_blosum_map,
    24,
    11,
    -6,
    NULL,
    PARASAIL_MATRIX_TYPE_SQUARE,
    24,
    "ARNDCQEGHILKMFPSTWYVBZX*"
};

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_BLOSUM85_H_ */

