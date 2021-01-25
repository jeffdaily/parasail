/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 *
 * This file was converted to C code from the raw file found at
 * ftp://ftp.cbi.pku.edu.cn/pub/software/blast/matrices/BLOSUMN, the
 * Center for Bioinformatics, Peking University, China.
 */
#ifndef _PARASAIL_BLOSUMN_H_
#define _PARASAIL_BLOSUMN_H_

#include "parasail.h"
#include "blosum_map.h"

#ifdef __cplusplus
extern "C" {
#endif

/* #  Matrix made by matblas from blosumn.iij */
/* #  * column uses minimum score */
/* #  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units */
/* #  Blocks Database = /data/blocks_5.0/blocks.dat */
/* #  Cluster Percentage: >= -2 */
/* #  Entropy =   1.5172, Expected =  -1.1484 */

static const int parasail_blosumn_[] = {
/*        A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   * */
/* A */   6, -2, -2, -3, -2, -1, -2, -1, -3, -3, -3, -2, -2, -4, -1,  1, -1, -4, -4, -1, -3, -2, -1, -7,
/* R */  -2,  7, -1, -3, -6,  0, -2, -4, -1, -5, -4,  2, -3, -4, -3, -2, -2, -5, -4, -4, -2, -1, -2, -7,
/* N */  -2, -1,  7,  1, -4, -1, -1, -2,  0, -5, -5, -1, -4, -5, -4,  0, -1, -6, -4, -4,  4, -1, -2, -7,
/* D */  -3, -3,  1,  7, -6, -2,  1, -3, -2, -6, -6, -2, -5, -5, -3, -2, -2, -7, -5, -5,  4,  0, -3, -7,
/* C */  -2, -6, -4, -6,  9, -5, -7, -5, -6, -2, -3, -5, -3, -3, -5, -2, -2, -5, -4, -2, -5, -6, -4, -7,
/* Q */  -1,  0, -1, -2, -5,  7,  1, -4,  0, -4, -3,  1, -1, -4, -2, -1, -2, -4, -3, -4, -1,  4, -2, -7,
/* E */  -2, -2, -1,  1, -7,  1,  6, -4, -1, -5, -5,  0, -4, -5, -3, -1, -2, -5, -4, -4,  0,  5, -2, -7,
/* G */  -1, -4, -2, -3, -5, -4, -4,  6, -4, -6, -6, -3, -5, -5, -4, -1, -3, -5, -6, -5, -2, -4, -3, -7,
/* H */  -3, -1,  0, -2, -6,  0, -1, -4,  9, -5, -4, -2, -3, -3, -4, -2, -3, -4,  1, -5, -1, -1, -3, -7,
/* I */  -3, -5, -5, -6, -2, -4, -5, -6, -5,  6,  1, -4,  1, -1, -5, -4, -2, -4, -3,  2, -5, -5, -2, -7,
/* L */  -3, -4, -5, -6, -3, -3, -5, -6, -4,  1,  5, -4,  2,  0, -5, -4, -3, -4, -3,  0, -5, -4, -2, -7,
/* K */  -2,  2, -1, -2, -5,  1,  0, -3, -2, -4, -4,  6, -2, -4, -2, -1, -2, -6, -4, -4, -1,  0, -2, -7,
/* M */  -2, -3, -4, -5, -3, -1, -4, -5, -3,  1,  2, -2,  8, -1, -4, -3, -2, -2, -3,  0, -5, -3, -2, -7,
/* F */  -4, -4, -5, -5, -3, -4, -5, -5, -3, -1,  0, -4, -1,  7, -5, -4, -3,  0,  3, -2, -5, -5, -3, -7,
/* P */  -1, -3, -4, -3, -5, -2, -3, -4, -4, -5, -5, -2, -4, -5,  8, -2, -3, -5, -5, -4, -4, -3, -3, -7,
/* S */   1, -2,  0, -2, -2, -1, -1, -1, -2, -4, -4, -1, -3, -4, -2,  6,  1, -4, -3, -3, -1, -1, -1, -7,
/* T */  -1, -2, -1, -2, -2, -2, -2, -3, -3, -2, -3, -2, -2, -3, -3,  1,  6, -5, -3, -1, -2, -2, -1, -7,
/* W */  -4, -5, -6, -7, -5, -4, -5, -5, -4, -4, -4, -6, -2,  0, -5, -4, -5, 11,  1, -3, -6, -4, -4, -7,
/* Y */  -4, -4, -4, -5, -4, -3, -4, -6,  1, -3, -3, -4, -3,  3, -5, -3, -3,  1,  8, -3, -4, -4, -3, -7,
/* V */  -1, -4, -4, -5, -2, -4, -4, -5, -5,  2,  0, -4,  0, -2, -4, -3, -1, -3, -3,  5, -5, -4, -2, -7,
/* B */  -3, -2,  4,  4, -5, -1,  0, -2, -1, -5, -5, -1, -5, -5, -4, -1, -2, -6, -4, -5,  4,  1, -2, -7,
/* Z */  -2, -1, -1,  0, -6,  4,  5, -4, -1, -5, -4,  0, -3, -5, -3, -1, -2, -4, -4, -4,  1,  4, -2, -7,
/* X */  -1, -2, -2, -3, -4, -2, -2, -3, -3, -2, -2, -2, -2, -3, -3, -1, -1, -4, -3, -2, -2, -2, -2, -7,
/* * */  -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7,  1
};

static const parasail_matrix_t parasail_blosumn = {
    "blosumn",
    parasail_blosumn_,
    parasail_blosum_map,
    24,
    11,
    -7,
    NULL,
    PARASAIL_MATRIX_TYPE_SQUARE,
    24,
    "ARNDCQEGHILKMFPSTWYVBZX*"
};

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_BLOSUMN_H_ */

