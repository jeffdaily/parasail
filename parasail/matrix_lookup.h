/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#ifndef _PARASAIL_MATRIX_LOOKUP_H_
#define _PARASAIL_MATRIX_LOOKUP_H_

#include "parasail/matrices/blosum100.h"
#include "parasail/matrices/blosum30.h"
#include "parasail/matrices/blosum35.h"
#include "parasail/matrices/blosum40.h"
#include "parasail/matrices/blosum45.h"
#include "parasail/matrices/blosum50.h"
#include "parasail/matrices/blosum55.h"
#include "parasail/matrices/blosum60.h"
#include "parasail/matrices/blosum62.h"
#include "parasail/matrices/blosum65.h"
#include "parasail/matrices/blosum70.h"
#include "parasail/matrices/blosum75.h"
#include "parasail/matrices/blosum80.h"
#include "parasail/matrices/blosum85.h"
#include "parasail/matrices/blosum90.h"
#include "parasail/matrices/pam10.h"
#include "parasail/matrices/pam100.h"
#include "parasail/matrices/pam110.h"
#include "parasail/matrices/pam120.h"
#include "parasail/matrices/pam130.h"
#include "parasail/matrices/pam140.h"
#include "parasail/matrices/pam150.h"
#include "parasail/matrices/pam160.h"
#include "parasail/matrices/pam170.h"
#include "parasail/matrices/pam180.h"
#include "parasail/matrices/pam190.h"
#include "parasail/matrices/pam20.h"
#include "parasail/matrices/pam200.h"
#include "parasail/matrices/pam210.h"
#include "parasail/matrices/pam220.h"
#include "parasail/matrices/pam230.h"
#include "parasail/matrices/pam240.h"
#include "parasail/matrices/pam250.h"
#include "parasail/matrices/pam260.h"
#include "parasail/matrices/pam270.h"
#include "parasail/matrices/pam280.h"
#include "parasail/matrices/pam290.h"
#include "parasail/matrices/pam30.h"
#include "parasail/matrices/pam300.h"
#include "parasail/matrices/pam310.h"
#include "parasail/matrices/pam320.h"
#include "parasail/matrices/pam330.h"
#include "parasail/matrices/pam340.h"
#include "parasail/matrices/pam350.h"
#include "parasail/matrices/pam360.h"
#include "parasail/matrices/pam370.h"
#include "parasail/matrices/pam380.h"
#include "parasail/matrices/pam390.h"
#include "parasail/matrices/pam40.h"
#include "parasail/matrices/pam400.h"
#include "parasail/matrices/pam410.h"
#include "parasail/matrices/pam420.h"
#include "parasail/matrices/pam430.h"
#include "parasail/matrices/pam440.h"
#include "parasail/matrices/pam450.h"
#include "parasail/matrices/pam460.h"
#include "parasail/matrices/pam470.h"
#include "parasail/matrices/pam480.h"
#include "parasail/matrices/pam490.h"
#include "parasail/matrices/pam50.h"
#include "parasail/matrices/pam500.h"
#include "parasail/matrices/pam60.h"
#include "parasail/matrices/pam70.h"
#include "parasail/matrices/pam80.h"
#include "parasail/matrices/pam90.h"
#include "parasail/matrices/blosum_map.h"
#include "parasail/matrices/pam_map.h"

parasail_matrix_t parasail_matrices[] = {
    {PARASAIL_MATRIX_BLOSUM100},
    {PARASAIL_MATRIX_BLOSUM30},
    {PARASAIL_MATRIX_BLOSUM35},
    {PARASAIL_MATRIX_BLOSUM40},
    {PARASAIL_MATRIX_BLOSUM45},
    {PARASAIL_MATRIX_BLOSUM50},
    {PARASAIL_MATRIX_BLOSUM55},
    {PARASAIL_MATRIX_BLOSUM60},
    {PARASAIL_MATRIX_BLOSUM62},
    {PARASAIL_MATRIX_BLOSUM65},
    {PARASAIL_MATRIX_BLOSUM70},
    {PARASAIL_MATRIX_BLOSUM75},
    {PARASAIL_MATRIX_BLOSUM80},
    {PARASAIL_MATRIX_BLOSUM85},
    {PARASAIL_MATRIX_BLOSUM90},
    {PARASAIL_MATRIX_PAM10},
    {PARASAIL_MATRIX_PAM100},
    {PARASAIL_MATRIX_PAM110},
    {PARASAIL_MATRIX_PAM120},
    {PARASAIL_MATRIX_PAM130},
    {PARASAIL_MATRIX_PAM140},
    {PARASAIL_MATRIX_PAM150},
    {PARASAIL_MATRIX_PAM160},
    {PARASAIL_MATRIX_PAM170},
    {PARASAIL_MATRIX_PAM180},
    {PARASAIL_MATRIX_PAM190},
    {PARASAIL_MATRIX_PAM20},
    {PARASAIL_MATRIX_PAM200},
    {PARASAIL_MATRIX_PAM210},
    {PARASAIL_MATRIX_PAM220},
    {PARASAIL_MATRIX_PAM230},
    {PARASAIL_MATRIX_PAM240},
    {PARASAIL_MATRIX_PAM250},
    {PARASAIL_MATRIX_PAM260},
    {PARASAIL_MATRIX_PAM270},
    {PARASAIL_MATRIX_PAM280},
    {PARASAIL_MATRIX_PAM290},
    {PARASAIL_MATRIX_PAM30},
    {PARASAIL_MATRIX_PAM300},
    {PARASAIL_MATRIX_PAM310},
    {PARASAIL_MATRIX_PAM320},
    {PARASAIL_MATRIX_PAM330},
    {PARASAIL_MATRIX_PAM340},
    {PARASAIL_MATRIX_PAM350},
    {PARASAIL_MATRIX_PAM360},
    {PARASAIL_MATRIX_PAM370},
    {PARASAIL_MATRIX_PAM380},
    {PARASAIL_MATRIX_PAM390},
    {PARASAIL_MATRIX_PAM40},
    {PARASAIL_MATRIX_PAM400},
    {PARASAIL_MATRIX_PAM410},
    {PARASAIL_MATRIX_PAM420},
    {PARASAIL_MATRIX_PAM430},
    {PARASAIL_MATRIX_PAM440},
    {PARASAIL_MATRIX_PAM450},
    {PARASAIL_MATRIX_PAM460},
    {PARASAIL_MATRIX_PAM470},
    {PARASAIL_MATRIX_PAM480},
    {PARASAIL_MATRIX_PAM490},
    {PARASAIL_MATRIX_PAM50},
    {PARASAIL_MATRIX_PAM500},
    {PARASAIL_MATRIX_PAM60},
    {PARASAIL_MATRIX_PAM70},
    {PARASAIL_MATRIX_PAM80},
    {PARASAIL_MATRIX_PAM90},
    {"NULL",NULL,NULL,NULL,0}
};

parasail_matrix_t* parasail_matrix_lookup(const char *matrixname)
{
    parasail_matrix_t *matrix = NULL;

    if (matrixname) {
        int index = 0;
        parasail_matrix_t *current;
        current = &parasail_matrices[index++];
        while (current->matrix) {
            if (0 == strcmp(matrixname, current->name)) {
                matrix = current;
                break;
            }
            current = &parasail_matrices[index++];
        }
    }

    return matrix;
}

#endif /* _PARASAIL_MATRIX_LOOKUP_H_ */

