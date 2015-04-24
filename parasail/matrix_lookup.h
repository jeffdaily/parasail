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

#include "parasail/matrices/BLOSUM100.h"
#include "parasail/matrices/BLOSUM30.h"
#include "parasail/matrices/BLOSUM35.h"
#include "parasail/matrices/BLOSUM40.h"
#include "parasail/matrices/BLOSUM45.h"
#include "parasail/matrices/BLOSUM50.h"
#include "parasail/matrices/BLOSUM55.h"
#include "parasail/matrices/BLOSUM60.h"
#include "parasail/matrices/BLOSUM62.h"
#include "parasail/matrices/BLOSUM65.h"
#include "parasail/matrices/BLOSUM70.h"
#include "parasail/matrices/BLOSUM75.h"
#include "parasail/matrices/BLOSUM80.h"
#include "parasail/matrices/BLOSUM85.h"
#include "parasail/matrices/BLOSUM90.h"
#include "parasail/matrices/PAM10.h"
#include "parasail/matrices/PAM100.h"
#include "parasail/matrices/PAM110.h"
#include "parasail/matrices/PAM120.h"
#include "parasail/matrices/PAM130.h"
#include "parasail/matrices/PAM140.h"
#include "parasail/matrices/PAM150.h"
#include "parasail/matrices/PAM160.h"
#include "parasail/matrices/PAM170.h"
#include "parasail/matrices/PAM180.h"
#include "parasail/matrices/PAM190.h"
#include "parasail/matrices/PAM20.h"
#include "parasail/matrices/PAM200.h"
#include "parasail/matrices/PAM210.h"
#include "parasail/matrices/PAM220.h"
#include "parasail/matrices/PAM230.h"
#include "parasail/matrices/PAM240.h"
#include "parasail/matrices/PAM250.h"
#include "parasail/matrices/PAM260.h"
#include "parasail/matrices/PAM270.h"
#include "parasail/matrices/PAM280.h"
#include "parasail/matrices/PAM290.h"
#include "parasail/matrices/PAM30.h"
#include "parasail/matrices/PAM300.h"
#include "parasail/matrices/PAM310.h"
#include "parasail/matrices/PAM320.h"
#include "parasail/matrices/PAM330.h"
#include "parasail/matrices/PAM340.h"
#include "parasail/matrices/PAM350.h"
#include "parasail/matrices/PAM360.h"
#include "parasail/matrices/PAM370.h"
#include "parasail/matrices/PAM380.h"
#include "parasail/matrices/PAM390.h"
#include "parasail/matrices/PAM40.h"
#include "parasail/matrices/PAM400.h"
#include "parasail/matrices/PAM410.h"
#include "parasail/matrices/PAM420.h"
#include "parasail/matrices/PAM430.h"
#include "parasail/matrices/PAM440.h"
#include "parasail/matrices/PAM450.h"
#include "parasail/matrices/PAM460.h"
#include "parasail/matrices/PAM470.h"
#include "parasail/matrices/PAM480.h"
#include "parasail/matrices/PAM490.h"
#include "parasail/matrices/PAM50.h"
#include "parasail/matrices/PAM500.h"
#include "parasail/matrices/PAM60.h"
#include "parasail/matrices/PAM70.h"
#include "parasail/matrices/PAM80.h"
#include "parasail/matrices/PAM90.h"
#include "parasail/matrices/blosum_map.h"
#include "parasail/matrices/pam_map.h"

typedef struct parasail_matrix {
    const char * name;
    const int8_t *matrix;
    const int (*matrix_)[24];
    const int *mapper;
    int size;
} parasail_matrix_t;

parasail_matrix_t parasail_matrices[] = {
    {"blosum100",parasail_blosum100,parasail_blosum100_,parasail_blosum_map,24},
    {"blosum30",parasail_blosum30,parasail_blosum30_,parasail_blosum_map,24},
    {"blosum35",parasail_blosum35,parasail_blosum35_,parasail_blosum_map,24},
    {"blosum40",parasail_blosum40,parasail_blosum40_,parasail_blosum_map,24},
    {"blosum45",parasail_blosum45,parasail_blosum45_,parasail_blosum_map,24},
    {"blosum50",parasail_blosum50,parasail_blosum50_,parasail_blosum_map,24},
    {"blosum55",parasail_blosum55,parasail_blosum55_,parasail_blosum_map,24},
    {"blosum60",parasail_blosum60,parasail_blosum60_,parasail_blosum_map,24},
    {"blosum62",parasail_blosum62,parasail_blosum62_,parasail_blosum_map,24},
    {"blosum65",parasail_blosum65,parasail_blosum65_,parasail_blosum_map,24},
    {"blosum70",parasail_blosum70,parasail_blosum70_,parasail_blosum_map,24},
    {"blosum75",parasail_blosum75,parasail_blosum75_,parasail_blosum_map,24},
    {"blosum80",parasail_blosum80,parasail_blosum80_,parasail_blosum_map,24},
    {"blosum85",parasail_blosum85,parasail_blosum85_,parasail_blosum_map,24},
    {"blosum90",parasail_blosum90,parasail_blosum90_,parasail_blosum_map,24},
    {"pam10",parasail_pam10,parasail_pam10_,parasail_pam_map,24},
    {"pam100",parasail_pam100,parasail_pam100_,parasail_pam_map,24},
    {"pam110",parasail_pam110,parasail_pam110_,parasail_pam_map,24},
    {"pam120",parasail_pam120,parasail_pam120_,parasail_pam_map,24},
    {"pam130",parasail_pam130,parasail_pam130_,parasail_pam_map,24},
    {"pam140",parasail_pam140,parasail_pam140_,parasail_pam_map,24},
    {"pam150",parasail_pam150,parasail_pam150_,parasail_pam_map,24},
    {"pam160",parasail_pam160,parasail_pam160_,parasail_pam_map,24},
    {"pam170",parasail_pam170,parasail_pam170_,parasail_pam_map,24},
    {"pam180",parasail_pam180,parasail_pam180_,parasail_pam_map,24},
    {"pam190",parasail_pam190,parasail_pam190_,parasail_pam_map,24},
    {"pam20",parasail_pam20,parasail_pam20_,parasail_pam_map,24},
    {"pam200",parasail_pam200,parasail_pam200_,parasail_pam_map,24},
    {"pam210",parasail_pam210,parasail_pam210_,parasail_pam_map,24},
    {"pam220",parasail_pam220,parasail_pam220_,parasail_pam_map,24},
    {"pam230",parasail_pam230,parasail_pam230_,parasail_pam_map,24},
    {"pam240",parasail_pam240,parasail_pam240_,parasail_pam_map,24},
    {"pam250",parasail_pam250,parasail_pam250_,parasail_pam_map,24},
    {"pam260",parasail_pam260,parasail_pam260_,parasail_pam_map,24},
    {"pam270",parasail_pam270,parasail_pam270_,parasail_pam_map,24},
    {"pam280",parasail_pam280,parasail_pam280_,parasail_pam_map,24},
    {"pam290",parasail_pam290,parasail_pam290_,parasail_pam_map,24},
    {"pam30",parasail_pam30,parasail_pam30_,parasail_pam_map,24},
    {"pam300",parasail_pam300,parasail_pam300_,parasail_pam_map,24},
    {"pam310",parasail_pam310,parasail_pam310_,parasail_pam_map,24},
    {"pam320",parasail_pam320,parasail_pam320_,parasail_pam_map,24},
    {"pam330",parasail_pam330,parasail_pam330_,parasail_pam_map,24},
    {"pam340",parasail_pam340,parasail_pam340_,parasail_pam_map,24},
    {"pam350",parasail_pam350,parasail_pam350_,parasail_pam_map,24},
    {"pam360",parasail_pam360,parasail_pam360_,parasail_pam_map,24},
    {"pam370",parasail_pam370,parasail_pam370_,parasail_pam_map,24},
    {"pam380",parasail_pam380,parasail_pam380_,parasail_pam_map,24},
    {"pam390",parasail_pam390,parasail_pam390_,parasail_pam_map,24},
    {"pam40",parasail_pam40,parasail_pam40_,parasail_pam_map,24},
    {"pam400",parasail_pam400,parasail_pam400_,parasail_pam_map,24},
    {"pam410",parasail_pam410,parasail_pam410_,parasail_pam_map,24},
    {"pam420",parasail_pam420,parasail_pam420_,parasail_pam_map,24},
    {"pam430",parasail_pam430,parasail_pam430_,parasail_pam_map,24},
    {"pam440",parasail_pam440,parasail_pam440_,parasail_pam_map,24},
    {"pam450",parasail_pam450,parasail_pam450_,parasail_pam_map,24},
    {"pam460",parasail_pam460,parasail_pam460_,parasail_pam_map,24},
    {"pam470",parasail_pam470,parasail_pam470_,parasail_pam_map,24},
    {"pam480",parasail_pam480,parasail_pam480_,parasail_pam_map,24},
    {"pam490",parasail_pam490,parasail_pam490_,parasail_pam_map,24},
    {"pam50",parasail_pam50,parasail_pam50_,parasail_pam_map,24},
    {"pam500",parasail_pam500,parasail_pam500_,parasail_pam_map,24},
    {"pam60",parasail_pam60,parasail_pam60_,parasail_pam_map,24},
    {"pam70",parasail_pam70,parasail_pam70_,parasail_pam_map,24},
    {"pam80",parasail_pam80,parasail_pam80_,parasail_pam_map,24},
    {"pam90",parasail_pam90,parasail_pam90_,parasail_pam_map,24},
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

