/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#ifndef _PARASAIL_BLOSUM_LOOKUP_H_
#define _PARASAIL_BLOSUM_LOOKUP_H_

typedef const int (*parasail_blosum_t)[24];

typedef struct blosum {
    const char * name;
    parasail_blosum_t pointer;
} blosum_t;

blosum_t blosums[] = {
    {"blosum40",blosum40},
    {"blosum45",blosum45},
    {"blosum50",blosum50},
    {"blosum62",blosum62},
    {"blosum75",blosum75},
    {"blosum80",blosum80},
    {"blosum90",blosum90},
    {"NULL",NULL},
};

parasail_blosum_t lookup_blosum(const char *blosumname)
{
    parasail_blosum_t blosum = NULL;

    if (blosumname) {
        int index = 0;
        blosum_t b;
        b = blosums[index++];
        while (b.pointer) {
            if (0 == strcmp(blosumname, b.name)) {
                blosum = b.pointer;
                break;
            }
            b = blosums[index++];
        }
    }

    return blosum;
}

#endif /* _PARASAIL_BLOSUM_LOOKUP_H_ */
