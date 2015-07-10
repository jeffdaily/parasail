/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_IO_H_
#define _PARASAIL_IO_H_

#include <sys/stat.h>

typedef struct parasail_file {
    int fd;
    off_t size;
    char *buf;
} parasail_file_t;

typedef struct parasail_file_stat {
    unsigned long sequences;
    unsigned long characters;
    unsigned long shortest;
    unsigned long longest;
    float mean;
    float stddev;
} parasail_file_stat_t;

parasail_file_t* parasail_open(const char *fname);

/** Closes file and frees file parameter. */
void parasail_close(parasail_file_t *file);

int parasail_is_fasta(const parasail_file_t *pf);

int parasail_is_fastq(const parasail_file_t *pf);

parasail_file_stat_t* parasail_stat(const parasail_file_t *pf);

parasail_file_stat_t* parasail_stat_fasta(const parasail_file_t *pf);

parasail_file_stat_t* parasail_stat_fastq(const parasail_file_t *pf);

#endif /* _PARASAIL_IO_H_ */
