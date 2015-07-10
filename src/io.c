/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "parasail/io.h"
#include "parasail/stats.h"

parasail_file_t* parasail_open(const char *fname)
{
    int fd = -1;
    struct stat fs;
    char *buf = NULL;
    parasail_file_t *pf = NULL;

    fd = open(fname, O_RDONLY);
    if (fd == -1) {
        fprintf(stderr, "Cannot open input file `%s': ", fname);
        perror("open");
        exit(EXIT_FAILURE);
    }

    if (-1 == fstat(fd, &fs)) {
        fprintf(stderr, "Cannont stat input file `%s': ", fname);
        perror("fstat");
        exit(EXIT_FAILURE);
    }

    buf = (char*)mmap(NULL, fs.st_size, PROT_READ, MAP_SHARED, fd, 0);
    if (MAP_FAILED == buf) {
        fprintf(stderr, "Cannont mmap input file `%s': ", fname);
        perror("mmap");
        exit(EXIT_FAILURE);
    }

    pf = (parasail_file_t*)malloc(sizeof(parasail_file_t));
    if (NULL == pf) {
        fprintf(stderr, "Cannont allocate parasail_file_t");
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    pf->fd = fd;
    pf->size = fs.st_size;
    pf->buf = buf;
    return pf;
}

void parasail_close(parasail_file_t *pf)
{
    munmap(pf->buf, pf->size);
    close(pf->fd);
    free(pf);
}

int parasail_is_fasta(const parasail_file_t *pf)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail_is_fasta given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return '>' == pf->buf[0];
}

int parasail_is_fastq(const parasail_file_t *pf)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail_is_fastq given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return '@' == pf->buf[0];
}

parasail_file_stat_t* parasail_stat(const parasail_file_t *pf)
{
    parasail_file_stat_t *stat = NULL;

    if (NULL == pf) {
        fprintf(stderr, "parasail_stat given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (parasail_is_fasta(pf)) {
        stat = parasail_stat_fasta(pf);
    }
    else if (parasail_is_fastq(pf)) {
        stat = parasail_stat_fastq(pf);
    }
    else {
        fprintf(stderr, "parasail_stat: cannot determine file format\n");
        exit(EXIT_FAILURE);
    }

    return stat;
}

parasail_file_stat_t* parasail_stat_fasta(const parasail_file_t *pf)
{
    char *T = pf->buf;
    off_t i = 0;
    unsigned long seq = 0;
    unsigned long c = 0;
    unsigned long c_tot = 0;
    stats_t stats;
    parasail_file_stat_t *pfs = NULL;

    stats_clear(&stats);

    if (NULL == pf) {
        fprintf(stderr, "parasail_stat_fasta given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    /* first line is always first sequence ID */
    if (T[i] != '>') {
        fprintf(stderr, "poorly formatted FASTA file\n");
        exit(EXIT_FAILURE);
    }

    /* read entire first line */
    while (T[i] != '\n' && T[i] != '\r') {
        ++i;
    }
    /* for the case of "\r\n"  or "\n\r" */
    if (T[i+1] == '\n' || T[i+1] == '\r') {
        ++i;
    }
    ++i;

    /* count that first sequence */
    ++seq;

    /* read rest of file */
    while (i<pf->size) {
        if (T[i] == '>') {
            /* encountered a new sequence */
            ++seq;
            stats_sample_value(&stats, c);
            c = 0;
            /* skip rest of this line */
            while (T[i] != '\n' && T[i] != '\r') {
                ++i;
            }
            /* for the case of "\r\n" or "\n\r" */
            if (T[i+1] == '\n' || T[i+1] == '\r') {
                ++i;
            }
        }
        else if (isalpha(T[i])) {
            ++c;
            ++c_tot;
        }
        else if (T[i] == '\n' || T[i] == '\r') {
            /* ignore newline */
            /* for the case of "\r\n" or "\n\r" */
            if (T[i+1] == '\n' || T[i+1] == '\r') {
                ++i;
            }
        }
        else if (isprint(T[i])) {
            fprintf(stderr, "error: non-alpha character ('%c')\n", T[i]);
            exit(EXIT_FAILURE);
        }
        else {
            fprintf(stderr, "error: non-printing character ('%d')\n", T[i]);
            exit(EXIT_FAILURE);
        }
        ++i;
    }

    /* still should have one sequence in the pipe */
    if (0 == c) {
        fprintf(stderr, "error: empty sequence at end of input\n");
        exit(EXIT_FAILURE);
    }
    stats_sample_value(&stats, c);

    pfs = (parasail_file_stat_t*)malloc(sizeof(parasail_file_stat_t));
    if (NULL == pfs) {
        fprintf(stderr, "Cannont allocate parasail_file_stat_t");
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    pfs->sequences = seq;
    pfs->characters = c_tot;
    pfs->shortest = (unsigned long)stats._min;
    pfs->longest = (unsigned long)stats._max;
    pfs->mean = (float)stats._mean;
    pfs->stddev = (float)stats_stddev(&stats);

    return pfs;
}

/* increments i until T[i] points to final newline character, returns index */
inline static off_t skip_line(char *T, off_t i)
{
    while (T[i] != '\n' && T[i] != '\r') {
        ++i;
    }

    /* for the case of "\r\n" or "\n\r" */
    if (T[i+1] == '\n' || T[i+1] == '\r') {
        ++i;
    }

    return i;
}

/*
 * Line 1 begins with a '@' character and is followed by a sequence
 * identifier and an optional description (like a FASTA title line).
 *
 * Line 2 is the raw sequence letters.
 *
 * Line 3 begins with a '+' character and is optionally followed by the
 * same sequence identifier (and any description) again.
 *
 * Line 4 encodes the quality values for the sequence in Line 2, and
 * must contain the same number of symbols as letters in the sequence.
 */
parasail_file_stat_t* parasail_stat_fastq(const parasail_file_t *pf)
{
    char *T = pf->buf;
    int first = 1;
    off_t i = 0;
    unsigned long seq = 0;
    unsigned long c = 0;
    unsigned long c_tot = 0;
    unsigned long line = 0;
    stats_t stats;
    parasail_file_stat_t *pfs = NULL;

    stats_clear(&stats);

    if (NULL == pf) {
        fprintf(stderr, "parasail_stat_fastq given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    /* read file */
    while (i<pf->size) {
        if (T[i] != '@') {
            fprintf(stderr, "poorly formatted FASTQ file\n");
            fprintf(stderr, "line %lu\n", line);
            exit(EXIT_FAILURE);
        }

        /* encountered a new sequence */
        ++seq;
        if (first) {
            first = 0;
        }
        else {
            stats_sample_value(&stats, c);
        }
        c = 0;

        i = skip_line(T, i);

        /* go to next line */
        ++i;
        ++line;

        /* rest of next line is the sequence */
        while (T[i] != '\n' && T[i] != '\r') {
            ++c;
            ++i;
        }

        /* for the case of "\r\n" or "\n\r" */
        if (T[i+1] == '\n' || T[i+1] == '\r') {
            ++i;
        }

        stats_sample_value(&stats, c);

        /* go to next line */
        ++i;
        ++line;

        if (T[i] != '+') {
            fprintf(stderr, "poorly formatted FASTQ file\n");
            fprintf(stderr, "line %lu\n", line);
            exit(EXIT_FAILURE);
        }

        i = skip_line(T, i);

        /* go to next line */
        ++i;
        ++line;

        /* rest of next line are the quality control values */
        i = skip_line(T, i);

        /* go to next line */
        ++i;
        ++line;
    }

    pfs = (parasail_file_stat_t*)malloc(sizeof(parasail_file_stat_t));
    if (NULL == pfs) {
        fprintf(stderr, "Cannont allocate parasail_file_stat_t");
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    pfs->sequences = seq;
    pfs->characters = c_tot;
    pfs->shortest = (unsigned long)stats._min;
    pfs->longest = (unsigned long)stats._max;
    pfs->mean = (float)stats._mean;
    pfs->stddev = (float)stats_stddev(&stats);

    return pfs;
}

