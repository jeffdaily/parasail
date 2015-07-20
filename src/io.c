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

#ifdef UNUSED
#elif defined(__GNUC__)
# define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#elif defined(__LCLINT__)
# define UNUSED(x) /*@unused@*/ x
#elif defined(__cplusplus)
# define UNUSED(x)
#else
# define UNUSED(x) x
#endif

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
    munmap((void*)pf->buf, pf->size);
    close(pf->fd);
    free(pf);
}

int parasail_is_fasta(const parasail_file_t *pf)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail_is_fasta given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return parasail_is_fasta_buffer(pf->buf, pf->size);
}

int parasail_is_fasta_buffer(const char *buf, off_t UNUSED(size))
{
    if (NULL == buf) {
        fprintf(stderr, "parasail_is_fasta_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return '>' == buf[0];
}

int parasail_is_fastq(const parasail_file_t *pf)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail_is_fastq given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return parasail_is_fastq_buffer(pf->buf, pf->size);
}

int parasail_is_fastq_buffer(const char *buf, off_t UNUSED(size))
{
    if (NULL == buf) {
        fprintf(stderr, "parasail_is_fastq_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return '@' == buf[0];
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

parasail_file_stat_t* parasail_stat_buffer(const char *buf, off_t size)
{
    parasail_file_stat_t *stat = NULL;

    if (NULL == buf) {
        fprintf(stderr, "parasail_stat_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (parasail_is_fasta_buffer(buf, size)) {
        stat = parasail_stat_fasta_buffer(buf, size);
    }
    else if (parasail_is_fastq_buffer(buf, size)) {
        stat = parasail_stat_fastq_buffer(buf, size);
    }
    else {
        fprintf(stderr, "parasail_stat: cannot determine file format\n");
        exit(EXIT_FAILURE);
    }

    return stat;
}

/* increments i until T[i] points to final newline character, returns index */
inline static off_t skip_line(const char *T, off_t i)
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

parasail_file_stat_t* parasail_stat_fasta(const parasail_file_t *pf)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail_stat_fasta given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return parasail_stat_fasta_buffer(pf->buf, pf->size);
}

parasail_file_stat_t* parasail_stat_fasta_buffer(const char *T, off_t size)
{
    off_t i = 0;
    unsigned long seq = 0;
    unsigned long c = 0;
    unsigned long c_tot = 0;
    stats_t stats;
    parasail_file_stat_t *pfs = NULL;

    stats_clear(&stats);

    if (NULL == T) {
        fprintf(stderr, "parasail_stat_fasta_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    /* first line is always first sequence ID */
    if (T[i] != '>') {
        fprintf(stderr, "poorly formatted FASTA file\n");
        exit(EXIT_FAILURE);
    }

    i = skip_line(T, i);
    ++i;

    /* count that first sequence */
    ++seq;

    /* read rest of file */
    while (i<size) {
        if (T[i] == '>') {
            /* encountered a new sequence */
            ++seq;
            stats_sample_value(&stats, c);
            c = 0;
            i = skip_line(T, i);
        }
        else if (isalpha(T[i])) {
            ++c;
            ++c_tot;
        }
        else if (T[i] == '\n' || T[i] == '\r') {
            /* ignore newline */
            /* for the case of "\r\n" or "\n\r" */
            if (i+1<size && (T[i+1] == '\n' || T[i+1] == '\r')) {
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
    if (NULL == pf) {
        fprintf(stderr, "parasail_stat_fastq given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    return parasail_stat_fastq_buffer(pf->buf, pf->size);
}

parasail_file_stat_t* parasail_stat_fastq_buffer(const char *T, off_t size)
{
    int first = 1;
    off_t i = 0;
    unsigned long seq = 0;
    unsigned long c = 0;
    unsigned long c_tot = 0;
    unsigned long line = 0;
    stats_t stats;
    parasail_file_stat_t *pfs = NULL;

    stats_clear(&stats);

    if (NULL == T) {
        fprintf(stderr, "parasail_stat_fastq_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    /* read file */
    while (i<size) {
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
        if (i+1<size && (T[i+1] == '\n' || T[i+1] == '\r')) {
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

char * parasail_pack(const parasail_file_t *pf, long * size)
{
    char *packed = NULL;

    if (NULL == pf) {
        fprintf(stderr, "parasail_pack given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (parasail_is_fasta(pf)) {
        packed = parasail_pack_fasta(pf, size);
    }
    else if (parasail_is_fastq(pf)) {
        packed = parasail_pack_fastq(pf, size);
    }
    else {
        fprintf(stderr, "parasail_pack: cannot determine file format\n");
        exit(EXIT_FAILURE);
    }

    return packed;
}

char * parasail_pack_buffer(const char *buf, off_t size, long * packed_size)
{
    char *packed = NULL;

    if (NULL == buf) {
        fprintf(stderr, "parasail_pack_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (parasail_is_fasta_buffer(buf, size)) {
        packed = parasail_pack_fasta_buffer(buf, size, packed_size);
    }
    else if (parasail_is_fastq_buffer(buf, size)) {
        packed = parasail_pack_fastq_buffer(buf, size, packed_size);
    }
    else {
        fprintf(stderr, "parasail_pack: cannot determine file format\n");
        exit(EXIT_FAILURE);
    }

    return packed;
}

char * parasail_pack_fasta(const parasail_file_t *pf, long * packed_size)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail_pack_fasta given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (NULL == packed_size) {
        fprintf(stderr, "parasail_pack_fasta given NULL size pointer\n");
        exit(EXIT_FAILURE);
    }

    return parasail_pack_fasta_buffer(pf->buf, pf->size, packed_size);
}

char * parasail_pack_fasta_buffer(const char *T, off_t size, long * packed_size)
{
    parasail_file_stat_t *pfs = NULL;
    off_t i = 0;
    off_t w = 0;
    char *P = NULL;

    if (NULL == T) {
        fprintf(stderr, "parasail_pack_fasta_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (NULL == packed_size) {
        fprintf(stderr, "parasail_pack_fasta_buffer given NULL size pointer\n");
        exit(EXIT_FAILURE);
    }

    pfs = parasail_stat_fasta_buffer(T, size);

    P = (char*)malloc(sizeof(char) * (pfs->characters+pfs->sequences+1));

    /* first line is always first sequence ID */
    if (T[i] != '>') {
        fprintf(stderr, "poorly formatted FASTA file\n");
        exit(EXIT_FAILURE);
    }

    i = skip_line(T, i);
    ++i;

    /* read rest of file */
    while (i<size) {
        if (T[i] == '>') {
            /* encountered a new sequence */
            P[w++] = '$';
            i = skip_line(T, i);
        }
        else if (isalpha(T[i])) {
            P[w++] = T[i];
        }
        else if (T[i] == '\n' || T[i] == '\r') {
            /* ignore newline */
            /* for the case of "\r\n" or "\n\r" */
            if (i+1<size && (T[i+1] == '\n' || T[i+1] == '\r')) {
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

    free(pfs);

    P[w++] = '$';
    P[w] = '\0';
    *packed_size = w;
    return P;
}

char * parasail_pack_fastq(const parasail_file_t *pf, long * size)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail_pack_fastq given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (NULL == size) {
        fprintf(stderr, "parasail_pack_fastq given NULL size pointer\n");
        exit(EXIT_FAILURE);
    }

    return parasail_pack_fastq_buffer(pf->buf, pf->size, size);
}

char * parasail_pack_fastq_buffer(const char *T, off_t size, long * packed_size)
{
    char *P = NULL;
    int first = 1;
    off_t i = 0;
    off_t w = 0;
    unsigned long line = 0;
    parasail_file_stat_t *pfs = NULL;

    if (NULL == T) {
        fprintf(stderr, "parasail_pack_fastq_buffer given NULL pointer\n");
        exit(EXIT_FAILURE);
    }

    if (NULL == packed_size) {
        fprintf(stderr, "parasail_pack_fastq_buffer given NULL size pointer\n");
        exit(EXIT_FAILURE);
    }

    pfs = parasail_stat_fastq_buffer(T, size);

    P = (char*)malloc(sizeof(char) * (pfs->characters+pfs->sequences+1));

    /* read file */
    while (i<size) {
        if (T[i] != '@') {
            fprintf(stderr, "poorly formatted FASTQ file\n");
            fprintf(stderr, "line %lu\n", line);
            exit(EXIT_FAILURE);
        }

        /* encountered a new sequence */
        if (first) {
            first = 0;
        }
        else {
            P[w++] = '$';
        }

        i = skip_line(T, i);

        /* go to next line */
        ++i;
        ++line;

        /* rest of next line is the sequence */
        while (T[i] != '\n' && T[i] != '\r') {
            P[w++] = T[i];
            ++i;
        }

        /* for the case of "\r\n" or "\n\r" */
        if (i+1<size && (T[i+1] == '\n' || T[i+1] == '\r')) {
            ++i;
        }

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

    free(pfs);

    P[w++] = '$';
    P[w] = '\0';
    *packed_size = w;
    return P;
}

