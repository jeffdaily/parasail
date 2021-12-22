/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef HAVE_SYS_MMAN_H
#include <sys/mman.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "parasail.h"
#include "parasail/io.h"
#include "parasail/memory.h"
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

/* special functions may exist for windows, msys, mingw */
#if defined(HAVE_STRUCT___STAT64) && defined(HAVE__STAT64) && defined(HAVE__FSTAT64)
#define STATBUF struct __stat64
#define STATFUNC _stat64
#define FSTATFUNC _fstat64
#else
#define STATBUF struct stat
#define STATFUNC stat
#define FSTATFUNC fstat
#endif

parasail_file_t* parasail_open(const char *fname)
{
    parasail_file_t *pf = NULL;
    char *buf = NULL;
#if defined(HAVE_SYS_MMAN_H)
    int fd = -1;
    STATBUF fs;
#else
    FILE *fd = NULL;
    STATBUF fs;
#endif

    if (NULL == fname) {
        fprintf(stderr, "parasail_open: NULL filename\n");
        return NULL;
    }

#if defined(HAVE_SYS_MMAN_H)
    fd = open(fname, O_RDONLY);
    if (fd == -1) {
        perror("open");
        fprintf(stderr, "parasail_open: "
                "cannot open input file `%s'\n", fname);
        return NULL;
    }

    if (-1 == FSTATFUNC(fd, &fs)) {
        perror("fstat");
        fprintf(stderr, "parasail_open: "
                "cannont stat input file `%s'\n", fname);
        return NULL;
    }

    buf = (char*)mmap(NULL, fs.st_size, PROT_READ, MAP_SHARED, fd, 0);
    if (MAP_FAILED == buf) {
        perror("mmap");
        fprintf(stderr, "parasail_open: "
                "cannont mmap input file `%s'\n", fname);
        return NULL;
    }
#else
    fd = fopen(fname, "rb");
    if (NULL == fd) {
        perror("fopen");
        fprintf(stderr, "parasail_open: "
                "cannot open input file `%s'\n", fname);
        return NULL;
    }

    if (0 != STATFUNC(fname, &fs)) {
        perror("_stat");
        fprintf(stderr, "parasail_open: "
                "cannont stat input file `%s'\n", fname);
        return NULL;
    }

    /* Allocate a buffer to hold the whole file */
    buf = (char*)malloc(fs.st_size + 1);
    if (NULL == buf) {
        perror("malloc");
        fprintf(stderr, "parasail_open: "
                "cannont malloc buffer for input file `%s'\n", fname);
        return NULL;
    }
    /* Slurp file into buffer */
    if (fs.st_size != fread(buf, 1, fs.st_size, fd)) {
        perror("fread");
        fprintf(stderr, "parasail_open: "
                "cannont read input file `%s'\n", fname);
        free(buf);
        return NULL;
    }
    /* Close the file early */
    fclose(fd);
#endif

    pf = (parasail_file_t*)malloc(sizeof(parasail_file_t));
    if (NULL == pf) {
        perror("malloc");
        fprintf(stderr, "parasail_open: "
                "cannont allocate parasail_file_t\n");
        free(buf);
        return NULL;
    }

#if defined(HAVE_SYS_MMAN_H)
    pf->fd = fd;
#else
    pf->fd = 0;
#endif
    pf->size = fs.st_size;
    pf->buf = buf;
    return pf;
}

void parasail_close(parasail_file_t *pf)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail_close: NULL file\n");
        return;
    }
#if defined(HAVE_SYS_MMAN_H)
    if (-1 == munmap((void*)pf->buf, pf->size)) {
        perror("munmap");
        fprintf(stderr, "parasail_close: cannot munmap file buffer\n");
    }
    if (-1 == close(pf->fd)) {
        perror("close");
        fprintf(stderr, "parasail_close: cannot close file descriptor\n");
    }
#else
    free((void*)pf->buf);
    /* file was already closed */
#endif
    free(pf);
}

int parasail_is_fasta(const parasail_file_t *pf)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail_is_fasta: NULL pointer\n");
        return -1;
    }

    return parasail_is_fasta_buffer(pf->buf, pf->size);
}

int parasail_is_fasta_buffer(const char *buf, off_t UNUSED(size))
{
    if (NULL == buf) {
        fprintf(stderr, "parasail_is_fasta_buffer: NULL pointer\n");
        return -1;
    }

    return '>' == buf[0];
}

int parasail_is_fastq(const parasail_file_t *pf)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail_is_fastq: NULL pointer\n");
        return -1;
    }

    return parasail_is_fastq_buffer(pf->buf, pf->size);
}

int parasail_is_fastq_buffer(const char *buf, off_t UNUSED(size))
{
    if (NULL == buf) {
        fprintf(stderr, "parasail_is_fastq_buffer: NULL pointer\n");
        return -1;
    }

    return '@' == buf[0];
}

parasail_file_stat_t* parasail_stat(const parasail_file_t *pf)
{
    parasail_file_stat_t *stat = NULL;

    if (NULL == pf) {
        fprintf(stderr, "parasail_stat: NULL pointer\n");
        return NULL;
    }

    if (1 == parasail_is_fasta(pf)) {
        stat = parasail_stat_fasta(pf);
    }
    else if (1 == parasail_is_fastq(pf)) {
        stat = parasail_stat_fastq(pf);
    }
    else {
        fprintf(stderr, "parasail_stat: cannot determine file format\n");
        return NULL;
    }

    if (NULL == stat) {
        fprintf(stderr, "parasail_stat: failed\n");
        return NULL;
    }

    return stat;
}

parasail_file_stat_t* parasail_stat_buffer(const char *buf, off_t size)
{
    parasail_file_stat_t *stat = NULL;

    if (NULL == buf) {
        fprintf(stderr, "parasail_stat_buffer: NULL pointer\n");
        return NULL;
    }

    if (1 == parasail_is_fasta_buffer(buf, size)) {
        stat = parasail_stat_fasta_buffer(buf, size);
    }
    else if (parasail_is_fastq_buffer(buf, size)) {
        stat = parasail_stat_fastq_buffer(buf, size);
    }
    else {
        fprintf(stderr, "parasail_stat: cannot determine file format\n");
        return NULL;
    }

    if (NULL == stat) {
        fprintf(stderr, "parasail_stat_buffer: failed\n");
        return NULL;
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
        fprintf(stderr, "parasail_stat_fasta: NULL pointer\n");
        return NULL;
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
        fprintf(stderr, "parasail_stat_fasta_buffer: NULL pointer\n");
        return NULL;
    }

    /* first line is always first sequence ID */
    if (T[i] != '>') {
        fprintf(stderr, "parasail_stat_fasta_buffer: "
                "poorly formatted FASTA file\n");
        return NULL;
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
            fprintf(stderr, "parasail_stat_fasta_buffer: "
                    "non-alpha character ('%c')\n", T[i]);
            return NULL;
        }
        else {
            fprintf(stderr, "parasail_stat_fasta_buffer: "
                    "non-printing character ('%d')\n", T[i]);
            return NULL;
        }
        ++i;
    }

    /* still should have one sequence in the pipe */
    if (0 == c) {
        fprintf(stderr, "parasail_stat_fasta_buffer: "
                "empty sequence at end of input\n");
        return NULL;
    }
    stats_sample_value(&stats, c);

    pfs = (parasail_file_stat_t*)malloc(sizeof(parasail_file_stat_t));
    if (NULL == pfs) {
        perror("malloc");
        fprintf(stderr, "parasail_stat_fasta_buffer: "
                "cannont allocate parasail_file_stat_t");
        return NULL;
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
        fprintf(stderr, "parasail_stat_fastq: NULL pointer\n");
        return NULL;
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
        fprintf(stderr, "parasail_stat_fastq_buffer: NULL pointer\n");
        return NULL;
    }

    /* read file */
    while (i<size) {
        if (T[i] != '@') {
            fprintf(stderr, "parasail_stat_fastq_buffer: "
                    "poorly formatted FASTQ file, line %lu\n", line);
            return NULL;
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
            fprintf(stderr, "parasail_stat_fastq_buffer: "
                    "poorly formatted FASTQ file, line %lu\n", line);
            return NULL;
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
        perror("malloc");
        fprintf(stderr, "parasail_stat_fastq_buffer: "
                "cannont allocate parasail_file_stat_t");
        return NULL;
    }

    pfs->sequences = seq;
    pfs->characters = c_tot;
    pfs->shortest = (unsigned long)stats._min;
    pfs->longest = (unsigned long)stats._max;
    pfs->mean = (float)stats._mean;
    pfs->stddev = (float)stats_stddev(&stats);

    return pfs;
}

char * parasail_read(const parasail_file_t *pf, long * size)
{
    char * buffer = NULL;

    if (NULL == pf) {
        fprintf(stderr, "parasail_read: NULL pointer\n");
        return NULL;
    }

    if (NULL == size) {
        fprintf(stderr, "parasail_read: NULL size pointer\n");
        return NULL;
    }

    buffer = (char*)malloc(sizeof(char) * (pf->size+1));
    if (NULL == buffer) {
        perror("malloc");
        fprintf(stderr, "parasail_read: "
                "cannont malloc buffer for input file");
        return NULL;
    }

    (void)memcpy(buffer, pf->buf, pf->size);
    buffer[pf->size] = '\0';
    *size = pf->size;
    return buffer;
}

char * parasail_pack(const parasail_file_t *pf, long * size)
{
    char *packed = NULL;

    if (NULL == pf) {
        fprintf(stderr, "parasail_pack: NULL pointer\n");
        return NULL;
    }

    if (NULL == size) {
        fprintf(stderr, "parasail_pack: NULL size pointer\n");
        return NULL;
    }

    if (1 == parasail_is_fasta(pf)) {
        packed = parasail_pack_fasta(pf, size);
    }
    else if (1 == parasail_is_fastq(pf)) {
        packed = parasail_pack_fastq(pf, size);
    }
    else {
        fprintf(stderr, "parasail_pack: cannot determine file format\n");
        return NULL;
    }

    if (NULL == packed) {
        fprintf(stderr, "parasail_pack: failed\n");
        return NULL;
    }

    return packed;
}

char * parasail_pack_buffer(const char *buf, off_t size, long * packed_size)
{
    char *packed = NULL;

    if (NULL == buf) {
        fprintf(stderr, "parasail_pack_buffer: NULL pointer\n");
        return NULL;
    }

    if (NULL == packed_size) {
        fprintf(stderr, "parasail_pack_buffer: NULL size pointer\n");
        return NULL;
    }

    if (1 == parasail_is_fasta_buffer(buf, size)) {
        packed = parasail_pack_fasta_buffer(buf, size, packed_size);
    }
    else if (1 == parasail_is_fastq_buffer(buf, size)) {
        packed = parasail_pack_fastq_buffer(buf, size, packed_size);
    }
    else {
        fprintf(stderr, "parasail_pack: cannot determine file format\n");
        return NULL;
    }

    if (NULL == packed) {
        fprintf(stderr, "parasail_pack_buffer: failed\n");
        return NULL;
    }

    return packed;
}

char * parasail_pack_fasta(const parasail_file_t *pf, long * packed_size)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail_pack_fasta: NULL pointer\n");
        return NULL;
    }

    if (NULL == packed_size) {
        fprintf(stderr, "parasail_pack_fasta: NULL size pointer\n");
        return NULL;
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
        fprintf(stderr, "parasail_pack_fasta_buffer: NULL pointer\n");
        return NULL;
    }

    if (NULL == packed_size) {
        fprintf(stderr, "parasail_pack_fasta_buffer: NULL size pointer\n");
        return NULL;
    }

    pfs = parasail_stat_fasta_buffer(T, size);
    if (NULL == pfs) {
        fprintf(stderr, "parasail_stat_fasta_buffer: fasta stat failed\n");
        return NULL;
    }

    P = (char*)malloc(sizeof(char) * (pfs->characters+pfs->sequences+1));
    if (NULL == P) {
        perror("malloc");
        fprintf(stderr, "parasail_pack_fasta_buffer: malloc failed\n");
        free(pfs);
        return NULL;
    }
    free(pfs);

    /* first line is always first sequence ID */
    if (T[i] != '>') {
        fprintf(stderr, "parasail_pack_fasta_buffer: "
                "poorly formatted FASTA file\n");
        free(P);
        return NULL;
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
            fprintf(stderr, "parasail_pack_fasta_buffer: "
                    "non-alpha character ('%c')\n", T[i]);
            free(P);
            return NULL;
        }
        else {
            fprintf(stderr, "parasail_pack_fasta_buffer: "
                    "non-printing character ('%d')\n", T[i]);
            free(P);
            return NULL;
        }
        ++i;
    }

    P[w++] = '$';
    P[w] = '\0';
    *packed_size = w;
    return P;
}

char * parasail_pack_fastq(const parasail_file_t *pf, long * size)
{
    if (NULL == pf) {
        fprintf(stderr, "parasail_pack_fastq: NULL pointer\n");
        return NULL;
    }

    if (NULL == size) {
        fprintf(stderr, "parasail_pack_fastq: NULL size pointer\n");
        return NULL;
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
        fprintf(stderr, "parasail_pack_fastq_buffer: NULL pointer\n");
        return NULL;
    }

    if (NULL == packed_size) {
        fprintf(stderr, "parasail_pack_fastq_buffer: NULL size pointer\n");
        return NULL;
    }

    pfs = parasail_stat_fastq_buffer(T, size);
    if (NULL == pfs) {
        fprintf(stderr, "parasail_stat_fastq_buffer: fastq stat failed\n");
        return NULL;
    }

    P = (char*)malloc(sizeof(char) * (pfs->characters+pfs->sequences+1));
    if (NULL == P) {
        perror("malloc");
        fprintf(stderr, "parasail_pack_fastq_buffer: malloc failed\n");
        free(pfs);
        return NULL;
    }
    free(pfs);

    /* read file */
    while (i<size) {
        if (T[i] != '@') {
            fprintf(stderr, "parasail_pack_fastq_buffer: "
                    "poorly formatted FASTQ file, line %lu\n", line);
            free(P);
            return NULL;
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
            fprintf(stderr, "parasail_pack_fastq_buffer: "
                    "poorly formatted FASTQ file, line %lu\n", line);
            free(P);
            return NULL;
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

    P[w++] = '$';
    P[w] = '\0';
    *packed_size = w;
    return P;
}

/* increments i until T[i] points non-number, returns number */
#define TOKEN_MAX 10
inline static int get_num(const char *T, off_t *i_, int *result)
{
    off_t i = *i_;
    int retval = 0;
    int p = 0;
    char token[TOKEN_MAX];

    /* skip any whitespace */
    while (T[i] == ' ' || T[i] == '\t') {
        ++i;
    }

    if (isdigit(T[i]) || T[i] == '-') {
        token[0] = T[i];
    }
    else {
        return -1;
    }

    ++i;
    for (p=1; p<TOKEN_MAX; ++p) {
        if (isdigit(T[i])) {
            token[p] = T[i];
            ++i;
        }
        else {
            break;
        }
    }
    if (TOKEN_MAX == p) {
        return -1;
    }
    token[p] = '\0';

    retval = sscanf(token, "%d", result);
    if (1 != retval) {
        return -1;
    }

    *i_ = i;

    return 1;
}

#define ALPHABET_MAX 256
inline static char*  get_alphabet(const char *T, off_t i, off_t size)
{
    char *alphabet = NULL;
    off_t _i = i;
    size_t count = 0;

    /* count number of letters first */
    while (i<size) {
        if (T[i] == '\n' || T[i] == '\r' || T[i] == '#') {
            break;
        }
        else if (T[i] == ' ' || T[i] == '\t') {
            /* ignore newline */
        }
        else if (isalpha(T[i]) || T[i] == '*') {
            ++count;
        }
        else {
            return NULL;
        }
        ++i;
    }

    if (0 == count) {
        return NULL;
    }

    alphabet = (char*)malloc(sizeof(char)*(count+1));
    if (NULL == alphabet) {
        perror("malloc");
        return NULL;
    }

    i = _i;
    count = 0;
    while (i<size) {
        if (T[i] == '\n' || T[i] == '\r' || T[i] == '#') {
            break;
        }
        else if (T[i] == ' ' || T[i] == '\t') {
            /* ignore newline */
        }
        else if (isalpha(T[i]) || T[i] == '*') {
            alphabet[count] = T[i];
            ++count;
        }
        else {
            return NULL;
        }
        ++i;
    }

    alphabet[count] = '\0';
    return alphabet;
}

/* expects a matrix in the form
 * # FREQS A 0.325 C 0.175 G 0.175 T 0.325
 *     A   R   G   C   Y   T   K   M   S   W   N   X
 * A   8   3  -7 -17 -19 -21 -14  -4 -12  -6  -1 -30
 * R   0   3   2 -16 -18 -19  -8  -8  -7 -10  -1 -30
 * G -10  12  12 -16 -17 -18  -2 -13  -1 -14  -1 -30
 * C -18 -17 -16  12  12 -10 -13  -2  -1 -14  -1 -30
 * Y -19 -18 -16   2   0   0  -8  -8  -7 -10  -1 -30
 * T -21 -19 -17  -7   3   8  -4 -14 -12  -6  -1 -30
 * K -15  -9  -2 -11  -8  -4  -3 -13  -7 -10  -1 -30
 * M  -4  -8 -11  -2  -9 -15 -13  -3  -7 -10  -1 -30
 * S -14  -8  -1  -1  -8 -14  -8  -8  -1 -14  -1 -30
 * W  -6  -9 -12 -12  -9  -6  -9  -9 -12  -6  -1 -30
 * N  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
 * X -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30
 */
static parasail_matrix_t* parasail_matrix_from_file_internal(const char *filename, int case_sensitive)
{
    parasail_matrix_t *retval = NULL;
    int *matrix = NULL;
    size_t matrix_i = 0;
    size_t matrix_capacity = 0;
    int *mapper = NULL;
    char *alphabet = NULL;
    char *alphabet_query = NULL;
    size_t alphabet_query_i = 0;
    size_t alphabet_query_capacity = 0;
    parasail_file_t *pf = NULL;
    const char *T = NULL;
    off_t i = 0;
    off_t size = 0;
    int first_alpha = 1;
    size_t count = 0;
    size_t asize = 0;
    int max = INT_MIN;
    int min = INT_MAX;
    size_t c = 0;
    int type = PARASAIL_MATRIX_TYPE_SQUARE;

    if (NULL == filename) {
        fprintf(stderr, "parasail_matrix_from_file: NULL pointer\n");
        return NULL;
    }

    pf = parasail_open(filename);
    if (NULL == pf) {
        fprintf(stderr, "parasail_matrix_from_file: "
                "parasail_open(%s) failed\n", filename);
        return NULL;
    }
    T = pf->buf;
    size = pf->size;

    while (i<size) {
        if (T[i] == '#') {
            /* ignore comments */
            i = skip_line(T, i);
        }
        else if (isalnum(T[i]) || T[i] == '*' || T[i] == '-') {
            if (first_alpha) {
                first_alpha = 0;
                alphabet = get_alphabet(T, i, size);
                if (NULL == alphabet) {
                    fprintf(stderr, "parasail_matrix_from_file: "
                            "poorly formed matrix file alphabet\n");
                    parasail_close(pf);
                    return NULL;
                }
                asize = strlen(alphabet);
                matrix_capacity = asize*asize;
                matrix = (int*)malloc(sizeof(int)*matrix_capacity);
                if (NULL == matrix) {
                    perror("malloc");
                    fprintf(stderr, "parasail_matrix_from_file: "
                            "cannont malloc buffer for matrix\n");
                    free(alphabet);
                    parasail_close(pf);
                    return NULL;
                }
                alphabet_query_capacity = asize;
                alphabet_query = (char*)malloc(sizeof(char)*(alphabet_query_capacity+1));
                if (NULL == alphabet_query) {
                    perror("malloc");
                    fprintf(stderr, "parasail_matrix_from_file: "
                            "cannont malloc buffer for matrix alphabet\n");
                    free(matrix);
                    free(alphabet);
                    parasail_close(pf);
                    return NULL;
                }
                i = skip_line(T, i);
            }
            else {
                size_t j=0;
                /* store the letter */
                if (alphabet_query_i >= alphabet_query_capacity) {
                    alphabet_query_capacity *= 2;
                    alphabet_query = realloc(alphabet_query, sizeof(char)*(alphabet_query_capacity+1));
                    if (NULL == alphabet_query) {
                        perror("realloc");
                        fprintf(stderr, "parasail_matrix_from_file: "
                                "couldn't grow query size\n");
                        if (alphabet) free(alphabet);
                        if (alphabet_query) free(alphabet_query);
                        if (matrix) free(matrix);
                        parasail_close(pf);
                        return NULL;
                    }
                }
                ++count;
                /* allow PSSM matrix to omit query */
                if (isalpha(T[i]) || T[i] == '*') {
                    alphabet_query[alphabet_query_i++] = T[i];
                    ++i; /* skip over the letter */
                }
                else {
                    alphabet_query[alphabet_query_i++] = '*';
                    /* don't increment read index i */
                }
                for (j=0; j<asize; ++j) {
                    int val = 0;
                    int retcode = get_num(T, &i, &val);
                    if (-1 == retcode) {
                        fprintf(stderr, "parasail_matrix_from_file: "
                                "poorly formed matrix file\n");
                        if (alphabet) free(alphabet);
                        if (alphabet_query) free(alphabet_query);
                        if (matrix) free(matrix);
                        parasail_close(pf);
                        return NULL;
                    }
                    if (matrix_i >= matrix_capacity) {
                        matrix_capacity *= 2;
                        matrix = realloc(matrix, sizeof(int)*matrix_capacity);
                        if (NULL == matrix) {
                            perror("realloc");
                            fprintf(stderr, "parasail_matrix_from_file: "
                                    "couldn't grow matrix size\n");
                            if (alphabet) free(alphabet);
                            if (alphabet_query) free(alphabet_query);
                            if (matrix) free(matrix);
                            parasail_close(pf);
                            return NULL;
                        }
                    }
                    matrix[matrix_i++] = val;
                    max = val > max ? val : max;
                    min = val < min ? val : min;
                }
            }
        }
        else if (T[i] == '\n' || T[i] == '\r') {
            /* ignore newline */
            /* for the case of "\r\n" or "\n\r" */
            if (i+1<size && (T[i+1] == '\n' || T[i+1] == '\r')) {
                ++i;
            }
        }
        else if (T[i] == ' ' || T[i] == '\t') {
            /* ignore spaces and tabs */
        }
        else if (isprint(T[i])) {
            fprintf(stderr, "parasail_matrix_from_file: "
                    "non-alpha character in matrix file ('%c')\n", T[i]);
            if (alphabet) free(alphabet);
            if (alphabet_query) free(alphabet_query);
            if (matrix) free(matrix);
            parasail_close(pf);
            return NULL;
        }
        else {
            fprintf(stderr, "parasail_matrix_from_file: "
                    "non-printing character in matrix file ('%d')\n", T[i]);
            if (alphabet) free(alphabet);
            if (alphabet_query) free(alphabet_query);
            if (matrix) free(matrix);
            parasail_close(pf);
            return NULL;
        }
        ++i;
    }

    parasail_close(pf);

    /* we consider this a square matrix if alphabet sizes match and are identical */
    alphabet_query[alphabet_query_i] = '\0';
    if (asize == alphabet_query_i && 0 == strcmp(alphabet, alphabet_query)) {
        /* square */
        if (matrix_i != asize*asize) {
            fprintf(stderr, "parasail_matrix_from_file: "
                    "matrix is missing values\n");
            free(alphabet);
            free(alphabet_query);
            free(matrix);
            return NULL;
        }
        if (count != asize) {
            fprintf(stderr, "parasail_matrix_from_file: "
                    "matrix is missing rows\n");
            free(alphabet);
            free(alphabet_query);
            free(matrix);
            return NULL;
        }
    }
    else {
        /* PSSM */
        type = PARASAIL_MATRIX_TYPE_PSSM;
        if (matrix_i != asize*count) {
            fprintf(stderr, "parasail_matrix_from_file: "
                    "matrix is missing values\n");
            free(alphabet);
            free(alphabet_query);
            free(matrix);
            return NULL;
        }
    }

    mapper = (int*)malloc(sizeof(int)*256);
    if (NULL == mapper) {
        perror("malloc");
        fprintf(stderr, "parasail_matrix_from_file: "
                "cannont malloc mapper buffer for matrix file `%s'\n",
                filename);
        free(alphabet);
        free(alphabet_query);
        free(matrix);
        return NULL;
    }
    parasail_memset_int(mapper, asize-1, 256);
    if (case_sensitive) {
        for (c=0; c<asize; ++c) {
            mapper[(unsigned char)alphabet[c]] = (int)c;
        }
    }
    else {
        for (c=0; c<asize; ++c) {
            mapper[toupper((unsigned char)alphabet[c])] = (int)c;
            mapper[tolower((unsigned char)alphabet[c])] = (int)c;
        }
    }

    retval = (parasail_matrix_t*)malloc(sizeof(parasail_matrix_t));
    if (NULL == retval) {
        perror("malloc");
        fprintf(stderr, "parasail_matrix_from_file: "
                "cannont malloc buffer for matrix file `%s'\n", filename);
        free(matrix);
        return NULL;
    }

    retval->name = filename;
    retval->matrix = matrix;
    retval->mapper = mapper;
    retval->size = (int)asize;
    retval->max = max;
    retval->min = min;
    retval->user_matrix = matrix;
    retval->type = type;
    retval->length = (int)count;
    retval->alphabet = alphabet;
    retval->query = alphabet_query;
    return retval;
}

parasail_matrix_t* parasail_matrix_from_file(const char *filename)
{
    return parasail_matrix_from_file_internal(filename, 0);
}

parasail_matrix_t* parasail_matrix_from_file_case_sensitive(const char *filename)
{
    return parasail_matrix_from_file_internal(filename, 1);
}

