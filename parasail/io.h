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
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct parasail_file {
    int fd;
    off_t size;
    const char *buf;
} parasail_file_t;

typedef struct parasail_file_stat {
    unsigned long sequences;
    unsigned long characters;
    unsigned long label_characters;
    unsigned long shortest;
    unsigned long longest;
    float mean;
    float stddev;
} parasail_file_stat_t;

typedef struct parasail_file_labels {
    unsigned long sequences;
    char ** labels;
    unsigned long * sizes;
    char * memory;
} parasail_file_labels_t;

extern parasail_file_t* parasail_open(const char *fname);

/** Closes file and frees file parameter. */
extern void parasail_close(parasail_file_t *file);

/** Frees the label structure. */
extern void parasail_free_file_labels(parasail_file_labels_t *labels);

extern int parasail_is_fasta(const parasail_file_t *pf);

extern int parasail_is_fastq(const parasail_file_t *pf);

extern parasail_file_stat_t* parasail_stat(const parasail_file_t *pf);

extern parasail_file_stat_t* parasail_stat_fasta(const parasail_file_t *pf);

extern parasail_file_stat_t* parasail_stat_fastq(const parasail_file_t *pf);

extern char * parasail_read(const parasail_file_t *pf, long * size);

extern char * parasail_pack(const parasail_file_t *pf, long * packed_size);

extern char * parasail_pack_fasta(const parasail_file_t *pf, long * packed_size);

extern char * parasail_pack_fastq(const parasail_file_t *pf, long * packed_size);

extern char * parasail_pack_with_labels(const parasail_file_t *pf, long * packed_size, parasail_file_labels_t **labels);

extern char * parasail_pack_fasta_with_labels(const parasail_file_t *pf, long * packed_size, parasail_file_labels_t **labels);

extern char * parasail_pack_fastq_with_labels(const parasail_file_t *pf, long * packed_size, parasail_file_labels_t **labels);


/* char buffer versions of io functions */

extern int parasail_is_fasta_buffer(const char *, off_t size);

extern int parasail_is_fastq_buffer(const char *pf, off_t size);

extern parasail_file_stat_t* parasail_stat_buffer(const char *buffer, off_t size);

extern parasail_file_stat_t* parasail_stat_fasta_buffer(const char *buffer, off_t size);

extern parasail_file_stat_t* parasail_stat_fastq_buffer(const char *buffer, off_t size);

extern char * parasail_pack_buffer(const char *buffer, off_t size, long * packed_size);

extern char * parasail_pack_fasta_buffer(const char *buffer, off_t size, long * packed_size);

extern char * parasail_pack_fastq_buffer(const char *buffer, off_t size, long * packed_size);

extern char * parasail_pack_buffer_with_labels(const char *buffer, off_t size, long * packed_size, parasail_file_labels_t **labels);

extern char * parasail_pack_fasta_buffer_with_labels(const char *buffer, off_t size, long * packed_size, parasail_file_labels_t **labels);

extern char * parasail_pack_fastq_buffer_with_labels(const char *buffer, off_t size, long * packed_size, parasail_file_labels_t **labels);

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_IO_H_ */
