/**
 * @file parasail_stats
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2015 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Reads fasta/fastq file of database sequences and computes statistics.
 */
#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "parasail.h"
#include "parasail/io.h"

static void print_help(const char *progname, int status) {
    fprintf(stderr, "\nusage: %s "
            "file "
            "\n\n",
            progname);
    exit(status);
}

int main(int argc, char **argv) {
    parasail_file_t *pf = NULL;
    parasail_file_stat_t *pfs = NULL;
    const char *progname = "parasail_aligner";
    const char *type = NULL;

    /* Check arguments. */
    if (argc > 2) {
        fprintf(stderr, "Too many arguments.\n");
        print_help(progname, EXIT_FAILURE);
    }
    else if (argc < 2) {
        fprintf(stderr, "Missing input file.\n");
        print_help(progname, EXIT_FAILURE);
    }

    /* open file */
    pf = parasail_open(argv[1]);

    /* check type */
    if (parasail_is_fasta(pf)) {
        type = "FASTA";
    }
    else if (parasail_is_fastq(pf)) {
        type = "FASTQ";
    }
    else {
        fprintf(stderr, "unrecognized file format\n");
        exit(EXIT_FAILURE);
    }

    /* compute stats */
    pfs = parasail_stat(pf);

    /* print the stats */
    fprintf(stdout,
            "%25s: %s\n"
            "%25s: %lu\n"
            "%25s: %lu\n"
            "%25s: %lu\n"
            "%25s: %lu\n"
            "%25s: %f\n"
            "%25s: %f\n",
            "file type", type,
            "sequence count", pfs->sequences,
            "character count", pfs->characters,
            "shortest sequence", pfs->shortest,
            "longest sequence", pfs->longest,
            "sequence length mean", pfs->mean,
            "sequence length stddev", pfs->stddev
            );

    parasail_close(pf);
    free(pfs);

    return 0;
}

