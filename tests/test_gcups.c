#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#define USE_ZLIB 0

#include "kseq.h"

#if USE_ZLIB
#include <zlib.h>
KSEQ_INIT(gzFile, gzread)
#else
KSEQ_INIT(int, read)
#endif

static inline void parse_sequences(
        const char *filename, char ***strings_, size_t **sizes_, size_t *count_)
{
#if USE_ZLIB
    gzFile fp;
#else
    FILE *fp;
#endif
    kseq_t *seq = NULL;
    int l = 0;
    char **strings = NULL;
    size_t *sizes = NULL;
    size_t count = 0;
    size_t memory = 1000;

#if USE_ZLIB
    fp = gzopen(filename, "r");
#else
    fp = fopen(filename, "r");
    if(fp == NULL) {
        perror("fopen");
        exit(1);
    }
#endif
    strings = (char**)malloc(sizeof(char*) * memory);
    sizes = (size_t*)malloc(sizeof(size_t) * memory);
#if USE_ZLIB
    seq = kseq_init(fp);
#else
    seq = kseq_init(fileno(fp));
#endif
    while ((l = kseq_read(seq)) >= 0) {
        strings[count] = strdup(seq->seq.s);
        if (NULL == strings[count]) {
            perror("strdup");
            exit(1);
        }
        sizes[count] = seq->seq.l;
        ++count;
        if (count >= memory) {
            char **new_strings = NULL;
            size_t *new_sizes = NULL;
            memory *= 2;
            new_strings = realloc(strings, sizeof(char*) * memory);
            if (NULL == new_strings) {
                perror("realloc");
                exit(1);
            }
            strings = new_strings;
            new_sizes = realloc(sizes, sizeof(size_t) * memory);
            if (NULL == new_sizes) {
                perror("realloc");
                exit(1);
            }
            sizes = new_sizes;
        }
    }
    kseq_destroy(seq);
#if USE_ZLIB
    gzclose(fp);
#else
    fclose(fp);
#endif

    *strings_ = strings;
    *sizes_ = sizes;
    *count_ = count;
}

static inline unsigned long binomial_coefficient(unsigned long n, unsigned long k)
{
    /* from http://blog.plover.com/math/choose.html */
    unsigned long r = 1;
    unsigned long d;
    if (k > n) {
        return 0;
    }
    for (d = 1; d <= k; d++) {
        r *= n--;
        r /= d;
    }
    return r;
}

static inline void k_combination2(unsigned long pos, unsigned long *a, unsigned long *b)
{
    double s;
    double i = floor(sqrt(2.0 * pos)) - 1.0;
    if (i <= 1.0) {
        i = 1.0;
    }
    s = i * (i - 1.0) / 2.0;
    while (pos - s >= i) {
        s += i;
        i += 1;
    }
    *a = (unsigned long)(pos - s);
    *b = (unsigned long)(i);
}

int main(int argc, char **argv)
{
    unsigned long i = 0;
    size_t seq_count = 10;
    size_t limit = 0;
    char **sequences = NULL;
    size_t *sizes = NULL;
    char *filename = NULL;
    int c = 0;
    int distribution = 0;

    while ((c = getopt(argc, argv, "f:d")) != -1) {
        switch (c) {
            case 'd':
                distribution = 1;
                break;
            case 'f':
                filename = optarg;
                break;
            case '?':
                if (optopt == 'f' || optopt == 'n') {
                    fprintf(stderr,
                            "Option -%c requires an argument.\n",
                            optopt);
                }
                else if (isprint(optopt)) {
                    fprintf(stderr, "Unknown option `-%c'.\n",
                            optopt);
                }
                else {
                    fprintf(stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                }
                exit(1);
            default:
                fprintf(stderr, "default case in getopt\n");
                exit(1);
        }
    }

    if (filename) {
        parse_sequences(filename, &sequences, &sizes, &seq_count);
    }
    else {
        fprintf(stderr, "missing filename\n");
        exit(1);
    }

    if (distribution) {
        for (i=0; i<seq_count; ++i) {
            printf("%lu\n", (unsigned long)sizes[i]);
        }
    }
    else {
        limit = binomial_coefficient(seq_count, 2);
        printf("%lu choose 2 is %lu\n",
                (unsigned long)seq_count, (unsigned long)limit);

        {
            unsigned long a=0;
            unsigned long b=1;
            unsigned long work=0;
            unsigned long columns=0;
            for (i=0; i<limit; ++i) {
                k_combination2(i, &a, &b);
                work += sizes[a]*sizes[b];
                columns += sizes[b];
            }
            printf("work=%lu\n", work);
            printf("columns=%lu\n", columns);
        }
    }

    return 0;
}

