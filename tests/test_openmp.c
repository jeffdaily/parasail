#include "config.h"

#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <omp.h>

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "parasail.h"
#include "blosum/blosum62.h"
#include "timer.h"
#include "timer_real.h"

static inline void parse_sequences(
        const char *filename, char ***strings_, size_t **sizes_, size_t *count_)
{
    gzFile fp;
    kseq_t *seq = NULL;
    int l = 0;
    char **strings = NULL;
    size_t *sizes = NULL;
    size_t count = 0;
    size_t memory = 1000;
    size_t i = 0;

    fp = gzopen(filename, "r");
    if(fp == Z_NULL) {
        perror("gzopen");
        exit(1);
    }
    strings = malloc(sizeof(char*) * memory);
    sizes = malloc(sizeof(size_t) * memory);
    seq = kseq_init(fp);
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
    gzclose(fp);

    *strings_ = strings;
    *sizes_ = sizes;
    *count_ = count;
}

static inline char* rand_string(size_t size)
{
    char *str = NULL;
    const char charset[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    if (size) {
        size_t n;
        --size;
        str = malloc(size + 1);
        for (n = 0; n < size; n++) {
            int key = rand() % (int) (sizeof charset - 1);
            str[n] = charset[key];
        }
        str[size] = '\0';
    }
    return str;
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
    int shortest = INT_MAX;
    int longest = 0;
    double timer_clock = 0.0;
    unsigned long long timer_rtdsc = 0U;
    unsigned long i = 0;
    size_t seq_count = 10;
    size_t limit = 0;
    char **sequences = NULL;
    size_t *sizes = NULL;
    char *endptr = NULL;
    char *filename = NULL;
    int c = 0;

    while ((c = getopt(argc, argv, "f:n:")) != -1) {
        switch (c) {
            case 'f':
                filename = optarg;
                break;
            case 'n':
                errno = 0;
                seq_count = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
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
        /* generate 'seq_count' number of random strings */
        sequences = (char**)malloc(sizeof(char*)*seq_count);
        sizes = (size_t*)malloc(sizeof(size_t)*seq_count);
        for (i=0; i<seq_count; ++i) {
            sizes[i] = (rand()%32767)+10;
            shortest = sizes[i] < shortest ? sizes[i] : shortest;
            longest = sizes[i] > longest ? sizes[i] : longest;
            sequences[i] = rand_string(sizes[i]);
        }
        printf("done generating %lu random srings, "
                "shortest is %d, longest is %d\n",
                seq_count, shortest, longest);
    }

    limit = binomial_coefficient(seq_count, 2);
    printf("%lu choose 2 is %lu\n", seq_count, limit);

    timer_init();
    printf("%s timer\n", timer_name());

#pragma omp parallel
    {
#pragma omp single
        {
            int N = omp_get_max_threads();
            printf("omp_get_max_threads()=%d\n", N);
        }
    }

    timer_clock = timer_real();
    timer_rtdsc = timer_start();
#pragma omp parallel
    {
        unsigned long a=0;
        unsigned long b=1;
        parasail_result_t *result = NULL;
#pragma omp for schedule(dynamic)
        for (i=0; i<limit; ++i) {
            k_combination2(i, &a, &b);
            result = nw(sequences[a], sizes[a], sequences[b], sizes[b],
                    10, 1, blosum62);
            parasail_result_free(result);
        }
    }
    timer_rtdsc = timer_end(timer_rtdsc);
    timer_clock = timer_real() - timer_clock;
    printf("nw\t%llu\t%f\n", timer_rtdsc, timer_clock);

    return 0;
}

