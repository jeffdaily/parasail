#include "config.h"

#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "kseq.h"
KSEQ_INIT(int, read)

#include "parasail.h"
#include "parasail_internal.h"
#include "blosum/blosum40.h"
#include "blosum/blosum45.h"
#include "blosum/blosum50.h"
#include "blosum/blosum62.h"
#include "blosum/blosum75.h"
#include "blosum/blosum80.h"
#include "blosum/blosum90.h"
#include "stats.h"
//#include "timer.h"
#include "timer_real.h"

#include "blosum_lookup.h"
#include "function_lookup.h"

static inline size_t parse_sequences(
        const char *filename, char ***strings_, size_t **sizes_, size_t *count_)
{
    FILE* fp;
    kseq_t *seq = NULL;
    int l = 0;
    char **strings = NULL;
    size_t *sizes = NULL;
    size_t count = 0;
    size_t memory = 1000;
    size_t biggest = 0;

    fp = fopen(filename, "r");
    if(fp == NULL) {
        perror("fopen");
        exit(1);
    }
    strings = malloc(sizeof(char*) * memory);
    sizes = malloc(sizeof(size_t) * memory);
    seq = kseq_init(fileno(fp));
    while ((l = kseq_read(seq)) >= 0) {
        strings[count] = strdup(seq->seq.s);
        if (NULL == strings[count]) {
            perror("strdup");
            exit(1);
        }
        sizes[count] = seq->seq.l;
        if (sizes[count] > biggest)
            biggest = sizes[count];
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
    fclose(fp);

    *strings_ = strings;
    *sizes_ = sizes;
    *count_ = count;
    return biggest;
}

int main(int argc, char **argv)
{
    unsigned long i = 0;
    unsigned long j = 0;
    char *filename_database = NULL;
    char **sequences_database = NULL;
    size_t *sizes_database = NULL;
    size_t seq_count_database = 0;
    char *filename_queries = NULL;
    char **sequences_queries = NULL;
    size_t *sizes_queries = NULL;
    size_t seq_count_queries = 0;
    char *endptr = NULL;
    char *funcname = NULL;
    parasail_function_t function = NULL;
    int c = 0;
    char *blosumname = NULL;
    parasail_blosum_t blosum = blosum62;
    int gap_open = 10;
    int gap_extend = 1;
    int truncate = 0;
    int exact_length = 0;
    stats_t stats_time;
#if ENABLE_CORRECTION_STATS
    unsigned long long corrections = 0;
#endif
    size_t biggest = 0;

    stats_clear(&stats_time);

    while ((c = getopt(argc, argv, "a:b:f:q:o:e:t:x:")) != -1) {
        switch (c) {
            case 'a':
                funcname = optarg;
                break;
            case 'b':
                blosumname = optarg;
                break;
            case 'f':
                filename_database = optarg;
                break;
            case 'q':
                filename_queries = optarg;
                break;
            case 'o':
                errno = 0;
                gap_open = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case 'e':
                errno = 0;
                gap_extend = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case 't':
                errno = 0;
                truncate = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case 'x':
                errno = 0;
                exact_length = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case '?':
                if (optopt == 'a'
                        || optopt == 'b'
                        || optopt == 'f'
                        || optopt == 'q'
                        || optopt == 'o'
                        || optopt == 'e'
                        || optopt == 't'
                        || optopt == 'x')
                {
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

    /* select the function */
    if (funcname) {
        function = lookup_function(funcname);
        if (NULL == function) {
            fprintf(stderr, "Specified function not found.\n");
            exit(1);
        }
    }
    else {
        fprintf(stderr, "No alignment function specified.\n");
        exit(1);
    }

    /* select the blosum matrix */
    if (blosumname) {
        blosum = lookup_blosum(blosumname);
        if (NULL == blosum) {
            fprintf(stderr, "Specified blosum matrix not found.\n");
            fprintf(stderr, "Choices are {"
                    "blosum40,"
                    "blosum45,"
                    "blosum50,"
                    "blosum62,"
                    "blosum75,"
                    "blosum80,"
                    "blosum90}\n");
            exit(1);
        }
    }

    if (filename_database) {
        (void)parse_sequences(filename_database, &sequences_database, &sizes_database, &seq_count_database);
    }
    else {
        fprintf(stderr, "missing database filename\n");
        exit(1);
    }

    if (filename_queries) {
        biggest = parse_sequences(filename_queries,
                &sequences_queries, &sizes_queries, &seq_count_queries);
    }
    else {
        fprintf(stderr, "missing query filename\n");
        exit(1);
    }

    char *seen = malloc(sizeof(char)*biggest);
    for (i=0; i<biggest; ++i) {
        seen[i] = 0;
    }

    for (i=0; i<seq_count_queries; ++i) {
        int saturated_query = 0;
#if ENABLE_CORRECTION_STATS
        unsigned long long corrections_query = 0;
#endif
        double local_timer = timer_real();
        if (truncate > 0 && sizes_queries[i] > (unsigned long)truncate) {
            continue;
        }
        if (exact_length > 0 && sizes_queries[i] != (unsigned long)exact_length) {
            continue;
        }
        if (seen[sizes_queries[i]]) {
            //printf("skipping %d\n", i);
            continue;
        }
        else {
            seen[sizes_queries[i]] = 1;
        }
        for (j=0; j<seq_count_database; ++j) {
            parasail_result_t *result = function(
                    sequences_queries[i], sizes_queries[i],
                    sequences_database[j], sizes_database[j],
                    gap_open, gap_extend, blosum);
            saturated_query += result->saturated;
#if ENABLE_CORRECTION_STATS
            corrections_query += result->corrections;
#endif
            parasail_result_free(result);
        }
        local_timer = timer_real() - local_timer;
        printf("%lu\t %lu\t %d\t %llu\t %f\n",
                i, sizes_queries[i],
                saturated_query,
#if ENABLE_CORRECTION_STATS
                corrections_query,
#else
                0ULL,
#endif
                local_timer);
        if (exact_length != 0) {
            /* if we got this far, we found our query, so break */
            break;
        }
    }

    return 0;
}

