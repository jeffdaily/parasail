#include "config.h"

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <errno.h>
#include <pwd.h>
#include <unistd.h>
#include <sys/types.h>

#include "kseq.h"
KSEQ_INIT(int, read)

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/cpuid.h"
#include "parasail/function_lookup.h"
#include "parasail/stats.h"
#include "timer.h"
#include "timer_real.h"

static double pctf(double orig, double new)
{
    return orig / new;
}

#ifdef __MIC__
static const char *get_user_name()
{
    uid_t uid = geteuid();
    struct passwd *pw = getpwuid(uid);
    if (pw) {
        return pw->pw_name;
    }
    return "";
}
#endif

static void print_array(
        const char * filename_,
        const int * const restrict array,
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len)
{
    int i;
    int j;
    FILE *f = NULL;
#ifdef __MIC__
    const char *username = get_user_name();
    char filename[4096] = {0};
    strcat(filename, "/tmp/");
    if (username[0] != '\0') {
        strcat(filename, username);
        strcat(filename, "/");
    }
    strcat(filename, filename_);
#else
    const char *filename = filename_;
#endif
    f = fopen(filename, "w");
    if (NULL == f) {
        printf("fopen(\"%s\") error: %s\n", filename, strerror(errno));
        exit(-1);
    }
    fprintf(f, " ");
    for (j=0; j<s2Len; ++j) {
        fprintf(f, "%4c", s2[j]);
    }
    fprintf(f, "\n");
    for (i=0; i<s1Len; ++i) {
        fprintf(f, "%c", s1[i]);
        for (j=0; j<s2Len; ++j) {
            fprintf(f, "%4d", array[i*s2Len + j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

static void print_rowcol(
        const char * filename_,
        const int * const restrict row,
        const int * const restrict col,
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len)
{
    int i;
    int j;
    FILE *f = NULL;
#ifdef __MIC__
    const char *username = get_user_name();
    char filename[4096] = {0};
    strcat(filename, "/tmp/");
    if (username[0] != '\0') {
        strcat(filename, username);
        strcat(filename, "/");
    }
    strcat(filename, filename_);
#else
    const char *filename = filename_;
#endif
    f = fopen(filename, "w");
    if (NULL == f) {
        printf("fopen(\"%s\") error: %s\n", filename, strerror(errno));
        exit(-1);
    }
    fprintf(f, "%c", s1[s1Len-1]);
    if (NULL == row) {
        for (j=0; j<s2Len; ++j) {
            fprintf(f, "%4c", '!');
        }
    }
    else {
        for (j=0; j<s2Len; ++j) {
            fprintf(f, "%4d", row[j]);
        }
    }
    fprintf(f, "\n");
    fprintf(f, "%c", s2[s2Len-1]);
    if (NULL == col) {
        for (i=0; i<s1Len; ++i) {
            fprintf(f, "%4c", '!');
        }
    }
    else {
        for (i=0; i<s1Len; ++i) {
            fprintf(f, "%4d", col[i]);
        }
    }
    fprintf(f, "\n");
    fclose(f);
}

static inline void parse_sequences(
        const char *filename,
        char ***strings_,
        unsigned long **sizes_,
        unsigned long *count_)
{
    FILE* fp;
    kseq_t *seq = NULL;
    int l = 0;
    char **strings = NULL;
    unsigned long *sizes = NULL;
    unsigned long count = 0;
    unsigned long memory = 1000;

    errno = 0;
    fp = fopen(filename, "r");
    if(fp == NULL) {
        perror("fopen");
        exit(1);
    }
    strings = malloc(sizeof(char*) * memory);
    sizes = malloc(sizeof(unsigned long) * memory);
    seq = kseq_init(fileno(fp));
    while ((l = kseq_read(seq)) >= 0) {
        errno = 0;
        strings[count] = strdup(seq->seq.s);
        if (NULL == strings[count]) {
            perror("strdup");
            exit(1);
        }
        sizes[count] = seq->seq.l;
        ++count;
        if (count >= memory) {
            char **new_strings = NULL;
            unsigned long *new_sizes = NULL;
            memory *= 2;
            errno = 0;
            new_strings = realloc(strings, sizeof(char*) * memory);
            if (NULL == new_strings) {
                perror("realloc");
                exit(1);
            }
            strings = new_strings;
            errno = 0;
            new_sizes = realloc(sizes, sizeof(unsigned long) * memory);
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
}


int main(int argc, char **argv)
{
    long seqA_index = LONG_MAX;
    long seqB_index = LONG_MAX;
    const char *seqA = NULL;
    const char *seqB = NULL;
    int lena = 0;
    int lenb = 0;
    int score = 0;
    int matches = 0;
    int similar = 0;
    int length = 0;
    int end_query = 0;
    int end_ref = 0;
    unsigned long long timer_rdtsc = 0;
    unsigned long long timer_rdtsc_single = 0;
    double timer_nsecs = 0.0;
    double timer_nsecs_single = 0.0;
    double timer_nsecs_ref_mean = 0.0;
    double timer_rdtsc_ref_mean = 0.0;
    int limit = 2;
    int i = 0;
    int index = 0;
    parasail_function_info_t f;
    parasail_result_t *result = NULL;
    stats_t stats_rdtsc;
    stats_t stats_nsecs;
    int c = 0;
    char *filename = NULL;
    char **sequences = NULL;
    unsigned long *sizes = NULL;
    unsigned long seq_count = 0;
    unsigned long s = 0;
    char *endptr = NULL;
    char *matrixname = NULL;
    const parasail_matrix_t *matrix = NULL;
    int open = 10;
    int extend = 1;
    int do_normal = 1;
    int do_stats = 1;
    int do_nonstats = 1;
    int do_table = 1;
    int do_rowcol = 1;
    int use_rdtsc = 0;

    while ((c = getopt(argc, argv, "a:b:f:m:n:o:e:rRTNSs")) != -1) {
        switch (c) {
            case 'a':
                errno = 0;
                seqA_index = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol seqA_index");
                    fprintf(stderr, "invalid seqA index\n");
                    exit(1);
                }
                break;
            case 'b':
                errno = 0;
                seqB_index = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol seqB_index");
                    fprintf(stderr, "invalid seqB index\n");
                    exit(1);
                }
                break;
            case 'f':
                filename = optarg;
                break;
            case 'm':
                matrixname = optarg;
                break;
            case 'n':
                errno = 0;
                limit = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol limit");
                    exit(1);
                }
                break;
            case 'o':
                errno = 0;
                open = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol open");
                    exit(1);
                }
                break;
            case 'e':
                errno = 0;
                extend = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol extend");
                    exit(1);
                }
                break;
            case 'r':
                use_rdtsc = 1;
            case 'R':
                do_rowcol = 0;
                break;
            case 'T':
                do_table = 0;
                break;
            case 'N':
                do_normal = 0;
                break;
            case 'S':
                do_stats = 0;
                break;
            case 's':
                do_nonstats = 0;
                break;
            case '?':
                if (optopt == 'a'
                        || optopt == 'b'
                        || optopt == 'f'
                        || optopt == 'n') {
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
        fprintf(stderr, "no filename specified\n");
        exit(1);
    }

    /* select the matrix */
    if (matrixname) {
        matrix = parasail_matrix_lookup(matrixname);
        if (NULL == matrix) {
            fprintf(stderr, "Specified substitution matrix not found.\n");
            exit(1);
        }
    }
    else {
        fprintf(stderr, "No substitution matrix specified.\n");
        exit(1);
    }

    if (seqA_index == LONG_MAX) {
        fprintf(stderr, "seqA index not specified\n");
        exit(1);
    }

    if (seqB_index == LONG_MAX) {
        fprintf(stderr, "seqB index not specified\n");
        exit(1);
    }

    if ((unsigned long)seqA_index >= seq_count) {
        fprintf(stderr, "seqA index out of bounds\n");
        exit(1);
    }

    if ((unsigned long)seqB_index >= seq_count) {
        fprintf(stderr, "seqB index out of bounds\n");
        exit(1);
    }

    seqA = sequences[seqA_index];
    seqB = sequences[seqB_index];
    lena = strlen(seqA);
    lenb = strlen(seqB);

    printf("file: %s\n", filename);
    printf("matrix: %s\n", matrixname);
    printf("gap open: %d\n", open);
    printf("gap extend: %d\n", extend);
    printf("seq pair %lu,%lu\n", seqA_index, seqB_index);

    printf("%-15s %8s %6s %4s %5s %5s "
           "%8s %8s %8s %8s %8s %8s "
           "%8s %5s %8s %8s %8s\n",
            "name", "type", "isa", "bits", "width", "elem",
            "score", "matches", "similar", "length", "end_query", "end_ref",
            "avg", "imp", "stddev", "min", "max");

    stats_clear(&stats_rdtsc);
    stats_clear(&stats_nsecs);
    index = 0;
    f = functions[index++];
    while (f.pointer) {
        char name[16] = {'\0'};
        int new_limit = f.is_table ? 1 : limit;
        int saturated = 0;
#if 0
        if (f.is_table && HAVE_KNC) {
            f = functions[index++];
            continue;
        }
#endif
        if ((0 == strncmp(f.isa, "sse2",  4) && 0 == parasail_can_use_sse2()) 
                || (0 == strncmp(f.isa, "sse41", 5) && 0 == parasail_can_use_sse41())
                || (0 == strncmp(f.isa, "avx2",  4) && 0 == parasail_can_use_avx2())) {
            f = functions[index++];
            continue;
        }
        if (f.is_stats) {
            if (!do_stats) {
                f = functions[index++];
                continue;
            }
        }
        else {
            if (!do_nonstats) {
                f = functions[index++];
                continue;
            }
        }
        if (f.is_table) {
            if (!do_table) {
                f = functions[index++];
                continue;
            }
        }
        else if (f.is_rowcol) {
            if (!do_rowcol) {
                f = functions[index++];
                continue;
            }
        }
        else {
            if (!do_normal) {
                f = functions[index++];
                continue;
            }
        }
        stats_clear(&stats_rdtsc);
        timer_rdtsc = timer_start();
        timer_nsecs = timer_real();
        for (i=0; i<new_limit; ++i) {
            timer_rdtsc_single = timer_start();
            timer_nsecs_single = timer_real();
            result = f.pointer(seqA, lena, seqB, lenb, open, extend, matrix);
            timer_rdtsc_single = timer_start()-(timer_rdtsc_single);
            timer_nsecs_single = timer_real() - timer_nsecs_single;
            stats_sample_value(&stats_rdtsc, timer_rdtsc_single);
            stats_sample_value(&stats_nsecs, timer_nsecs_single);
            score = result->score;
            similar = result->similar;
            matches = result->matches;
            length = result->length;
            end_query = result->end_query;
            end_ref = result->end_ref;
            saturated = result->saturated;
            parasail_result_free(result);
        }
        timer_rdtsc = timer_start()-(timer_rdtsc);
        timer_nsecs = timer_real() - timer_nsecs;
        if (f.is_ref) {
            timer_nsecs_ref_mean = stats_nsecs._mean;
            timer_rdtsc_ref_mean = stats_rdtsc._mean;
        }
        strcpy(name, f.alg);
        /* xeon phi was unable to perform I/O running natively */
        if (f.is_table) {
            char suffix[256] = {0};
            if (strlen(f.type)) {
                strcat(suffix, "_");
                strcat(suffix, f.type);
            }
            if (strlen(f.isa)) {
                strcat(suffix, "_");
                strcat(suffix, f.isa);
            }
            if (strlen(f.bits)) {
                strcat(suffix, "_");
                strcat(suffix, f.bits);
            }
            if (strlen(f.width)) {
                strcat(suffix, "_");
                strcat(suffix, f.width);
            }
            strcat(suffix, ".txt");
            result = f.pointer(seqA, lena, seqB, lenb, open, extend, matrix);
            {
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_scr");
                strcat(filename, suffix);
                print_array(filename, result->score_table, seqA, lena, seqB, lenb);
            }
            if (f.is_stats) {
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_mch");
                strcat(filename, suffix);
                print_array(filename, result->matches_table, seqA, lena, seqB, lenb);
            }
            if (f.is_stats) {
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_sim");
                strcat(filename, suffix);
                print_array(filename, result->similar_table, seqA, lena, seqB, lenb);
            }
            if (f.is_stats) {
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_len");
                strcat(filename, suffix);
                print_array(filename, result->length_table, seqA, lena, seqB, lenb);
            }
            parasail_result_free(result);
        }
        else if (f.is_rowcol) {
            char suffix[256] = {0};
            if (strlen(f.type)) {
                strcat(suffix, "_");
                strcat(suffix, f.type);
            }
            if (strlen(f.isa)) {
                strcat(suffix, "_");
                strcat(suffix, f.isa);
            }
            if (strlen(f.bits)) {
                strcat(suffix, "_");
                strcat(suffix, f.bits);
            }
            if (strlen(f.width)) {
                strcat(suffix, "_");
                strcat(suffix, f.width);
            }
            strcat(suffix, ".txt");
            result = f.pointer(seqA, lena, seqB, lenb, open, extend, matrix);
            {
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_rowcol_scr");
                strcat(filename, suffix);
                print_rowcol(filename, result->score_row, result->score_col, seqA, lena, seqB, lenb);
            }
            if (f.is_stats) {
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_rowcol_mch");
                strcat(filename, suffix);
                print_rowcol(filename, result->matches_row, result->matches_col, seqA, lena, seqB, lenb);
            }
            if (f.is_stats) {
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_rowcol_sim");
                strcat(filename, suffix);
                print_rowcol(filename, result->similar_row, result->similar_col, seqA, lena, seqB, lenb);
            }
            if (f.is_stats) {
                char filename[256] = {'\0'};
                strcpy(filename, f.alg);
                strcat(filename, "_rowcol_len");
                strcat(filename, suffix);
                print_rowcol(filename, result->length_row, result->length_col, seqA, lena, seqB, lenb);
            }
            parasail_result_free(result);
        }
        if (f.is_table) {
            strcat(name, "_table");
        }
        else if (f.is_rowcol) {
            strcat(name, "_rowcol");
        }
        if (use_rdtsc) {
            printf(
                "%-15s %8s %6s %4s %5s %5d "
                "%8d %8d %8d %8d "
                "%8d %8d "
                "%8.1f %5.1f %8.1f %8.0f %8.0f\n",
                name, f.type, f.isa, f.bits, f.width, f.lanes,
                score, matches, similar, length,
                end_query, end_ref,
                saturated ? 0 : stats_rdtsc._mean,
                saturated ? 0 : pctf(timer_rdtsc_ref_mean, stats_rdtsc._mean),
                saturated ? 0 : stats_stddev(&stats_rdtsc),
                saturated ? 0 : stats_rdtsc._min,
                saturated ? 0 : stats_rdtsc._max);
        }
        else {
            printf(
                "%-15s %8s %6s %4s %5s %5d "
                "%8d %8d %8d %8d "
                "%8d %8d "
                "%8.3f %5.2f %8.3f %8.3f %8.3f\n",
                name, f.type, f.isa, f.bits, f.width, f.lanes,
                score, matches, similar, length,
                end_query, end_ref,
                saturated ? 0 : stats_nsecs._mean,
                saturated ? 0 : pctf(timer_nsecs_ref_mean, stats_nsecs._mean),
                saturated ? 0 : stats_stddev(&stats_nsecs),
                saturated ? 0 : stats_nsecs._min,
                saturated ? 0 : stats_nsecs._max);
        }
        f = functions[index++];
    }

    if (filename) {
        for (s=0; s<seq_count; ++s) {
            free(sequences[s]);
        }
        free(sequences);
        free(sizes);
    }

    return 0;
}
