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

#include "kseq.h"
KSEQ_INIT(int, read)

#include "parasail.h"
#include "parasail/cpuid.h"
#include "parasail/memory.h"
#include "parasail/matrix_lookup.h"

#include "func_verify.h"

static int verbose = 0;

typedef struct gap_score {
    int open;
    int extend;
} gap_score_t;

gap_score_t gap_scores[] = {
    {10,1},
    {10,2},
    {14,2},
    {40,2},
    {INT_MIN,INT_MIN}
};

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

    fp = fopen(filename, "r");
    if(fp == NULL) {
        perror("fopen");
        exit(1);
    }
    strings = malloc(sizeof(char*) * memory);
    sizes = malloc(sizeof(unsigned long) * memory);
    seq = kseq_init(fileno(fp));
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
            unsigned long *new_sizes = NULL;
            memory *= 2;
            new_strings = realloc(strings, sizeof(char*) * memory);
            if (NULL == new_strings) {
                perror("realloc");
                exit(1);
            }
            strings = new_strings;
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

static inline unsigned long binomial_coefficient(
        unsigned long n,
        unsigned long k)
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

static inline void k_combination2(
        unsigned long pos,
        unsigned long *a,
        unsigned long *b)
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

static void check_functions(
        parasail_function_group_t f,
        char **sequences,
        unsigned long *sizes,
        unsigned long pair_limit,
        const parasail_matrix_t *matrix_,
        gap_score_t gap)
{
    const parasail_function_info_t *functions = f.fs;
    unsigned long matrix_index = 0;
    unsigned long gap_index = 0;
    unsigned long function_index = 0;
    unsigned long pair_index = 0;
    parasail_function_t *reference_function = NULL;
    const parasail_matrix_t ** matrices = parasail_matrices;
    const parasail_matrix_t * single_matrix[] = {
        matrix_,
        NULL
    };

    if (NULL != matrix_) {
        matrices = single_matrix;
    }

    printf("checking %s functions\n", f.name);
    for (matrix_index=0; NULL!=matrices[matrix_index]; ++matrix_index) {
        const parasail_matrix_t *matrix = matrices[matrix_index];
        const char *matrixname = matrix->name;
        if (verbose) printf("\t%s\n", matrixname);
        for (gap_index=0; INT_MIN!=gap_scores[gap_index].open; ++gap_index) {
            int open = gap_scores[gap_index].open;
            int extend = gap_scores[gap_index].extend;
            if (gap.open != INT_MIN && gap.extend != INT_MIN) {
                open = gap.open;
                extend = gap.extend;
            }
            if (verbose) printf("\t\topen=%d extend=%d\n", open, extend);
            reference_function = functions[0].pointer;
            for (function_index=1;
                    NULL!=functions[function_index].pointer;
                    ++function_index) {
                if (verbose) printf("\t\t\t%s\n", functions[function_index].name);
                unsigned long saturated = 0;
#pragma omp parallel for
                for (pair_index=0; pair_index<pair_limit; ++pair_index) {
                    parasail_result_t *reference_result = NULL;
                    parasail_result_t *result = NULL;
                    unsigned long a = 0;
                    unsigned long b = 1;
                    k_combination2(pair_index, &a, &b);
                    //printf("\t\t\t\tpair=%lu (%lu,%lu)\n", pair_index, a, b);
                    reference_result = reference_function(
                            sequences[a], sizes[a],
                            sequences[b], sizes[b],
                            open, extend,
                            matrix);
                    result = functions[function_index].pointer(
                            sequences[a], sizes[a],
                            sequences[b], sizes[b],
                            open, extend,
                            matrix);
                    if (result->saturated) {
                        /* no point in comparing a result that saturated */
                        parasail_result_free(reference_result);
                        parasail_result_free(result);
#pragma omp atomic
                        saturated += 1;
                        continue;
                    }
                    if (reference_result->score != result->score) {
#pragma omp critical(printer)
                        {
                            printf("%s(%lu,%lu,%d,%d,%s) wrong score (%d!=%d)\n",
                                    functions[function_index].name,
                                    a, b, open, extend,
                                    matrixname,
                                    reference_result->score, result->score);
                        }
                    }
                    if (reference_result->end_query != result->end_query) {
#pragma omp critical(printer)
                        {
                            printf("%s(%lu,%lu,%d,%d,%s) wrong end_query (%d!=%d)\n",
                                    functions[function_index].name,
                                    a, b, open, extend,
                                    matrixname,
                                    reference_result->end_query, result->end_query);
                        }
                    }
                    if (reference_result->end_ref != result->end_ref) {
#pragma omp critical(printer)
                        {
                            printf("%s(%lu,%lu,%d,%d,%s) wrong end_ref (%d!=%d)\n",
                                    functions[function_index].name,
                                    a, b, open, extend,
                                    matrixname,
                                    reference_result->end_ref, result->end_ref);
                        }
                    }
                    if (reference_result->matches != result->matches) {
#pragma omp critical(printer)
                        {
                            printf("%s(%lu,%lu,%d,%d,%s) wrong matches (%d!=%d)\n",
                                    functions[function_index].name,
                                    a, b, open, extend,
                                    matrixname,
                                    reference_result->matches, result->matches);
                        }
                    }
                    if (reference_result->similar != result->similar) {
#pragma omp critical(printer)
                        {
                            printf("%s(%lu,%lu,%d,%d,%s) wrong similar (%d!=%d)\n",
                                    functions[function_index].name,
                                    a, b, open, extend,
                                    matrixname,
                                    reference_result->similar, result->similar);
                        }
                    }
                    if (reference_result->length != result->length) {
#pragma omp critical(printer)
                        {
                            printf("%s(%lu,%lu,%d,%d,%s) wrong length (%d!=%d)\n",
                                    functions[function_index].name,
                                    a, b, open, extend,
                                    matrixname,
                                    reference_result->length, result->length);
                        }
                    }
                    parasail_result_free(reference_result);
                    parasail_result_free(result);
                }
                if (verbose && saturated) {
                    printf("%s %d %d %s saturated %lu times\n",
                            functions[function_index].name,
                            open, extend,
                            matrixname,
                            saturated);
                }
            }
            if (gap.open != INT_MIN && gap.extend != INT_MIN) {
                /* user-specified gap, don't loop */
                break;
            }
        }
    }
}

int main(int argc, char **argv)
{
    unsigned long i = 0;
    unsigned long seq_count = 0;
    unsigned long limit = 0;
    char **sequences = NULL;
    unsigned long *sizes = NULL;
    char *endptr = NULL;
    char *filename = NULL;
    int c = 0;
    int test_scores = 1;
    int test_stats = 0;
    char *matrixname = NULL;
    const parasail_matrix_t *matrix = NULL;
    gap_score_t gap = {INT_MIN,INT_MIN};
    int do_sse2 = 1;
    int do_sse41 = 1;
    int do_avx2 = 1;
    int do_disp = 1;
    int do_nw = 1;
    int do_sg = 1;
    int do_sw = 1;

    while ((c = getopt(argc, argv, "f:m:n:o:e:vsSi:")) != -1) {
        switch (c) {
            case 'f':
                filename = optarg;
                break;
            case 'm':
                matrixname = optarg;
                break;
            case 'n':
                errno = 0;
                seq_count = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case 'o':
                errno = 0;
                gap.open = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol gap.open");
                    exit(1);
                }
                break;
            case 'e':
                errno = 0;
                gap.extend = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol gap.extend");
                    exit(1);
                }
                break;
            case 'v':
                verbose = 1;
                break;
            case 's':
                test_stats = 1;
                break;
            case 'S':
                test_scores = 0;
                break;
            case 'i':
                do_sse2 = (NULL == strstr(optarg, "sse2"));
                do_sse41 = (NULL == strstr(optarg, "sse41"));
                do_avx2 = (NULL == strstr(optarg, "avx2"));
                do_disp = (NULL == strstr(optarg, "disp"));
                do_nw = (NULL == strstr(optarg, "nw"));
                do_sg = (NULL == strstr(optarg, "sg"));
                do_sw = (NULL == strstr(optarg, "sw"));
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

    limit = binomial_coefficient(seq_count, 2);
    printf("%lu choose 2 is %lu\n", seq_count, limit);


#if HAVE_SSE2
    if (do_sse2 && parasail_can_use_sse2()) {
        if (test_scores) {
            if (do_nw) check_functions(parasail_nw_sse2, sequences, sizes, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_sse2, sequences, sizes, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_sse2, sequences, sizes, limit, matrix, gap);
        }
        if (test_stats) {
            if (do_nw) check_functions(parasail_nw_stats_sse2, sequences, sizes, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_stats_sse2, sequences, sizes, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_stats_sse2, sequences, sizes, limit, matrix, gap);
        }
    }
#endif

#if HAVE_SSE41
    if (do_sse41 && parasail_can_use_sse41()) {
        if (test_scores) {
            if (do_nw) check_functions(parasail_nw_sse41, sequences, sizes, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_sse41, sequences, sizes, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_sse41, sequences, sizes, limit, matrix, gap);
        }
        if (test_stats) {
            if (do_nw) check_functions(parasail_nw_stats_sse41, sequences, sizes, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_stats_sse41, sequences, sizes, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_stats_sse41, sequences, sizes, limit, matrix, gap);
        }
    }
#endif

#if HAVE_AVX2
    if (do_avx2 && parasail_can_use_avx2()) {
        if (test_scores) {
            if (do_nw) check_functions(parasail_nw_avx2, sequences, sizes, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_avx2, sequences, sizes, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_avx2, sequences, sizes, limit, matrix, gap);
        }
        if (test_stats) {
            if (do_nw) check_functions(parasail_nw_stats_avx2, sequences, sizes, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_stats_avx2, sequences, sizes, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_stats_avx2, sequences, sizes, limit, matrix, gap);
        }
    }
#endif

#if HAVE_KNC
    {
        if (test_scores) {
            if (do_nw) check_functions(parasail_nw_knc, sequences, sizes, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_knc, sequences, sizes, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_knc, sequences, sizes, limit, matrix, gap);
        }
        if (test_stats) {
            if (do_nw) check_functions(parasail_nw_stats_knc, sequences, sizes, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_stats_knc, sequences, sizes, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_stats_knc, sequences, sizes, limit, matrix, gap);
        }
    }
#endif

    if (do_disp) {
        if (test_scores) {
            if (do_nw) check_functions(parasail_nw_disp, sequences, sizes, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_disp, sequences, sizes, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_disp, sequences, sizes, limit, matrix, gap);
        }
        if (test_stats) {
            if (do_nw) check_functions(parasail_nw_stats_disp, sequences, sizes, limit, matrix, gap);
            if (do_sg) check_functions(parasail_sg_stats_disp, sequences, sizes, limit, matrix, gap);
            if (do_sw) check_functions(parasail_sw_stats_disp, sequences, sizes, limit, matrix, gap);
        }
    }

    for (i=0; i<seq_count; ++i) {
        free(sequences[i]);
    }
    free(sequences);
    free(sizes);

    return 0;
}

