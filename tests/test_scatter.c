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

#if HAVE_SSE2
#include "ssw.h"
#endif

#include "parasail.h"
#include "parasail/function_lookup.h"
//#include "timer.h"
#include "timer_real.h"

/* This table is used to transform amino acid letters into numbers. */
static const int8_t table[128] = {
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
    14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
    23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
    14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
};

static inline void parse_sequences(
        const char *filename, char ***strings_, size_t **sizes_, size_t *count_)
{
    FILE* fp;
    kseq_t *seq = NULL;
    int l = 0;
    char **strings = NULL;
    size_t *sizes = NULL;
    size_t count = 0;
    size_t memory = 1000;

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
    unsigned long shortest = INT_MAX;
    unsigned long longest = 0;
    double timer_clock = 0.0;
    unsigned long i = 0;
    size_t seq_count = 10;
    size_t limit = 0;
    char **sequences = NULL;
    size_t *sizes = NULL;
    char *endptr = NULL;
    char *funcname = NULL;
    parasail_function_t *function = NULL;
    int lanes = 1;
    char *filename = NULL;
    int c = 0;
    const char *matrixname = "blosum62";
    const parasail_matrix_t *matrix = NULL;
    int gap_open = 10;
    int gap_extend = 1;
    int saturated = 0;

    while ((c = getopt(argc, argv, "a:b:f:n:o:e:")) != -1) {
        switch (c) {
            case 'a':
                funcname = optarg;
                break;
            case 'b':
                matrixname = optarg;
                break;
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

    /* select the function */
    if (funcname) {
        int index = 0;
        parasail_function_info_t f;
        f = functions[index++];
        while (f.pointer) {
            if (0 == strcmp(funcname, f.name)) {
                function = f.pointer;
                lanes = f.lanes;
                break;
            }
            f = functions[index++];
        }
        if (NULL == function) {
            fprintf(stderr, "Specified function not found.\n");
            exit(1);
        }
    }
    else {
        fprintf(stderr, "No alignment function specified.\n");
        exit(1);
    }

    /* select the substitution matrix */
    if (matrixname) {
        matrix = parasail_matrix_lookup(matrixname);
        if (NULL == matrix) {
            fprintf(stderr, "Specified substitution matrix not found.\n");
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
    }

    limit = binomial_coefficient(seq_count, 2);

    printf("size_A,segLen,size_B,score,matches,similar,length,corrections,cells,time,");
    printf("A_a,R_a,N_a,D_a,C_a,Q_a,E_a,G_a,H_a,I_a,L_a,K_a,M_a,F_a,P_a,S_a,T_a,W_a,Y_a,V_a,B_a,Z_a,X_a,NA_a,");
    printf("A_b,R_b,N_b,D_b,C_b,Q_b,E_b,G_b,H_b,I_b,L_b,K_b,M_b,F_b,P_b,S_b,T_b,W_b,Y_b,V_b,B_b,Z_b,X_b,NA_b,");
    printf("CUPS\n");

    timer_clock = timer_real();
#pragma omp parallel
    {
        unsigned long a=0;
        unsigned long b=1;
        double timer_local = 0.0;
        unsigned long a_counts[24];
        unsigned long b_counts[24];
        unsigned long j;
#pragma omp for schedule(dynamic)
        for (i=0; i<limit; ++i) {
            parasail_result_t *result = NULL;
            k_combination2(i, &a, &b);
            timer_local = timer_real();
            result = function(sequences[a], sizes[a], sequences[b], sizes[b],
                    gap_open, gap_extend, matrix);
            timer_local = timer_real() - timer_local;
            for (j=0; j<24; ++j) {
                a_counts[j] = 0;
                b_counts[j] = 0;
            }
            for (j=0; j<sizes[a]; ++j) {
                a_counts[table[(unsigned)sequences[a][j]]] += 1;
            }
            for (j=0; j<sizes[b]; ++j) {
                b_counts[table[(unsigned)sequences[b][j]]] += 1;
            }
#pragma omp critical
            printf("%lu,%lu,%lu,%d,%d,%d,%d,%lu,%f",
                    (unsigned long)sizes[a],
                    (unsigned long)(sizes[a]+lanes-1)/lanes,
                    (unsigned long)sizes[b],
                    result->score, result->matches,
                    result->similar, result->length,
                    (unsigned long)(sizes[a]*sizes[b]),
                    timer_local);
            for (j=0; j<24; ++j) {
                //printf(",%lu", a_counts[j]);
                printf(",%f", (double)(a_counts[j])/sizes[a]);
            }
            for (j=0; j<24; ++j) {
                //printf(",%lu", b_counts[j]);
                printf(",%f", (double)(b_counts[j])/sizes[b]);
            }
            printf(",%f\n", sizes[a]*sizes[b]/timer_local);
#pragma omp atomic
            saturated += result->saturated;
            parasail_result_free(result);
        }
    }
    timer_clock = timer_real() - timer_clock;

    return 0;
}

