#include "config.h"

#include <math.h>
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parasail.h"
#include "blosum/blosum62.h"
#include "timer.h"
#include "timer_real.h"

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
    int longest = 0;
    double timer_clock = 0.0;
    unsigned long long timer_rtdsc = 0U;
    unsigned long i = 0;
    unsigned long seq_count = 10;
    unsigned long limit = binomial_coefficient(seq_count, 2);
    parasail_result_t **result = NULL;
    parasail_workspace_t **workspace = NULL;
    const char **sequences = NULL;
    int *sizes = NULL;

    timer_init();
    printf("%s timer\n", timer_name());

    /* generate 'seq_count' number of random strings */
    sequences = (const char**)malloc(sizeof(char*)*seq_count);
    sizes = (int*)malloc(sizeof(int)*seq_count);
    for (i=0; i<seq_count; ++i) {
        sizes[i] = (rand()%32767)+10;
        longest = sizes[i] > longest ? sizes[i] : longest;
        sequences[i] = rand_string(sizes[i]);
    }
    printf("done generating %lu random srings, longest is %d\n",
            seq_count, longest);
    printf("%lu choose 2 is %lu\n", seq_count, limit);

#pragma omp parallel
    {
#pragma omp single
        {
            int N = omp_get_max_threads();
            result = (parasail_result_t**)malloc(sizeof(parasail_result_t*)*N);
            workspace = (parasail_workspace_t**)malloc(sizeof(parasail_workspace_t*)*N);
            for (i=0; i<N; ++i) {
                result[i] = parasail_result_allocate(longest);
                workspace[i] = parasail_workspace_allocate(longest);
            }
            printf("omp_get_max_threads()=%d\n", N);
        }
    }

    timer_clock = timer_real();
    timer_rtdsc = timer_start();
#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        unsigned long a=0, b=1;
#pragma omp for schedule(dynamic) private(a,b)
        for (i=0; i<limit; ++i) {
            k_combination2(i, &a, &b);
            //printf("%lu=(%lu,%lu) sizes %d %d\n", i, a, b, sizes[a], sizes[b]);
            nw(sequences[a], sizes[a], sequences[b], sizes[b],
                    10, 1, blosum62, result[thread_num]);
        }
    }
    timer_rtdsc = timer_end(timer_rtdsc);
    timer_clock = timer_real() - timer_clock;
    printf("nw\t%llu\t%f\n", timer_rtdsc, timer_clock);

    timer_clock = timer_real();
    timer_rtdsc = timer_start();
#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        unsigned long a=0, b=1;
#pragma omp for schedule(dynamic)
        for (i=0; i<limit; ++i) {
            k_combination2(i, &a, &b);
            nw_ext(sequences[a], sizes[a], sequences[b], sizes[b],
                    10, 1, blosum62,
                    result[thread_num], workspace[thread_num]);
        }
    }
    timer_rtdsc = timer_end(timer_rtdsc);
    timer_clock = timer_real() - timer_clock;
    printf("nw_ext\t%llu\t%f\n", timer_rtdsc, timer_clock);

#pragma omp parallel
    {
#pragma omp single
        {
            int N = omp_get_max_threads();
            for (i=0; i<N; ++i) {
                parasail_result_free(result[i]);
                parasail_workspace_free(workspace[i]);
            }
            free(result);
            free(workspace);
        }
    }

    return 0;
}

