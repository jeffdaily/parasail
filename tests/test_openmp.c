#include "config.h"

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

int main(int argc, char **argv)
{
    const char *seqA = "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR";
    const char *seqB = "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE";
    //const char *seqA = "MEFYDVAVTV";
    //const char *seqB = "AALGVAARAGFLAAGFASSS";
    //const char *seqA = rand_string(32000);
    //const char *seqB = rand_string(16000);
    const int lena = strlen(seqA);
    const int lenb = strlen(seqB);
    const int longest = (lena>lenb?lena:lenb);
    int score;
    double timer_clock;
    unsigned long long timer_rtdsc;
    int i = 0;
    int limit = 100;
    parasail_result_t **result = NULL;
    parasail_workspace_t **workspace = NULL;

    timer_init();
    printf("%s timer\n", timer_name());

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
#pragma omp for schedule(dynamic)
        for (i=0; i<limit; ++i) {
            nw(seqA, lena, seqB, lenb, 10, 1, blosum62,
                    result[thread_num]);
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
#pragma omp for schedule(dynamic)
        for (i=0; i<limit; ++i) {
            nw_ext(seqA, lena, seqB, lenb, 10, 1, blosum62,
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

