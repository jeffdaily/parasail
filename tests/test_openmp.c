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

int main(int argc, char **argv)
{
    const char *seqA = "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR";
    const char *seqB = "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE";
    //const char *seqA = "MEFYDVAVTV";
    //const char *seqB = "AALGVAARAGFLAAGFASSS";
    const int lena = strlen(seqA);
    const int lenb = strlen(seqB);
    const int longest = (lena>lenb?lena:lenb);
    int score;
    unsigned long long timer;
    unsigned long long timer_ref;
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

    timer_ref = timer_start();
#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
#pragma omp for schedule(dynamic)
        for (i=0; i<limit; ++i) {
            nw(seqA, lena, seqB, lenb, 10, 1, blosum62,
                    result[thread_num]);
        }
    }
    timer_ref = timer_end(timer_ref);
    printf("nw\t%llu\n", timer_ref);

    timer_ref = timer_start();
#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
#pragma omp for schedule(dynamic)
        for (i=0; i<limit; ++i) {
            nw_ext(seqA, lena, seqB, lenb, 10, 1, blosum62,
                    result[thread_num], workspace[thread_num]);
        }
    }
    timer_ref = timer_end(timer_ref);
    printf("nw_ext\t%llu\n", timer_ref);

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

