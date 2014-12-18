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

#define USE_PERCENT_IMPROVED 0

static float pct(unsigned long long orig_, unsigned long long new_)
{
    float orig = (float)orig_;
    float new = (float)new_;
#if USE_PERCENT_IMPROVED
    return 100.0*(orig - new)/orig;
#else
    return orig / new;
#endif
}


int main(int argc, char **argv)
{
    const char *seqA = "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR";
    const char *seqB = "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE";
    //const char *seqA = "MEFYDVAVTV";
    //const char *seqB = "AALGVAARAGFLAAGFASSS";
    const int lena = strlen(seqA);
    const int lenb = strlen(seqB);
    const int longest = (lena>lenb?lena:lenb) + 32 /* +32 for woz padding */;
    int score;
    int matches;
    int length;
    unsigned long long timer;
    unsigned long long timer_ref;
    size_t limit = 1000;
    size_t i;
    parasail_result_t *result = NULL;

    timer_init();
    printf("%s timer\n", timer_name());
#if USE_PERCENT_IMPROVED
    printf("alg\ttype\tvec_w\tval_w\ttime\t%%imp\tscore\tmatches\tlength\n");
#else
    printf("alg\ttype\tvec_w\tval_w\ttime\tx_imp\tscore\tmatches\tlength\n");
#endif

    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        result = nw_stats(seqA, lena, seqB, lenb, 10, 1, blosum62);
        score = result->score;
        matches = result->matches;
        length = result->length;
        parasail_result_free(result);
    }
    timer_ref = timer_end(timer_ref);
    printf("nw\tref\t\t\t%llu\t\t%d\t%d\t%d\n", timer_ref/limit, score, matches, length);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result = nw_stats_scan(seqA, lena, seqB, lenb, 10, 1, blosum62);
        score = result->score;
        matches = result->matches;
        length = result->length;
        parasail_result_free(result);
    }
    timer = timer_end(timer);
    printf("nw\tscan\t\t\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);

#if 0
#if HAVE_SSE2
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats_scan_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length);
    }
    timer = timer_end(timer);
    printf("nw\tscan\t128\t16\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if HAVE_SSE41
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats_scan_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length);
    }
    timer = timer_end(timer);
    printf("nw\tscan\t128\t8\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if HAVE_SSE2
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats_wozniak_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer = timer_end(timer);
    printf("nw\twozniak\t128\t16\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if HAVE_SSE41
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats_wozniak_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer = timer_end(timer);
    printf("nw\twozniak\t128\t8\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if HAVE_SSE2
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats_striped_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length);
    }
    timer = timer_end(timer);
    printf("nw\tstriped\t128\t16\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if HAVE_SSE41
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = nw_stats_striped_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length);
    }
    timer = timer_end(timer);
    printf("nw\tstriped\t128\t8\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif
#endif

    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        result = sg_stats(seqA, lena, seqB, lenb, 10, 1, blosum62);
        score = result->score;
        matches = result->matches;
        length = result->length;
        parasail_result_free(result);
    }
    timer_ref = timer_end(timer_ref);
    printf("sg\tref\t\t\t%llu\t\t%d\t%d\t%d\n", timer_ref/limit, score, matches, length);

    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result = sg_stats_scan(seqA, lena, seqB, lenb, 10, 1, blosum62);
        score = result->score;
        matches = result->matches;
        length = result->length;
        parasail_result_free(result);
    }
    timer = timer_end(timer);
    printf("sg\tscan\t\t\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);

#if 0
#if HAVE_SSE2
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_stats_scan_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__,
                &matches, &length);
    }
    timer = timer_end(timer);
    printf("sg\tscan\t128\t16\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if HAVE_SSE41
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_stats_scan_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62__,
                &matches, &length);
    }
    timer = timer_end(timer);
    printf("sg\tscan\t128\t8\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if HAVE_SSE2
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_stats_wozniak_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer = timer_end(timer);
    printf("sg\twozniak\t128\t16\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if HAVE_SSE41
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_stats_wozniak_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer = timer_end(timer);
    printf("sg\twozniak\t128\t8\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if HAVE_SSE2
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sg_stats_striped_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__,
                &matches, &length);
    }
    timer = timer_end(timer);
    printf("sg\tstriped\t128\t16\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif
#endif

    timer_ref = timer_start();
    for (i=0; i<limit; ++i) {
        result = sw_stats(seqA, lena, seqB, lenb, 10, 1, blosum62);
        score = result->score;
        matches = result->matches;
        length = result->length;
        parasail_result_free(result);
    }
    timer_ref = timer_end(timer_ref);
    printf("sw\tref\t\t\t%llu\t\t%d\t%d\t%d\n", timer_ref/limit, score, matches, length);


    timer = timer_start();
    for (i=0; i<limit; ++i) {
        result = sw_stats_scan(seqA, lena, seqB, lenb, 10, 1, blosum62);
        score = result->score;
        matches = result->matches;
        length = result->length;
        parasail_result_free(result);
    }
    timer = timer_end(timer);
    printf("sw\tscan\t\t\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);

#if 0
#if HAVE_SSE2
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_stats_scan_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__,
                &matches, &length);
    }
    timer = timer_end(timer);
    printf("sw\tscan\t128\t16\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if HAVE_SSE41
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_stats_scan_128_8(seqA, lena, seqB, lenb, 10, 1, blosum62__,
                &matches, &length);
    }
    timer = timer_end(timer);
    printf("sw\tscan\t128\t8\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if HAVE_SSE2
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_stats_wozniak_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62,
                &matches, &length, tbl_pr, del_pr, mch_pr, len_pr);
    }
    timer = timer_end(timer);
    printf("sw\twozniak\t128\t16\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif

#if HAVE_SSE2
    timer = timer_start();
    for (i=0; i<limit; ++i) {
        score = sw_stats_striped_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62__,
                &matches, &length);
    }
    timer = timer_end(timer);
    printf("sw\tstriped\t128\t16\t%llu\t%4.1f\t%d\t%d\t%d\n", timer/limit, pct(timer_ref,timer), score, matches, length);
#endif
#endif

    return 0;
}