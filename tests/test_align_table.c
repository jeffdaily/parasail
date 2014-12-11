#include "config.h"

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parasail.h"
#include "blosum/blosum62.h"

static void print_array(
        const char * filename,
        const int * const restrict array,
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len)
{
    int i;
    int j;
    FILE *f = fopen(filename, "w");
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



int main(int argc, char **argv)
{
    const char *seqA = "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR";
    const char *seqB = "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE";
    //const char *seqA = "MEFYDVAVTV"
    //                   "MEFYDVAVTV";
    //const char *seqB = "AALGVAARAGFLAAGFASSS"
    //                   "AALGVAARAGFLAAGFASSS";
    const int lena = strlen(seqA);
    const int lenb = strlen(seqB);
    parasail_result_t *result = NULL;

    printf("alg\t\ttype\tvec_w\tval_w\tscore\tmatches\tlength\n");

    result = nw_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    printf("nw\t\t\t\t\t%d\t%d\t%d\n", result->score, result->matches, result->length);
    parasail_result_free(result);

    result = nw_scan_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_scan_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    printf("nw_scan\t\t\t\t\t%d\t%d\t%d\n", result->score, result->matches, result->length);
    parasail_result_free(result);

#if HAVE_SSE2
    result = nw_table_scan_sse2_128_16(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_scr_scan_128_16.txt", result->score_table, seqA, lena, seqB, lenb);
    printf("nw\t\tscan\t128\t16\t%d\t%d\t%d\n", result->score, result->matches, result->length);
    parasail_result_free(result);
#endif

#if 0
#if HAVE_SSE41
    score = nw_scan_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, result->score_table);
    print_array("nw_scan_128_8_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_AVX_512
    score = nw_scan_512_32_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, result->score_table);
    print_array("nw_scan_512_32_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE2
    score = nw_wozniak_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_wozniak_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE41
    score = nw_wozniak_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_wozniak_128_8_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE2
    score = nw_striped_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, result->score_table);
    print_array("nw_striped_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE41
    score = nw_striped_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, result->score_table);
    print_array("nw_striped_128_8_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif
#endif

    result = nw_stats_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_stats_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_len.txt", result->length_table, seqA, lena, seqB, lenb);
    printf("nw_stats\t\t\t\t%d\t%d\t%d\n", result->score, result->matches, result->length);
    parasail_result_free(result);

    result = nw_scan_stats_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_scan_stats_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("nw_scan_stats_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("nw_scan_stats_len.txt", result->length_table, seqA, lena, seqB, lenb);
    printf("nw_scan_stats\t\t\t\t%d\t%d\t%d\n", result->score, result->matches, result->length);
    parasail_result_free(result);

#if 0
#if HAVE_SSE2
    score = nw_stats_scan_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, result->score_table, result->matches_table, result->length_table);
    print_array("nw_stats_scan_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_scan_128_16_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_scan_128_16_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE41
    score = nw_stats_scan_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, result->score_table, result->matches_table, result->length_table);
    print_array("nw_stats_scan_128_8_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_scan_128_8_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_scan_128_8_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE2
    score = nw_stats_wozniak_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_stats_wozniak_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_wozniak_128_16_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_wozniak_128_16_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE41
    score = nw_stats_wozniak_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_stats_wozniak_128_8_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_wozniak_128_8_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_wozniak_128_8_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE2
    score = nw_stats_striped_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, result->score_table, result->matches_table, result->length_table);
    print_array("nw_stats_striped_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_striped_128_16_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_striped_128_16_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE41
    score = nw_stats_striped_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, result->score_table, result->matches_table, result->length_table);
    print_array("nw_stats_striped_128_8_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_striped_128_8_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_striped_128_8_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif
#endif

    result = sg_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    printf("sg\t\t\t\t\t%d\t%d\t%d\n", result->score, result->matches, result->length);
    parasail_result_free(result);

    result = sg_scan_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_scan_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    printf("sg_scan\t\t\t\t\t%d\t%d\t%d\n", result->score, result->matches, result->length);
    parasail_result_free(result);


#if 0
#if HAVE_SSE2
    score = sg_scan_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, result->score_table);
    print_array("sg_scan_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE41
    score = sg_scan_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, result->score_table);
    print_array("sg_scan_128_8_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE2
    score = sg_wozniak_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_wozniak_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE41
    score = sg_wozniak_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_wozniak_128_8_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE2
    score = sg_striped_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, result->score_table);
    print_array("sg_striped_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE41
    score = sg_striped_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, result->score_table);
    print_array("sg_striped_128_8_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif
#endif

    result = sg_stats_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_stats_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("sg_stats_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("sg_stats_len.txt", result->length_table, seqA, lena, seqB, lenb);
    printf("sg_stats\t\t\t\t%d\t%d\t%d\n", result->score, result->matches, result->length);
    parasail_result_free(result);

    result = sg_scan_stats_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_scan_stats_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("sg_scan_stats_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("sg_scan_stats_len.txt", result->length_table, seqA, lena, seqB, lenb);
    printf("sg_scan_stats\t\t\t\t%d\t%d\t%d\n", result->score, result->matches, result->length);
    parasail_result_free(result);

#if 0
#if HAVE_SSE2
    score = sg_stats_scan_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, result->score_table, result->matches_table, result->length_table);
    print_array("sg_stats_scan_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("sg_stats_scan_128_16_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("sg_stats_scan_128_16_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE41
    score = sg_stats_scan_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, result->score_table, result->matches_table, result->length_table);
    print_array("sg_stats_scan_128_8_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("sg_stats_scan_128_8_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("sg_stats_scan_128_8_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE2
    score = sg_stats_wozniak_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_stats_wozniak_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("sg_stats_wozniak_128_16_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("sg_stats_wozniak_128_16_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE41
    score = sg_stats_wozniak_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_stats_wozniak_128_8_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("sg_stats_wozniak_128_8_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("sg_stats_wozniak_128_8_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE2
    score = sg_stats_striped_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, result->score_table, result->matches_table, result->length_table);
    print_array("sg_stats_striped_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("sg_stats_striped_128_16_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("sg_stats_striped_128_16_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif
#endif

    result = sw_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sw_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    printf("sw\t\t\t\t\t%d\t%d\t%d\n", result->score, result->matches, result->length);
    parasail_result_free(result);

    result = sw_scan_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sw_scan_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    printf("sw_scan\t\t\t\t\t%d\t%d\t%d\n", result->score, result->matches, result->length);
    parasail_result_free(result);

#if 0
#if HAVE_SSE2
    score = sw_scan_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, result->score_table);
    print_array("sw_scan_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE41
    score = sw_scan_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, result->score_table);
    print_array("sw_scan_128_8_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE2
    score = sw_wozniak_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sw_wozniak_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE2
    score = sw_striped_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, result->score_table);
    print_array("sw_striped_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
#endif
#endif

    result = sw_stats_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sw_stats_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("sw_stats_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("sw_stats_len.txt", result->length_table, seqA, lena, seqB, lenb);
    printf("sw_stats\t\t\t\t%d\t%d\t%d\n", result->score, result->matches, result->length);
    parasail_result_free(result);

    result = sw_scan_stats_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sw_scan_stats_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("sw_scan_stats_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("sw_scan_stats_len.txt", result->length_table, seqA, lena, seqB, lenb);
    printf("sw_scan_stats\t\t\t\t%d\t%d\t%d\n", result->score, result->matches, result->length);
    parasail_result_free(result);

#if 0
#if HAVE_SSE2
    score = sw_stats_scan_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, result->score_table, result->matches_table, result->length_table);
    print_array("sw_stats_scan_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("sw_stats_scan_128_16_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("sw_stats_scan_128_16_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE41
    score = sw_stats_scan_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, result->score_table, result->matches_table, result->length_table);
    print_array("sw_stats_scan_128_8_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("sw_stats_scan_128_8_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("sw_stats_scan_128_8_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE2
    score = sw_stats_wozniak_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sw_stats_wozniak_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("sw_stats_wozniak_128_16_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("sw_stats_wozniak_128_16_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif

#if HAVE_SSE2
    score = sw_stats_striped_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, result->score_table, result->matches_table, result->length_table);
    print_array("sw_stats_striped_128_16_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("sw_stats_striped_128_16_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("sw_stats_striped_128_16_len.txt", result->length_table, seqA, lena, seqB, lenb);
#endif
#endif

    return 0;
}
