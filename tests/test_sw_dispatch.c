#include "config.h"

#include <stdint.h>
#include <stdio.h>
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
    const char *seqA = "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR"
                       "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR";
    const char *seqB = "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE"
                       "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE";
    size_t lena = strlen(seqA);
    size_t lenb = strlen(seqB);
    parasail_result_t * result = NULL;

    result = parasail_sw(seqA, lena, seqB, lenb, 10, 1, blosum62);
    printf("score=%d\n", result->score);
    parasail_result_free(result);

    result = parasail_sw_32(seqA, lena, seqB, lenb, 10, 1, blosum62);
    printf("score32=%d\n", result->score);
    parasail_result_free(result);

    result = parasail_sw_16(seqA, lena, seqB, lenb, 10, 1, blosum62);
    printf("score16=%d\n", result->score);
    parasail_result_free(result);

    result = parasail_sw_8(seqA, lena, seqB, lenb, 10, 1, blosum62);
    printf("score8=%d\n", result->score);
    parasail_result_free(result);

#ifdef __MIC__
    result = sw_table_diag_knc_512_32(seqA, lena, seqB, lenb, 10, 1, blosum62);
    printf("table score32=%d\n", result->score);
    print_array("sw_scr_diag_knc_512_32.txt", result->score_table, seqA, lena, seqB, lenb);
    parasail_result_free(result);
#endif

#ifdef __MIC__
    result = sw_stats_table_scan_knc_512_32(seqA, lena, seqB, lenb, 10, 1, blosum62);
    printf("table score32=%d\n", result->score);
    print_array("sw_stats_scr_diag_knc_512_32.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("sw_stats_mch_diag_knc_512_32.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("sw_stats_len_diag_knc_512_32.txt", result->length_table, seqA, lena, seqB, lenb);
    parasail_result_free(result);
#endif

    return 0;
}

