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

    result = nw_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_ref_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    parasail_result_free(result);

#if 0
    score = nw_scan_row_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_scan_row_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

    score = nw_scan_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_scan_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

#if HAVE_SSE2
    score = nw_scan_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("nw_scan_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE41
    score = nw_scan_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("nw_scan_128_8_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_AVX_512
    score = nw_scan_512_32_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("nw_scan_512_32_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE2
    score = nw_wozniak_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_wozniak_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE41
    score = nw_wozniak_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_wozniak_128_8_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE2
    score = nw_striped_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("nw_striped_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE41
    score = nw_striped_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("nw_striped_128_8_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif
#endif

    result = nw_stats_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_stats_ref_scr.txt", result->score_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_ref_mch.txt", result->matches_table, seqA, lena, seqB, lenb);
    print_array("nw_stats_ref_len.txt", result->length_table, seqA, lena, seqB, lenb);
    parasail_result_free(result);

#if 0
    score = nw_stats_scan_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_stats_scan_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_scan_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_scan_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);

#if HAVE_SSE2
    score = nw_stats_scan_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, scr_tbl, mch_tbl, len_tbl);
    print_array("nw_stats_scan_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_scan_128_16_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_scan_128_16_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE41
    score = nw_stats_scan_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, scr_tbl, mch_tbl, len_tbl);
    print_array("nw_stats_scan_128_8_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_scan_128_8_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_scan_128_8_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE2
    score = nw_stats_wozniak_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_stats_wozniak_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_wozniak_128_16_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_wozniak_128_16_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE41
    score = nw_stats_wozniak_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("nw_stats_wozniak_128_8_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_wozniak_128_8_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_wozniak_128_8_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE2
    score = nw_stats_striped_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, scr_tbl, mch_tbl, len_tbl);
    print_array("nw_stats_striped_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_striped_128_16_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_striped_128_16_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE41
    score = nw_stats_striped_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, scr_tbl, mch_tbl, len_tbl);
    print_array("nw_stats_striped_128_8_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_striped_128_8_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("nw_stats_striped_128_8_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif

    score = sg_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_ref_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

    score = sg_scan_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_scan_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

#if HAVE_SSE2
    score = sg_scan_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("sg_scan_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE41
    score = sg_scan_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("sg_scan_128_8_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE2
    score = sg_wozniak_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_wozniak_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE41
    score = sg_wozniak_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_wozniak_128_8_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE2
    score = sg_striped_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("sg_striped_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE41
    score = sg_striped_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("sg_striped_128_8_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

    score = sg_stats_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_stats_ref_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_ref_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_ref_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);

    score = sg_stats_scan_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_stats_scan_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_scan_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_scan_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);

#if HAVE_SSE2
    score = sg_stats_scan_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, scr_tbl, mch_tbl, len_tbl);
    print_array("sg_stats_scan_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_scan_128_16_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_scan_128_16_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE41
    score = sg_stats_scan_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, scr_tbl, mch_tbl, len_tbl);
    print_array("sg_stats_scan_128_8_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_scan_128_8_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_scan_128_8_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE2
    score = sg_stats_wozniak_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_stats_wozniak_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_wozniak_128_16_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_wozniak_128_16_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE41
    score = sg_stats_wozniak_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sg_stats_wozniak_128_8_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_wozniak_128_8_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_wozniak_128_8_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE2
    score = sg_stats_striped_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, scr_tbl, mch_tbl, len_tbl);
    print_array("sg_stats_striped_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_striped_128_16_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sg_stats_striped_128_16_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif

    score = sw_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sw_ref_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

    score = sw_scan_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sw_scan_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);

#if HAVE_SSE2
    score = sw_scan_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("sw_scan_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE41
    score = sw_scan_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("sw_scan_128_8_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE2
    score = sw_wozniak_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sw_wozniak_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE2
    score = sw_striped_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, scr_tbl);
    print_array("sw_striped_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
#endif

    score = sw_stats_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sw_stats_ref_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_ref_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_ref_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);

    score = sw_stats_scan_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sw_stats_scan_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_scan_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_scan_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);

#if HAVE_SSE2
    score = sw_stats_scan_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, scr_tbl, mch_tbl, len_tbl);
    print_array("sw_stats_scan_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_scan_128_16_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_scan_128_16_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE41
    score = sw_stats_scan_128_8_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, scr_tbl, mch_tbl, len_tbl);
    print_array("sw_stats_scan_128_8_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_scan_128_8_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_scan_128_8_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE2
    score = sw_stats_wozniak_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62);
    print_array("sw_stats_wozniak_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_wozniak_128_16_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_wozniak_128_16_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif

#if HAVE_SSE2
    score = sw_stats_striped_128_16_table(seqA, lena, seqB, lenb, 10, 1, blosum62__, &matches, &length, scr_tbl, mch_tbl, len_tbl);
    print_array("sw_stats_striped_128_16_scr.txt", scr_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_striped_128_16_mch.txt", mch_tbl, seqA, lena, seqB, lenb);
    print_array("sw_stats_striped_128_16_len.txt", len_tbl, seqA, lena, seqB, lenb);
    memset(scr_tbl, 0, sizeof(int)*tbl_size);
    memset(mch_tbl, 0, sizeof(int)*tbl_size);
    memset(len_tbl, 0, sizeof(int)*tbl_size);
#endif
#endif

    return 0;
}
