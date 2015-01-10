#include "config.h"

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parasail.h"
#include "parasail_cpuid.h"
#include "blosum/blosum62.h"
#include "stats.h"
#include "timer.h"
#include "timer_real.h"

#define USE_TIMER_REAL 0

static double pctull(unsigned long long orig_, unsigned long long new_)
{
    double orig = (double)orig_;
    double new = (double)new_;
    return orig / new;
}

static double pctf(double orig, double new)
{
    return orig / new;
}

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

typedef struct func {
    parasail_result_t* (*f)(const char * const restrict s1, const int s1Len,
                            const char * const restrict s2, const int s2Len,
                            const int open, const int gap,
                            const int matrix[24][24]);
    const char * alg;
    const char * type;
    const char * isa;
    const char * bits;
    const char * width;
    char is_table;
    char is_stats;
    char is_ref;
} func_t;


static inline int elem(func_t f) {
    int i_bits = atoi(f.bits);
    int i_width = atoi(f.width);
    return i_bits / i_width;
}


int main(int argc, char **argv)
{
#if 0
    const char *seqA = "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR"
                       "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR"
                       "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR"
                       "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR";
    const char *seqB = "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE"
                       "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE"
                       "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE"
                       "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE";
#endif
#if 1
    const char *seqA = "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR"
                       "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR";
    const char *seqB = "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE"
                       "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE";
#endif
#if 0
    const char *seqA = "MEFYDVAVTVGMLCIIIYLLLVRQFRYWTERNVPQLNPHLLFGDVRDVNKTHHIGEKFRQLYNELKGKHPFGGIYMFTKPVALVTDLELVKNVFVKDFQYFHDRGTYYDEKHDPLSAHLFNLEGYKWKSLRNKITPTFTSGKMKMMFPTVAAAGKQFKDYLEDAIGEQEEFELKELLARYTTDVIGTCAFGIECNSMRNPNAEFRVMGKKIFGRSRSNLQLLLMNAFPSVAKLVGIKLILPEVSDFFMNAVRDTIKYRVENNVQRNDFMDILIRMRSDKETKSDDGTLTFHEIAAQAFVFFVAGFETSSSLMAFTLYELALDQDMQDKARKCVTDVLERHNGELTYEAAMEMDYLDCVLKGWVR";
    const char *seqB = "AALGVAARAGFLAAGFASSSELSSELSSEDSAAFLAAAAGVAAFAGVFTIAAFGVAATADLLAAGLHSSSELSSELSSEDSAAFFAATAGVAALAGVLAAAAAFGVAATADFFAAGLESSSELSSELSSDDSAVFFAAAAGVATFAGVLAAAATFGVAACAGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVATFTGVLAAAAACAAAACVGFFAAGLDSSSELSSELSSEDSAAFFAAAAGVAALAGVLAAAAACAGFFAAGLESSSELSSE";
#endif
#if 0
    const char *seqA = "MEFYDVAVTV"
                       "MEFYDVAVTV"
                       "MEFYDVAVTV"
                       "MEFYDVAVTV";
    const char *seqB = "AALGVAARAGFLAAGFASSS"
                       "AALGVAARAGFLAAGFASSS"
                       "AALGVAARAGFLAAGFASSS"
                       "AALGVAARAGFLAAGFASSS";
#endif
    const int lena = strlen(seqA);
    const int lenb = strlen(seqB);
    int score;
    int matches;
    int length;
    unsigned long long timer_rdtsc;
    unsigned long long timer_rdtsc_single;
    double timer_rdtsc_ref_mean;
    double timer_nsecs;
    double timer_nsecs_single;
    double timer_nsecs_ref_mean;
    int limit = 100;
    int i;
    int index;
    func_t f;
    parasail_result_t *result = NULL;
    stats_t stats_rdtsc;
    stats_t stats_nsecs;

    if (argc > 1) { limit = atoi(argv[1]); }

    func_t functions[] = {
        {nw,                        "nw", "orig",    "NA",    "32",  "32", 0, 0, 1},
        {nw_scan,                   "nw", "scan",    "NA",    "32",  "32", 0, 0, 0},
#if HAVE_SSE2
        {nw_scan_sse2_128_32,       "nw", "scan",    "sse2",  "128", "32", 0, 0, 0},
        {nw_scan_sse2_128_16,       "nw", "scan",    "sse2",  "128", "16", 0, 0, 0},
        {nw_scan_sse2_128_8,        "nw", "scan",    "sse2",  "128", "8",  0, 0, 0},
        {nw_diag_sse2_128_32,       "nw", "diag",    "sse2",  "128", "32", 0, 0, 0},
        {nw_diag_sse2_128_16,       "nw", "diag",    "sse2",  "128", "16", 0, 0, 0},
        {nw_diag_sse2_128_8,        "nw", "diag",    "sse2",  "128", "8",  0, 0, 0},
        {nw_striped_sse2_128_32,    "nw", "striped", "sse2",  "128", "32", 0, 0, 0},
        {nw_striped_sse2_128_16,    "nw", "striped", "sse2",  "128", "16", 0, 0, 0},
        {nw_striped_sse2_128_8,     "nw", "striped", "sse2",  "128", "8",  0, 0, 0},
#endif
#if HAVE_SSE41
        {nw_scan_sse41_128_32,      "nw", "scan",    "sse41", "128", "32", 0, 0, 0},
        {nw_scan_sse41_128_16,      "nw", "scan",    "sse41", "128", "16", 0, 0, 0},
        {nw_scan_sse41_128_8,       "nw", "scan",    "sse41", "128", "8",  0, 0, 0},
        {nw_diag_sse41_128_32,      "nw", "diag",    "sse41", "128", "32", 0, 0, 0},
        {nw_diag_sse41_128_16,      "nw", "diag",    "sse41", "128", "16", 0, 0, 0},
        {nw_diag_sse41_128_8,       "nw", "diag",    "sse41", "128", "8",  0, 0, 0},
        {nw_striped_sse41_128_32,   "nw", "striped", "sse41", "128", "32", 0, 0, 0},
        {nw_striped_sse41_128_16,   "nw", "striped", "sse41", "128", "16", 0, 0, 0},
        {nw_striped_sse41_128_8,    "nw", "striped", "sse41", "128", "8",  0, 0, 0},
#endif
#if HAVE_AVX2
        {nw_scan_avx2_256_32,       "nw", "scan",    "avx2",  "256", "32", 0, 0, 0},
        {nw_scan_avx2_256_16,       "nw", "scan",    "avx2",  "256", "16", 0, 0, 0},
        {nw_scan_avx2_256_8,        "nw", "scan",    "avx2",  "256", "8",  0, 0, 0},
        {nw_diag_avx2_256_32,       "nw", "diag",    "avx2",  "256", "32", 0, 0, 0},
        {nw_diag_avx2_256_16,       "nw", "diag",    "avx2",  "256", "16", 0, 0, 0},
        //{nw_diag_avx2_256_8,        "nw", "diag",    "avx2",  "256", "8",  0, 0, 0},
        {nw_striped_avx2_256_32,    "nw", "striped", "avx2",  "256", "32", 0, 0, 0},
        {nw_striped_avx2_256_16,    "nw", "striped", "avx2",  "256", "16", 0, 0, 0},
        //{nw_striped_avx2_256_8,     "nw", "striped", "avx2",  "256", "8",  0, 0, 0},
#endif
#if HAVE_KNC
        {nw_scan_knc_512_32,        "nw", "scan",    "knc",   "512", "32", 0, 0, 0},
        //{nw_diag_knc_512_32,        "nw", "diag",    "knc",   "512", "32", 0, 0, 0},
        //{nw_striped_knc_512_32,     "nw", "striped", "knc",   "512", "32", 0, 0, 0},
#endif

        {sg,                        "sg", "orig",    "NA",    "32",  "32", 0, 0, 1},
        {sg_scan,                   "sg", "scan",    "NA",    "32",  "32", 0, 0, 0},
#if HAVE_SSE2
        {sg_scan_sse2_128_32,       "sg", "scan",    "sse2",  "128", "32", 0, 0, 0},
        {sg_scan_sse2_128_16,       "sg", "scan",    "sse2",  "128", "16", 0, 0, 0},
        {sg_scan_sse2_128_8,        "sg", "scan",    "sse2",  "128", "8",  0, 0, 0},
        {sg_diag_sse2_128_32,       "sg", "diag",    "sse2",  "128", "32", 0, 0, 0},
        {sg_diag_sse2_128_16,       "sg", "diag",    "sse2",  "128", "16", 0, 0, 0},
        {sg_diag_sse2_128_8,        "sg", "diag",    "sse2",  "128", "8",  0, 0, 0},
        {sg_striped_sse2_128_32,    "sg", "striped", "sse2",  "128", "32", 0, 0, 0},
        {sg_striped_sse2_128_16,    "sg", "striped", "sse2",  "128", "16", 0, 0, 0},
        {sg_striped_sse2_128_8,     "sg", "striped", "sse2",  "128", "8",  0, 0, 0},
#endif
#if HAVE_SSE41
        {sg_scan_sse41_128_32,      "sg", "scan",    "sse41", "128", "32", 0, 0, 0},
        {sg_scan_sse41_128_16,      "sg", "scan",    "sse41", "128", "16", 0, 0, 0},
        {sg_scan_sse41_128_8,       "sg", "scan",    "sse41", "128", "8",  0, 0, 0},
        {sg_diag_sse41_128_32,      "sg", "diag",    "sse41", "128", "32", 0, 0, 0},
        {sg_diag_sse41_128_16,      "sg", "diag",    "sse41", "128", "16", 0, 0, 0},
        {sg_diag_sse41_128_8,       "sg", "diag",    "sse41", "128", "8",  0, 0, 0},
        {sg_striped_sse41_128_32,   "sg", "striped", "sse41", "128", "32", 0, 0, 0},
        {sg_striped_sse41_128_16,   "sg", "striped", "sse41", "128", "16", 0, 0, 0},
        {sg_striped_sse41_128_8,    "sg", "striped", "sse41", "128", "8",  0, 0, 0},
#endif
#if HAVE_AVX2
        {sg_scan_avx2_256_32,       "sg", "scan",    "avx2",  "256", "32", 0, 0, 0},
        {sg_scan_avx2_256_16,       "sg", "scan",    "avx2",  "256", "16", 0, 0, 0},
        {sg_scan_avx2_256_8,        "sg", "scan",    "avx2",  "256", "8",  0, 0, 0},
        {sg_diag_avx2_256_32,       "sg", "diag",    "avx2",  "256", "32", 0, 0, 0},
        {sg_diag_avx2_256_16,       "sg", "diag",    "avx2",  "256", "16", 0, 0, 0},
        //{sg_diag_avx2_256_8,        "sg", "diag",    "avx2",  "256", "8",  0, 0, 0},
        {sg_striped_avx2_256_32,    "sg", "striped", "avx2",  "256", "32", 0, 0, 0},
        {sg_striped_avx2_256_16,    "sg", "striped", "avx2",  "256", "16", 0, 0, 0},
        //{sg_striped_avx2_256_8,     "sg", "striped", "avx2",  "256", "8",  0, 0, 0},
#endif
#if HAVE_KNC
        {sg_scan_knc_512_32,        "sg", "scan",    "knc",   "512", "32", 0, 0, 0},
        //{sg_diag_knc_512_32,        "sg", "diag",    "knc",   "512", "32", 0, 0, 0},
        //{sg_striped_knc_512_32,     "sg", "striped", "knc",   "512", "32", 0, 0, 0},
#endif

        {sw,                        "sw", "orig",    "NA",    "32",  "32", 0, 0, 1},
        {sw_scan,                   "sw", "scan",    "NA",    "32",  "32", 0, 0, 0},
#if HAVE_SSE2
        {sw_scan_sse2_128_32,       "sw", "scan",    "sse2",  "128", "32", 0, 0, 0},
        {sw_scan_sse2_128_16,       "sw", "scan",    "sse2",  "128", "16", 0, 0, 0},
        {sw_scan_sse2_128_8,        "sw", "scan",    "sse2",  "128", "8",  0, 0, 0},
        {sw_diag_sse2_128_32,       "sw", "diag",    "sse2",  "128", "32", 0, 0, 0},
        {sw_diag_sse2_128_16,       "sw", "diag",    "sse2",  "128", "16", 0, 0, 0},
        {sw_diag_sse2_128_8,        "sw", "diag",    "sse2",  "128", "8",  0, 0, 0},
        {sw_striped_sse2_128_32,    "sw", "striped", "sse2",  "128", "32", 0, 0, 0},
        {sw_striped_sse2_128_16,    "sw", "striped", "sse2",  "128", "16", 0, 0, 0},
        {sw_striped_sse2_128_8,     "sw", "striped", "sse2",  "128", "8",  0, 0, 0},
#endif
#if HAVE_SSE41
        {sw_scan_sse41_128_32,      "sw", "scan",    "sse41", "128", "32", 0, 0, 0},
        {sw_scan_sse41_128_16,      "sw", "scan",    "sse41", "128", "16", 0, 0, 0},
        {sw_scan_sse41_128_8,       "sw", "scan",    "sse41", "128", "8",  0, 0, 0},
        {sw_diag_sse41_128_32,      "sw", "diag",    "sse41", "128", "32", 0, 0, 0},
        {sw_diag_sse41_128_16,      "sw", "diag",    "sse41", "128", "16", 0, 0, 0},
        {sw_diag_sse41_128_8,       "sw", "diag",    "sse41", "128", "8",  0, 0, 0},
        {sw_striped_sse41_128_32,   "sw", "striped", "sse41", "128", "32", 0, 0, 0},
        {sw_striped_sse41_128_16,   "sw", "striped", "sse41", "128", "16", 0, 0, 0},
        {sw_striped_sse41_128_8,    "sw", "striped", "sse41", "128", "8",  0, 0, 0},
#endif
#if HAVE_AVX2
        {sw_scan_avx2_256_32,       "sw", "scan",    "avx2",  "256", "32", 0, 0, 0},
        {sw_scan_avx2_256_16,       "sw", "scan",    "avx2",  "256", "16", 0, 0, 0},
        {sw_scan_avx2_256_8,        "sw", "scan",    "avx2",  "256", "8",  0, 0, 0},
        {sw_diag_avx2_256_32,       "sw", "diag",    "avx2",  "256", "32", 0, 0, 0},
        {sw_diag_avx2_256_16,       "sw", "diag",    "avx2",  "256", "16", 0, 0, 0},
        //{sw_diag_avx2_256_8,        "sw", "diag",    "avx2",  "256", "8",  0, 0, 0},
        {sw_striped_avx2_256_32,    "sw", "striped", "avx2",  "256", "32", 0, 0, 0},
        {sw_striped_avx2_256_16,    "sw", "striped", "avx2",  "256", "16", 0, 0, 0},
        //{sw_striped_avx2_256_8,     "sw", "striped", "avx2",  "256", "8",  0, 0, 0},
#endif
#if HAVE_KNC
        {sw_scan_knc_512_32,        "sw", "scan",    "knc",   "512", "32", 0, 0, 0},
        //{sw_diag_knc_512_32,        "sw", "diag",    "knc",   "512", "32", 0, 0, 0},
        //{sw_striped_knc_512_32,     "sw", "striped", "knc",   "512", "32", 0, 0, 0},
#endif
                                   
        {nw_table,                     "nw", "orig",    "NA",    "32",  "32", 1, 0, 1},
        {nw_table_scan,                "nw", "scan",    "NA",    "32",  "32", 1, 0, 0},
#if HAVE_SSE2
        {nw_table_scan_sse2_128_32,    "nw", "scan",    "sse2",  "128", "32", 1, 0, 0},
        {nw_table_scan_sse2_128_16,    "nw", "scan",    "sse2",  "128", "16", 1, 0, 0},
        {nw_table_scan_sse2_128_8,     "nw", "scan",    "sse2",  "128", "8",  1, 0, 0},
        {nw_table_diag_sse2_128_32,    "nw", "diag",    "sse2",  "128", "32", 1, 0, 0},
        {nw_table_diag_sse2_128_16,    "nw", "diag",    "sse2",  "128", "16", 1, 0, 0},
        {nw_table_diag_sse2_128_8,     "nw", "diag",    "sse2",  "128", "8",  1, 0, 0},
        {nw_table_striped_sse2_128_32, "nw", "striped", "sse2",  "128", "32", 1, 0, 0},
        {nw_table_striped_sse2_128_16, "nw", "striped", "sse2",  "128", "16", 1, 0, 0},
        {nw_table_striped_sse2_128_8,  "nw", "striped", "sse2",  "128", "8",  1, 0, 0},
#endif
#if HAVE_SSE41
        {nw_table_scan_sse41_128_32,    "nw", "scan",    "sse41", "128", "32", 1, 0, 0},
        {nw_table_scan_sse41_128_16,    "nw", "scan",    "sse41", "128", "16", 1, 0, 0},
        {nw_table_scan_sse41_128_8,     "nw", "scan",    "sse41", "128", "8",  1, 0, 0},
        {nw_table_diag_sse41_128_32,    "nw", "diag",    "sse41", "128", "32", 1, 0, 0},
        {nw_table_diag_sse41_128_16,    "nw", "diag",    "sse41", "128", "16", 1, 0, 0},
        {nw_table_diag_sse41_128_8,     "nw", "diag",    "sse41", "128", "8",  1, 0, 0},
        {nw_table_striped_sse41_128_32, "nw", "striped", "sse41", "128", "32", 1, 0, 0},
        {nw_table_striped_sse41_128_16, "nw", "striped", "sse41", "128", "16", 1, 0, 0},
        {nw_table_striped_sse41_128_8,  "nw", "striped", "sse41", "128", "8",  1, 0, 0},
#endif
#if HAVE_AVX2
        {nw_table_scan_avx2_256_32,       "nw", "scan",    "avx2",  "256", "32", 1, 0, 0},
        {nw_table_scan_avx2_256_16,       "nw", "scan",    "avx2",  "256", "16", 1, 0, 0},
        {nw_table_scan_avx2_256_8,        "nw", "scan",    "avx2",  "256", "8",  1, 0, 0},
        {nw_table_diag_avx2_256_32,       "nw", "diag",    "avx2",  "256", "32", 1, 0, 0},
        {nw_table_diag_avx2_256_16,       "nw", "diag",    "avx2",  "256", "16", 1, 0, 0},
        //{nw_table_diag_avx2_256_8,        "nw", "diag",    "avx2",  "256", "8",  1, 0, 0},
        {nw_table_striped_avx2_256_32,    "nw", "striped", "avx2",  "256", "32", 1, 0, 0},
        {nw_table_striped_avx2_256_16,    "nw", "striped", "avx2",  "256", "16", 1, 0, 0},
        //{nw_table_striped_avx2_256_8,     "nw", "striped", "avx2",  "256", "8",  1, 0, 0},
#endif
#if HAVE_KNC
        {nw_table_scan_knc_512_32,        "nw", "scan",    "knc",   "512", "32", 1, 0, 0},
        //{nw_table_diag_knc_512_32,        "nw", "diag",    "knc",   "512", "32", 1, 0, 0},
        //{nw_table_striped_knc_512_32,     "nw", "striped", "knc",   "512", "32", 1, 0, 0},
#endif

        {sg_table,                     "sg", "orig",    "NA",   "32",  "32", 1, 0, 1},
        {sg_table_scan,                "sg", "scan",    "NA",   "32",  "32", 1, 0, 0},
#if HAVE_SSE2
        {sg_table_scan_sse2_128_32,    "sg", "scan",    "sse2", "128", "32", 1, 0, 0},
        {sg_table_scan_sse2_128_16,    "sg", "scan",    "sse2", "128", "16", 1, 0, 0},
        {sg_table_scan_sse2_128_8,     "sg", "scan",    "sse2", "128", "8",  1, 0, 0},
        {sg_table_diag_sse2_128_32,    "sg", "diag",    "sse2", "128", "32", 1, 0, 0},
        {sg_table_diag_sse2_128_16,    "sg", "diag",    "sse2", "128", "16", 1, 0, 0},
        {sg_table_diag_sse2_128_8,     "sg", "diag",    "sse2", "128", "8",  1, 0, 0},
        {sg_table_striped_sse2_128_32, "sg", "striped", "sse2", "128", "32", 1, 0, 0},
        {sg_table_striped_sse2_128_16, "sg", "striped", "sse2", "128", "16", 1, 0, 0},
        {sg_table_striped_sse2_128_8,  "sg", "striped", "sse2", "128", "8",  1, 0, 0},
#endif
#if HAVE_SSE41
        {sg_table_scan_sse41_128_32,   "sg", "scan",    "sse41", "128", "32", 1, 0, 0},
        {sg_table_scan_sse41_128_16,   "sg", "scan",    "sse41", "128", "16", 1, 0, 0},
        {sg_table_scan_sse41_128_8,    "sg", "scan",    "sse41", "128", "8",  1, 0, 0},
        {sg_table_diag_sse41_128_32,   "sg", "diag",    "sse41", "128", "32", 1, 0, 0},
        {sg_table_diag_sse41_128_16,   "sg", "diag",    "sse41", "128", "16", 1, 0, 0},
        {sg_table_diag_sse41_128_8,    "sg", "diag",    "sse41", "128", "8",  1, 0, 0},
        {sg_table_striped_sse41_128_32,"sg", "striped", "sse41", "128", "32", 1, 0, 0},
        {sg_table_striped_sse41_128_16,"sg", "striped", "sse41", "128", "16", 1, 0, 0},
        {sg_table_striped_sse41_128_8, "sg", "striped", "sse41", "128", "8",  1, 0, 0},
#endif
#if HAVE_AVX2
        {sg_table_scan_avx2_256_32,       "sg", "scan",    "avx2",  "256", "32", 1, 0, 0},
        {sg_table_scan_avx2_256_16,       "sg", "scan",    "avx2",  "256", "16", 1, 0, 0},
        {sg_table_scan_avx2_256_8,        "sg", "scan",    "avx2",  "256", "8",  1, 0, 0},
        {sg_table_diag_avx2_256_32,       "sg", "diag",    "avx2",  "256", "32", 1, 0, 0},
        {sg_table_diag_avx2_256_16,       "sg", "diag",    "avx2",  "256", "16", 1, 0, 0},
        //{sg_table_diag_avx2_256_8,        "sg", "diag",    "avx2",  "256", "8",  1, 0, 0},
        {sg_table_striped_avx2_256_32,    "sg", "striped", "avx2",  "256", "32", 1, 0, 0},
        {sg_table_striped_avx2_256_16,    "sg", "striped", "avx2",  "256", "16", 1, 0, 0},
        //{sg_table_striped_avx2_256_8,     "sg", "striped", "avx2",  "256", "8",  1, 0, 0},
#endif
#if HAVE_KNC
        {sg_table_scan_knc_512_32,        "sg", "scan",    "knc",   "512", "32", 1, 0, 0},
        //{sg_table_diag_knc_512_32,        "sg", "diag",    "knc",   "512", "32", 1, 0, 0},
        //{sg_table_striped_knc_512_32,     "sg", "striped", "knc",   "512", "32", 1, 0, 0},
#endif

        {sw_table,                     "sw", "orig",    "NA",   "32",  "32", 1, 0, 1},
        {sw_table_scan,                "sw", "scan",    "NA",   "32",  "32", 1, 0, 0},
#if HAVE_SSE2
        {sw_table_scan_sse2_128_32,    "sw", "scan",    "sse2", "128", "32", 1, 0, 0},
        {sw_table_scan_sse2_128_16,    "sw", "scan",    "sse2", "128", "16", 1, 0, 0},
        {sw_table_scan_sse2_128_8,     "sw", "scan",    "sse2", "128", "8",  1, 0, 0},
        {sw_table_diag_sse2_128_32,    "sw", "diag",    "sse2", "128", "32", 1, 0, 0},
        {sw_table_diag_sse2_128_16,    "sw", "diag",    "sse2", "128", "16", 1, 0, 0},
        {sw_table_diag_sse2_128_8,     "sw", "diag",    "sse2", "128", "8",  1, 0, 0},
        {sw_table_striped_sse2_128_32, "sw", "striped", "sse2", "128", "32", 1, 0, 0},
        {sw_table_striped_sse2_128_16, "sw", "striped", "sse2", "128", "16", 1, 0, 0},
        {sw_table_striped_sse2_128_8,  "sw", "striped", "sse2", "128", "8",  1, 0, 0},
#endif
#if HAVE_SSE41
        {sw_table_scan_sse41_128_32,    "sw", "scan",    "sse41", "128", "32", 1, 0, 0},
        {sw_table_scan_sse41_128_16,    "sw", "scan",    "sse41", "128", "16", 1, 0, 0},
        {sw_table_scan_sse41_128_8,     "sw", "scan",    "sse41", "128", "8",  1, 0, 0},
        {sw_table_diag_sse41_128_32,    "sw", "diag",    "sse41", "128", "32", 1, 0, 0},
        {sw_table_diag_sse41_128_16,    "sw", "diag",    "sse41", "128", "16", 1, 0, 0},
        {sw_table_diag_sse41_128_8,     "sw", "diag",    "sse41", "128", "8",  1, 0, 0},
        {sw_table_striped_sse41_128_32, "sw", "striped", "sse41", "128", "32", 1, 0, 0},
        {sw_table_striped_sse41_128_16, "sw", "striped", "sse41", "128", "16", 1, 0, 0},
        {sw_table_striped_sse41_128_8,  "sw", "striped", "sse41", "128", "8",  1, 0, 0},
#endif
#if HAVE_AVX2
        {sw_table_scan_avx2_256_32,     "sw", "scan",    "avx2",  "256", "32", 1, 0, 0},
        {sw_table_scan_avx2_256_16,     "sw", "scan",    "avx2",  "256", "16", 1, 0, 0},
        {sw_table_scan_avx2_256_8,      "sw", "scan",    "avx2",  "256", "8",  1, 0, 0},
        {sw_table_diag_avx2_256_32,     "sw", "diag",    "avx2",  "256", "32", 1, 0, 0},
        {sw_table_diag_avx2_256_16,     "sw", "diag",    "avx2",  "256", "16", 1, 0, 0},
        //{sw_table_diag_avx2_256_8,      "sw", "diag",    "avx2",  "256", "8",  1, 0, 0},
        {sw_table_striped_avx2_256_32,  "sw", "striped", "avx2",  "256", "32", 1, 0, 0},
        {sw_table_striped_avx2_256_16,  "sw", "striped", "avx2",  "256", "16", 1, 0, 0},
        //{sw_table_striped_avx2_256_8,   "sw", "striped", "avx2",  "256", "8",  1, 0, 0},
#endif
#if HAVE_KNC
        {sw_table_scan_knc_512_32,        "sw", "scan",    "knc",   "512", "32", 1, 0, 0},
        //{sw_table_diag_knc_512_32,        "sw", "diag",    "knc",   "512", "32", 1, 0, 0},
        //{sw_table_striped_knc_512_32,     "sw", "striped", "knc",   "512", "32", 1, 0, 0},
#endif

        {nw_stats,                        "nw_stats", "orig",    "NA",    "32",  "32", 0, 1, 1},
        {nw_stats_scan,                   "nw_stats", "scan",    "NA",    "32",  "32", 0, 1, 0},
#if HAVE_SSE2
        {nw_stats_scan_sse2_128_32,       "nw_stats", "scan",    "sse2",  "128", "32", 0, 1, 0},
        {nw_stats_scan_sse2_128_16,       "nw_stats", "scan",    "sse2",  "128", "16", 0, 1, 0},
        {nw_stats_scan_sse2_128_8,        "nw_stats", "scan",    "sse2",  "128", "8",  0, 1, 0},
        {nw_stats_diag_sse2_128_32,       "nw_stats", "diag",    "sse2",  "128", "32", 0, 1, 0},
        {nw_stats_diag_sse2_128_16,       "nw_stats", "diag",    "sse2",  "128", "16", 0, 1, 0},
        {nw_stats_diag_sse2_128_8,        "nw_stats", "diag",    "sse2",  "128", "8",  0, 1, 0},
        {nw_stats_striped_sse2_128_32,    "nw_stats", "striped", "sse2",  "128", "32", 0, 1, 0},
        {nw_stats_striped_sse2_128_16,    "nw_stats", "striped", "sse2",  "128", "16", 0, 1, 0},
        {nw_stats_striped_sse2_128_8,     "nw_stats", "striped", "sse2",  "128", "8",  0, 1, 0},
#endif
#if HAVE_SSE41
        {nw_stats_scan_sse41_128_32,      "nw_stats", "scan",    "sse41", "128", "32", 0, 1, 0},
        {nw_stats_scan_sse41_128_16,      "nw_stats", "scan",    "sse41", "128", "16", 0, 1, 0},
        {nw_stats_scan_sse41_128_8,       "nw_stats", "scan",    "sse41", "128", "8",  0, 1, 0},
        {nw_stats_diag_sse41_128_32,      "nw_stats", "diag",    "sse41", "128", "32", 0, 1, 0},
        {nw_stats_diag_sse41_128_16,      "nw_stats", "diag",    "sse41", "128", "16", 0, 1, 0},
        {nw_stats_diag_sse41_128_8,       "nw_stats", "diag",    "sse41", "128", "8",  0, 1, 0},
        {nw_stats_striped_sse41_128_32,   "nw_stats", "striped", "sse41", "128", "32", 0, 1, 0},
        {nw_stats_striped_sse41_128_16,   "nw_stats", "striped", "sse41", "128", "16", 0, 1, 0},
        {nw_stats_striped_sse41_128_8,    "nw_stats", "striped", "sse41", "128", "8",  0, 1, 0},
#endif
#if HAVE_AVX2
        {nw_stats_scan_avx2_256_32,      "nw_stats", "scan",    "avx2", "256", "32", 0, 1, 0},
        {nw_stats_scan_avx2_256_16,      "nw_stats", "scan",    "avx2", "256", "16", 0, 1, 0},
        {nw_stats_scan_avx2_256_8,       "nw_stats", "scan",    "avx2", "256", "8",  0, 1, 0},
        //{nw_stats_diag_avx2_256_32,      "nw_stats", "diag",    "avx2", "256", "32", 0, 1, 0},
        //{nw_stats_diag_avx2_256_16,      "nw_stats", "diag",    "avx2", "256", "16", 0, 1, 0},
        //{nw_stats_diag_avx2_256_8,       "nw_stats", "diag",    "avx2", "256", "8",  0, 1, 0},
        //{nw_stats_striped_avx2_256_32,   "nw_stats", "striped", "avx2", "256", "32", 0, 1, 0},
        //{nw_stats_striped_avx2_256_16,   "nw_stats", "striped", "avx2", "256", "16", 0, 1, 0},
        //{nw_stats_striped_avx2_256_8,    "nw_stats", "striped", "avx2", "256", "8",  0, 1, 0},
#endif
                                   
        {sg_stats,                        "sg_stats", "orig",    "NA",    "32",  "32", 0, 1, 1},
        {sg_stats_scan,                   "sg_stats", "scan",    "NA",    "32",  "32", 0, 1, 0},
#if HAVE_SSE2
        {sg_stats_scan_sse2_128_32,       "sg_stats", "scan",    "sse2",  "128", "32", 0, 1, 0},
        {sg_stats_scan_sse2_128_16,       "sg_stats", "scan",    "sse2",  "128", "16", 0, 1, 0},
        {sg_stats_scan_sse2_128_8,        "sg_stats", "scan",    "sse2",  "128", "8",  0, 1, 0},
        {sg_stats_diag_sse2_128_32,       "sg_stats", "diag",    "sse2",  "128", "32", 0, 1, 0},
        {sg_stats_diag_sse2_128_16,       "sg_stats", "diag",    "sse2",  "128", "16", 0, 1, 0},
        {sg_stats_diag_sse2_128_8,        "sg_stats", "diag",    "sse2",  "128", "8",  0, 1, 0},
        {sg_stats_striped_sse2_128_32,    "sg_stats", "striped", "sse2",  "128", "32", 0, 1, 0},
        {sg_stats_striped_sse2_128_16,    "sg_stats", "striped", "sse2",  "128", "16", 0, 1, 0},
        {sg_stats_striped_sse2_128_8,     "sg_stats", "striped", "sse2",  "128", "8",  0, 1, 0},
#endif
#if HAVE_SSE41
        {sg_stats_scan_sse41_128_32,      "sg_stats", "scan",    "sse41", "128", "32", 0, 1, 0},
        {sg_stats_scan_sse41_128_16,      "sg_stats", "scan",    "sse41", "128", "16", 0, 1, 0},
        {sg_stats_scan_sse41_128_8,       "sg_stats", "scan",    "sse41", "128", "8",  0, 1, 0},
        {sg_stats_diag_sse41_128_32,      "sg_stats", "diag",    "sse41", "128", "32", 0, 1, 0},
        {sg_stats_diag_sse41_128_16,      "sg_stats", "diag",    "sse41", "128", "16", 0, 1, 0},
        {sg_stats_diag_sse41_128_8,       "sg_stats", "diag",    "sse41", "128", "8",  0, 1, 0},
        {sg_stats_striped_sse41_128_32,   "sg_stats", "striped", "sse41", "128", "32", 0, 1, 0},
        {sg_stats_striped_sse41_128_16,   "sg_stats", "striped", "sse41", "128", "16", 0, 1, 0},
        {sg_stats_striped_sse41_128_8,    "sg_stats", "striped", "sse41", "128", "8",  0, 1, 0},
#endif
#if HAVE_AVX2
        {sg_stats_scan_avx2_256_32,       "sg_stats", "scan",    "avx2",  "256", "32", 0, 1, 0},
        {sg_stats_scan_avx2_256_16,       "sg_stats", "scan",    "avx2",  "256", "16", 0, 1, 0},
        {sg_stats_scan_avx2_256_8,        "sg_stats", "scan",    "avx2",  "256", "8",  0, 1, 0},
        //{sg_stats_diag_avx2_256_32,       "sg_stats", "diag",    "avx2",  "256", "32", 0, 1, 0},
        //{sg_stats_diag_avx2_256_16,       "sg_stats", "diag",    "avx2",  "256", "16", 0, 1, 0},
        //{sg_stats_diag_avx2_256_8,        "sg_stats", "diag",    "avx2",  "256", "8",  0, 1, 0},
        //{sg_stats_striped_avx2_256_32,    "sg_stats", "striped", "avx2",  "256", "32", 0, 1, 0},
        //{sg_stats_striped_avx2_256_16,    "sg_stats", "striped", "avx2",  "256", "16", 0, 1, 0},
        //{sg_stats_striped_avx2_256_8,     "sg_stats", "striped", "avx2",  "256", "8",  0, 1, 0},
#endif
                                   
        {sw_stats,                        "sw_stats", "orig",    "NA",    "32",  "32", 0, 1, 1},
        {sw_stats_scan,                   "sw_stats", "scan",    "NA",    "32",  "32", 0, 1, 0},
#if HAVE_SSE2
        {sw_stats_scan_sse2_128_32,       "sw_stats", "scan",    "sse2",  "128", "32", 0, 1, 0},
        {sw_stats_scan_sse2_128_16,       "sw_stats", "scan",    "sse2",  "128", "16", 0, 1, 0},
        {sw_stats_scan_sse2_128_8,        "sw_stats", "scan",    "sse2",  "128", "8",  0, 1, 0},
        {sw_stats_diag_sse2_128_32,       "sw_stats", "diag",    "sse2",  "128", "32", 0, 1, 0},
        {sw_stats_diag_sse2_128_16,       "sw_stats", "diag",    "sse2",  "128", "16", 0, 1, 0},
        {sw_stats_diag_sse2_128_8,        "sw_stats", "diag",    "sse2",  "128", "8",  0, 1, 0},
        {sw_stats_striped_sse2_128_32,    "sw_stats", "striped", "sse2",  "128", "32", 0, 1, 0},
        {sw_stats_striped_sse2_128_16,    "sw_stats", "striped", "sse2",  "128", "16", 0, 1, 0},
        {sw_stats_striped_sse2_128_8,     "sw_stats", "striped", "sse2",  "128", "8",  0, 1, 0},
#endif
#if HAVE_SSE41
        {sw_stats_scan_sse41_128_32,      "sw_stats", "scan",    "sse41", "128", "32", 0, 1, 0},
        {sw_stats_scan_sse41_128_16,      "sw_stats", "scan",    "sse41", "128", "16", 0, 1, 0},
        {sw_stats_scan_sse41_128_8,       "sw_stats", "scan",    "sse41", "128", "8",  0, 1, 0},
        {sw_stats_diag_sse41_128_32,      "sw_stats", "diag",    "sse41", "128", "32", 0, 1, 0},
        {sw_stats_diag_sse41_128_16,      "sw_stats", "diag",    "sse41", "128", "16", 0, 1, 0},
        {sw_stats_diag_sse41_128_8,       "sw_stats", "diag",    "sse41", "128", "8",  0, 1, 0},
        {sw_stats_striped_sse41_128_32,   "sw_stats", "striped", "sse41", "128", "32", 0, 1, 0},
        {sw_stats_striped_sse41_128_16,   "sw_stats", "striped", "sse41", "128", "16", 0, 1, 0},
        {sw_stats_striped_sse41_128_8,    "sw_stats", "striped", "sse41", "128", "8",  0, 1, 0},
#endif
#if HAVE_AVX2
        {sw_stats_scan_avx2_256_32,      "sw_stats", "scan",    "avx2", "256", "32", 0, 1, 0},
        {sw_stats_scan_avx2_256_16,      "sw_stats", "scan",    "avx2", "256", "16", 0, 1, 0},
        {sw_stats_scan_avx2_256_8,       "sw_stats", "scan",    "avx2", "256", "8",  0, 1, 0},
        //{sw_stats_diag_avx2_256_32,      "sw_stats", "diag",    "avx2", "256", "32", 0, 1, 0},
        //{sw_stats_diag_avx2_256_16,      "sw_stats", "diag",    "avx2", "256", "16", 0, 1, 0},
        //{sw_stats_diag_avx2_256_8,       "sw_stats", "diag",    "avx2", "256", "8",  0, 1, 0},
        //{sw_stats_striped_avx2_256_32,   "sw_stats", "striped", "avx2", "256", "32", 0, 1, 0},
        //{sw_stats_striped_avx2_256_16,   "sw_stats", "striped", "avx2", "256", "16", 0, 1, 0},
        //{sw_stats_striped_avx2_256_8,    "sw_stats", "striped", "avx2", "256", "8",  0, 1, 0},
#endif
                                   
        {nw_stats_table,                        "nw_stats", "orig",    "NA",    "32",  "32", 1, 1, 1},
        {nw_stats_table_scan,                   "nw_stats", "scan",    "NA",    "32",  "32", 1, 1, 0},
#if HAVE_SSE2
        {nw_stats_table_scan_sse2_128_32,       "nw_stats", "scan",    "sse2",  "128", "32", 1, 1, 0},
        {nw_stats_table_scan_sse2_128_16,       "nw_stats", "scan",    "sse2",  "128", "16", 1, 1, 0},
        {nw_stats_table_scan_sse2_128_8,        "nw_stats", "scan",    "sse2",  "128", "8",  1, 1, 0},
        {nw_stats_table_diag_sse2_128_32,       "nw_stats", "diag",    "sse2",  "128", "32", 1, 1, 0},
        {nw_stats_table_diag_sse2_128_16,       "nw_stats", "diag",    "sse2",  "128", "16", 1, 1, 0},
        {nw_stats_table_diag_sse2_128_8,        "nw_stats", "diag",    "sse2",  "128", "8",  1, 1, 0},
        {nw_stats_table_striped_sse2_128_32,    "nw_stats", "striped", "sse2",  "128", "32", 1, 1, 0},
        {nw_stats_table_striped_sse2_128_16,    "nw_stats", "striped", "sse2",  "128", "16", 1, 1, 0},
        {nw_stats_table_striped_sse2_128_8,     "nw_stats", "striped", "sse2",  "128", "8",  1, 1, 0},
#endif
#if HAVE_SSE41
        {nw_stats_table_scan_sse41_128_32,      "nw_stats", "scan",    "sse41", "128", "32", 1, 1, 0},
        {nw_stats_table_scan_sse41_128_16,      "nw_stats", "scan",    "sse41", "128", "16", 1, 1, 0},
        {nw_stats_table_scan_sse41_128_8,       "nw_stats", "scan",    "sse41", "128", "8",  1, 1, 0},
        {nw_stats_table_diag_sse41_128_32,      "nw_stats", "diag",    "sse41", "128", "32", 1, 1, 0},
        {nw_stats_table_diag_sse41_128_16,      "nw_stats", "diag",    "sse41", "128", "16", 1, 1, 0},
        {nw_stats_table_diag_sse41_128_8,       "nw_stats", "diag",    "sse41", "128", "8",  1, 1, 0},
        {nw_stats_table_striped_sse41_128_32,   "nw_stats", "striped", "sse41", "128", "32", 1, 1, 0},
        {nw_stats_table_striped_sse41_128_16,   "nw_stats", "striped", "sse41", "128", "16", 1, 1, 0},
        {nw_stats_table_striped_sse41_128_8,    "nw_stats", "striped", "sse41", "128", "8",  1, 1, 0},
#endif
#if HAVE_AVX2
        {nw_stats_table_scan_avx2_256_32,       "nw_stats", "scan",    "avx2",  "256", "32", 1, 1, 0},
        {nw_stats_table_scan_avx2_256_16,       "nw_stats", "scan",    "avx2",  "256", "16", 1, 1, 0},
        {nw_stats_table_scan_avx2_256_8,        "nw_stats", "scan",    "avx2",  "256", "8",  1, 1, 0},
        //{nw_stats_table_diag_avx2_256_32,       "nw_stats", "diag",    "avx2",  "256", "32", 1, 1, 0},
        //{nw_stats_table_diag_avx2_256_16,       "nw_stats", "diag",    "avx2",  "256", "16", 1, 1, 0},
        //{nw_stats_table_diag_avx2_256_8,        "nw_stats", "diag",    "avx2",  "256", "8",  1, 1, 0},
        //{nw_stats_table_striped_avx2_256_32,    "nw_stats", "striped", "avx2",  "256", "32", 1, 1, 0},
        //{nw_stats_table_striped_avx2_256_16,    "nw_stats", "striped", "avx2",  "256", "16", 1, 1, 0},
        //{nw_stats_table_striped_avx2_256_8,     "nw_stats", "striped", "avx2",  "256", "8",  1, 1, 0},
#endif

        {sg_stats_table,                        "sg_stats", "orig",    "NA",    "32",  "32", 1, 1, 1},
        {sg_stats_table_scan,                   "sg_stats", "scan",    "NA",    "32",  "32", 1, 1, 0},
#if HAVE_SSE2
        {sg_stats_table_scan_sse2_128_32,       "sg_stats", "scan",    "sse2",  "128", "32", 1, 1, 0},
        {sg_stats_table_scan_sse2_128_16,       "sg_stats", "scan",    "sse2",  "128", "16", 1, 1, 0},
        {sg_stats_table_scan_sse2_128_8,        "sg_stats", "scan",    "sse2",  "128", "8",  1, 1, 0},
        {sg_stats_table_diag_sse2_128_32,       "sg_stats", "diag",    "sse2",  "128", "32", 1, 1, 0},
        {sg_stats_table_diag_sse2_128_16,       "sg_stats", "diag",    "sse2",  "128", "16", 1, 1, 0},
        {sg_stats_table_diag_sse2_128_8,        "sg_stats", "diag",    "sse2",  "128", "8",  1, 1, 0},
        {sg_stats_table_striped_sse2_128_32,    "sg_stats", "striped", "sse2",  "128", "32", 1, 1, 0},
        {sg_stats_table_striped_sse2_128_16,    "sg_stats", "striped", "sse2",  "128", "16", 1, 1, 0},
        {sg_stats_table_striped_sse2_128_8,     "sg_stats", "striped", "sse2",  "128", "8",  1, 1, 0},
#endif
#if HAVE_SSE41
        {sg_stats_table_scan_sse41_128_32,      "sg_stats", "scan",    "sse41", "128", "32", 1, 1, 0},
        {sg_stats_table_scan_sse41_128_16,      "sg_stats", "scan",    "sse41", "128", "16", 1, 1, 0},
        {sg_stats_table_scan_sse41_128_8,       "sg_stats", "scan",    "sse41", "128", "8",  1, 1, 0},
        {sg_stats_table_diag_sse41_128_32,      "sg_stats", "diag",    "sse41", "128", "32", 1, 1, 0},
        {sg_stats_table_diag_sse41_128_16,      "sg_stats", "diag",    "sse41", "128", "16", 1, 1, 0},
        {sg_stats_table_diag_sse41_128_8,       "sg_stats", "diag",    "sse41", "128", "8",  1, 1, 0},
        {sg_stats_table_striped_sse41_128_32,   "sg_stats", "striped", "sse41", "128", "32", 1, 1, 0},
        {sg_stats_table_striped_sse41_128_16,   "sg_stats", "striped", "sse41", "128", "16", 1, 1, 0},
        {sg_stats_table_striped_sse41_128_8,    "sg_stats", "striped", "sse41", "128", "8",  1, 1, 0},
#endif
#if HAVE_AVX2
        {sg_stats_table_scan_avx2_256_32,       "sg_stats", "scan",    "avx2",  "256", "32", 1, 1, 0},
        {sg_stats_table_scan_avx2_256_16,       "sg_stats", "scan",    "avx2",  "256", "16", 1, 1, 0},
        {sg_stats_table_scan_avx2_256_8,        "sg_stats", "scan",    "avx2",  "256", "8",  1, 1, 0},
        //{sg_stats_table_diag_avx2_256_32,       "sg_stats", "diag",    "avx2",  "256", "32", 1, 1, 0},
        //{sg_stats_table_diag_avx2_256_16,       "sg_stats", "diag",    "avx2",  "256", "16", 1, 1, 0},
        //{sg_stats_table_diag_avx2_256_8,        "sg_stats", "diag",    "avx2",  "256", "8",  1, 1, 0},
        //{sg_stats_table_striped_avx2_256_32,    "sg_stats", "striped", "avx2",  "256", "32", 1, 1, 0},
        //{sg_stats_table_striped_avx2_256_16,    "sg_stats", "striped", "avx2",  "256", "16", 1, 1, 0},
        //{sg_stats_table_striped_avx2_256_8,     "sg_stats", "striped", "avx2",  "256", "8",  1, 1, 0},
#endif

        {sw_stats_table,                        "sw_stats", "orig",    "NA",    "32",  "32", 1, 1, 1},
        {sw_stats_table_scan,                   "sw_stats", "scan",    "NA",    "32",  "32", 1, 1, 0},
#if HAVE_SSE2
        {sw_stats_table_scan_sse2_128_32,       "sw_stats", "scan",    "sse2",  "128", "32", 1, 1, 0},
        {sw_stats_table_scan_sse2_128_16,       "sw_stats", "scan",    "sse2",  "128", "16", 1, 1, 0},
        {sw_stats_table_scan_sse2_128_8,        "sw_stats", "scan",    "sse2",  "128", "8",  1, 1, 0},
        {sw_stats_table_diag_sse2_128_32,       "sw_stats", "diag",    "sse2",  "128", "32", 1, 1, 0},
        {sw_stats_table_diag_sse2_128_16,       "sw_stats", "diag",    "sse2",  "128", "16", 1, 1, 0},
        {sw_stats_table_diag_sse2_128_8,        "sw_stats", "diag",    "sse2",  "128", "8",  1, 1, 0},
        {sw_stats_table_striped_sse2_128_32,    "sw_stats", "striped", "sse2",  "128", "32", 1, 1, 0},
        {sw_stats_table_striped_sse2_128_16,    "sw_stats", "striped", "sse2",  "128", "16", 1, 1, 0},
        {sw_stats_table_striped_sse2_128_8,     "sw_stats", "striped", "sse2",  "128", "8",  1, 1, 0},
#endif
#if HAVE_SSE41
        {sw_stats_table_scan_sse41_128_32,      "sw_stats", "scan",    "sse41", "128", "32", 1, 1, 0},
        {sw_stats_table_scan_sse41_128_16,      "sw_stats", "scan",    "sse41", "128", "16", 1, 1, 0},
        {sw_stats_table_scan_sse41_128_8,       "sw_stats", "scan",    "sse41", "128", "8",  1, 1, 0},
        {sw_stats_table_diag_sse41_128_32,      "sw_stats", "diag",    "sse41", "128", "32", 1, 1, 0},
        {sw_stats_table_diag_sse41_128_16,      "sw_stats", "diag",    "sse41", "128", "16", 1, 1, 0},
        {sw_stats_table_diag_sse41_128_8,       "sw_stats", "diag",    "sse41", "128", "8",  1, 1, 0},
        {sw_stats_table_striped_sse41_128_32,   "sw_stats", "striped", "sse41", "128", "32", 1, 1, 0},
        {sw_stats_table_striped_sse41_128_16,   "sw_stats", "striped", "sse41", "128", "16", 1, 1, 0},
        {sw_stats_table_striped_sse41_128_8,    "sw_stats", "striped", "sse41", "128", "8",  1, 1, 0},
#endif
#if HAVE_AVX2
        {sw_stats_table_scan_avx2_256_32,       "sw_stats", "scan",    "avx2",  "256", "32", 1, 1, 0},
        {sw_stats_table_scan_avx2_256_16,       "sw_stats", "scan",    "avx2",  "256", "16", 1, 1, 0},
        {sw_stats_table_scan_avx2_256_8,        "sw_stats", "scan",    "avx2",  "256", "8",  1, 1, 0},
        //{sw_stats_table_diag_avx2_256_32,       "sw_stats", "diag",    "avx2",  "256", "32", 1, 1, 0},
        //{sw_stats_table_diag_avx2_256_16,       "sw_stats", "diag",    "avx2",  "256", "16", 1, 1, 0},
        //{sw_stats_table_diag_avx2_256_8,        "sw_stats", "diag",    "avx2",  "256", "8",  1, 1, 0},
        //{sw_stats_table_striped_avx2_256_32,    "sw_stats", "striped", "avx2",  "256", "32", 1, 1, 0},
        //{sw_stats_table_striped_avx2_256_16,    "sw_stats", "striped", "avx2",  "256", "16", 1, 1, 0},
        //{sw_stats_table_striped_avx2_256_8,     "sw_stats", "striped", "avx2",  "256", "8",  1, 1, 0},
#endif

        {NULL, "", "", "", "", "", 0, 0, 0}
    };

    timer_init();

    printf("%-15s %8s %6s %4s %5s %5s %8s %8s %8s %8s %5s %8s %8s %8s\n",
            "name", "type", "isa", "bits", "width", "elem",
            "score", "matches", "length",
            "avg", "imp", "stddev", "min", "max");

    stats_clear(&stats_rdtsc);
    stats_clear(&stats_nsecs);
    index = 0;
    f = functions[index++];
    while (f.f) {
        char name[16] = {'\0'};
        int new_limit = f.is_table ? 1 : limit;
        if ((0 == strncmp(f.isa, "sse2",  4) && 0 == parasail_can_use_sse2()) 
                || (0 == strncmp(f.isa, "sse41", 5) && 0 == parasail_can_use_sse41())
                || (0 == strncmp(f.isa, "avx2",  4) && 0 == parasail_can_use_avx2())) {
            f = functions[index++];
            continue;
        }
        stats_clear(&stats_rdtsc);
        timer_rdtsc = timer_start();
        timer_nsecs = timer_real();
        for (i=0; i<new_limit; ++i) {
            timer_rdtsc_single = timer_start();
            timer_nsecs_single = timer_real();
            result = f.f(seqA, lena, seqB, lenb, 10, 1, blosum62);
            timer_rdtsc_single = timer_end(timer_rdtsc_single);
            timer_nsecs_single = timer_real() - timer_nsecs_single;
            stats_sample_value(&stats_rdtsc, timer_rdtsc_single);
            stats_sample_value(&stats_nsecs, timer_nsecs_single);
            score = result->score;
            matches = result->matches;
            length = result->length;
            parasail_result_free(result);
        }
        timer_rdtsc = timer_end(timer_rdtsc);
        timer_nsecs = timer_real() - timer_nsecs;
        if (f.is_ref) {
            timer_rdtsc_ref_mean = stats_rdtsc._mean;
            timer_nsecs_ref_mean = stats_nsecs._mean;
        }
        strcpy(name, f.alg);
        /* xeon phi was unable to perform I/O running natively */
        if (f.is_table && 0 == HAVE_KNC) {
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
            result = f.f(seqA, lena, seqB, lenb, 10, 1, blosum62);
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
                strcat(filename, "_len");
                strcat(filename, suffix);
                print_array(filename, result->length_table, seqA, lena, seqB, lenb);
            }
            parasail_result_free(result);
        }
        if (f.is_table) {
            strcat(name, "_table");
        }
#if USE_TIMER_REAL
        printf("%-15s %8s %6s %4s %5s %5d %8d %8d %8d %8.2f %5.2f %8.2f %8.2f %8.2f %8.7f %5.7f %8.7f %8.7f %8.7f\n",
                name, f.type, f.isa, f.bits, f.width, elem(f),
                score, matches, length,
                stats_rdtsc._mean, pctf(timer_rdtsc_ref_mean, stats_rdtsc._mean),
                stats_stddev(&stats_rdtsc), stats_rdtsc._min, stats_rdtsc._max,
                stats_nsecs._mean, pctf(timer_nsecs_ref_mean, stats_nsecs._mean),
                stats_stddev(&stats_nsecs), stats_nsecs._min, stats_nsecs._max);
#else
        printf("%-15s %8s %6s %4s %5s %5d %8d %8d %8d %8.2f %5.2f %8.2f %8.2f %8.2f\n",
                name, f.type, f.isa, f.bits, f.width, elem(f),
                score, matches, length,
                stats_rdtsc._mean, pctf(timer_rdtsc_ref_mean, stats_rdtsc._mean),
                stats_stddev(&stats_rdtsc), stats_rdtsc._min, stats_rdtsc._max);
#endif
        f = functions[index++];
    }

    return 0;
}
