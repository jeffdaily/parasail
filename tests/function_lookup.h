/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#ifndef _PARASAIL_FUNCTION_TYPE_H_
#define _PARASAIL_FUNCTION_TYPE_H_

typedef parasail_result_t* (*parasail_function_t)(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24]);

typedef struct func {
    parasail_function_t pointer;
    const char * name;
    const char * alg;
    const char * type;
    const char * isa;
    const char * bits;
    const char * width;
    int lanes;
    char is_table;
    char is_stats;
    char is_ref;
} func_t;

typedef struct funcs {
    const char * name;
    func_t *fs;
} funcs_t;


func_t functions[] = {
{nw,                                  "nw",                                  "nw",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{nw_scan,                             "nw_scan",                             "nw",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
#if HAVE_SSE2
{nw_scan_sse2_128_64,                 "nw_scan_sse2_128_64",                 "nw",    "scan", "sse2",  "128", "64",  2, 0, 0, 0},
{nw_scan_sse2_128_32,                 "nw_scan_sse2_128_32",                 "nw",    "scan", "sse2",  "128", "32",  4, 0, 0, 0},
{nw_scan_sse2_128_16,                 "nw_scan_sse2_128_16",                 "nw",    "scan", "sse2",  "128", "16",  8, 0, 0, 0},
{nw_scan_sse2_128_8,                  "nw_scan_sse2_128_8",                  "nw",    "scan", "sse2",  "128",  "8", 16, 0, 0, 0},
{nw_striped_sse2_128_64,              "nw_striped_sse2_128_64",              "nw", "striped", "sse2",  "128", "64",  2, 0, 0, 0},
{nw_striped_sse2_128_32,              "nw_striped_sse2_128_32",              "nw", "striped", "sse2",  "128", "32",  4, 0, 0, 0},
{nw_striped_sse2_128_16,              "nw_striped_sse2_128_16",              "nw", "striped", "sse2",  "128", "16",  8, 0, 0, 0},
{nw_striped_sse2_128_8,               "nw_striped_sse2_128_8",               "nw", "striped", "sse2",  "128",  "8", 16, 0, 0, 0},
{nw_diag_sse2_128_64,                 "nw_diag_sse2_128_64",                 "nw",    "diag", "sse2",  "128", "64",  2, 0, 0, 0},
{nw_diag_sse2_128_32,                 "nw_diag_sse2_128_32",                 "nw",    "diag", "sse2",  "128", "32",  4, 0, 0, 0},
{nw_diag_sse2_128_16,                 "nw_diag_sse2_128_16",                 "nw",    "diag", "sse2",  "128", "16",  8, 0, 0, 0},
{nw_diag_sse2_128_8,                  "nw_diag_sse2_128_8",                  "nw",    "diag", "sse2",  "128",  "8", 16, 0, 0, 0},
#endif
#if HAVE_SSE41
{nw_scan_sse41_128_64,                "nw_scan_sse41_128_64",                "nw",    "scan", "sse41", "128", "64",  2, 0, 0, 0},
{nw_scan_sse41_128_32,                "nw_scan_sse41_128_32",                "nw",    "scan", "sse41", "128", "32",  4, 0, 0, 0},
{nw_scan_sse41_128_16,                "nw_scan_sse41_128_16",                "nw",    "scan", "sse41", "128", "16",  8, 0, 0, 0},
{nw_scan_sse41_128_8,                 "nw_scan_sse41_128_8",                 "nw",    "scan", "sse41", "128",  "8", 16, 0, 0, 0},
{nw_striped_sse41_128_64,             "nw_striped_sse41_128_64",             "nw", "striped", "sse41", "128", "64",  2, 0, 0, 0},
{nw_striped_sse41_128_32,             "nw_striped_sse41_128_32",             "nw", "striped", "sse41", "128", "32",  4, 0, 0, 0},
{nw_striped_sse41_128_16,             "nw_striped_sse41_128_16",             "nw", "striped", "sse41", "128", "16",  8, 0, 0, 0},
{nw_striped_sse41_128_8,              "nw_striped_sse41_128_8",              "nw", "striped", "sse41", "128",  "8", 16, 0, 0, 0},
{nw_diag_sse41_128_64,                "nw_diag_sse41_128_64",                "nw",    "diag", "sse41", "128", "64",  2, 0, 0, 0},
{nw_diag_sse41_128_32,                "nw_diag_sse41_128_32",                "nw",    "diag", "sse41", "128", "32",  4, 0, 0, 0},
{nw_diag_sse41_128_16,                "nw_diag_sse41_128_16",                "nw",    "diag", "sse41", "128", "16",  8, 0, 0, 0},
{nw_diag_sse41_128_8,                 "nw_diag_sse41_128_8",                 "nw",    "diag", "sse41", "128",  "8", 16, 0, 0, 0},
#endif
#if HAVE_AVX2
{nw_scan_avx2_256_64,                 "nw_scan_avx2_256_64",                 "nw",    "scan", "avx2",  "256", "64",  4, 0, 0, 0},
{nw_scan_avx2_256_32,                 "nw_scan_avx2_256_32",                 "nw",    "scan", "avx2",  "256", "32",  8, 0, 0, 0},
{nw_scan_avx2_256_16,                 "nw_scan_avx2_256_16",                 "nw",    "scan", "avx2",  "256", "16", 16, 0, 0, 0},
{nw_scan_avx2_256_8,                  "nw_scan_avx2_256_8",                  "nw",    "scan", "avx2",  "256",  "8", 32, 0, 0, 0},
{nw_striped_avx2_256_64,              "nw_striped_avx2_256_64",              "nw", "striped", "avx2",  "256", "64",  4, 0, 0, 0},
{nw_striped_avx2_256_32,              "nw_striped_avx2_256_32",              "nw", "striped", "avx2",  "256", "32",  8, 0, 0, 0},
{nw_striped_avx2_256_16,              "nw_striped_avx2_256_16",              "nw", "striped", "avx2",  "256", "16", 16, 0, 0, 0},
{nw_striped_avx2_256_8,               "nw_striped_avx2_256_8",               "nw", "striped", "avx2",  "256",  "8", 32, 0, 0, 0},
{nw_diag_avx2_256_64,                 "nw_diag_avx2_256_64",                 "nw",    "diag", "avx2",  "256", "64",  4, 0, 0, 0},
{nw_diag_avx2_256_32,                 "nw_diag_avx2_256_32",                 "nw",    "diag", "avx2",  "256", "32",  8, 0, 0, 0},
{nw_diag_avx2_256_16,                 "nw_diag_avx2_256_16",                 "nw",    "diag", "avx2",  "256", "16", 16, 0, 0, 0},
{nw_diag_avx2_256_8,                  "nw_diag_avx2_256_8",                  "nw",    "diag", "avx2",  "256",  "8", 32, 0, 0, 0},
#endif
#if HAVE_KNC
{nw_scan_knc_512_32,                  "nw_scan_knc_512_32",                  "nw",    "scan", "knc",   "512", "32", 16, 0, 0, 0},
{nw_striped_knc_512_32,               "nw_striped_knc_512_32",               "nw", "striped", "knc",   "512", "32", 16, 0, 0, 0},
{nw_diag_knc_512_32,                  "nw_diag_knc_512_32",                  "nw",    "diag", "knc",   "512", "32", 16, 0, 0, 0},
#endif
{sg,                                  "sg",                                  "sg",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{sg_scan,                             "sg_scan",                             "sg",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
#if HAVE_SSE2
{sg_scan_sse2_128_64,                 "sg_scan_sse2_128_64",                 "sg",    "scan", "sse2",  "128", "64",  2, 0, 0, 0},
{sg_scan_sse2_128_32,                 "sg_scan_sse2_128_32",                 "sg",    "scan", "sse2",  "128", "32",  4, 0, 0, 0},
{sg_scan_sse2_128_16,                 "sg_scan_sse2_128_16",                 "sg",    "scan", "sse2",  "128", "16",  8, 0, 0, 0},
{sg_scan_sse2_128_8,                  "sg_scan_sse2_128_8",                  "sg",    "scan", "sse2",  "128",  "8", 16, 0, 0, 0},
{sg_striped_sse2_128_64,              "sg_striped_sse2_128_64",              "sg", "striped", "sse2",  "128", "64",  2, 0, 0, 0},
{sg_striped_sse2_128_32,              "sg_striped_sse2_128_32",              "sg", "striped", "sse2",  "128", "32",  4, 0, 0, 0},
{sg_striped_sse2_128_16,              "sg_striped_sse2_128_16",              "sg", "striped", "sse2",  "128", "16",  8, 0, 0, 0},
{sg_striped_sse2_128_8,               "sg_striped_sse2_128_8",               "sg", "striped", "sse2",  "128",  "8", 16, 0, 0, 0},
{sg_diag_sse2_128_64,                 "sg_diag_sse2_128_64",                 "sg",    "diag", "sse2",  "128", "64",  2, 0, 0, 0},
{sg_diag_sse2_128_32,                 "sg_diag_sse2_128_32",                 "sg",    "diag", "sse2",  "128", "32",  4, 0, 0, 0},
{sg_diag_sse2_128_16,                 "sg_diag_sse2_128_16",                 "sg",    "diag", "sse2",  "128", "16",  8, 0, 0, 0},
{sg_diag_sse2_128_8,                  "sg_diag_sse2_128_8",                  "sg",    "diag", "sse2",  "128",  "8", 16, 0, 0, 0},
#endif
#if HAVE_SSE41
{sg_scan_sse41_128_64,                "sg_scan_sse41_128_64",                "sg",    "scan", "sse41", "128", "64",  2, 0, 0, 0},
{sg_scan_sse41_128_32,                "sg_scan_sse41_128_32",                "sg",    "scan", "sse41", "128", "32",  4, 0, 0, 0},
{sg_scan_sse41_128_16,                "sg_scan_sse41_128_16",                "sg",    "scan", "sse41", "128", "16",  8, 0, 0, 0},
{sg_scan_sse41_128_8,                 "sg_scan_sse41_128_8",                 "sg",    "scan", "sse41", "128",  "8", 16, 0, 0, 0},
{sg_striped_sse41_128_64,             "sg_striped_sse41_128_64",             "sg", "striped", "sse41", "128", "64",  2, 0, 0, 0},
{sg_striped_sse41_128_32,             "sg_striped_sse41_128_32",             "sg", "striped", "sse41", "128", "32",  4, 0, 0, 0},
{sg_striped_sse41_128_16,             "sg_striped_sse41_128_16",             "sg", "striped", "sse41", "128", "16",  8, 0, 0, 0},
{sg_striped_sse41_128_8,              "sg_striped_sse41_128_8",              "sg", "striped", "sse41", "128",  "8", 16, 0, 0, 0},
{sg_diag_sse41_128_64,                "sg_diag_sse41_128_64",                "sg",    "diag", "sse41", "128", "64",  2, 0, 0, 0},
{sg_diag_sse41_128_32,                "sg_diag_sse41_128_32",                "sg",    "diag", "sse41", "128", "32",  4, 0, 0, 0},
{sg_diag_sse41_128_16,                "sg_diag_sse41_128_16",                "sg",    "diag", "sse41", "128", "16",  8, 0, 0, 0},
{sg_diag_sse41_128_8,                 "sg_diag_sse41_128_8",                 "sg",    "diag", "sse41", "128",  "8", 16, 0, 0, 0},
#endif
#if HAVE_AVX2
{sg_scan_avx2_256_64,                 "sg_scan_avx2_256_64",                 "sg",    "scan", "avx2",  "256", "64",  4, 0, 0, 0},
{sg_scan_avx2_256_32,                 "sg_scan_avx2_256_32",                 "sg",    "scan", "avx2",  "256", "32",  8, 0, 0, 0},
{sg_scan_avx2_256_16,                 "sg_scan_avx2_256_16",                 "sg",    "scan", "avx2",  "256", "16", 16, 0, 0, 0},
{sg_scan_avx2_256_8,                  "sg_scan_avx2_256_8",                  "sg",    "scan", "avx2",  "256",  "8", 32, 0, 0, 0},
{sg_striped_avx2_256_64,              "sg_striped_avx2_256_64",              "sg", "striped", "avx2",  "256", "64",  4, 0, 0, 0},
{sg_striped_avx2_256_32,              "sg_striped_avx2_256_32",              "sg", "striped", "avx2",  "256", "32",  8, 0, 0, 0},
{sg_striped_avx2_256_16,              "sg_striped_avx2_256_16",              "sg", "striped", "avx2",  "256", "16", 16, 0, 0, 0},
{sg_striped_avx2_256_8,               "sg_striped_avx2_256_8",               "sg", "striped", "avx2",  "256",  "8", 32, 0, 0, 0},
{sg_diag_avx2_256_64,                 "sg_diag_avx2_256_64",                 "sg",    "diag", "avx2",  "256", "64",  4, 0, 0, 0},
{sg_diag_avx2_256_32,                 "sg_diag_avx2_256_32",                 "sg",    "diag", "avx2",  "256", "32",  8, 0, 0, 0},
{sg_diag_avx2_256_16,                 "sg_diag_avx2_256_16",                 "sg",    "diag", "avx2",  "256", "16", 16, 0, 0, 0},
{sg_diag_avx2_256_8,                  "sg_diag_avx2_256_8",                  "sg",    "diag", "avx2",  "256",  "8", 32, 0, 0, 0},
#endif
#if HAVE_KNC
{sg_scan_knc_512_32,                  "sg_scan_knc_512_32",                  "sg",    "scan", "knc",   "512", "32", 16, 0, 0, 0},
{sg_striped_knc_512_32,               "sg_striped_knc_512_32",               "sg", "striped", "knc",   "512", "32", 16, 0, 0, 0},
{sg_diag_knc_512_32,                  "sg_diag_knc_512_32",                  "sg",    "diag", "knc",   "512", "32", 16, 0, 0, 0},
#endif
{sw,                                  "sw",                                  "sw",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{sw_scan,                             "sw_scan",                             "sw",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
#if HAVE_SSE2
{sw_scan_sse2_128_64,                 "sw_scan_sse2_128_64",                 "sw",    "scan", "sse2",  "128", "64",  2, 0, 0, 0},
{sw_scan_sse2_128_32,                 "sw_scan_sse2_128_32",                 "sw",    "scan", "sse2",  "128", "32",  4, 0, 0, 0},
{sw_scan_sse2_128_16,                 "sw_scan_sse2_128_16",                 "sw",    "scan", "sse2",  "128", "16",  8, 0, 0, 0},
{sw_scan_sse2_128_8,                  "sw_scan_sse2_128_8",                  "sw",    "scan", "sse2",  "128",  "8", 16, 0, 0, 0},
{sw_striped_sse2_128_64,              "sw_striped_sse2_128_64",              "sw", "striped", "sse2",  "128", "64",  2, 0, 0, 0},
{sw_striped_sse2_128_32,              "sw_striped_sse2_128_32",              "sw", "striped", "sse2",  "128", "32",  4, 0, 0, 0},
{sw_striped_sse2_128_16,              "sw_striped_sse2_128_16",              "sw", "striped", "sse2",  "128", "16",  8, 0, 0, 0},
{sw_striped_sse2_128_8,               "sw_striped_sse2_128_8",               "sw", "striped", "sse2",  "128",  "8", 16, 0, 0, 0},
{sw_diag_sse2_128_64,                 "sw_diag_sse2_128_64",                 "sw",    "diag", "sse2",  "128", "64",  2, 0, 0, 0},
{sw_diag_sse2_128_32,                 "sw_diag_sse2_128_32",                 "sw",    "diag", "sse2",  "128", "32",  4, 0, 0, 0},
{sw_diag_sse2_128_16,                 "sw_diag_sse2_128_16",                 "sw",    "diag", "sse2",  "128", "16",  8, 0, 0, 0},
{sw_diag_sse2_128_8,                  "sw_diag_sse2_128_8",                  "sw",    "diag", "sse2",  "128",  "8", 16, 0, 0, 0},
#endif
#if HAVE_SSE41
{sw_scan_sse41_128_64,                "sw_scan_sse41_128_64",                "sw",    "scan", "sse41", "128", "64",  2, 0, 0, 0},
{sw_scan_sse41_128_32,                "sw_scan_sse41_128_32",                "sw",    "scan", "sse41", "128", "32",  4, 0, 0, 0},
{sw_scan_sse41_128_16,                "sw_scan_sse41_128_16",                "sw",    "scan", "sse41", "128", "16",  8, 0, 0, 0},
{sw_scan_sse41_128_8,                 "sw_scan_sse41_128_8",                 "sw",    "scan", "sse41", "128",  "8", 16, 0, 0, 0},
{sw_striped_sse41_128_64,             "sw_striped_sse41_128_64",             "sw", "striped", "sse41", "128", "64",  2, 0, 0, 0},
{sw_striped_sse41_128_32,             "sw_striped_sse41_128_32",             "sw", "striped", "sse41", "128", "32",  4, 0, 0, 0},
{sw_striped_sse41_128_16,             "sw_striped_sse41_128_16",             "sw", "striped", "sse41", "128", "16",  8, 0, 0, 0},
{sw_striped_sse41_128_8,              "sw_striped_sse41_128_8",              "sw", "striped", "sse41", "128",  "8", 16, 0, 0, 0},
{sw_diag_sse41_128_64,                "sw_diag_sse41_128_64",                "sw",    "diag", "sse41", "128", "64",  2, 0, 0, 0},
{sw_diag_sse41_128_32,                "sw_diag_sse41_128_32",                "sw",    "diag", "sse41", "128", "32",  4, 0, 0, 0},
{sw_diag_sse41_128_16,                "sw_diag_sse41_128_16",                "sw",    "diag", "sse41", "128", "16",  8, 0, 0, 0},
{sw_diag_sse41_128_8,                 "sw_diag_sse41_128_8",                 "sw",    "diag", "sse41", "128",  "8", 16, 0, 0, 0},
#endif
#if HAVE_AVX2
{sw_scan_avx2_256_64,                 "sw_scan_avx2_256_64",                 "sw",    "scan", "avx2",  "256", "64",  4, 0, 0, 0},
{sw_scan_avx2_256_32,                 "sw_scan_avx2_256_32",                 "sw",    "scan", "avx2",  "256", "32",  8, 0, 0, 0},
{sw_scan_avx2_256_16,                 "sw_scan_avx2_256_16",                 "sw",    "scan", "avx2",  "256", "16", 16, 0, 0, 0},
{sw_scan_avx2_256_8,                  "sw_scan_avx2_256_8",                  "sw",    "scan", "avx2",  "256",  "8", 32, 0, 0, 0},
{sw_striped_avx2_256_64,              "sw_striped_avx2_256_64",              "sw", "striped", "avx2",  "256", "64",  4, 0, 0, 0},
{sw_striped_avx2_256_32,              "sw_striped_avx2_256_32",              "sw", "striped", "avx2",  "256", "32",  8, 0, 0, 0},
{sw_striped_avx2_256_16,              "sw_striped_avx2_256_16",              "sw", "striped", "avx2",  "256", "16", 16, 0, 0, 0},
{sw_striped_avx2_256_8,               "sw_striped_avx2_256_8",               "sw", "striped", "avx2",  "256",  "8", 32, 0, 0, 0},
{sw_diag_avx2_256_64,                 "sw_diag_avx2_256_64",                 "sw",    "diag", "avx2",  "256", "64",  4, 0, 0, 0},
{sw_diag_avx2_256_32,                 "sw_diag_avx2_256_32",                 "sw",    "diag", "avx2",  "256", "32",  8, 0, 0, 0},
{sw_diag_avx2_256_16,                 "sw_diag_avx2_256_16",                 "sw",    "diag", "avx2",  "256", "16", 16, 0, 0, 0},
{sw_diag_avx2_256_8,                  "sw_diag_avx2_256_8",                  "sw",    "diag", "avx2",  "256",  "8", 32, 0, 0, 0},
#endif
#if HAVE_KNC
{sw_scan_knc_512_32,                  "sw_scan_knc_512_32",                  "sw",    "scan", "knc",   "512", "32", 16, 0, 0, 0},
{sw_striped_knc_512_32,               "sw_striped_knc_512_32",               "sw", "striped", "knc",   "512", "32", 16, 0, 0, 0},
{sw_diag_knc_512_32,                  "sw_diag_knc_512_32",                  "sw",    "diag", "knc",   "512", "32", 16, 0, 0, 0},
#endif
{nw_stats,                            "nw_stats",                            "nw_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{nw_stats_scan,                       "nw_stats_scan",                       "nw_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
#if HAVE_SSE2
{nw_stats_scan_sse2_128_32,           "nw_stats_scan_sse2_128_32",           "nw_stats",    "scan", "sse2",  "128", "32",  4, 0, 1, 0},
{nw_stats_scan_sse2_128_16,           "nw_stats_scan_sse2_128_16",           "nw_stats",    "scan", "sse2",  "128", "16",  8, 0, 1, 0},
{nw_stats_scan_sse2_128_8,            "nw_stats_scan_sse2_128_8",            "nw_stats",    "scan", "sse2",  "128",  "8", 16, 0, 1, 0},
{nw_stats_striped_sse2_128_32,        "nw_stats_striped_sse2_128_32",        "nw_stats", "striped", "sse2",  "128", "32",  4, 0, 1, 0},
{nw_stats_striped_sse2_128_16,        "nw_stats_striped_sse2_128_16",        "nw_stats", "striped", "sse2",  "128", "16",  8, 0, 1, 0},
{nw_stats_striped_sse2_128_8,         "nw_stats_striped_sse2_128_8",         "nw_stats", "striped", "sse2",  "128",  "8", 16, 0, 1, 0},
{nw_stats_diag_sse2_128_32,           "nw_stats_diag_sse2_128_32",           "nw_stats",    "diag", "sse2",  "128", "32",  4, 0, 1, 0},
{nw_stats_diag_sse2_128_16,           "nw_stats_diag_sse2_128_16",           "nw_stats",    "diag", "sse2",  "128", "16",  8, 0, 1, 0},
{nw_stats_diag_sse2_128_8,            "nw_stats_diag_sse2_128_8",            "nw_stats",    "diag", "sse2",  "128",  "8", 16, 0, 1, 0},
#endif
#if HAVE_SSE41
{nw_stats_scan_sse41_128_32,          "nw_stats_scan_sse41_128_32",          "nw_stats",    "scan", "sse41", "128", "32",  4, 0, 1, 0},
{nw_stats_scan_sse41_128_16,          "nw_stats_scan_sse41_128_16",          "nw_stats",    "scan", "sse41", "128", "16",  8, 0, 1, 0},
{nw_stats_scan_sse41_128_8,           "nw_stats_scan_sse41_128_8",           "nw_stats",    "scan", "sse41", "128",  "8", 16, 0, 1, 0},
{nw_stats_striped_sse41_128_32,       "nw_stats_striped_sse41_128_32",       "nw_stats", "striped", "sse41", "128", "32",  4, 0, 1, 0},
{nw_stats_striped_sse41_128_16,       "nw_stats_striped_sse41_128_16",       "nw_stats", "striped", "sse41", "128", "16",  8, 0, 1, 0},
{nw_stats_striped_sse41_128_8,        "nw_stats_striped_sse41_128_8",        "nw_stats", "striped", "sse41", "128",  "8", 16, 0, 1, 0},
{nw_stats_diag_sse41_128_32,          "nw_stats_diag_sse41_128_32",          "nw_stats",    "diag", "sse41", "128", "32",  4, 0, 1, 0},
{nw_stats_diag_sse41_128_16,          "nw_stats_diag_sse41_128_16",          "nw_stats",    "diag", "sse41", "128", "16",  8, 0, 1, 0},
{nw_stats_diag_sse41_128_8,           "nw_stats_diag_sse41_128_8",           "nw_stats",    "diag", "sse41", "128",  "8", 16, 0, 1, 0},
#endif
#if HAVE_AVX2
{nw_stats_scan_avx2_256_32,           "nw_stats_scan_avx2_256_32",           "nw_stats",    "scan", "avx2",  "256", "32",  8, 0, 1, 0},
{nw_stats_scan_avx2_256_16,           "nw_stats_scan_avx2_256_16",           "nw_stats",    "scan", "avx2",  "256", "16", 16, 0, 1, 0},
{nw_stats_scan_avx2_256_8,            "nw_stats_scan_avx2_256_8",            "nw_stats",    "scan", "avx2",  "256",  "8", 32, 0, 1, 0},
{nw_stats_striped_avx2_256_32,        "nw_stats_striped_avx2_256_32",        "nw_stats", "striped", "avx2",  "256", "32",  8, 0, 1, 0},
{nw_stats_striped_avx2_256_16,        "nw_stats_striped_avx2_256_16",        "nw_stats", "striped", "avx2",  "256", "16", 16, 0, 1, 0},
{nw_stats_striped_avx2_256_8,         "nw_stats_striped_avx2_256_8",         "nw_stats", "striped", "avx2",  "256",  "8", 32, 0, 1, 0},
{nw_stats_diag_avx2_256_32,           "nw_stats_diag_avx2_256_32",           "nw_stats",    "diag", "avx2",  "256", "32",  8, 0, 1, 0},
{nw_stats_diag_avx2_256_16,           "nw_stats_diag_avx2_256_16",           "nw_stats",    "diag", "avx2",  "256", "16", 16, 0, 1, 0},
{nw_stats_diag_avx2_256_8,            "nw_stats_diag_avx2_256_8",            "nw_stats",    "diag", "avx2",  "256",  "8", 32, 0, 1, 0},
#endif
#if HAVE_KNC
{nw_stats_scan_knc_512_32,            "nw_stats_scan_knc_512_32",            "nw_stats",    "scan", "knc",   "512", "32", 16, 0, 1, 0},
{nw_stats_striped_knc_512_32,         "nw_stats_striped_knc_512_32",         "nw_stats", "striped", "knc",   "512", "32", 16, 0, 1, 0},
{nw_stats_diag_knc_512_32,            "nw_stats_diag_knc_512_32",            "nw_stats",    "diag", "knc",   "512", "32", 16, 0, 1, 0},
#endif
{sg_stats,                            "sg_stats",                            "sg_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{sg_stats_scan,                       "sg_stats_scan",                       "sg_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
#if HAVE_SSE2
{sg_stats_scan_sse2_128_32,           "sg_stats_scan_sse2_128_32",           "sg_stats",    "scan", "sse2",  "128", "32",  4, 0, 1, 0},
{sg_stats_scan_sse2_128_16,           "sg_stats_scan_sse2_128_16",           "sg_stats",    "scan", "sse2",  "128", "16",  8, 0, 1, 0},
{sg_stats_scan_sse2_128_8,            "sg_stats_scan_sse2_128_8",            "sg_stats",    "scan", "sse2",  "128",  "8", 16, 0, 1, 0},
{sg_stats_striped_sse2_128_32,        "sg_stats_striped_sse2_128_32",        "sg_stats", "striped", "sse2",  "128", "32",  4, 0, 1, 0},
{sg_stats_striped_sse2_128_16,        "sg_stats_striped_sse2_128_16",        "sg_stats", "striped", "sse2",  "128", "16",  8, 0, 1, 0},
{sg_stats_striped_sse2_128_8,         "sg_stats_striped_sse2_128_8",         "sg_stats", "striped", "sse2",  "128",  "8", 16, 0, 1, 0},
{sg_stats_diag_sse2_128_32,           "sg_stats_diag_sse2_128_32",           "sg_stats",    "diag", "sse2",  "128", "32",  4, 0, 1, 0},
{sg_stats_diag_sse2_128_16,           "sg_stats_diag_sse2_128_16",           "sg_stats",    "diag", "sse2",  "128", "16",  8, 0, 1, 0},
{sg_stats_diag_sse2_128_8,            "sg_stats_diag_sse2_128_8",            "sg_stats",    "diag", "sse2",  "128",  "8", 16, 0, 1, 0},
#endif
#if HAVE_SSE41
{sg_stats_scan_sse41_128_32,          "sg_stats_scan_sse41_128_32",          "sg_stats",    "scan", "sse41", "128", "32",  4, 0, 1, 0},
{sg_stats_scan_sse41_128_16,          "sg_stats_scan_sse41_128_16",          "sg_stats",    "scan", "sse41", "128", "16",  8, 0, 1, 0},
{sg_stats_scan_sse41_128_8,           "sg_stats_scan_sse41_128_8",           "sg_stats",    "scan", "sse41", "128",  "8", 16, 0, 1, 0},
{sg_stats_striped_sse41_128_32,       "sg_stats_striped_sse41_128_32",       "sg_stats", "striped", "sse41", "128", "32",  4, 0, 1, 0},
{sg_stats_striped_sse41_128_16,       "sg_stats_striped_sse41_128_16",       "sg_stats", "striped", "sse41", "128", "16",  8, 0, 1, 0},
{sg_stats_striped_sse41_128_8,        "sg_stats_striped_sse41_128_8",        "sg_stats", "striped", "sse41", "128",  "8", 16, 0, 1, 0},
{sg_stats_diag_sse41_128_32,          "sg_stats_diag_sse41_128_32",          "sg_stats",    "diag", "sse41", "128", "32",  4, 0, 1, 0},
{sg_stats_diag_sse41_128_16,          "sg_stats_diag_sse41_128_16",          "sg_stats",    "diag", "sse41", "128", "16",  8, 0, 1, 0},
{sg_stats_diag_sse41_128_8,           "sg_stats_diag_sse41_128_8",           "sg_stats",    "diag", "sse41", "128",  "8", 16, 0, 1, 0},
#endif
#if HAVE_AVX2
{sg_stats_scan_avx2_256_32,           "sg_stats_scan_avx2_256_32",           "sg_stats",    "scan", "avx2",  "256", "32",  8, 0, 1, 0},
{sg_stats_scan_avx2_256_16,           "sg_stats_scan_avx2_256_16",           "sg_stats",    "scan", "avx2",  "256", "16", 16, 0, 1, 0},
{sg_stats_scan_avx2_256_8,            "sg_stats_scan_avx2_256_8",            "sg_stats",    "scan", "avx2",  "256",  "8", 32, 0, 1, 0},
{sg_stats_striped_avx2_256_32,        "sg_stats_striped_avx2_256_32",        "sg_stats", "striped", "avx2",  "256", "32",  8, 0, 1, 0},
{sg_stats_striped_avx2_256_16,        "sg_stats_striped_avx2_256_16",        "sg_stats", "striped", "avx2",  "256", "16", 16, 0, 1, 0},
{sg_stats_striped_avx2_256_8,         "sg_stats_striped_avx2_256_8",         "sg_stats", "striped", "avx2",  "256",  "8", 32, 0, 1, 0},
{sg_stats_diag_avx2_256_32,           "sg_stats_diag_avx2_256_32",           "sg_stats",    "diag", "avx2",  "256", "32",  8, 0, 1, 0},
{sg_stats_diag_avx2_256_16,           "sg_stats_diag_avx2_256_16",           "sg_stats",    "diag", "avx2",  "256", "16", 16, 0, 1, 0},
{sg_stats_diag_avx2_256_8,            "sg_stats_diag_avx2_256_8",            "sg_stats",    "diag", "avx2",  "256",  "8", 32, 0, 1, 0},
#endif
#if HAVE_KNC
{sg_stats_scan_knc_512_32,            "sg_stats_scan_knc_512_32",            "sg_stats",    "scan", "knc",   "512", "32", 16, 0, 1, 0},
{sg_stats_striped_knc_512_32,         "sg_stats_striped_knc_512_32",         "sg_stats", "striped", "knc",   "512", "32", 16, 0, 1, 0},
{sg_stats_diag_knc_512_32,            "sg_stats_diag_knc_512_32",            "sg_stats",    "diag", "knc",   "512", "32", 16, 0, 1, 0},
#endif
{sw_stats,                            "sw_stats",                            "sw_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{sw_stats_scan,                       "sw_stats_scan",                       "sw_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
#if HAVE_SSE2
{sw_stats_scan_sse2_128_32,           "sw_stats_scan_sse2_128_32",           "sw_stats",    "scan", "sse2",  "128", "32",  4, 0, 1, 0},
{sw_stats_scan_sse2_128_16,           "sw_stats_scan_sse2_128_16",           "sw_stats",    "scan", "sse2",  "128", "16",  8, 0, 1, 0},
{sw_stats_scan_sse2_128_8,            "sw_stats_scan_sse2_128_8",            "sw_stats",    "scan", "sse2",  "128",  "8", 16, 0, 1, 0},
{sw_stats_striped_sse2_128_32,        "sw_stats_striped_sse2_128_32",        "sw_stats", "striped", "sse2",  "128", "32",  4, 0, 1, 0},
{sw_stats_striped_sse2_128_16,        "sw_stats_striped_sse2_128_16",        "sw_stats", "striped", "sse2",  "128", "16",  8, 0, 1, 0},
{sw_stats_striped_sse2_128_8,         "sw_stats_striped_sse2_128_8",         "sw_stats", "striped", "sse2",  "128",  "8", 16, 0, 1, 0},
{sw_stats_diag_sse2_128_32,           "sw_stats_diag_sse2_128_32",           "sw_stats",    "diag", "sse2",  "128", "32",  4, 0, 1, 0},
{sw_stats_diag_sse2_128_16,           "sw_stats_diag_sse2_128_16",           "sw_stats",    "diag", "sse2",  "128", "16",  8, 0, 1, 0},
{sw_stats_diag_sse2_128_8,            "sw_stats_diag_sse2_128_8",            "sw_stats",    "diag", "sse2",  "128",  "8", 16, 0, 1, 0},
#endif
#if HAVE_SSE41
{sw_stats_scan_sse41_128_32,          "sw_stats_scan_sse41_128_32",          "sw_stats",    "scan", "sse41", "128", "32",  4, 0, 1, 0},
{sw_stats_scan_sse41_128_16,          "sw_stats_scan_sse41_128_16",          "sw_stats",    "scan", "sse41", "128", "16",  8, 0, 1, 0},
{sw_stats_scan_sse41_128_8,           "sw_stats_scan_sse41_128_8",           "sw_stats",    "scan", "sse41", "128",  "8", 16, 0, 1, 0},
{sw_stats_striped_sse41_128_32,       "sw_stats_striped_sse41_128_32",       "sw_stats", "striped", "sse41", "128", "32",  4, 0, 1, 0},
{sw_stats_striped_sse41_128_16,       "sw_stats_striped_sse41_128_16",       "sw_stats", "striped", "sse41", "128", "16",  8, 0, 1, 0},
{sw_stats_striped_sse41_128_8,        "sw_stats_striped_sse41_128_8",        "sw_stats", "striped", "sse41", "128",  "8", 16, 0, 1, 0},
{sw_stats_diag_sse41_128_32,          "sw_stats_diag_sse41_128_32",          "sw_stats",    "diag", "sse41", "128", "32",  4, 0, 1, 0},
{sw_stats_diag_sse41_128_16,          "sw_stats_diag_sse41_128_16",          "sw_stats",    "diag", "sse41", "128", "16",  8, 0, 1, 0},
{sw_stats_diag_sse41_128_8,           "sw_stats_diag_sse41_128_8",           "sw_stats",    "diag", "sse41", "128",  "8", 16, 0, 1, 0},
#endif
#if HAVE_AVX2
{sw_stats_scan_avx2_256_32,           "sw_stats_scan_avx2_256_32",           "sw_stats",    "scan", "avx2",  "256", "32",  8, 0, 1, 0},
{sw_stats_scan_avx2_256_16,           "sw_stats_scan_avx2_256_16",           "sw_stats",    "scan", "avx2",  "256", "16", 16, 0, 1, 0},
{sw_stats_scan_avx2_256_8,            "sw_stats_scan_avx2_256_8",            "sw_stats",    "scan", "avx2",  "256",  "8", 32, 0, 1, 0},
{sw_stats_striped_avx2_256_32,        "sw_stats_striped_avx2_256_32",        "sw_stats", "striped", "avx2",  "256", "32",  8, 0, 1, 0},
{sw_stats_striped_avx2_256_16,        "sw_stats_striped_avx2_256_16",        "sw_stats", "striped", "avx2",  "256", "16", 16, 0, 1, 0},
{sw_stats_striped_avx2_256_8,         "sw_stats_striped_avx2_256_8",         "sw_stats", "striped", "avx2",  "256",  "8", 32, 0, 1, 0},
{sw_stats_diag_avx2_256_32,           "sw_stats_diag_avx2_256_32",           "sw_stats",    "diag", "avx2",  "256", "32",  8, 0, 1, 0},
{sw_stats_diag_avx2_256_16,           "sw_stats_diag_avx2_256_16",           "sw_stats",    "diag", "avx2",  "256", "16", 16, 0, 1, 0},
{sw_stats_diag_avx2_256_8,            "sw_stats_diag_avx2_256_8",            "sw_stats",    "diag", "avx2",  "256",  "8", 32, 0, 1, 0},
#endif
#if HAVE_KNC
{sw_stats_scan_knc_512_32,            "sw_stats_scan_knc_512_32",            "sw_stats",    "scan", "knc",   "512", "32", 16, 0, 1, 0},
{sw_stats_striped_knc_512_32,         "sw_stats_striped_knc_512_32",         "sw_stats", "striped", "knc",   "512", "32", 16, 0, 1, 0},
{sw_stats_diag_knc_512_32,            "sw_stats_diag_knc_512_32",            "sw_stats",    "diag", "knc",   "512", "32", 16, 0, 1, 0},
#endif
{nw_table,                            "nw_table",                            "nw",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{nw_table_scan,                       "nw_table_scan",                       "nw",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
#if HAVE_SSE2
{nw_table_scan_sse2_128_64,           "nw_table_scan_sse2_128_64",           "nw",    "scan", "sse2",  "128", "64",  2, 1, 0, 0},
{nw_table_scan_sse2_128_32,           "nw_table_scan_sse2_128_32",           "nw",    "scan", "sse2",  "128", "32",  4, 1, 0, 0},
{nw_table_scan_sse2_128_16,           "nw_table_scan_sse2_128_16",           "nw",    "scan", "sse2",  "128", "16",  8, 1, 0, 0},
{nw_table_scan_sse2_128_8,            "nw_table_scan_sse2_128_8",            "nw",    "scan", "sse2",  "128",  "8", 16, 1, 0, 0},
{nw_table_striped_sse2_128_64,        "nw_table_striped_sse2_128_64",        "nw", "striped", "sse2",  "128", "64",  2, 1, 0, 0},
{nw_table_striped_sse2_128_32,        "nw_table_striped_sse2_128_32",        "nw", "striped", "sse2",  "128", "32",  4, 1, 0, 0},
{nw_table_striped_sse2_128_16,        "nw_table_striped_sse2_128_16",        "nw", "striped", "sse2",  "128", "16",  8, 1, 0, 0},
{nw_table_striped_sse2_128_8,         "nw_table_striped_sse2_128_8",         "nw", "striped", "sse2",  "128",  "8", 16, 1, 0, 0},
{nw_table_diag_sse2_128_64,           "nw_table_diag_sse2_128_64",           "nw",    "diag", "sse2",  "128", "64",  2, 1, 0, 0},
{nw_table_diag_sse2_128_32,           "nw_table_diag_sse2_128_32",           "nw",    "diag", "sse2",  "128", "32",  4, 1, 0, 0},
{nw_table_diag_sse2_128_16,           "nw_table_diag_sse2_128_16",           "nw",    "diag", "sse2",  "128", "16",  8, 1, 0, 0},
{nw_table_diag_sse2_128_8,            "nw_table_diag_sse2_128_8",            "nw",    "diag", "sse2",  "128",  "8", 16, 1, 0, 0},
#endif
#if HAVE_SSE41
{nw_table_scan_sse41_128_64,          "nw_table_scan_sse41_128_64",          "nw",    "scan", "sse41", "128", "64",  2, 1, 0, 0},
{nw_table_scan_sse41_128_32,          "nw_table_scan_sse41_128_32",          "nw",    "scan", "sse41", "128", "32",  4, 1, 0, 0},
{nw_table_scan_sse41_128_16,          "nw_table_scan_sse41_128_16",          "nw",    "scan", "sse41", "128", "16",  8, 1, 0, 0},
{nw_table_scan_sse41_128_8,           "nw_table_scan_sse41_128_8",           "nw",    "scan", "sse41", "128",  "8", 16, 1, 0, 0},
{nw_table_striped_sse41_128_64,       "nw_table_striped_sse41_128_64",       "nw", "striped", "sse41", "128", "64",  2, 1, 0, 0},
{nw_table_striped_sse41_128_32,       "nw_table_striped_sse41_128_32",       "nw", "striped", "sse41", "128", "32",  4, 1, 0, 0},
{nw_table_striped_sse41_128_16,       "nw_table_striped_sse41_128_16",       "nw", "striped", "sse41", "128", "16",  8, 1, 0, 0},
{nw_table_striped_sse41_128_8,        "nw_table_striped_sse41_128_8",        "nw", "striped", "sse41", "128",  "8", 16, 1, 0, 0},
{nw_table_diag_sse41_128_64,          "nw_table_diag_sse41_128_64",          "nw",    "diag", "sse41", "128", "64",  2, 1, 0, 0},
{nw_table_diag_sse41_128_32,          "nw_table_diag_sse41_128_32",          "nw",    "diag", "sse41", "128", "32",  4, 1, 0, 0},
{nw_table_diag_sse41_128_16,          "nw_table_diag_sse41_128_16",          "nw",    "diag", "sse41", "128", "16",  8, 1, 0, 0},
{nw_table_diag_sse41_128_8,           "nw_table_diag_sse41_128_8",           "nw",    "diag", "sse41", "128",  "8", 16, 1, 0, 0},
#endif
#if HAVE_AVX2
{nw_table_scan_avx2_256_64,           "nw_table_scan_avx2_256_64",           "nw",    "scan", "avx2",  "256", "64",  4, 1, 0, 0},
{nw_table_scan_avx2_256_32,           "nw_table_scan_avx2_256_32",           "nw",    "scan", "avx2",  "256", "32",  8, 1, 0, 0},
{nw_table_scan_avx2_256_16,           "nw_table_scan_avx2_256_16",           "nw",    "scan", "avx2",  "256", "16", 16, 1, 0, 0},
{nw_table_scan_avx2_256_8,            "nw_table_scan_avx2_256_8",            "nw",    "scan", "avx2",  "256",  "8", 32, 1, 0, 0},
{nw_table_striped_avx2_256_64,        "nw_table_striped_avx2_256_64",        "nw", "striped", "avx2",  "256", "64",  4, 1, 0, 0},
{nw_table_striped_avx2_256_32,        "nw_table_striped_avx2_256_32",        "nw", "striped", "avx2",  "256", "32",  8, 1, 0, 0},
{nw_table_striped_avx2_256_16,        "nw_table_striped_avx2_256_16",        "nw", "striped", "avx2",  "256", "16", 16, 1, 0, 0},
{nw_table_striped_avx2_256_8,         "nw_table_striped_avx2_256_8",         "nw", "striped", "avx2",  "256",  "8", 32, 1, 0, 0},
{nw_table_diag_avx2_256_64,           "nw_table_diag_avx2_256_64",           "nw",    "diag", "avx2",  "256", "64",  4, 1, 0, 0},
{nw_table_diag_avx2_256_32,           "nw_table_diag_avx2_256_32",           "nw",    "diag", "avx2",  "256", "32",  8, 1, 0, 0},
{nw_table_diag_avx2_256_16,           "nw_table_diag_avx2_256_16",           "nw",    "diag", "avx2",  "256", "16", 16, 1, 0, 0},
{nw_table_diag_avx2_256_8,            "nw_table_diag_avx2_256_8",            "nw",    "diag", "avx2",  "256",  "8", 32, 1, 0, 0},
#endif
#if HAVE_KNC
{nw_table_scan_knc_512_32,            "nw_table_scan_knc_512_32",            "nw",    "scan", "knc",   "512", "32", 16, 1, 0, 0},
{nw_table_striped_knc_512_32,         "nw_table_striped_knc_512_32",         "nw", "striped", "knc",   "512", "32", 16, 1, 0, 0},
{nw_table_diag_knc_512_32,            "nw_table_diag_knc_512_32",            "nw",    "diag", "knc",   "512", "32", 16, 1, 0, 0},
#endif
{sg_table,                            "sg_table",                            "sg",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{sg_table_scan,                       "sg_table_scan",                       "sg",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
#if HAVE_SSE2
{sg_table_scan_sse2_128_64,           "sg_table_scan_sse2_128_64",           "sg",    "scan", "sse2",  "128", "64",  2, 1, 0, 0},
{sg_table_scan_sse2_128_32,           "sg_table_scan_sse2_128_32",           "sg",    "scan", "sse2",  "128", "32",  4, 1, 0, 0},
{sg_table_scan_sse2_128_16,           "sg_table_scan_sse2_128_16",           "sg",    "scan", "sse2",  "128", "16",  8, 1, 0, 0},
{sg_table_scan_sse2_128_8,            "sg_table_scan_sse2_128_8",            "sg",    "scan", "sse2",  "128",  "8", 16, 1, 0, 0},
{sg_table_striped_sse2_128_64,        "sg_table_striped_sse2_128_64",        "sg", "striped", "sse2",  "128", "64",  2, 1, 0, 0},
{sg_table_striped_sse2_128_32,        "sg_table_striped_sse2_128_32",        "sg", "striped", "sse2",  "128", "32",  4, 1, 0, 0},
{sg_table_striped_sse2_128_16,        "sg_table_striped_sse2_128_16",        "sg", "striped", "sse2",  "128", "16",  8, 1, 0, 0},
{sg_table_striped_sse2_128_8,         "sg_table_striped_sse2_128_8",         "sg", "striped", "sse2",  "128",  "8", 16, 1, 0, 0},
{sg_table_diag_sse2_128_64,           "sg_table_diag_sse2_128_64",           "sg",    "diag", "sse2",  "128", "64",  2, 1, 0, 0},
{sg_table_diag_sse2_128_32,           "sg_table_diag_sse2_128_32",           "sg",    "diag", "sse2",  "128", "32",  4, 1, 0, 0},
{sg_table_diag_sse2_128_16,           "sg_table_diag_sse2_128_16",           "sg",    "diag", "sse2",  "128", "16",  8, 1, 0, 0},
{sg_table_diag_sse2_128_8,            "sg_table_diag_sse2_128_8",            "sg",    "diag", "sse2",  "128",  "8", 16, 1, 0, 0},
#endif
#if HAVE_SSE41
{sg_table_scan_sse41_128_64,          "sg_table_scan_sse41_128_64",          "sg",    "scan", "sse41", "128", "64",  2, 1, 0, 0},
{sg_table_scan_sse41_128_32,          "sg_table_scan_sse41_128_32",          "sg",    "scan", "sse41", "128", "32",  4, 1, 0, 0},
{sg_table_scan_sse41_128_16,          "sg_table_scan_sse41_128_16",          "sg",    "scan", "sse41", "128", "16",  8, 1, 0, 0},
{sg_table_scan_sse41_128_8,           "sg_table_scan_sse41_128_8",           "sg",    "scan", "sse41", "128",  "8", 16, 1, 0, 0},
{sg_table_striped_sse41_128_64,       "sg_table_striped_sse41_128_64",       "sg", "striped", "sse41", "128", "64",  2, 1, 0, 0},
{sg_table_striped_sse41_128_32,       "sg_table_striped_sse41_128_32",       "sg", "striped", "sse41", "128", "32",  4, 1, 0, 0},
{sg_table_striped_sse41_128_16,       "sg_table_striped_sse41_128_16",       "sg", "striped", "sse41", "128", "16",  8, 1, 0, 0},
{sg_table_striped_sse41_128_8,        "sg_table_striped_sse41_128_8",        "sg", "striped", "sse41", "128",  "8", 16, 1, 0, 0},
{sg_table_diag_sse41_128_64,          "sg_table_diag_sse41_128_64",          "sg",    "diag", "sse41", "128", "64",  2, 1, 0, 0},
{sg_table_diag_sse41_128_32,          "sg_table_diag_sse41_128_32",          "sg",    "diag", "sse41", "128", "32",  4, 1, 0, 0},
{sg_table_diag_sse41_128_16,          "sg_table_diag_sse41_128_16",          "sg",    "diag", "sse41", "128", "16",  8, 1, 0, 0},
{sg_table_diag_sse41_128_8,           "sg_table_diag_sse41_128_8",           "sg",    "diag", "sse41", "128",  "8", 16, 1, 0, 0},
#endif
#if HAVE_AVX2
{sg_table_scan_avx2_256_64,           "sg_table_scan_avx2_256_64",           "sg",    "scan", "avx2",  "256", "64",  4, 1, 0, 0},
{sg_table_scan_avx2_256_32,           "sg_table_scan_avx2_256_32",           "sg",    "scan", "avx2",  "256", "32",  8, 1, 0, 0},
{sg_table_scan_avx2_256_16,           "sg_table_scan_avx2_256_16",           "sg",    "scan", "avx2",  "256", "16", 16, 1, 0, 0},
{sg_table_scan_avx2_256_8,            "sg_table_scan_avx2_256_8",            "sg",    "scan", "avx2",  "256",  "8", 32, 1, 0, 0},
{sg_table_striped_avx2_256_64,        "sg_table_striped_avx2_256_64",        "sg", "striped", "avx2",  "256", "64",  4, 1, 0, 0},
{sg_table_striped_avx2_256_32,        "sg_table_striped_avx2_256_32",        "sg", "striped", "avx2",  "256", "32",  8, 1, 0, 0},
{sg_table_striped_avx2_256_16,        "sg_table_striped_avx2_256_16",        "sg", "striped", "avx2",  "256", "16", 16, 1, 0, 0},
{sg_table_striped_avx2_256_8,         "sg_table_striped_avx2_256_8",         "sg", "striped", "avx2",  "256",  "8", 32, 1, 0, 0},
{sg_table_diag_avx2_256_64,           "sg_table_diag_avx2_256_64",           "sg",    "diag", "avx2",  "256", "64",  4, 1, 0, 0},
{sg_table_diag_avx2_256_32,           "sg_table_diag_avx2_256_32",           "sg",    "diag", "avx2",  "256", "32",  8, 1, 0, 0},
{sg_table_diag_avx2_256_16,           "sg_table_diag_avx2_256_16",           "sg",    "diag", "avx2",  "256", "16", 16, 1, 0, 0},
{sg_table_diag_avx2_256_8,            "sg_table_diag_avx2_256_8",            "sg",    "diag", "avx2",  "256",  "8", 32, 1, 0, 0},
#endif
#if HAVE_KNC
{sg_table_scan_knc_512_32,            "sg_table_scan_knc_512_32",            "sg",    "scan", "knc",   "512", "32", 16, 1, 0, 0},
{sg_table_striped_knc_512_32,         "sg_table_striped_knc_512_32",         "sg", "striped", "knc",   "512", "32", 16, 1, 0, 0},
{sg_table_diag_knc_512_32,            "sg_table_diag_knc_512_32",            "sg",    "diag", "knc",   "512", "32", 16, 1, 0, 0},
#endif
{sw_table,                            "sw_table",                            "sw",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{sw_table_scan,                       "sw_table_scan",                       "sw",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
#if HAVE_SSE2
{sw_table_scan_sse2_128_64,           "sw_table_scan_sse2_128_64",           "sw",    "scan", "sse2",  "128", "64",  2, 1, 0, 0},
{sw_table_scan_sse2_128_32,           "sw_table_scan_sse2_128_32",           "sw",    "scan", "sse2",  "128", "32",  4, 1, 0, 0},
{sw_table_scan_sse2_128_16,           "sw_table_scan_sse2_128_16",           "sw",    "scan", "sse2",  "128", "16",  8, 1, 0, 0},
{sw_table_scan_sse2_128_8,            "sw_table_scan_sse2_128_8",            "sw",    "scan", "sse2",  "128",  "8", 16, 1, 0, 0},
{sw_table_striped_sse2_128_64,        "sw_table_striped_sse2_128_64",        "sw", "striped", "sse2",  "128", "64",  2, 1, 0, 0},
{sw_table_striped_sse2_128_32,        "sw_table_striped_sse2_128_32",        "sw", "striped", "sse2",  "128", "32",  4, 1, 0, 0},
{sw_table_striped_sse2_128_16,        "sw_table_striped_sse2_128_16",        "sw", "striped", "sse2",  "128", "16",  8, 1, 0, 0},
{sw_table_striped_sse2_128_8,         "sw_table_striped_sse2_128_8",         "sw", "striped", "sse2",  "128",  "8", 16, 1, 0, 0},
{sw_table_diag_sse2_128_64,           "sw_table_diag_sse2_128_64",           "sw",    "diag", "sse2",  "128", "64",  2, 1, 0, 0},
{sw_table_diag_sse2_128_32,           "sw_table_diag_sse2_128_32",           "sw",    "diag", "sse2",  "128", "32",  4, 1, 0, 0},
{sw_table_diag_sse2_128_16,           "sw_table_diag_sse2_128_16",           "sw",    "diag", "sse2",  "128", "16",  8, 1, 0, 0},
{sw_table_diag_sse2_128_8,            "sw_table_diag_sse2_128_8",            "sw",    "diag", "sse2",  "128",  "8", 16, 1, 0, 0},
#endif
#if HAVE_SSE41
{sw_table_scan_sse41_128_64,          "sw_table_scan_sse41_128_64",          "sw",    "scan", "sse41", "128", "64",  2, 1, 0, 0},
{sw_table_scan_sse41_128_32,          "sw_table_scan_sse41_128_32",          "sw",    "scan", "sse41", "128", "32",  4, 1, 0, 0},
{sw_table_scan_sse41_128_16,          "sw_table_scan_sse41_128_16",          "sw",    "scan", "sse41", "128", "16",  8, 1, 0, 0},
{sw_table_scan_sse41_128_8,           "sw_table_scan_sse41_128_8",           "sw",    "scan", "sse41", "128",  "8", 16, 1, 0, 0},
{sw_table_striped_sse41_128_64,       "sw_table_striped_sse41_128_64",       "sw", "striped", "sse41", "128", "64",  2, 1, 0, 0},
{sw_table_striped_sse41_128_32,       "sw_table_striped_sse41_128_32",       "sw", "striped", "sse41", "128", "32",  4, 1, 0, 0},
{sw_table_striped_sse41_128_16,       "sw_table_striped_sse41_128_16",       "sw", "striped", "sse41", "128", "16",  8, 1, 0, 0},
{sw_table_striped_sse41_128_8,        "sw_table_striped_sse41_128_8",        "sw", "striped", "sse41", "128",  "8", 16, 1, 0, 0},
{sw_table_diag_sse41_128_64,          "sw_table_diag_sse41_128_64",          "sw",    "diag", "sse41", "128", "64",  2, 1, 0, 0},
{sw_table_diag_sse41_128_32,          "sw_table_diag_sse41_128_32",          "sw",    "diag", "sse41", "128", "32",  4, 1, 0, 0},
{sw_table_diag_sse41_128_16,          "sw_table_diag_sse41_128_16",          "sw",    "diag", "sse41", "128", "16",  8, 1, 0, 0},
{sw_table_diag_sse41_128_8,           "sw_table_diag_sse41_128_8",           "sw",    "diag", "sse41", "128",  "8", 16, 1, 0, 0},
#endif
#if HAVE_AVX2
{sw_table_scan_avx2_256_64,           "sw_table_scan_avx2_256_64",           "sw",    "scan", "avx2",  "256", "64",  4, 1, 0, 0},
{sw_table_scan_avx2_256_32,           "sw_table_scan_avx2_256_32",           "sw",    "scan", "avx2",  "256", "32",  8, 1, 0, 0},
{sw_table_scan_avx2_256_16,           "sw_table_scan_avx2_256_16",           "sw",    "scan", "avx2",  "256", "16", 16, 1, 0, 0},
{sw_table_scan_avx2_256_8,            "sw_table_scan_avx2_256_8",            "sw",    "scan", "avx2",  "256",  "8", 32, 1, 0, 0},
{sw_table_striped_avx2_256_64,        "sw_table_striped_avx2_256_64",        "sw", "striped", "avx2",  "256", "64",  4, 1, 0, 0},
{sw_table_striped_avx2_256_32,        "sw_table_striped_avx2_256_32",        "sw", "striped", "avx2",  "256", "32",  8, 1, 0, 0},
{sw_table_striped_avx2_256_16,        "sw_table_striped_avx2_256_16",        "sw", "striped", "avx2",  "256", "16", 16, 1, 0, 0},
{sw_table_striped_avx2_256_8,         "sw_table_striped_avx2_256_8",         "sw", "striped", "avx2",  "256",  "8", 32, 1, 0, 0},
{sw_table_diag_avx2_256_64,           "sw_table_diag_avx2_256_64",           "sw",    "diag", "avx2",  "256", "64",  4, 1, 0, 0},
{sw_table_diag_avx2_256_32,           "sw_table_diag_avx2_256_32",           "sw",    "diag", "avx2",  "256", "32",  8, 1, 0, 0},
{sw_table_diag_avx2_256_16,           "sw_table_diag_avx2_256_16",           "sw",    "diag", "avx2",  "256", "16", 16, 1, 0, 0},
{sw_table_diag_avx2_256_8,            "sw_table_diag_avx2_256_8",            "sw",    "diag", "avx2",  "256",  "8", 32, 1, 0, 0},
#endif
#if HAVE_KNC
{sw_table_scan_knc_512_32,            "sw_table_scan_knc_512_32",            "sw",    "scan", "knc",   "512", "32", 16, 1, 0, 0},
{sw_table_striped_knc_512_32,         "sw_table_striped_knc_512_32",         "sw", "striped", "knc",   "512", "32", 16, 1, 0, 0},
{sw_table_diag_knc_512_32,            "sw_table_diag_knc_512_32",            "sw",    "diag", "knc",   "512", "32", 16, 1, 0, 0},
#endif
{nw_stats_table,                      "nw_stats_table",                      "nw_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{nw_stats_table_scan,                 "nw_stats_table_scan",                 "nw_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
#if HAVE_SSE2
{nw_stats_table_scan_sse2_128_32,     "nw_stats_table_scan_sse2_128_32",     "nw_stats",    "scan", "sse2",  "128", "32",  4, 1, 1, 0},
{nw_stats_table_scan_sse2_128_16,     "nw_stats_table_scan_sse2_128_16",     "nw_stats",    "scan", "sse2",  "128", "16",  8, 1, 1, 0},
{nw_stats_table_scan_sse2_128_8,      "nw_stats_table_scan_sse2_128_8",      "nw_stats",    "scan", "sse2",  "128",  "8", 16, 1, 1, 0},
{nw_stats_table_striped_sse2_128_32,  "nw_stats_table_striped_sse2_128_32",  "nw_stats", "striped", "sse2",  "128", "32",  4, 1, 1, 0},
{nw_stats_table_striped_sse2_128_16,  "nw_stats_table_striped_sse2_128_16",  "nw_stats", "striped", "sse2",  "128", "16",  8, 1, 1, 0},
{nw_stats_table_striped_sse2_128_8,   "nw_stats_table_striped_sse2_128_8",   "nw_stats", "striped", "sse2",  "128",  "8", 16, 1, 1, 0},
{nw_stats_table_diag_sse2_128_32,     "nw_stats_table_diag_sse2_128_32",     "nw_stats",    "diag", "sse2",  "128", "32",  4, 1, 1, 0},
{nw_stats_table_diag_sse2_128_16,     "nw_stats_table_diag_sse2_128_16",     "nw_stats",    "diag", "sse2",  "128", "16",  8, 1, 1, 0},
{nw_stats_table_diag_sse2_128_8,      "nw_stats_table_diag_sse2_128_8",      "nw_stats",    "diag", "sse2",  "128",  "8", 16, 1, 1, 0},
#endif
#if HAVE_SSE41
{nw_stats_table_scan_sse41_128_32,    "nw_stats_table_scan_sse41_128_32",    "nw_stats",    "scan", "sse41", "128", "32",  4, 1, 1, 0},
{nw_stats_table_scan_sse41_128_16,    "nw_stats_table_scan_sse41_128_16",    "nw_stats",    "scan", "sse41", "128", "16",  8, 1, 1, 0},
{nw_stats_table_scan_sse41_128_8,     "nw_stats_table_scan_sse41_128_8",     "nw_stats",    "scan", "sse41", "128",  "8", 16, 1, 1, 0},
{nw_stats_table_striped_sse41_128_32, "nw_stats_table_striped_sse41_128_32", "nw_stats", "striped", "sse41", "128", "32",  4, 1, 1, 0},
{nw_stats_table_striped_sse41_128_16, "nw_stats_table_striped_sse41_128_16", "nw_stats", "striped", "sse41", "128", "16",  8, 1, 1, 0},
{nw_stats_table_striped_sse41_128_8,  "nw_stats_table_striped_sse41_128_8",  "nw_stats", "striped", "sse41", "128",  "8", 16, 1, 1, 0},
{nw_stats_table_diag_sse41_128_32,    "nw_stats_table_diag_sse41_128_32",    "nw_stats",    "diag", "sse41", "128", "32",  4, 1, 1, 0},
{nw_stats_table_diag_sse41_128_16,    "nw_stats_table_diag_sse41_128_16",    "nw_stats",    "diag", "sse41", "128", "16",  8, 1, 1, 0},
{nw_stats_table_diag_sse41_128_8,     "nw_stats_table_diag_sse41_128_8",     "nw_stats",    "diag", "sse41", "128",  "8", 16, 1, 1, 0},
#endif
#if HAVE_AVX2
{nw_stats_table_scan_avx2_256_32,     "nw_stats_table_scan_avx2_256_32",     "nw_stats",    "scan", "avx2",  "256", "32",  8, 1, 1, 0},
{nw_stats_table_scan_avx2_256_16,     "nw_stats_table_scan_avx2_256_16",     "nw_stats",    "scan", "avx2",  "256", "16", 16, 1, 1, 0},
{nw_stats_table_scan_avx2_256_8,      "nw_stats_table_scan_avx2_256_8",      "nw_stats",    "scan", "avx2",  "256",  "8", 32, 1, 1, 0},
{nw_stats_table_striped_avx2_256_32,  "nw_stats_table_striped_avx2_256_32",  "nw_stats", "striped", "avx2",  "256", "32",  8, 1, 1, 0},
{nw_stats_table_striped_avx2_256_16,  "nw_stats_table_striped_avx2_256_16",  "nw_stats", "striped", "avx2",  "256", "16", 16, 1, 1, 0},
{nw_stats_table_striped_avx2_256_8,   "nw_stats_table_striped_avx2_256_8",   "nw_stats", "striped", "avx2",  "256",  "8", 32, 1, 1, 0},
{nw_stats_table_diag_avx2_256_32,     "nw_stats_table_diag_avx2_256_32",     "nw_stats",    "diag", "avx2",  "256", "32",  8, 1, 1, 0},
{nw_stats_table_diag_avx2_256_16,     "nw_stats_table_diag_avx2_256_16",     "nw_stats",    "diag", "avx2",  "256", "16", 16, 1, 1, 0},
{nw_stats_table_diag_avx2_256_8,      "nw_stats_table_diag_avx2_256_8",      "nw_stats",    "diag", "avx2",  "256",  "8", 32, 1, 1, 0},
#endif
#if HAVE_KNC
{nw_stats_table_scan_knc_512_32,      "nw_stats_table_scan_knc_512_32",      "nw_stats",    "scan", "knc",   "512", "32", 16, 1, 1, 0},
{nw_stats_table_striped_knc_512_32,   "nw_stats_table_striped_knc_512_32",   "nw_stats", "striped", "knc",   "512", "32", 16, 1, 1, 0},
{nw_stats_table_diag_knc_512_32,      "nw_stats_table_diag_knc_512_32",      "nw_stats",    "diag", "knc",   "512", "32", 16, 1, 1, 0},
#endif
{sg_stats_table,                      "sg_stats_table",                      "sg_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{sg_stats_table_scan,                 "sg_stats_table_scan",                 "sg_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
#if HAVE_SSE2
{sg_stats_table_scan_sse2_128_32,     "sg_stats_table_scan_sse2_128_32",     "sg_stats",    "scan", "sse2",  "128", "32",  4, 1, 1, 0},
{sg_stats_table_scan_sse2_128_16,     "sg_stats_table_scan_sse2_128_16",     "sg_stats",    "scan", "sse2",  "128", "16",  8, 1, 1, 0},
{sg_stats_table_scan_sse2_128_8,      "sg_stats_table_scan_sse2_128_8",      "sg_stats",    "scan", "sse2",  "128",  "8", 16, 1, 1, 0},
{sg_stats_table_striped_sse2_128_32,  "sg_stats_table_striped_sse2_128_32",  "sg_stats", "striped", "sse2",  "128", "32",  4, 1, 1, 0},
{sg_stats_table_striped_sse2_128_16,  "sg_stats_table_striped_sse2_128_16",  "sg_stats", "striped", "sse2",  "128", "16",  8, 1, 1, 0},
{sg_stats_table_striped_sse2_128_8,   "sg_stats_table_striped_sse2_128_8",   "sg_stats", "striped", "sse2",  "128",  "8", 16, 1, 1, 0},
{sg_stats_table_diag_sse2_128_32,     "sg_stats_table_diag_sse2_128_32",     "sg_stats",    "diag", "sse2",  "128", "32",  4, 1, 1, 0},
{sg_stats_table_diag_sse2_128_16,     "sg_stats_table_diag_sse2_128_16",     "sg_stats",    "diag", "sse2",  "128", "16",  8, 1, 1, 0},
{sg_stats_table_diag_sse2_128_8,      "sg_stats_table_diag_sse2_128_8",      "sg_stats",    "diag", "sse2",  "128",  "8", 16, 1, 1, 0},
#endif
#if HAVE_SSE41
{sg_stats_table_scan_sse41_128_32,    "sg_stats_table_scan_sse41_128_32",    "sg_stats",    "scan", "sse41", "128", "32",  4, 1, 1, 0},
{sg_stats_table_scan_sse41_128_16,    "sg_stats_table_scan_sse41_128_16",    "sg_stats",    "scan", "sse41", "128", "16",  8, 1, 1, 0},
{sg_stats_table_scan_sse41_128_8,     "sg_stats_table_scan_sse41_128_8",     "sg_stats",    "scan", "sse41", "128",  "8", 16, 1, 1, 0},
{sg_stats_table_striped_sse41_128_32, "sg_stats_table_striped_sse41_128_32", "sg_stats", "striped", "sse41", "128", "32",  4, 1, 1, 0},
{sg_stats_table_striped_sse41_128_16, "sg_stats_table_striped_sse41_128_16", "sg_stats", "striped", "sse41", "128", "16",  8, 1, 1, 0},
{sg_stats_table_striped_sse41_128_8,  "sg_stats_table_striped_sse41_128_8",  "sg_stats", "striped", "sse41", "128",  "8", 16, 1, 1, 0},
{sg_stats_table_diag_sse41_128_32,    "sg_stats_table_diag_sse41_128_32",    "sg_stats",    "diag", "sse41", "128", "32",  4, 1, 1, 0},
{sg_stats_table_diag_sse41_128_16,    "sg_stats_table_diag_sse41_128_16",    "sg_stats",    "diag", "sse41", "128", "16",  8, 1, 1, 0},
{sg_stats_table_diag_sse41_128_8,     "sg_stats_table_diag_sse41_128_8",     "sg_stats",    "diag", "sse41", "128",  "8", 16, 1, 1, 0},
#endif
#if HAVE_AVX2
{sg_stats_table_scan_avx2_256_32,     "sg_stats_table_scan_avx2_256_32",     "sg_stats",    "scan", "avx2",  "256", "32",  8, 1, 1, 0},
{sg_stats_table_scan_avx2_256_16,     "sg_stats_table_scan_avx2_256_16",     "sg_stats",    "scan", "avx2",  "256", "16", 16, 1, 1, 0},
{sg_stats_table_scan_avx2_256_8,      "sg_stats_table_scan_avx2_256_8",      "sg_stats",    "scan", "avx2",  "256",  "8", 32, 1, 1, 0},
{sg_stats_table_striped_avx2_256_32,  "sg_stats_table_striped_avx2_256_32",  "sg_stats", "striped", "avx2",  "256", "32",  8, 1, 1, 0},
{sg_stats_table_striped_avx2_256_16,  "sg_stats_table_striped_avx2_256_16",  "sg_stats", "striped", "avx2",  "256", "16", 16, 1, 1, 0},
{sg_stats_table_striped_avx2_256_8,   "sg_stats_table_striped_avx2_256_8",   "sg_stats", "striped", "avx2",  "256",  "8", 32, 1, 1, 0},
{sg_stats_table_diag_avx2_256_32,     "sg_stats_table_diag_avx2_256_32",     "sg_stats",    "diag", "avx2",  "256", "32",  8, 1, 1, 0},
{sg_stats_table_diag_avx2_256_16,     "sg_stats_table_diag_avx2_256_16",     "sg_stats",    "diag", "avx2",  "256", "16", 16, 1, 1, 0},
{sg_stats_table_diag_avx2_256_8,      "sg_stats_table_diag_avx2_256_8",      "sg_stats",    "diag", "avx2",  "256",  "8", 32, 1, 1, 0},
#endif
#if HAVE_KNC
{sg_stats_table_scan_knc_512_32,      "sg_stats_table_scan_knc_512_32",      "sg_stats",    "scan", "knc",   "512", "32", 16, 1, 1, 0},
{sg_stats_table_striped_knc_512_32,   "sg_stats_table_striped_knc_512_32",   "sg_stats", "striped", "knc",   "512", "32", 16, 1, 1, 0},
{sg_stats_table_diag_knc_512_32,      "sg_stats_table_diag_knc_512_32",      "sg_stats",    "diag", "knc",   "512", "32", 16, 1, 1, 0},
#endif
{sw_stats_table,                      "sw_stats_table",                      "sw_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{sw_stats_table_scan,                 "sw_stats_table_scan",                 "sw_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
#if HAVE_SSE2
{sw_stats_table_scan_sse2_128_32,     "sw_stats_table_scan_sse2_128_32",     "sw_stats",    "scan", "sse2",  "128", "32",  4, 1, 1, 0},
{sw_stats_table_scan_sse2_128_16,     "sw_stats_table_scan_sse2_128_16",     "sw_stats",    "scan", "sse2",  "128", "16",  8, 1, 1, 0},
{sw_stats_table_scan_sse2_128_8,      "sw_stats_table_scan_sse2_128_8",      "sw_stats",    "scan", "sse2",  "128",  "8", 16, 1, 1, 0},
{sw_stats_table_striped_sse2_128_32,  "sw_stats_table_striped_sse2_128_32",  "sw_stats", "striped", "sse2",  "128", "32",  4, 1, 1, 0},
{sw_stats_table_striped_sse2_128_16,  "sw_stats_table_striped_sse2_128_16",  "sw_stats", "striped", "sse2",  "128", "16",  8, 1, 1, 0},
{sw_stats_table_striped_sse2_128_8,   "sw_stats_table_striped_sse2_128_8",   "sw_stats", "striped", "sse2",  "128",  "8", 16, 1, 1, 0},
{sw_stats_table_diag_sse2_128_32,     "sw_stats_table_diag_sse2_128_32",     "sw_stats",    "diag", "sse2",  "128", "32",  4, 1, 1, 0},
{sw_stats_table_diag_sse2_128_16,     "sw_stats_table_diag_sse2_128_16",     "sw_stats",    "diag", "sse2",  "128", "16",  8, 1, 1, 0},
{sw_stats_table_diag_sse2_128_8,      "sw_stats_table_diag_sse2_128_8",      "sw_stats",    "diag", "sse2",  "128",  "8", 16, 1, 1, 0},
#endif
#if HAVE_SSE41
{sw_stats_table_scan_sse41_128_32,    "sw_stats_table_scan_sse41_128_32",    "sw_stats",    "scan", "sse41", "128", "32",  4, 1, 1, 0},
{sw_stats_table_scan_sse41_128_16,    "sw_stats_table_scan_sse41_128_16",    "sw_stats",    "scan", "sse41", "128", "16",  8, 1, 1, 0},
{sw_stats_table_scan_sse41_128_8,     "sw_stats_table_scan_sse41_128_8",     "sw_stats",    "scan", "sse41", "128",  "8", 16, 1, 1, 0},
{sw_stats_table_striped_sse41_128_32, "sw_stats_table_striped_sse41_128_32", "sw_stats", "striped", "sse41", "128", "32",  4, 1, 1, 0},
{sw_stats_table_striped_sse41_128_16, "sw_stats_table_striped_sse41_128_16", "sw_stats", "striped", "sse41", "128", "16",  8, 1, 1, 0},
{sw_stats_table_striped_sse41_128_8,  "sw_stats_table_striped_sse41_128_8",  "sw_stats", "striped", "sse41", "128",  "8", 16, 1, 1, 0},
{sw_stats_table_diag_sse41_128_32,    "sw_stats_table_diag_sse41_128_32",    "sw_stats",    "diag", "sse41", "128", "32",  4, 1, 1, 0},
{sw_stats_table_diag_sse41_128_16,    "sw_stats_table_diag_sse41_128_16",    "sw_stats",    "diag", "sse41", "128", "16",  8, 1, 1, 0},
{sw_stats_table_diag_sse41_128_8,     "sw_stats_table_diag_sse41_128_8",     "sw_stats",    "diag", "sse41", "128",  "8", 16, 1, 1, 0},
#endif
#if HAVE_AVX2
{sw_stats_table_scan_avx2_256_32,     "sw_stats_table_scan_avx2_256_32",     "sw_stats",    "scan", "avx2",  "256", "32",  8, 1, 1, 0},
{sw_stats_table_scan_avx2_256_16,     "sw_stats_table_scan_avx2_256_16",     "sw_stats",    "scan", "avx2",  "256", "16", 16, 1, 1, 0},
{sw_stats_table_scan_avx2_256_8,      "sw_stats_table_scan_avx2_256_8",      "sw_stats",    "scan", "avx2",  "256",  "8", 32, 1, 1, 0},
{sw_stats_table_striped_avx2_256_32,  "sw_stats_table_striped_avx2_256_32",  "sw_stats", "striped", "avx2",  "256", "32",  8, 1, 1, 0},
{sw_stats_table_striped_avx2_256_16,  "sw_stats_table_striped_avx2_256_16",  "sw_stats", "striped", "avx2",  "256", "16", 16, 1, 1, 0},
{sw_stats_table_striped_avx2_256_8,   "sw_stats_table_striped_avx2_256_8",   "sw_stats", "striped", "avx2",  "256",  "8", 32, 1, 1, 0},
{sw_stats_table_diag_avx2_256_32,     "sw_stats_table_diag_avx2_256_32",     "sw_stats",    "diag", "avx2",  "256", "32",  8, 1, 1, 0},
{sw_stats_table_diag_avx2_256_16,     "sw_stats_table_diag_avx2_256_16",     "sw_stats",    "diag", "avx2",  "256", "16", 16, 1, 1, 0},
{sw_stats_table_diag_avx2_256_8,      "sw_stats_table_diag_avx2_256_8",      "sw_stats",    "diag", "avx2",  "256",  "8", 32, 1, 1, 0},
#endif
#if HAVE_KNC
{sw_stats_table_scan_knc_512_32,      "sw_stats_table_scan_knc_512_32",      "sw_stats",    "scan", "knc",   "512", "32", 16, 1, 1, 0},
{sw_stats_table_striped_knc_512_32,   "sw_stats_table_striped_knc_512_32",   "sw_stats", "striped", "knc",   "512", "32", 16, 1, 1, 0},
{sw_stats_table_diag_knc_512_32,      "sw_stats_table_diag_knc_512_32",      "sw_stats",    "diag", "knc",   "512", "32", 16, 1, 1, 0},
#endif
{sw_blocked_sse41_128_32,             "sw_blocked_sse41_128_32",             "sw", "blocked", "sse41", "128", "32",  4, 0, 0, 0},
{sw_blocked_sse41_128_16,             "sw_blocked_sse41_128_16",             "sw", "blocked", "sse41", "128", "16",  8, 0, 0, 0},
{sw_table_blocked_sse41_128_32,       "sw_blocked_sse41_128_32",             "sw", "blocked", "sse41", "128", "32",  4, 1, 0, 0},
{sw_table_blocked_sse41_128_16,       "sw_blocked_sse41_128_16",             "sw", "blocked", "sse41", "128", "16",  8, 1, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
#if HAVE_SSE2
func_t nw_sse2_functions[] = {
{nw,                                  "nw",                                  "nw",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{nw_scan,                             "nw_scan",                             "nw",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
{nw_scan_sse2_128_64,                 "nw_scan_sse2_128_64",                 "nw",    "scan", "sse2",  "128", "64",  2, 0, 0, 0},
{nw_scan_sse2_128_32,                 "nw_scan_sse2_128_32",                 "nw",    "scan", "sse2",  "128", "32",  4, 0, 0, 0},
{nw_scan_sse2_128_16,                 "nw_scan_sse2_128_16",                 "nw",    "scan", "sse2",  "128", "16",  8, 0, 0, 0},
{nw_scan_sse2_128_8,                  "nw_scan_sse2_128_8",                  "nw",    "scan", "sse2",  "128",  "8", 16, 0, 0, 0},
{nw_striped_sse2_128_64,              "nw_striped_sse2_128_64",              "nw", "striped", "sse2",  "128", "64",  2, 0, 0, 0},
{nw_striped_sse2_128_32,              "nw_striped_sse2_128_32",              "nw", "striped", "sse2",  "128", "32",  4, 0, 0, 0},
{nw_striped_sse2_128_16,              "nw_striped_sse2_128_16",              "nw", "striped", "sse2",  "128", "16",  8, 0, 0, 0},
{nw_striped_sse2_128_8,               "nw_striped_sse2_128_8",               "nw", "striped", "sse2",  "128",  "8", 16, 0, 0, 0},
{nw_diag_sse2_128_64,                 "nw_diag_sse2_128_64",                 "nw",    "diag", "sse2",  "128", "64",  2, 0, 0, 0},
{nw_diag_sse2_128_32,                 "nw_diag_sse2_128_32",                 "nw",    "diag", "sse2",  "128", "32",  4, 0, 0, 0},
{nw_diag_sse2_128_16,                 "nw_diag_sse2_128_16",                 "nw",    "diag", "sse2",  "128", "16",  8, 0, 0, 0},
{nw_diag_sse2_128_8,                  "nw_diag_sse2_128_8",                  "nw",    "diag", "sse2",  "128",  "8", 16, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_sse2 = {"nw_sse2", nw_sse2_functions};
#endif
#if HAVE_SSE41
func_t nw_sse41_functions[] = {
{nw,                                  "nw",                                  "nw",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{nw_scan,                             "nw_scan",                             "nw",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
{nw_scan_sse41_128_64,                "nw_scan_sse41_128_64",                "nw",    "scan", "sse41", "128", "64",  2, 0, 0, 0},
{nw_scan_sse41_128_32,                "nw_scan_sse41_128_32",                "nw",    "scan", "sse41", "128", "32",  4, 0, 0, 0},
{nw_scan_sse41_128_16,                "nw_scan_sse41_128_16",                "nw",    "scan", "sse41", "128", "16",  8, 0, 0, 0},
{nw_scan_sse41_128_8,                 "nw_scan_sse41_128_8",                 "nw",    "scan", "sse41", "128",  "8", 16, 0, 0, 0},
{nw_striped_sse41_128_64,             "nw_striped_sse41_128_64",             "nw", "striped", "sse41", "128", "64",  2, 0, 0, 0},
{nw_striped_sse41_128_32,             "nw_striped_sse41_128_32",             "nw", "striped", "sse41", "128", "32",  4, 0, 0, 0},
{nw_striped_sse41_128_16,             "nw_striped_sse41_128_16",             "nw", "striped", "sse41", "128", "16",  8, 0, 0, 0},
{nw_striped_sse41_128_8,              "nw_striped_sse41_128_8",              "nw", "striped", "sse41", "128",  "8", 16, 0, 0, 0},
{nw_diag_sse41_128_64,                "nw_diag_sse41_128_64",                "nw",    "diag", "sse41", "128", "64",  2, 0, 0, 0},
{nw_diag_sse41_128_32,                "nw_diag_sse41_128_32",                "nw",    "diag", "sse41", "128", "32",  4, 0, 0, 0},
{nw_diag_sse41_128_16,                "nw_diag_sse41_128_16",                "nw",    "diag", "sse41", "128", "16",  8, 0, 0, 0},
{nw_diag_sse41_128_8,                 "nw_diag_sse41_128_8",                 "nw",    "diag", "sse41", "128",  "8", 16, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_sse41 = {"nw_sse41", nw_sse41_functions};
#endif
#if HAVE_AVX2
func_t nw_avx2_functions[] = {
{nw,                                  "nw",                                  "nw",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{nw_scan,                             "nw_scan",                             "nw",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
{nw_scan_avx2_256_64,                 "nw_scan_avx2_256_64",                 "nw",    "scan", "avx2",  "256", "64",  4, 0, 0, 0},
{nw_scan_avx2_256_32,                 "nw_scan_avx2_256_32",                 "nw",    "scan", "avx2",  "256", "32",  8, 0, 0, 0},
{nw_scan_avx2_256_16,                 "nw_scan_avx2_256_16",                 "nw",    "scan", "avx2",  "256", "16", 16, 0, 0, 0},
{nw_scan_avx2_256_8,                  "nw_scan_avx2_256_8",                  "nw",    "scan", "avx2",  "256",  "8", 32, 0, 0, 0},
{nw_striped_avx2_256_64,              "nw_striped_avx2_256_64",              "nw", "striped", "avx2",  "256", "64",  4, 0, 0, 0},
{nw_striped_avx2_256_32,              "nw_striped_avx2_256_32",              "nw", "striped", "avx2",  "256", "32",  8, 0, 0, 0},
{nw_striped_avx2_256_16,              "nw_striped_avx2_256_16",              "nw", "striped", "avx2",  "256", "16", 16, 0, 0, 0},
{nw_striped_avx2_256_8,               "nw_striped_avx2_256_8",               "nw", "striped", "avx2",  "256",  "8", 32, 0, 0, 0},
{nw_diag_avx2_256_64,                 "nw_diag_avx2_256_64",                 "nw",    "diag", "avx2",  "256", "64",  4, 0, 0, 0},
{nw_diag_avx2_256_32,                 "nw_diag_avx2_256_32",                 "nw",    "diag", "avx2",  "256", "32",  8, 0, 0, 0},
{nw_diag_avx2_256_16,                 "nw_diag_avx2_256_16",                 "nw",    "diag", "avx2",  "256", "16", 16, 0, 0, 0},
{nw_diag_avx2_256_8,                  "nw_diag_avx2_256_8",                  "nw",    "diag", "avx2",  "256",  "8", 32, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_avx2 = {"nw_avx2", nw_avx2_functions};
#endif
#if HAVE_KNC
func_t nw_knc_functions[] = {
{nw,                                  "nw",                                  "nw",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{nw_scan,                             "nw_scan",                             "nw",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
{nw_scan_knc_512_32,                  "nw_scan_knc_512_32",                  "nw",    "scan", "knc",   "512", "32", 16, 0, 0, 0},
{nw_striped_knc_512_32,               "nw_striped_knc_512_32",               "nw", "striped", "knc",   "512", "32", 16, 0, 0, 0},
{nw_diag_knc_512_32,                  "nw_diag_knc_512_32",                  "nw",    "diag", "knc",   "512", "32", 16, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_knc = {"nw_knc", nw_knc_functions};
#endif
#if HAVE_SSE2
func_t sg_sse2_functions[] = {
{sg,                                  "sg",                                  "sg",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{sg_scan,                             "sg_scan",                             "sg",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
{sg_scan_sse2_128_64,                 "sg_scan_sse2_128_64",                 "sg",    "scan", "sse2",  "128", "64",  2, 0, 0, 0},
{sg_scan_sse2_128_32,                 "sg_scan_sse2_128_32",                 "sg",    "scan", "sse2",  "128", "32",  4, 0, 0, 0},
{sg_scan_sse2_128_16,                 "sg_scan_sse2_128_16",                 "sg",    "scan", "sse2",  "128", "16",  8, 0, 0, 0},
{sg_scan_sse2_128_8,                  "sg_scan_sse2_128_8",                  "sg",    "scan", "sse2",  "128",  "8", 16, 0, 0, 0},
{sg_striped_sse2_128_64,              "sg_striped_sse2_128_64",              "sg", "striped", "sse2",  "128", "64",  2, 0, 0, 0},
{sg_striped_sse2_128_32,              "sg_striped_sse2_128_32",              "sg", "striped", "sse2",  "128", "32",  4, 0, 0, 0},
{sg_striped_sse2_128_16,              "sg_striped_sse2_128_16",              "sg", "striped", "sse2",  "128", "16",  8, 0, 0, 0},
{sg_striped_sse2_128_8,               "sg_striped_sse2_128_8",               "sg", "striped", "sse2",  "128",  "8", 16, 0, 0, 0},
{sg_diag_sse2_128_64,                 "sg_diag_sse2_128_64",                 "sg",    "diag", "sse2",  "128", "64",  2, 0, 0, 0},
{sg_diag_sse2_128_32,                 "sg_diag_sse2_128_32",                 "sg",    "diag", "sse2",  "128", "32",  4, 0, 0, 0},
{sg_diag_sse2_128_16,                 "sg_diag_sse2_128_16",                 "sg",    "diag", "sse2",  "128", "16",  8, 0, 0, 0},
{sg_diag_sse2_128_8,                  "sg_diag_sse2_128_8",                  "sg",    "diag", "sse2",  "128",  "8", 16, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_sse2 = {"sg_sse2", sg_sse2_functions};
#endif
#if HAVE_SSE41
func_t sg_sse41_functions[] = {
{sg,                                  "sg",                                  "sg",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{sg_scan,                             "sg_scan",                             "sg",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
{sg_scan_sse41_128_64,                "sg_scan_sse41_128_64",                "sg",    "scan", "sse41", "128", "64",  2, 0, 0, 0},
{sg_scan_sse41_128_32,                "sg_scan_sse41_128_32",                "sg",    "scan", "sse41", "128", "32",  4, 0, 0, 0},
{sg_scan_sse41_128_16,                "sg_scan_sse41_128_16",                "sg",    "scan", "sse41", "128", "16",  8, 0, 0, 0},
{sg_scan_sse41_128_8,                 "sg_scan_sse41_128_8",                 "sg",    "scan", "sse41", "128",  "8", 16, 0, 0, 0},
{sg_striped_sse41_128_64,             "sg_striped_sse41_128_64",             "sg", "striped", "sse41", "128", "64",  2, 0, 0, 0},
{sg_striped_sse41_128_32,             "sg_striped_sse41_128_32",             "sg", "striped", "sse41", "128", "32",  4, 0, 0, 0},
{sg_striped_sse41_128_16,             "sg_striped_sse41_128_16",             "sg", "striped", "sse41", "128", "16",  8, 0, 0, 0},
{sg_striped_sse41_128_8,              "sg_striped_sse41_128_8",              "sg", "striped", "sse41", "128",  "8", 16, 0, 0, 0},
{sg_diag_sse41_128_64,                "sg_diag_sse41_128_64",                "sg",    "diag", "sse41", "128", "64",  2, 0, 0, 0},
{sg_diag_sse41_128_32,                "sg_diag_sse41_128_32",                "sg",    "diag", "sse41", "128", "32",  4, 0, 0, 0},
{sg_diag_sse41_128_16,                "sg_diag_sse41_128_16",                "sg",    "diag", "sse41", "128", "16",  8, 0, 0, 0},
{sg_diag_sse41_128_8,                 "sg_diag_sse41_128_8",                 "sg",    "diag", "sse41", "128",  "8", 16, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_sse41 = {"sg_sse41", sg_sse41_functions};
#endif
#if HAVE_AVX2
func_t sg_avx2_functions[] = {
{sg,                                  "sg",                                  "sg",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{sg_scan,                             "sg_scan",                             "sg",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
{sg_scan_avx2_256_64,                 "sg_scan_avx2_256_64",                 "sg",    "scan", "avx2",  "256", "64",  4, 0, 0, 0},
{sg_scan_avx2_256_32,                 "sg_scan_avx2_256_32",                 "sg",    "scan", "avx2",  "256", "32",  8, 0, 0, 0},
{sg_scan_avx2_256_16,                 "sg_scan_avx2_256_16",                 "sg",    "scan", "avx2",  "256", "16", 16, 0, 0, 0},
{sg_scan_avx2_256_8,                  "sg_scan_avx2_256_8",                  "sg",    "scan", "avx2",  "256",  "8", 32, 0, 0, 0},
{sg_striped_avx2_256_64,              "sg_striped_avx2_256_64",              "sg", "striped", "avx2",  "256", "64",  4, 0, 0, 0},
{sg_striped_avx2_256_32,              "sg_striped_avx2_256_32",              "sg", "striped", "avx2",  "256", "32",  8, 0, 0, 0},
{sg_striped_avx2_256_16,              "sg_striped_avx2_256_16",              "sg", "striped", "avx2",  "256", "16", 16, 0, 0, 0},
{sg_striped_avx2_256_8,               "sg_striped_avx2_256_8",               "sg", "striped", "avx2",  "256",  "8", 32, 0, 0, 0},
{sg_diag_avx2_256_64,                 "sg_diag_avx2_256_64",                 "sg",    "diag", "avx2",  "256", "64",  4, 0, 0, 0},
{sg_diag_avx2_256_32,                 "sg_diag_avx2_256_32",                 "sg",    "diag", "avx2",  "256", "32",  8, 0, 0, 0},
{sg_diag_avx2_256_16,                 "sg_diag_avx2_256_16",                 "sg",    "diag", "avx2",  "256", "16", 16, 0, 0, 0},
{sg_diag_avx2_256_8,                  "sg_diag_avx2_256_8",                  "sg",    "diag", "avx2",  "256",  "8", 32, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_avx2 = {"sg_avx2", sg_avx2_functions};
#endif
#if HAVE_KNC
func_t sg_knc_functions[] = {
{sg,                                  "sg",                                  "sg",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{sg_scan,                             "sg_scan",                             "sg",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
{sg_scan_knc_512_32,                  "sg_scan_knc_512_32",                  "sg",    "scan", "knc",   "512", "32", 16, 0, 0, 0},
{sg_striped_knc_512_32,               "sg_striped_knc_512_32",               "sg", "striped", "knc",   "512", "32", 16, 0, 0, 0},
{sg_diag_knc_512_32,                  "sg_diag_knc_512_32",                  "sg",    "diag", "knc",   "512", "32", 16, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_knc = {"sg_knc", sg_knc_functions};
#endif
#if HAVE_SSE2
func_t sw_sse2_functions[] = {
{sw,                                  "sw",                                  "sw",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{sw_scan,                             "sw_scan",                             "sw",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
{sw_scan_sse2_128_64,                 "sw_scan_sse2_128_64",                 "sw",    "scan", "sse2",  "128", "64",  2, 0, 0, 0},
{sw_scan_sse2_128_32,                 "sw_scan_sse2_128_32",                 "sw",    "scan", "sse2",  "128", "32",  4, 0, 0, 0},
{sw_scan_sse2_128_16,                 "sw_scan_sse2_128_16",                 "sw",    "scan", "sse2",  "128", "16",  8, 0, 0, 0},
{sw_scan_sse2_128_8,                  "sw_scan_sse2_128_8",                  "sw",    "scan", "sse2",  "128",  "8", 16, 0, 0, 0},
{sw_striped_sse2_128_64,              "sw_striped_sse2_128_64",              "sw", "striped", "sse2",  "128", "64",  2, 0, 0, 0},
{sw_striped_sse2_128_32,              "sw_striped_sse2_128_32",              "sw", "striped", "sse2",  "128", "32",  4, 0, 0, 0},
{sw_striped_sse2_128_16,              "sw_striped_sse2_128_16",              "sw", "striped", "sse2",  "128", "16",  8, 0, 0, 0},
{sw_striped_sse2_128_8,               "sw_striped_sse2_128_8",               "sw", "striped", "sse2",  "128",  "8", 16, 0, 0, 0},
{sw_diag_sse2_128_64,                 "sw_diag_sse2_128_64",                 "sw",    "diag", "sse2",  "128", "64",  2, 0, 0, 0},
{sw_diag_sse2_128_32,                 "sw_diag_sse2_128_32",                 "sw",    "diag", "sse2",  "128", "32",  4, 0, 0, 0},
{sw_diag_sse2_128_16,                 "sw_diag_sse2_128_16",                 "sw",    "diag", "sse2",  "128", "16",  8, 0, 0, 0},
{sw_diag_sse2_128_8,                  "sw_diag_sse2_128_8",                  "sw",    "diag", "sse2",  "128",  "8", 16, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_sse2 = {"sw_sse2", sw_sse2_functions};
#endif
#if HAVE_SSE41
func_t sw_sse41_functions[] = {
{sw,                                  "sw",                                  "sw",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{sw_scan,                             "sw_scan",                             "sw",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
{sw_scan_sse41_128_64,                "sw_scan_sse41_128_64",                "sw",    "scan", "sse41", "128", "64",  2, 0, 0, 0},
{sw_scan_sse41_128_32,                "sw_scan_sse41_128_32",                "sw",    "scan", "sse41", "128", "32",  4, 0, 0, 0},
{sw_scan_sse41_128_16,                "sw_scan_sse41_128_16",                "sw",    "scan", "sse41", "128", "16",  8, 0, 0, 0},
{sw_scan_sse41_128_8,                 "sw_scan_sse41_128_8",                 "sw",    "scan", "sse41", "128",  "8", 16, 0, 0, 0},
{sw_striped_sse41_128_64,             "sw_striped_sse41_128_64",             "sw", "striped", "sse41", "128", "64",  2, 0, 0, 0},
{sw_striped_sse41_128_32,             "sw_striped_sse41_128_32",             "sw", "striped", "sse41", "128", "32",  4, 0, 0, 0},
{sw_striped_sse41_128_16,             "sw_striped_sse41_128_16",             "sw", "striped", "sse41", "128", "16",  8, 0, 0, 0},
{sw_striped_sse41_128_8,              "sw_striped_sse41_128_8",              "sw", "striped", "sse41", "128",  "8", 16, 0, 0, 0},
{sw_diag_sse41_128_64,                "sw_diag_sse41_128_64",                "sw",    "diag", "sse41", "128", "64",  2, 0, 0, 0},
{sw_diag_sse41_128_32,                "sw_diag_sse41_128_32",                "sw",    "diag", "sse41", "128", "32",  4, 0, 0, 0},
{sw_diag_sse41_128_16,                "sw_diag_sse41_128_16",                "sw",    "diag", "sse41", "128", "16",  8, 0, 0, 0},
{sw_diag_sse41_128_8,                 "sw_diag_sse41_128_8",                 "sw",    "diag", "sse41", "128",  "8", 16, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_sse41 = {"sw_sse41", sw_sse41_functions};
#endif
#if HAVE_AVX2
func_t sw_avx2_functions[] = {
{sw,                                  "sw",                                  "sw",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{sw_scan,                             "sw_scan",                             "sw",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
{sw_scan_avx2_256_64,                 "sw_scan_avx2_256_64",                 "sw",    "scan", "avx2",  "256", "64",  4, 0, 0, 0},
{sw_scan_avx2_256_32,                 "sw_scan_avx2_256_32",                 "sw",    "scan", "avx2",  "256", "32",  8, 0, 0, 0},
{sw_scan_avx2_256_16,                 "sw_scan_avx2_256_16",                 "sw",    "scan", "avx2",  "256", "16", 16, 0, 0, 0},
{sw_scan_avx2_256_8,                  "sw_scan_avx2_256_8",                  "sw",    "scan", "avx2",  "256",  "8", 32, 0, 0, 0},
{sw_striped_avx2_256_64,              "sw_striped_avx2_256_64",              "sw", "striped", "avx2",  "256", "64",  4, 0, 0, 0},
{sw_striped_avx2_256_32,              "sw_striped_avx2_256_32",              "sw", "striped", "avx2",  "256", "32",  8, 0, 0, 0},
{sw_striped_avx2_256_16,              "sw_striped_avx2_256_16",              "sw", "striped", "avx2",  "256", "16", 16, 0, 0, 0},
{sw_striped_avx2_256_8,               "sw_striped_avx2_256_8",               "sw", "striped", "avx2",  "256",  "8", 32, 0, 0, 0},
{sw_diag_avx2_256_64,                 "sw_diag_avx2_256_64",                 "sw",    "diag", "avx2",  "256", "64",  4, 0, 0, 0},
{sw_diag_avx2_256_32,                 "sw_diag_avx2_256_32",                 "sw",    "diag", "avx2",  "256", "32",  8, 0, 0, 0},
{sw_diag_avx2_256_16,                 "sw_diag_avx2_256_16",                 "sw",    "diag", "avx2",  "256", "16", 16, 0, 0, 0},
{sw_diag_avx2_256_8,                  "sw_diag_avx2_256_8",                  "sw",    "diag", "avx2",  "256",  "8", 32, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_avx2 = {"sw_avx2", sw_avx2_functions};
#endif
#if HAVE_KNC
func_t sw_knc_functions[] = {
{sw,                                  "sw",                                  "sw",    "orig", "NA",     "32", "32",  1, 0, 0, 1},
{sw_scan,                             "sw_scan",                             "sw",    "scan", "NA",     "32", "32",  1, 0, 0, 0},
{sw_scan_knc_512_32,                  "sw_scan_knc_512_32",                  "sw",    "scan", "knc",   "512", "32", 16, 0, 0, 0},
{sw_striped_knc_512_32,               "sw_striped_knc_512_32",               "sw", "striped", "knc",   "512", "32", 16, 0, 0, 0},
{sw_diag_knc_512_32,                  "sw_diag_knc_512_32",                  "sw",    "diag", "knc",   "512", "32", 16, 0, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_knc = {"sw_knc", sw_knc_functions};
#endif
#if HAVE_SSE2
func_t nw_stats_sse2_functions[] = {
{nw_stats,                            "nw_stats",                            "nw_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{nw_stats_scan,                       "nw_stats_scan",                       "nw_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
{nw_stats_scan_sse2_128_32,           "nw_stats_scan_sse2_128_32",           "nw_stats",    "scan", "sse2",  "128", "32",  4, 0, 1, 0},
{nw_stats_scan_sse2_128_16,           "nw_stats_scan_sse2_128_16",           "nw_stats",    "scan", "sse2",  "128", "16",  8, 0, 1, 0},
{nw_stats_scan_sse2_128_8,            "nw_stats_scan_sse2_128_8",            "nw_stats",    "scan", "sse2",  "128",  "8", 16, 0, 1, 0},
{nw_stats_striped_sse2_128_32,        "nw_stats_striped_sse2_128_32",        "nw_stats", "striped", "sse2",  "128", "32",  4, 0, 1, 0},
{nw_stats_striped_sse2_128_16,        "nw_stats_striped_sse2_128_16",        "nw_stats", "striped", "sse2",  "128", "16",  8, 0, 1, 0},
{nw_stats_striped_sse2_128_8,         "nw_stats_striped_sse2_128_8",         "nw_stats", "striped", "sse2",  "128",  "8", 16, 0, 1, 0},
{nw_stats_diag_sse2_128_32,           "nw_stats_diag_sse2_128_32",           "nw_stats",    "diag", "sse2",  "128", "32",  4, 0, 1, 0},
{nw_stats_diag_sse2_128_16,           "nw_stats_diag_sse2_128_16",           "nw_stats",    "diag", "sse2",  "128", "16",  8, 0, 1, 0},
{nw_stats_diag_sse2_128_8,            "nw_stats_diag_sse2_128_8",            "nw_stats",    "diag", "sse2",  "128",  "8", 16, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_stats_sse2 = {"nw_stats_sse2", nw_stats_sse2_functions};
#endif
#if HAVE_SSE41
func_t nw_stats_sse41_functions[] = {
{nw_stats,                            "nw_stats",                            "nw_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{nw_stats_scan,                       "nw_stats_scan",                       "nw_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
{nw_stats_scan_sse41_128_32,          "nw_stats_scan_sse41_128_32",          "nw_stats",    "scan", "sse41", "128", "32",  4, 0, 1, 0},
{nw_stats_scan_sse41_128_16,          "nw_stats_scan_sse41_128_16",          "nw_stats",    "scan", "sse41", "128", "16",  8, 0, 1, 0},
{nw_stats_scan_sse41_128_8,           "nw_stats_scan_sse41_128_8",           "nw_stats",    "scan", "sse41", "128",  "8", 16, 0, 1, 0},
{nw_stats_striped_sse41_128_32,       "nw_stats_striped_sse41_128_32",       "nw_stats", "striped", "sse41", "128", "32",  4, 0, 1, 0},
{nw_stats_striped_sse41_128_16,       "nw_stats_striped_sse41_128_16",       "nw_stats", "striped", "sse41", "128", "16",  8, 0, 1, 0},
{nw_stats_striped_sse41_128_8,        "nw_stats_striped_sse41_128_8",        "nw_stats", "striped", "sse41", "128",  "8", 16, 0, 1, 0},
{nw_stats_diag_sse41_128_32,          "nw_stats_diag_sse41_128_32",          "nw_stats",    "diag", "sse41", "128", "32",  4, 0, 1, 0},
{nw_stats_diag_sse41_128_16,          "nw_stats_diag_sse41_128_16",          "nw_stats",    "diag", "sse41", "128", "16",  8, 0, 1, 0},
{nw_stats_diag_sse41_128_8,           "nw_stats_diag_sse41_128_8",           "nw_stats",    "diag", "sse41", "128",  "8", 16, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_stats_sse41 = {"nw_stats_sse41", nw_stats_sse41_functions};
#endif
#if HAVE_AVX2
func_t nw_stats_avx2_functions[] = {
{nw_stats,                            "nw_stats",                            "nw_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{nw_stats_scan,                       "nw_stats_scan",                       "nw_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
{nw_stats_scan_avx2_256_32,           "nw_stats_scan_avx2_256_32",           "nw_stats",    "scan", "avx2",  "256", "32",  8, 0, 1, 0},
{nw_stats_scan_avx2_256_16,           "nw_stats_scan_avx2_256_16",           "nw_stats",    "scan", "avx2",  "256", "16", 16, 0, 1, 0},
{nw_stats_scan_avx2_256_8,            "nw_stats_scan_avx2_256_8",            "nw_stats",    "scan", "avx2",  "256",  "8", 32, 0, 1, 0},
{nw_stats_striped_avx2_256_32,        "nw_stats_striped_avx2_256_32",        "nw_stats", "striped", "avx2",  "256", "32",  8, 0, 1, 0},
{nw_stats_striped_avx2_256_16,        "nw_stats_striped_avx2_256_16",        "nw_stats", "striped", "avx2",  "256", "16", 16, 0, 1, 0},
{nw_stats_striped_avx2_256_8,         "nw_stats_striped_avx2_256_8",         "nw_stats", "striped", "avx2",  "256",  "8", 32, 0, 1, 0},
{nw_stats_diag_avx2_256_32,           "nw_stats_diag_avx2_256_32",           "nw_stats",    "diag", "avx2",  "256", "32",  8, 0, 1, 0},
{nw_stats_diag_avx2_256_16,           "nw_stats_diag_avx2_256_16",           "nw_stats",    "diag", "avx2",  "256", "16", 16, 0, 1, 0},
{nw_stats_diag_avx2_256_8,            "nw_stats_diag_avx2_256_8",            "nw_stats",    "diag", "avx2",  "256",  "8", 32, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_stats_avx2 = {"nw_stats_avx2", nw_stats_avx2_functions};
#endif
#if HAVE_KNC
func_t nw_stats_knc_functions[] = {
{nw_stats,                            "nw_stats",                            "nw_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{nw_stats_scan,                       "nw_stats_scan",                       "nw_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
{nw_stats_scan_knc_512_32,            "nw_stats_scan_knc_512_32",            "nw_stats",    "scan", "knc",   "512", "32", 16, 0, 1, 0},
{nw_stats_striped_knc_512_32,         "nw_stats_striped_knc_512_32",         "nw_stats", "striped", "knc",   "512", "32", 16, 0, 1, 0},
{nw_stats_diag_knc_512_32,            "nw_stats_diag_knc_512_32",            "nw_stats",    "diag", "knc",   "512", "32", 16, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_stats_knc = {"nw_stats_knc", nw_stats_knc_functions};
#endif
#if HAVE_SSE2
func_t sg_stats_sse2_functions[] = {
{sg_stats,                            "sg_stats",                            "sg_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{sg_stats_scan,                       "sg_stats_scan",                       "sg_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
{sg_stats_scan_sse2_128_32,           "sg_stats_scan_sse2_128_32",           "sg_stats",    "scan", "sse2",  "128", "32",  4, 0, 1, 0},
{sg_stats_scan_sse2_128_16,           "sg_stats_scan_sse2_128_16",           "sg_stats",    "scan", "sse2",  "128", "16",  8, 0, 1, 0},
{sg_stats_scan_sse2_128_8,            "sg_stats_scan_sse2_128_8",            "sg_stats",    "scan", "sse2",  "128",  "8", 16, 0, 1, 0},
{sg_stats_striped_sse2_128_32,        "sg_stats_striped_sse2_128_32",        "sg_stats", "striped", "sse2",  "128", "32",  4, 0, 1, 0},
{sg_stats_striped_sse2_128_16,        "sg_stats_striped_sse2_128_16",        "sg_stats", "striped", "sse2",  "128", "16",  8, 0, 1, 0},
{sg_stats_striped_sse2_128_8,         "sg_stats_striped_sse2_128_8",         "sg_stats", "striped", "sse2",  "128",  "8", 16, 0, 1, 0},
{sg_stats_diag_sse2_128_32,           "sg_stats_diag_sse2_128_32",           "sg_stats",    "diag", "sse2",  "128", "32",  4, 0, 1, 0},
{sg_stats_diag_sse2_128_16,           "sg_stats_diag_sse2_128_16",           "sg_stats",    "diag", "sse2",  "128", "16",  8, 0, 1, 0},
{sg_stats_diag_sse2_128_8,            "sg_stats_diag_sse2_128_8",            "sg_stats",    "diag", "sse2",  "128",  "8", 16, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_stats_sse2 = {"sg_stats_sse2", sg_stats_sse2_functions};
#endif
#if HAVE_SSE41
func_t sg_stats_sse41_functions[] = {
{sg_stats,                            "sg_stats",                            "sg_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{sg_stats_scan,                       "sg_stats_scan",                       "sg_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
{sg_stats_scan_sse41_128_32,          "sg_stats_scan_sse41_128_32",          "sg_stats",    "scan", "sse41", "128", "32",  4, 0, 1, 0},
{sg_stats_scan_sse41_128_16,          "sg_stats_scan_sse41_128_16",          "sg_stats",    "scan", "sse41", "128", "16",  8, 0, 1, 0},
{sg_stats_scan_sse41_128_8,           "sg_stats_scan_sse41_128_8",           "sg_stats",    "scan", "sse41", "128",  "8", 16, 0, 1, 0},
{sg_stats_striped_sse41_128_32,       "sg_stats_striped_sse41_128_32",       "sg_stats", "striped", "sse41", "128", "32",  4, 0, 1, 0},
{sg_stats_striped_sse41_128_16,       "sg_stats_striped_sse41_128_16",       "sg_stats", "striped", "sse41", "128", "16",  8, 0, 1, 0},
{sg_stats_striped_sse41_128_8,        "sg_stats_striped_sse41_128_8",        "sg_stats", "striped", "sse41", "128",  "8", 16, 0, 1, 0},
{sg_stats_diag_sse41_128_32,          "sg_stats_diag_sse41_128_32",          "sg_stats",    "diag", "sse41", "128", "32",  4, 0, 1, 0},
{sg_stats_diag_sse41_128_16,          "sg_stats_diag_sse41_128_16",          "sg_stats",    "diag", "sse41", "128", "16",  8, 0, 1, 0},
{sg_stats_diag_sse41_128_8,           "sg_stats_diag_sse41_128_8",           "sg_stats",    "diag", "sse41", "128",  "8", 16, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_stats_sse41 = {"sg_stats_sse41", sg_stats_sse41_functions};
#endif
#if HAVE_AVX2
func_t sg_stats_avx2_functions[] = {
{sg_stats,                            "sg_stats",                            "sg_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{sg_stats_scan,                       "sg_stats_scan",                       "sg_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
{sg_stats_scan_avx2_256_32,           "sg_stats_scan_avx2_256_32",           "sg_stats",    "scan", "avx2",  "256", "32",  8, 0, 1, 0},
{sg_stats_scan_avx2_256_16,           "sg_stats_scan_avx2_256_16",           "sg_stats",    "scan", "avx2",  "256", "16", 16, 0, 1, 0},
{sg_stats_scan_avx2_256_8,            "sg_stats_scan_avx2_256_8",            "sg_stats",    "scan", "avx2",  "256",  "8", 32, 0, 1, 0},
{sg_stats_striped_avx2_256_32,        "sg_stats_striped_avx2_256_32",        "sg_stats", "striped", "avx2",  "256", "32",  8, 0, 1, 0},
{sg_stats_striped_avx2_256_16,        "sg_stats_striped_avx2_256_16",        "sg_stats", "striped", "avx2",  "256", "16", 16, 0, 1, 0},
{sg_stats_striped_avx2_256_8,         "sg_stats_striped_avx2_256_8",         "sg_stats", "striped", "avx2",  "256",  "8", 32, 0, 1, 0},
{sg_stats_diag_avx2_256_32,           "sg_stats_diag_avx2_256_32",           "sg_stats",    "diag", "avx2",  "256", "32",  8, 0, 1, 0},
{sg_stats_diag_avx2_256_16,           "sg_stats_diag_avx2_256_16",           "sg_stats",    "diag", "avx2",  "256", "16", 16, 0, 1, 0},
{sg_stats_diag_avx2_256_8,            "sg_stats_diag_avx2_256_8",            "sg_stats",    "diag", "avx2",  "256",  "8", 32, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_stats_avx2 = {"sg_stats_avx2", sg_stats_avx2_functions};
#endif
#if HAVE_KNC
func_t sg_stats_knc_functions[] = {
{sg_stats,                            "sg_stats",                            "sg_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{sg_stats_scan,                       "sg_stats_scan",                       "sg_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
{sg_stats_scan_knc_512_32,            "sg_stats_scan_knc_512_32",            "sg_stats",    "scan", "knc",   "512", "32", 16, 0, 1, 0},
{sg_stats_striped_knc_512_32,         "sg_stats_striped_knc_512_32",         "sg_stats", "striped", "knc",   "512", "32", 16, 0, 1, 0},
{sg_stats_diag_knc_512_32,            "sg_stats_diag_knc_512_32",            "sg_stats",    "diag", "knc",   "512", "32", 16, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_stats_knc = {"sg_stats_knc", sg_stats_knc_functions};
#endif
#if HAVE_SSE2
func_t sw_stats_sse2_functions[] = {
{sw_stats,                            "sw_stats",                            "sw_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{sw_stats_scan,                       "sw_stats_scan",                       "sw_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
{sw_stats_scan_sse2_128_32,           "sw_stats_scan_sse2_128_32",           "sw_stats",    "scan", "sse2",  "128", "32",  4, 0, 1, 0},
{sw_stats_scan_sse2_128_16,           "sw_stats_scan_sse2_128_16",           "sw_stats",    "scan", "sse2",  "128", "16",  8, 0, 1, 0},
{sw_stats_scan_sse2_128_8,            "sw_stats_scan_sse2_128_8",            "sw_stats",    "scan", "sse2",  "128",  "8", 16, 0, 1, 0},
{sw_stats_striped_sse2_128_32,        "sw_stats_striped_sse2_128_32",        "sw_stats", "striped", "sse2",  "128", "32",  4, 0, 1, 0},
{sw_stats_striped_sse2_128_16,        "sw_stats_striped_sse2_128_16",        "sw_stats", "striped", "sse2",  "128", "16",  8, 0, 1, 0},
{sw_stats_striped_sse2_128_8,         "sw_stats_striped_sse2_128_8",         "sw_stats", "striped", "sse2",  "128",  "8", 16, 0, 1, 0},
{sw_stats_diag_sse2_128_32,           "sw_stats_diag_sse2_128_32",           "sw_stats",    "diag", "sse2",  "128", "32",  4, 0, 1, 0},
{sw_stats_diag_sse2_128_16,           "sw_stats_diag_sse2_128_16",           "sw_stats",    "diag", "sse2",  "128", "16",  8, 0, 1, 0},
{sw_stats_diag_sse2_128_8,            "sw_stats_diag_sse2_128_8",            "sw_stats",    "diag", "sse2",  "128",  "8", 16, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_stats_sse2 = {"sw_stats_sse2", sw_stats_sse2_functions};
#endif
#if HAVE_SSE41
func_t sw_stats_sse41_functions[] = {
{sw_stats,                            "sw_stats",                            "sw_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{sw_stats_scan,                       "sw_stats_scan",                       "sw_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
{sw_stats_scan_sse41_128_32,          "sw_stats_scan_sse41_128_32",          "sw_stats",    "scan", "sse41", "128", "32",  4, 0, 1, 0},
{sw_stats_scan_sse41_128_16,          "sw_stats_scan_sse41_128_16",          "sw_stats",    "scan", "sse41", "128", "16",  8, 0, 1, 0},
{sw_stats_scan_sse41_128_8,           "sw_stats_scan_sse41_128_8",           "sw_stats",    "scan", "sse41", "128",  "8", 16, 0, 1, 0},
{sw_stats_striped_sse41_128_32,       "sw_stats_striped_sse41_128_32",       "sw_stats", "striped", "sse41", "128", "32",  4, 0, 1, 0},
{sw_stats_striped_sse41_128_16,       "sw_stats_striped_sse41_128_16",       "sw_stats", "striped", "sse41", "128", "16",  8, 0, 1, 0},
{sw_stats_striped_sse41_128_8,        "sw_stats_striped_sse41_128_8",        "sw_stats", "striped", "sse41", "128",  "8", 16, 0, 1, 0},
{sw_stats_diag_sse41_128_32,          "sw_stats_diag_sse41_128_32",          "sw_stats",    "diag", "sse41", "128", "32",  4, 0, 1, 0},
{sw_stats_diag_sse41_128_16,          "sw_stats_diag_sse41_128_16",          "sw_stats",    "diag", "sse41", "128", "16",  8, 0, 1, 0},
{sw_stats_diag_sse41_128_8,           "sw_stats_diag_sse41_128_8",           "sw_stats",    "diag", "sse41", "128",  "8", 16, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_stats_sse41 = {"sw_stats_sse41", sw_stats_sse41_functions};
#endif
#if HAVE_AVX2
func_t sw_stats_avx2_functions[] = {
{sw_stats,                            "sw_stats",                            "sw_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{sw_stats_scan,                       "sw_stats_scan",                       "sw_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
{sw_stats_scan_avx2_256_32,           "sw_stats_scan_avx2_256_32",           "sw_stats",    "scan", "avx2",  "256", "32",  8, 0, 1, 0},
{sw_stats_scan_avx2_256_16,           "sw_stats_scan_avx2_256_16",           "sw_stats",    "scan", "avx2",  "256", "16", 16, 0, 1, 0},
{sw_stats_scan_avx2_256_8,            "sw_stats_scan_avx2_256_8",            "sw_stats",    "scan", "avx2",  "256",  "8", 32, 0, 1, 0},
{sw_stats_striped_avx2_256_32,        "sw_stats_striped_avx2_256_32",        "sw_stats", "striped", "avx2",  "256", "32",  8, 0, 1, 0},
{sw_stats_striped_avx2_256_16,        "sw_stats_striped_avx2_256_16",        "sw_stats", "striped", "avx2",  "256", "16", 16, 0, 1, 0},
{sw_stats_striped_avx2_256_8,         "sw_stats_striped_avx2_256_8",         "sw_stats", "striped", "avx2",  "256",  "8", 32, 0, 1, 0},
{sw_stats_diag_avx2_256_32,           "sw_stats_diag_avx2_256_32",           "sw_stats",    "diag", "avx2",  "256", "32",  8, 0, 1, 0},
{sw_stats_diag_avx2_256_16,           "sw_stats_diag_avx2_256_16",           "sw_stats",    "diag", "avx2",  "256", "16", 16, 0, 1, 0},
{sw_stats_diag_avx2_256_8,            "sw_stats_diag_avx2_256_8",            "sw_stats",    "diag", "avx2",  "256",  "8", 32, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_stats_avx2 = {"sw_stats_avx2", sw_stats_avx2_functions};
#endif
#if HAVE_KNC
func_t sw_stats_knc_functions[] = {
{sw_stats,                            "sw_stats",                            "sw_stats",    "orig", "NA",     "32", "32",  1, 0, 1, 1},
{sw_stats_scan,                       "sw_stats_scan",                       "sw_stats",    "scan", "NA",     "32", "32",  1, 0, 1, 0},
{sw_stats_scan_knc_512_32,            "sw_stats_scan_knc_512_32",            "sw_stats",    "scan", "knc",   "512", "32", 16, 0, 1, 0},
{sw_stats_striped_knc_512_32,         "sw_stats_striped_knc_512_32",         "sw_stats", "striped", "knc",   "512", "32", 16, 0, 1, 0},
{sw_stats_diag_knc_512_32,            "sw_stats_diag_knc_512_32",            "sw_stats",    "diag", "knc",   "512", "32", 16, 0, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_stats_knc = {"sw_stats_knc", sw_stats_knc_functions};
#endif
#if HAVE_SSE2
func_t nw_table_sse2_functions[] = {
{nw_table,                            "nw_table",                            "nw",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{nw_table_scan,                       "nw_table_scan",                       "nw",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
{nw_table_scan_sse2_128_64,           "nw_table_scan_sse2_128_64",           "nw",    "scan", "sse2",  "128", "64",  2, 1, 0, 0},
{nw_table_scan_sse2_128_32,           "nw_table_scan_sse2_128_32",           "nw",    "scan", "sse2",  "128", "32",  4, 1, 0, 0},
{nw_table_scan_sse2_128_16,           "nw_table_scan_sse2_128_16",           "nw",    "scan", "sse2",  "128", "16",  8, 1, 0, 0},
{nw_table_scan_sse2_128_8,            "nw_table_scan_sse2_128_8",            "nw",    "scan", "sse2",  "128",  "8", 16, 1, 0, 0},
{nw_table_striped_sse2_128_64,        "nw_table_striped_sse2_128_64",        "nw", "striped", "sse2",  "128", "64",  2, 1, 0, 0},
{nw_table_striped_sse2_128_32,        "nw_table_striped_sse2_128_32",        "nw", "striped", "sse2",  "128", "32",  4, 1, 0, 0},
{nw_table_striped_sse2_128_16,        "nw_table_striped_sse2_128_16",        "nw", "striped", "sse2",  "128", "16",  8, 1, 0, 0},
{nw_table_striped_sse2_128_8,         "nw_table_striped_sse2_128_8",         "nw", "striped", "sse2",  "128",  "8", 16, 1, 0, 0},
{nw_table_diag_sse2_128_64,           "nw_table_diag_sse2_128_64",           "nw",    "diag", "sse2",  "128", "64",  2, 1, 0, 0},
{nw_table_diag_sse2_128_32,           "nw_table_diag_sse2_128_32",           "nw",    "diag", "sse2",  "128", "32",  4, 1, 0, 0},
{nw_table_diag_sse2_128_16,           "nw_table_diag_sse2_128_16",           "nw",    "diag", "sse2",  "128", "16",  8, 1, 0, 0},
{nw_table_diag_sse2_128_8,            "nw_table_diag_sse2_128_8",            "nw",    "diag", "sse2",  "128",  "8", 16, 1, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_table_sse2 = {"nw_table_sse2", nw_table_sse2_functions};
#endif
#if HAVE_SSE41
func_t nw_table_sse41_functions[] = {
{nw_table,                            "nw_table",                            "nw",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{nw_table_scan,                       "nw_table_scan",                       "nw",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
{nw_table_scan_sse41_128_64,          "nw_table_scan_sse41_128_64",          "nw",    "scan", "sse41", "128", "64",  2, 1, 0, 0},
{nw_table_scan_sse41_128_32,          "nw_table_scan_sse41_128_32",          "nw",    "scan", "sse41", "128", "32",  4, 1, 0, 0},
{nw_table_scan_sse41_128_16,          "nw_table_scan_sse41_128_16",          "nw",    "scan", "sse41", "128", "16",  8, 1, 0, 0},
{nw_table_scan_sse41_128_8,           "nw_table_scan_sse41_128_8",           "nw",    "scan", "sse41", "128",  "8", 16, 1, 0, 0},
{nw_table_striped_sse41_128_64,       "nw_table_striped_sse41_128_64",       "nw", "striped", "sse41", "128", "64",  2, 1, 0, 0},
{nw_table_striped_sse41_128_32,       "nw_table_striped_sse41_128_32",       "nw", "striped", "sse41", "128", "32",  4, 1, 0, 0},
{nw_table_striped_sse41_128_16,       "nw_table_striped_sse41_128_16",       "nw", "striped", "sse41", "128", "16",  8, 1, 0, 0},
{nw_table_striped_sse41_128_8,        "nw_table_striped_sse41_128_8",        "nw", "striped", "sse41", "128",  "8", 16, 1, 0, 0},
{nw_table_diag_sse41_128_64,          "nw_table_diag_sse41_128_64",          "nw",    "diag", "sse41", "128", "64",  2, 1, 0, 0},
{nw_table_diag_sse41_128_32,          "nw_table_diag_sse41_128_32",          "nw",    "diag", "sse41", "128", "32",  4, 1, 0, 0},
{nw_table_diag_sse41_128_16,          "nw_table_diag_sse41_128_16",          "nw",    "diag", "sse41", "128", "16",  8, 1, 0, 0},
{nw_table_diag_sse41_128_8,           "nw_table_diag_sse41_128_8",           "nw",    "diag", "sse41", "128",  "8", 16, 1, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_table_sse41 = {"nw_table_sse41", nw_table_sse41_functions};
#endif
#if HAVE_AVX2
func_t nw_table_avx2_functions[] = {
{nw_table,                            "nw_table",                            "nw",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{nw_table_scan,                       "nw_table_scan",                       "nw",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
{nw_table_scan_avx2_256_64,           "nw_table_scan_avx2_256_64",           "nw",    "scan", "avx2",  "256", "64",  4, 1, 0, 0},
{nw_table_scan_avx2_256_32,           "nw_table_scan_avx2_256_32",           "nw",    "scan", "avx2",  "256", "32",  8, 1, 0, 0},
{nw_table_scan_avx2_256_16,           "nw_table_scan_avx2_256_16",           "nw",    "scan", "avx2",  "256", "16", 16, 1, 0, 0},
{nw_table_scan_avx2_256_8,            "nw_table_scan_avx2_256_8",            "nw",    "scan", "avx2",  "256",  "8", 32, 1, 0, 0},
{nw_table_striped_avx2_256_64,        "nw_table_striped_avx2_256_64",        "nw", "striped", "avx2",  "256", "64",  4, 1, 0, 0},
{nw_table_striped_avx2_256_32,        "nw_table_striped_avx2_256_32",        "nw", "striped", "avx2",  "256", "32",  8, 1, 0, 0},
{nw_table_striped_avx2_256_16,        "nw_table_striped_avx2_256_16",        "nw", "striped", "avx2",  "256", "16", 16, 1, 0, 0},
{nw_table_striped_avx2_256_8,         "nw_table_striped_avx2_256_8",         "nw", "striped", "avx2",  "256",  "8", 32, 1, 0, 0},
{nw_table_diag_avx2_256_64,           "nw_table_diag_avx2_256_64",           "nw",    "diag", "avx2",  "256", "64",  4, 1, 0, 0},
{nw_table_diag_avx2_256_32,           "nw_table_diag_avx2_256_32",           "nw",    "diag", "avx2",  "256", "32",  8, 1, 0, 0},
{nw_table_diag_avx2_256_16,           "nw_table_diag_avx2_256_16",           "nw",    "diag", "avx2",  "256", "16", 16, 1, 0, 0},
{nw_table_diag_avx2_256_8,            "nw_table_diag_avx2_256_8",            "nw",    "diag", "avx2",  "256",  "8", 32, 1, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_table_avx2 = {"nw_table_avx2", nw_table_avx2_functions};
#endif
#if HAVE_KNC
func_t nw_table_knc_functions[] = {
{nw_table,                            "nw_table",                            "nw",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{nw_table_scan,                       "nw_table_scan",                       "nw",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
{nw_table_scan_knc_512_32,            "nw_table_scan_knc_512_32",            "nw",    "scan", "knc",   "512", "32", 16, 1, 0, 0},
{nw_table_striped_knc_512_32,         "nw_table_striped_knc_512_32",         "nw", "striped", "knc",   "512", "32", 16, 1, 0, 0},
{nw_table_diag_knc_512_32,            "nw_table_diag_knc_512_32",            "nw",    "diag", "knc",   "512", "32", 16, 1, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_table_knc = {"nw_table_knc", nw_table_knc_functions};
#endif
#if HAVE_SSE2
func_t sg_table_sse2_functions[] = {
{sg_table,                            "sg_table",                            "sg",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{sg_table_scan,                       "sg_table_scan",                       "sg",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
{sg_table_scan_sse2_128_64,           "sg_table_scan_sse2_128_64",           "sg",    "scan", "sse2",  "128", "64",  2, 1, 0, 0},
{sg_table_scan_sse2_128_32,           "sg_table_scan_sse2_128_32",           "sg",    "scan", "sse2",  "128", "32",  4, 1, 0, 0},
{sg_table_scan_sse2_128_16,           "sg_table_scan_sse2_128_16",           "sg",    "scan", "sse2",  "128", "16",  8, 1, 0, 0},
{sg_table_scan_sse2_128_8,            "sg_table_scan_sse2_128_8",            "sg",    "scan", "sse2",  "128",  "8", 16, 1, 0, 0},
{sg_table_striped_sse2_128_64,        "sg_table_striped_sse2_128_64",        "sg", "striped", "sse2",  "128", "64",  2, 1, 0, 0},
{sg_table_striped_sse2_128_32,        "sg_table_striped_sse2_128_32",        "sg", "striped", "sse2",  "128", "32",  4, 1, 0, 0},
{sg_table_striped_sse2_128_16,        "sg_table_striped_sse2_128_16",        "sg", "striped", "sse2",  "128", "16",  8, 1, 0, 0},
{sg_table_striped_sse2_128_8,         "sg_table_striped_sse2_128_8",         "sg", "striped", "sse2",  "128",  "8", 16, 1, 0, 0},
{sg_table_diag_sse2_128_64,           "sg_table_diag_sse2_128_64",           "sg",    "diag", "sse2",  "128", "64",  2, 1, 0, 0},
{sg_table_diag_sse2_128_32,           "sg_table_diag_sse2_128_32",           "sg",    "diag", "sse2",  "128", "32",  4, 1, 0, 0},
{sg_table_diag_sse2_128_16,           "sg_table_diag_sse2_128_16",           "sg",    "diag", "sse2",  "128", "16",  8, 1, 0, 0},
{sg_table_diag_sse2_128_8,            "sg_table_diag_sse2_128_8",            "sg",    "diag", "sse2",  "128",  "8", 16, 1, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_table_sse2 = {"sg_table_sse2", sg_table_sse2_functions};
#endif
#if HAVE_SSE41
func_t sg_table_sse41_functions[] = {
{sg_table,                            "sg_table",                            "sg",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{sg_table_scan,                       "sg_table_scan",                       "sg",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
{sg_table_scan_sse41_128_64,          "sg_table_scan_sse41_128_64",          "sg",    "scan", "sse41", "128", "64",  2, 1, 0, 0},
{sg_table_scan_sse41_128_32,          "sg_table_scan_sse41_128_32",          "sg",    "scan", "sse41", "128", "32",  4, 1, 0, 0},
{sg_table_scan_sse41_128_16,          "sg_table_scan_sse41_128_16",          "sg",    "scan", "sse41", "128", "16",  8, 1, 0, 0},
{sg_table_scan_sse41_128_8,           "sg_table_scan_sse41_128_8",           "sg",    "scan", "sse41", "128",  "8", 16, 1, 0, 0},
{sg_table_striped_sse41_128_64,       "sg_table_striped_sse41_128_64",       "sg", "striped", "sse41", "128", "64",  2, 1, 0, 0},
{sg_table_striped_sse41_128_32,       "sg_table_striped_sse41_128_32",       "sg", "striped", "sse41", "128", "32",  4, 1, 0, 0},
{sg_table_striped_sse41_128_16,       "sg_table_striped_sse41_128_16",       "sg", "striped", "sse41", "128", "16",  8, 1, 0, 0},
{sg_table_striped_sse41_128_8,        "sg_table_striped_sse41_128_8",        "sg", "striped", "sse41", "128",  "8", 16, 1, 0, 0},
{sg_table_diag_sse41_128_64,          "sg_table_diag_sse41_128_64",          "sg",    "diag", "sse41", "128", "64",  2, 1, 0, 0},
{sg_table_diag_sse41_128_32,          "sg_table_diag_sse41_128_32",          "sg",    "diag", "sse41", "128", "32",  4, 1, 0, 0},
{sg_table_diag_sse41_128_16,          "sg_table_diag_sse41_128_16",          "sg",    "diag", "sse41", "128", "16",  8, 1, 0, 0},
{sg_table_diag_sse41_128_8,           "sg_table_diag_sse41_128_8",           "sg",    "diag", "sse41", "128",  "8", 16, 1, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_table_sse41 = {"sg_table_sse41", sg_table_sse41_functions};
#endif
#if HAVE_AVX2
func_t sg_table_avx2_functions[] = {
{sg_table,                            "sg_table",                            "sg",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{sg_table_scan,                       "sg_table_scan",                       "sg",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
{sg_table_scan_avx2_256_64,           "sg_table_scan_avx2_256_64",           "sg",    "scan", "avx2",  "256", "64",  4, 1, 0, 0},
{sg_table_scan_avx2_256_32,           "sg_table_scan_avx2_256_32",           "sg",    "scan", "avx2",  "256", "32",  8, 1, 0, 0},
{sg_table_scan_avx2_256_16,           "sg_table_scan_avx2_256_16",           "sg",    "scan", "avx2",  "256", "16", 16, 1, 0, 0},
{sg_table_scan_avx2_256_8,            "sg_table_scan_avx2_256_8",            "sg",    "scan", "avx2",  "256",  "8", 32, 1, 0, 0},
{sg_table_striped_avx2_256_64,        "sg_table_striped_avx2_256_64",        "sg", "striped", "avx2",  "256", "64",  4, 1, 0, 0},
{sg_table_striped_avx2_256_32,        "sg_table_striped_avx2_256_32",        "sg", "striped", "avx2",  "256", "32",  8, 1, 0, 0},
{sg_table_striped_avx2_256_16,        "sg_table_striped_avx2_256_16",        "sg", "striped", "avx2",  "256", "16", 16, 1, 0, 0},
{sg_table_striped_avx2_256_8,         "sg_table_striped_avx2_256_8",         "sg", "striped", "avx2",  "256",  "8", 32, 1, 0, 0},
{sg_table_diag_avx2_256_64,           "sg_table_diag_avx2_256_64",           "sg",    "diag", "avx2",  "256", "64",  4, 1, 0, 0},
{sg_table_diag_avx2_256_32,           "sg_table_diag_avx2_256_32",           "sg",    "diag", "avx2",  "256", "32",  8, 1, 0, 0},
{sg_table_diag_avx2_256_16,           "sg_table_diag_avx2_256_16",           "sg",    "diag", "avx2",  "256", "16", 16, 1, 0, 0},
{sg_table_diag_avx2_256_8,            "sg_table_diag_avx2_256_8",            "sg",    "diag", "avx2",  "256",  "8", 32, 1, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_table_avx2 = {"sg_table_avx2", sg_table_avx2_functions};
#endif
#if HAVE_KNC
func_t sg_table_knc_functions[] = {
{sg_table,                            "sg_table",                            "sg",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{sg_table_scan,                       "sg_table_scan",                       "sg",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
{sg_table_scan_knc_512_32,            "sg_table_scan_knc_512_32",            "sg",    "scan", "knc",   "512", "32", 16, 1, 0, 0},
{sg_table_striped_knc_512_32,         "sg_table_striped_knc_512_32",         "sg", "striped", "knc",   "512", "32", 16, 1, 0, 0},
{sg_table_diag_knc_512_32,            "sg_table_diag_knc_512_32",            "sg",    "diag", "knc",   "512", "32", 16, 1, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_table_knc = {"sg_table_knc", sg_table_knc_functions};
#endif
#if HAVE_SSE2
func_t sw_table_sse2_functions[] = {
{sw_table,                            "sw_table",                            "sw",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{sw_table_scan,                       "sw_table_scan",                       "sw",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
{sw_table_scan_sse2_128_64,           "sw_table_scan_sse2_128_64",           "sw",    "scan", "sse2",  "128", "64",  2, 1, 0, 0},
{sw_table_scan_sse2_128_32,           "sw_table_scan_sse2_128_32",           "sw",    "scan", "sse2",  "128", "32",  4, 1, 0, 0},
{sw_table_scan_sse2_128_16,           "sw_table_scan_sse2_128_16",           "sw",    "scan", "sse2",  "128", "16",  8, 1, 0, 0},
{sw_table_scan_sse2_128_8,            "sw_table_scan_sse2_128_8",            "sw",    "scan", "sse2",  "128",  "8", 16, 1, 0, 0},
{sw_table_striped_sse2_128_64,        "sw_table_striped_sse2_128_64",        "sw", "striped", "sse2",  "128", "64",  2, 1, 0, 0},
{sw_table_striped_sse2_128_32,        "sw_table_striped_sse2_128_32",        "sw", "striped", "sse2",  "128", "32",  4, 1, 0, 0},
{sw_table_striped_sse2_128_16,        "sw_table_striped_sse2_128_16",        "sw", "striped", "sse2",  "128", "16",  8, 1, 0, 0},
{sw_table_striped_sse2_128_8,         "sw_table_striped_sse2_128_8",         "sw", "striped", "sse2",  "128",  "8", 16, 1, 0, 0},
{sw_table_diag_sse2_128_64,           "sw_table_diag_sse2_128_64",           "sw",    "diag", "sse2",  "128", "64",  2, 1, 0, 0},
{sw_table_diag_sse2_128_32,           "sw_table_diag_sse2_128_32",           "sw",    "diag", "sse2",  "128", "32",  4, 1, 0, 0},
{sw_table_diag_sse2_128_16,           "sw_table_diag_sse2_128_16",           "sw",    "diag", "sse2",  "128", "16",  8, 1, 0, 0},
{sw_table_diag_sse2_128_8,            "sw_table_diag_sse2_128_8",            "sw",    "diag", "sse2",  "128",  "8", 16, 1, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_table_sse2 = {"sw_table_sse2", sw_table_sse2_functions};
#endif
#if HAVE_SSE41
func_t sw_table_sse41_functions[] = {
{sw_table,                            "sw_table",                            "sw",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{sw_table_scan,                       "sw_table_scan",                       "sw",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
{sw_table_scan_sse41_128_64,          "sw_table_scan_sse41_128_64",          "sw",    "scan", "sse41", "128", "64",  2, 1, 0, 0},
{sw_table_scan_sse41_128_32,          "sw_table_scan_sse41_128_32",          "sw",    "scan", "sse41", "128", "32",  4, 1, 0, 0},
{sw_table_scan_sse41_128_16,          "sw_table_scan_sse41_128_16",          "sw",    "scan", "sse41", "128", "16",  8, 1, 0, 0},
{sw_table_scan_sse41_128_8,           "sw_table_scan_sse41_128_8",           "sw",    "scan", "sse41", "128",  "8", 16, 1, 0, 0},
{sw_table_striped_sse41_128_64,       "sw_table_striped_sse41_128_64",       "sw", "striped", "sse41", "128", "64",  2, 1, 0, 0},
{sw_table_striped_sse41_128_32,       "sw_table_striped_sse41_128_32",       "sw", "striped", "sse41", "128", "32",  4, 1, 0, 0},
{sw_table_striped_sse41_128_16,       "sw_table_striped_sse41_128_16",       "sw", "striped", "sse41", "128", "16",  8, 1, 0, 0},
{sw_table_striped_sse41_128_8,        "sw_table_striped_sse41_128_8",        "sw", "striped", "sse41", "128",  "8", 16, 1, 0, 0},
{sw_table_diag_sse41_128_64,          "sw_table_diag_sse41_128_64",          "sw",    "diag", "sse41", "128", "64",  2, 1, 0, 0},
{sw_table_diag_sse41_128_32,          "sw_table_diag_sse41_128_32",          "sw",    "diag", "sse41", "128", "32",  4, 1, 0, 0},
{sw_table_diag_sse41_128_16,          "sw_table_diag_sse41_128_16",          "sw",    "diag", "sse41", "128", "16",  8, 1, 0, 0},
{sw_table_diag_sse41_128_8,           "sw_table_diag_sse41_128_8",           "sw",    "diag", "sse41", "128",  "8", 16, 1, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_table_sse41 = {"sw_table_sse41", sw_table_sse41_functions};
#endif
#if HAVE_AVX2
func_t sw_table_avx2_functions[] = {
{sw_table,                            "sw_table",                            "sw",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{sw_table_scan,                       "sw_table_scan",                       "sw",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
{sw_table_scan_avx2_256_64,           "sw_table_scan_avx2_256_64",           "sw",    "scan", "avx2",  "256", "64",  4, 1, 0, 0},
{sw_table_scan_avx2_256_32,           "sw_table_scan_avx2_256_32",           "sw",    "scan", "avx2",  "256", "32",  8, 1, 0, 0},
{sw_table_scan_avx2_256_16,           "sw_table_scan_avx2_256_16",           "sw",    "scan", "avx2",  "256", "16", 16, 1, 0, 0},
{sw_table_scan_avx2_256_8,            "sw_table_scan_avx2_256_8",            "sw",    "scan", "avx2",  "256",  "8", 32, 1, 0, 0},
{sw_table_striped_avx2_256_64,        "sw_table_striped_avx2_256_64",        "sw", "striped", "avx2",  "256", "64",  4, 1, 0, 0},
{sw_table_striped_avx2_256_32,        "sw_table_striped_avx2_256_32",        "sw", "striped", "avx2",  "256", "32",  8, 1, 0, 0},
{sw_table_striped_avx2_256_16,        "sw_table_striped_avx2_256_16",        "sw", "striped", "avx2",  "256", "16", 16, 1, 0, 0},
{sw_table_striped_avx2_256_8,         "sw_table_striped_avx2_256_8",         "sw", "striped", "avx2",  "256",  "8", 32, 1, 0, 0},
{sw_table_diag_avx2_256_64,           "sw_table_diag_avx2_256_64",           "sw",    "diag", "avx2",  "256", "64",  4, 1, 0, 0},
{sw_table_diag_avx2_256_32,           "sw_table_diag_avx2_256_32",           "sw",    "diag", "avx2",  "256", "32",  8, 1, 0, 0},
{sw_table_diag_avx2_256_16,           "sw_table_diag_avx2_256_16",           "sw",    "diag", "avx2",  "256", "16", 16, 1, 0, 0},
{sw_table_diag_avx2_256_8,            "sw_table_diag_avx2_256_8",            "sw",    "diag", "avx2",  "256",  "8", 32, 1, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_table_avx2 = {"sw_table_avx2", sw_table_avx2_functions};
#endif
#if HAVE_KNC
func_t sw_table_knc_functions[] = {
{sw_table,                            "sw_table",                            "sw",    "orig", "NA",     "32", "32",  1, 1, 0, 1},
{sw_table_scan,                       "sw_table_scan",                       "sw",    "scan", "NA",     "32", "32",  1, 1, 0, 0},
{sw_table_scan_knc_512_32,            "sw_table_scan_knc_512_32",            "sw",    "scan", "knc",   "512", "32", 16, 1, 0, 0},
{sw_table_striped_knc_512_32,         "sw_table_striped_knc_512_32",         "sw", "striped", "knc",   "512", "32", 16, 1, 0, 0},
{sw_table_diag_knc_512_32,            "sw_table_diag_knc_512_32",            "sw",    "diag", "knc",   "512", "32", 16, 1, 0, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_table_knc = {"sw_table_knc", sw_table_knc_functions};
#endif
#if HAVE_SSE2
func_t nw_stats_table_sse2_functions[] = {
{nw_stats_table,                      "nw_stats_table",                      "nw_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{nw_stats_table_scan,                 "nw_stats_table_scan",                 "nw_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
{nw_stats_table_scan_sse2_128_32,     "nw_stats_table_scan_sse2_128_32",     "nw_stats",    "scan", "sse2",  "128", "32",  4, 1, 1, 0},
{nw_stats_table_scan_sse2_128_16,     "nw_stats_table_scan_sse2_128_16",     "nw_stats",    "scan", "sse2",  "128", "16",  8, 1, 1, 0},
{nw_stats_table_scan_sse2_128_8,      "nw_stats_table_scan_sse2_128_8",      "nw_stats",    "scan", "sse2",  "128",  "8", 16, 1, 1, 0},
{nw_stats_table_striped_sse2_128_32,  "nw_stats_table_striped_sse2_128_32",  "nw_stats", "striped", "sse2",  "128", "32",  4, 1, 1, 0},
{nw_stats_table_striped_sse2_128_16,  "nw_stats_table_striped_sse2_128_16",  "nw_stats", "striped", "sse2",  "128", "16",  8, 1, 1, 0},
{nw_stats_table_striped_sse2_128_8,   "nw_stats_table_striped_sse2_128_8",   "nw_stats", "striped", "sse2",  "128",  "8", 16, 1, 1, 0},
{nw_stats_table_diag_sse2_128_32,     "nw_stats_table_diag_sse2_128_32",     "nw_stats",    "diag", "sse2",  "128", "32",  4, 1, 1, 0},
{nw_stats_table_diag_sse2_128_16,     "nw_stats_table_diag_sse2_128_16",     "nw_stats",    "diag", "sse2",  "128", "16",  8, 1, 1, 0},
{nw_stats_table_diag_sse2_128_8,      "nw_stats_table_diag_sse2_128_8",      "nw_stats",    "diag", "sse2",  "128",  "8", 16, 1, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_stats_table_sse2 = {"nw_stats_table_sse2", nw_stats_table_sse2_functions};
#endif
#if HAVE_SSE41
func_t nw_stats_table_sse41_functions[] = {
{nw_stats_table,                      "nw_stats_table",                      "nw_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{nw_stats_table_scan,                 "nw_stats_table_scan",                 "nw_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
{nw_stats_table_scan_sse41_128_32,    "nw_stats_table_scan_sse41_128_32",    "nw_stats",    "scan", "sse41", "128", "32",  4, 1, 1, 0},
{nw_stats_table_scan_sse41_128_16,    "nw_stats_table_scan_sse41_128_16",    "nw_stats",    "scan", "sse41", "128", "16",  8, 1, 1, 0},
{nw_stats_table_scan_sse41_128_8,     "nw_stats_table_scan_sse41_128_8",     "nw_stats",    "scan", "sse41", "128",  "8", 16, 1, 1, 0},
{nw_stats_table_striped_sse41_128_32, "nw_stats_table_striped_sse41_128_32", "nw_stats", "striped", "sse41", "128", "32",  4, 1, 1, 0},
{nw_stats_table_striped_sse41_128_16, "nw_stats_table_striped_sse41_128_16", "nw_stats", "striped", "sse41", "128", "16",  8, 1, 1, 0},
{nw_stats_table_striped_sse41_128_8,  "nw_stats_table_striped_sse41_128_8",  "nw_stats", "striped", "sse41", "128",  "8", 16, 1, 1, 0},
{nw_stats_table_diag_sse41_128_32,    "nw_stats_table_diag_sse41_128_32",    "nw_stats",    "diag", "sse41", "128", "32",  4, 1, 1, 0},
{nw_stats_table_diag_sse41_128_16,    "nw_stats_table_diag_sse41_128_16",    "nw_stats",    "diag", "sse41", "128", "16",  8, 1, 1, 0},
{nw_stats_table_diag_sse41_128_8,     "nw_stats_table_diag_sse41_128_8",     "nw_stats",    "diag", "sse41", "128",  "8", 16, 1, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_stats_table_sse41 = {"nw_stats_table_sse41", nw_stats_table_sse41_functions};
#endif
#if HAVE_AVX2
func_t nw_stats_table_avx2_functions[] = {
{nw_stats_table,                      "nw_stats_table",                      "nw_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{nw_stats_table_scan,                 "nw_stats_table_scan",                 "nw_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
{nw_stats_table_scan_avx2_256_32,     "nw_stats_table_scan_avx2_256_32",     "nw_stats",    "scan", "avx2",  "256", "32",  8, 1, 1, 0},
{nw_stats_table_scan_avx2_256_16,     "nw_stats_table_scan_avx2_256_16",     "nw_stats",    "scan", "avx2",  "256", "16", 16, 1, 1, 0},
{nw_stats_table_scan_avx2_256_8,      "nw_stats_table_scan_avx2_256_8",      "nw_stats",    "scan", "avx2",  "256",  "8", 32, 1, 1, 0},
{nw_stats_table_striped_avx2_256_32,  "nw_stats_table_striped_avx2_256_32",  "nw_stats", "striped", "avx2",  "256", "32",  8, 1, 1, 0},
{nw_stats_table_striped_avx2_256_16,  "nw_stats_table_striped_avx2_256_16",  "nw_stats", "striped", "avx2",  "256", "16", 16, 1, 1, 0},
{nw_stats_table_striped_avx2_256_8,   "nw_stats_table_striped_avx2_256_8",   "nw_stats", "striped", "avx2",  "256",  "8", 32, 1, 1, 0},
{nw_stats_table_diag_avx2_256_32,     "nw_stats_table_diag_avx2_256_32",     "nw_stats",    "diag", "avx2",  "256", "32",  8, 1, 1, 0},
{nw_stats_table_diag_avx2_256_16,     "nw_stats_table_diag_avx2_256_16",     "nw_stats",    "diag", "avx2",  "256", "16", 16, 1, 1, 0},
{nw_stats_table_diag_avx2_256_8,      "nw_stats_table_diag_avx2_256_8",      "nw_stats",    "diag", "avx2",  "256",  "8", 32, 1, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_stats_table_avx2 = {"nw_stats_table_avx2", nw_stats_table_avx2_functions};
#endif
#if HAVE_KNC
func_t nw_stats_table_knc_functions[] = {
{nw_stats_table,                      "nw_stats_table",                      "nw_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{nw_stats_table_scan,                 "nw_stats_table_scan",                 "nw_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
{nw_stats_table_scan_knc_512_32,      "nw_stats_table_scan_knc_512_32",      "nw_stats",    "scan", "knc",   "512", "32", 16, 1, 1, 0},
{nw_stats_table_striped_knc_512_32,   "nw_stats_table_striped_knc_512_32",   "nw_stats", "striped", "knc",   "512", "32", 16, 1, 1, 0},
{nw_stats_table_diag_knc_512_32,      "nw_stats_table_diag_knc_512_32",      "nw_stats",    "diag", "knc",   "512", "32", 16, 1, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t nw_stats_table_knc = {"nw_stats_table_knc", nw_stats_table_knc_functions};
#endif
#if HAVE_SSE2
func_t sg_stats_table_sse2_functions[] = {
{sg_stats_table,                      "sg_stats_table",                      "sg_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{sg_stats_table_scan,                 "sg_stats_table_scan",                 "sg_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
{sg_stats_table_scan_sse2_128_32,     "sg_stats_table_scan_sse2_128_32",     "sg_stats",    "scan", "sse2",  "128", "32",  4, 1, 1, 0},
{sg_stats_table_scan_sse2_128_16,     "sg_stats_table_scan_sse2_128_16",     "sg_stats",    "scan", "sse2",  "128", "16",  8, 1, 1, 0},
{sg_stats_table_scan_sse2_128_8,      "sg_stats_table_scan_sse2_128_8",      "sg_stats",    "scan", "sse2",  "128",  "8", 16, 1, 1, 0},
{sg_stats_table_striped_sse2_128_32,  "sg_stats_table_striped_sse2_128_32",  "sg_stats", "striped", "sse2",  "128", "32",  4, 1, 1, 0},
{sg_stats_table_striped_sse2_128_16,  "sg_stats_table_striped_sse2_128_16",  "sg_stats", "striped", "sse2",  "128", "16",  8, 1, 1, 0},
{sg_stats_table_striped_sse2_128_8,   "sg_stats_table_striped_sse2_128_8",   "sg_stats", "striped", "sse2",  "128",  "8", 16, 1, 1, 0},
{sg_stats_table_diag_sse2_128_32,     "sg_stats_table_diag_sse2_128_32",     "sg_stats",    "diag", "sse2",  "128", "32",  4, 1, 1, 0},
{sg_stats_table_diag_sse2_128_16,     "sg_stats_table_diag_sse2_128_16",     "sg_stats",    "diag", "sse2",  "128", "16",  8, 1, 1, 0},
{sg_stats_table_diag_sse2_128_8,      "sg_stats_table_diag_sse2_128_8",      "sg_stats",    "diag", "sse2",  "128",  "8", 16, 1, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_stats_table_sse2 = {"sg_stats_table_sse2", sg_stats_table_sse2_functions};
#endif
#if HAVE_SSE41
func_t sg_stats_table_sse41_functions[] = {
{sg_stats_table,                      "sg_stats_table",                      "sg_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{sg_stats_table_scan,                 "sg_stats_table_scan",                 "sg_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
{sg_stats_table_scan_sse41_128_32,    "sg_stats_table_scan_sse41_128_32",    "sg_stats",    "scan", "sse41", "128", "32",  4, 1, 1, 0},
{sg_stats_table_scan_sse41_128_16,    "sg_stats_table_scan_sse41_128_16",    "sg_stats",    "scan", "sse41", "128", "16",  8, 1, 1, 0},
{sg_stats_table_scan_sse41_128_8,     "sg_stats_table_scan_sse41_128_8",     "sg_stats",    "scan", "sse41", "128",  "8", 16, 1, 1, 0},
{sg_stats_table_striped_sse41_128_32, "sg_stats_table_striped_sse41_128_32", "sg_stats", "striped", "sse41", "128", "32",  4, 1, 1, 0},
{sg_stats_table_striped_sse41_128_16, "sg_stats_table_striped_sse41_128_16", "sg_stats", "striped", "sse41", "128", "16",  8, 1, 1, 0},
{sg_stats_table_striped_sse41_128_8,  "sg_stats_table_striped_sse41_128_8",  "sg_stats", "striped", "sse41", "128",  "8", 16, 1, 1, 0},
{sg_stats_table_diag_sse41_128_32,    "sg_stats_table_diag_sse41_128_32",    "sg_stats",    "diag", "sse41", "128", "32",  4, 1, 1, 0},
{sg_stats_table_diag_sse41_128_16,    "sg_stats_table_diag_sse41_128_16",    "sg_stats",    "diag", "sse41", "128", "16",  8, 1, 1, 0},
{sg_stats_table_diag_sse41_128_8,     "sg_stats_table_diag_sse41_128_8",     "sg_stats",    "diag", "sse41", "128",  "8", 16, 1, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_stats_table_sse41 = {"sg_stats_table_sse41", sg_stats_table_sse41_functions};
#endif
#if HAVE_AVX2
func_t sg_stats_table_avx2_functions[] = {
{sg_stats_table,                      "sg_stats_table",                      "sg_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{sg_stats_table_scan,                 "sg_stats_table_scan",                 "sg_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
{sg_stats_table_scan_avx2_256_32,     "sg_stats_table_scan_avx2_256_32",     "sg_stats",    "scan", "avx2",  "256", "32",  8, 1, 1, 0},
{sg_stats_table_scan_avx2_256_16,     "sg_stats_table_scan_avx2_256_16",     "sg_stats",    "scan", "avx2",  "256", "16", 16, 1, 1, 0},
{sg_stats_table_scan_avx2_256_8,      "sg_stats_table_scan_avx2_256_8",      "sg_stats",    "scan", "avx2",  "256",  "8", 32, 1, 1, 0},
{sg_stats_table_striped_avx2_256_32,  "sg_stats_table_striped_avx2_256_32",  "sg_stats", "striped", "avx2",  "256", "32",  8, 1, 1, 0},
{sg_stats_table_striped_avx2_256_16,  "sg_stats_table_striped_avx2_256_16",  "sg_stats", "striped", "avx2",  "256", "16", 16, 1, 1, 0},
{sg_stats_table_striped_avx2_256_8,   "sg_stats_table_striped_avx2_256_8",   "sg_stats", "striped", "avx2",  "256",  "8", 32, 1, 1, 0},
{sg_stats_table_diag_avx2_256_32,     "sg_stats_table_diag_avx2_256_32",     "sg_stats",    "diag", "avx2",  "256", "32",  8, 1, 1, 0},
{sg_stats_table_diag_avx2_256_16,     "sg_stats_table_diag_avx2_256_16",     "sg_stats",    "diag", "avx2",  "256", "16", 16, 1, 1, 0},
{sg_stats_table_diag_avx2_256_8,      "sg_stats_table_diag_avx2_256_8",      "sg_stats",    "diag", "avx2",  "256",  "8", 32, 1, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_stats_table_avx2 = {"sg_stats_table_avx2", sg_stats_table_avx2_functions};
#endif
#if HAVE_KNC
func_t sg_stats_table_knc_functions[] = {
{sg_stats_table,                      "sg_stats_table",                      "sg_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{sg_stats_table_scan,                 "sg_stats_table_scan",                 "sg_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
{sg_stats_table_scan_knc_512_32,      "sg_stats_table_scan_knc_512_32",      "sg_stats",    "scan", "knc",   "512", "32", 16, 1, 1, 0},
{sg_stats_table_striped_knc_512_32,   "sg_stats_table_striped_knc_512_32",   "sg_stats", "striped", "knc",   "512", "32", 16, 1, 1, 0},
{sg_stats_table_diag_knc_512_32,      "sg_stats_table_diag_knc_512_32",      "sg_stats",    "diag", "knc",   "512", "32", 16, 1, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sg_stats_table_knc = {"sg_stats_table_knc", sg_stats_table_knc_functions};
#endif
#if HAVE_SSE2
func_t sw_stats_table_sse2_functions[] = {
{sw_stats_table,                      "sw_stats_table",                      "sw_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{sw_stats_table_scan,                 "sw_stats_table_scan",                 "sw_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
{sw_stats_table_scan_sse2_128_32,     "sw_stats_table_scan_sse2_128_32",     "sw_stats",    "scan", "sse2",  "128", "32",  4, 1, 1, 0},
{sw_stats_table_scan_sse2_128_16,     "sw_stats_table_scan_sse2_128_16",     "sw_stats",    "scan", "sse2",  "128", "16",  8, 1, 1, 0},
{sw_stats_table_scan_sse2_128_8,      "sw_stats_table_scan_sse2_128_8",      "sw_stats",    "scan", "sse2",  "128",  "8", 16, 1, 1, 0},
{sw_stats_table_striped_sse2_128_32,  "sw_stats_table_striped_sse2_128_32",  "sw_stats", "striped", "sse2",  "128", "32",  4, 1, 1, 0},
{sw_stats_table_striped_sse2_128_16,  "sw_stats_table_striped_sse2_128_16",  "sw_stats", "striped", "sse2",  "128", "16",  8, 1, 1, 0},
{sw_stats_table_striped_sse2_128_8,   "sw_stats_table_striped_sse2_128_8",   "sw_stats", "striped", "sse2",  "128",  "8", 16, 1, 1, 0},
{sw_stats_table_diag_sse2_128_32,     "sw_stats_table_diag_sse2_128_32",     "sw_stats",    "diag", "sse2",  "128", "32",  4, 1, 1, 0},
{sw_stats_table_diag_sse2_128_16,     "sw_stats_table_diag_sse2_128_16",     "sw_stats",    "diag", "sse2",  "128", "16",  8, 1, 1, 0},
{sw_stats_table_diag_sse2_128_8,      "sw_stats_table_diag_sse2_128_8",      "sw_stats",    "diag", "sse2",  "128",  "8", 16, 1, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_stats_table_sse2 = {"sw_stats_table_sse2", sw_stats_table_sse2_functions};
#endif
#if HAVE_SSE41
func_t sw_stats_table_sse41_functions[] = {
{sw_stats_table,                      "sw_stats_table",                      "sw_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{sw_stats_table_scan,                 "sw_stats_table_scan",                 "sw_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
{sw_stats_table_scan_sse41_128_32,    "sw_stats_table_scan_sse41_128_32",    "sw_stats",    "scan", "sse41", "128", "32",  4, 1, 1, 0},
{sw_stats_table_scan_sse41_128_16,    "sw_stats_table_scan_sse41_128_16",    "sw_stats",    "scan", "sse41", "128", "16",  8, 1, 1, 0},
{sw_stats_table_scan_sse41_128_8,     "sw_stats_table_scan_sse41_128_8",     "sw_stats",    "scan", "sse41", "128",  "8", 16, 1, 1, 0},
{sw_stats_table_striped_sse41_128_32, "sw_stats_table_striped_sse41_128_32", "sw_stats", "striped", "sse41", "128", "32",  4, 1, 1, 0},
{sw_stats_table_striped_sse41_128_16, "sw_stats_table_striped_sse41_128_16", "sw_stats", "striped", "sse41", "128", "16",  8, 1, 1, 0},
{sw_stats_table_striped_sse41_128_8,  "sw_stats_table_striped_sse41_128_8",  "sw_stats", "striped", "sse41", "128",  "8", 16, 1, 1, 0},
{sw_stats_table_diag_sse41_128_32,    "sw_stats_table_diag_sse41_128_32",    "sw_stats",    "diag", "sse41", "128", "32",  4, 1, 1, 0},
{sw_stats_table_diag_sse41_128_16,    "sw_stats_table_diag_sse41_128_16",    "sw_stats",    "diag", "sse41", "128", "16",  8, 1, 1, 0},
{sw_stats_table_diag_sse41_128_8,     "sw_stats_table_diag_sse41_128_8",     "sw_stats",    "diag", "sse41", "128",  "8", 16, 1, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_stats_table_sse41 = {"sw_stats_table_sse41", sw_stats_table_sse41_functions};
#endif
#if HAVE_AVX2
func_t sw_stats_table_avx2_functions[] = {
{sw_stats_table,                      "sw_stats_table",                      "sw_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{sw_stats_table_scan,                 "sw_stats_table_scan",                 "sw_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
{sw_stats_table_scan_avx2_256_32,     "sw_stats_table_scan_avx2_256_32",     "sw_stats",    "scan", "avx2",  "256", "32",  8, 1, 1, 0},
{sw_stats_table_scan_avx2_256_16,     "sw_stats_table_scan_avx2_256_16",     "sw_stats",    "scan", "avx2",  "256", "16", 16, 1, 1, 0},
{sw_stats_table_scan_avx2_256_8,      "sw_stats_table_scan_avx2_256_8",      "sw_stats",    "scan", "avx2",  "256",  "8", 32, 1, 1, 0},
{sw_stats_table_striped_avx2_256_32,  "sw_stats_table_striped_avx2_256_32",  "sw_stats", "striped", "avx2",  "256", "32",  8, 1, 1, 0},
{sw_stats_table_striped_avx2_256_16,  "sw_stats_table_striped_avx2_256_16",  "sw_stats", "striped", "avx2",  "256", "16", 16, 1, 1, 0},
{sw_stats_table_striped_avx2_256_8,   "sw_stats_table_striped_avx2_256_8",   "sw_stats", "striped", "avx2",  "256",  "8", 32, 1, 1, 0},
{sw_stats_table_diag_avx2_256_32,     "sw_stats_table_diag_avx2_256_32",     "sw_stats",    "diag", "avx2",  "256", "32",  8, 1, 1, 0},
{sw_stats_table_diag_avx2_256_16,     "sw_stats_table_diag_avx2_256_16",     "sw_stats",    "diag", "avx2",  "256", "16", 16, 1, 1, 0},
{sw_stats_table_diag_avx2_256_8,      "sw_stats_table_diag_avx2_256_8",      "sw_stats",    "diag", "avx2",  "256",  "8", 32, 1, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_stats_table_avx2 = {"sw_stats_table_avx2", sw_stats_table_avx2_functions};
#endif
#if HAVE_KNC
func_t sw_stats_table_knc_functions[] = {
{sw_stats_table,                      "sw_stats_table",                      "sw_stats",    "orig", "NA",     "32", "32",  1, 1, 1, 1},
{sw_stats_table_scan,                 "sw_stats_table_scan",                 "sw_stats",    "scan", "NA",     "32", "32",  1, 1, 1, 0},
{sw_stats_table_scan_knc_512_32,      "sw_stats_table_scan_knc_512_32",      "sw_stats",    "scan", "knc",   "512", "32", 16, 1, 1, 0},
{sw_stats_table_striped_knc_512_32,   "sw_stats_table_striped_knc_512_32",   "sw_stats", "striped", "knc",   "512", "32", 16, 1, 1, 0},
{sw_stats_table_diag_knc_512_32,      "sw_stats_table_diag_knc_512_32",      "sw_stats",    "diag", "knc",   "512", "32", 16, 1, 1, 0},
{NULL, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0}
};
funcs_t sw_stats_table_knc = {"sw_stats_table_knc", sw_stats_table_knc_functions};
#endif

parasail_function_t lookup_function(const char *funcname)
{
    parasail_function_t function = NULL;

    if (funcname) {
        int index = 0;
        func_t f;
        f = functions[index++];
        while (f.pointer) {
            if (0 == strcmp(funcname, f.name)) {
                function = f.pointer;
                break;
            }
            f = functions[index++];
        }
    }

    return function;
}

#endif /* _PARASAIL_FUNCTION_TYPE_H_ */
