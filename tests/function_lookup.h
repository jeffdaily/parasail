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
    const char * name;
    parasail_function_t pointer;
    int lanes;
} func_t;

func_t functions[] = {
    {"nw",                           nw,                         1},
    {"nw_scan",                      nw_scan,                    1},
#if HAVE_SSE2
    {"nw_scan_sse2_128_32",          nw_scan_sse2_128_32,        4},
    {"nw_scan_sse2_128_16",          nw_scan_sse2_128_16,        8},
    {"nw_scan_sse2_128_8",           nw_scan_sse2_128_8,        16},
    {"nw_diag_sse2_128_32",          nw_diag_sse2_128_32,        4},
    {"nw_diag_sse2_128_16",          nw_diag_sse2_128_16,        8},
    {"nw_diag_sse2_128_8",           nw_diag_sse2_128_8,        16},
    {"nw_striped_sse2_128_32",       nw_striped_sse2_128_32,     4},
    {"nw_striped_sse2_128_16",       nw_striped_sse2_128_16,     8},
    {"nw_striped_sse2_128_8",        nw_striped_sse2_128_8,     16},
#endif
#if HAVE_SSE41
    {"nw_scan_sse41_128_32",         nw_scan_sse41_128_32,       4},
    {"nw_scan_sse41_128_16",         nw_scan_sse41_128_16,       8},
    {"nw_scan_sse41_128_8",          nw_scan_sse41_128_8,       16},
    {"nw_diag_sse41_128_32",         nw_diag_sse41_128_32,       4},
    {"nw_diag_sse41_128_16",         nw_diag_sse41_128_16,       8},
    {"nw_diag_sse41_128_8",          nw_diag_sse41_128_8,       16},
    {"nw_striped_sse41_128_32",      nw_striped_sse41_128_32,    4},
    {"nw_striped_sse41_128_16",      nw_striped_sse41_128_16,    8},
    {"nw_striped_sse41_128_8",       nw_striped_sse41_128_8,    16},
#endif
#if HAVE_AVX2
    {"nw_scan_avx2_256_32",          nw_scan_avx2_256_32,        8},
    {"nw_scan_avx2_256_16",          nw_scan_avx2_256_16,       16},
    {"nw_scan_avx2_256_8",           nw_scan_avx2_256_8,        32},
    {"nw_diag_avx2_256_32",          nw_diag_avx2_256_32,        8},
    {"nw_diag_avx2_256_16",          nw_diag_avx2_256_16,       16},
    {"nw_diag_avx2_256_8",           nw_diag_avx2_256_8,        32},
    {"nw_striped_avx2_256_32",       nw_striped_avx2_256_32,     8},
    {"nw_striped_avx2_256_16",       nw_striped_avx2_256_16,    16},
    {"nw_striped_avx2_256_8",        nw_striped_avx2_256_8,     32},
#endif
#if HAVE_KNC
    {"nw_scan_knc_512_32",           nw_scan_knc_512_32,        16},
    {"nw_diag_knc_512_32",           nw_diag_knc_512_32,        16},
    {"nw_striped_knc_512_32",        nw_striped_knc_512_32,     16},
#endif

    {"sg",                           sg,                         1},
    {"sg_scan",                      sg_scan,                    1},
#if HAVE_SSE2
    {"sg_scan_sse2_128_32",          sg_scan_sse2_128_32,        4},
    {"sg_scan_sse2_128_16",          sg_scan_sse2_128_16,        8},
    {"sg_scan_sse2_128_8",           sg_scan_sse2_128_8,        16},
    {"sg_diag_sse2_128_32",          sg_diag_sse2_128_32,        4},
    {"sg_diag_sse2_128_16",          sg_diag_sse2_128_16,        8},
    {"sg_diag_sse2_128_8",           sg_diag_sse2_128_8,        16},
    {"sg_striped_sse2_128_32",       sg_striped_sse2_128_32,     4},
    {"sg_striped_sse2_128_16",       sg_striped_sse2_128_16,     8},
    {"sg_striped_sse2_128_8",        sg_striped_sse2_128_8,     16},
#endif
#if HAVE_SSE41
    {"sg_scan_sse41_128_32",         sg_scan_sse41_128_32,       4},
    {"sg_scan_sse41_128_16",         sg_scan_sse41_128_16,       8},
    {"sg_scan_sse41_128_8",          sg_scan_sse41_128_8,       16},
    {"sg_diag_sse41_128_32",         sg_diag_sse41_128_32,       4},
    {"sg_diag_sse41_128_16",         sg_diag_sse41_128_16,       8},
    {"sg_diag_sse41_128_8",          sg_diag_sse41_128_8,       16},
    {"sg_striped_sse41_128_32",      sg_striped_sse41_128_32,    4},
    {"sg_striped_sse41_128_16",      sg_striped_sse41_128_16,    8},
    {"sg_striped_sse41_128_8",       sg_striped_sse41_128_8,    16},
#endif
#if HAVE_AVX2
    {"sg_scan_avx2_256_32",          sg_scan_avx2_256_32,        8},
    {"sg_scan_avx2_256_16",          sg_scan_avx2_256_16,       16},
    {"sg_scan_avx2_256_8",           sg_scan_avx2_256_8,        32},
    {"sg_diag_avx2_256_32",          sg_diag_avx2_256_32,        8},
    {"sg_diag_avx2_256_16",          sg_diag_avx2_256_16,       16},
    {"sg_diag_avx2_256_8",           sg_diag_avx2_256_8,        32},
    {"sg_striped_avx2_256_32",       sg_striped_avx2_256_32,     8},
    {"sg_striped_avx2_256_16",       sg_striped_avx2_256_16,    16},
    {"sg_striped_avx2_256_8",        sg_striped_avx2_256_8,     32},
#endif
#if HAVE_KNC
    {"sg_scan_knc_512_32",           sg_scan_knc_512_32,        16},
    {"sg_diag_knc_512_32",           sg_diag_knc_512_32,        16},
    {"sg_striped_knc_512_32",        sg_striped_knc_512_32,     16},
#endif

    {"sw",                           sw,                         1},
    {"sw_scan",                      sw_scan,                    1},
#if HAVE_SSE2
    {"sw_scan_sse2_128_32",          sw_scan_sse2_128_32,        4},
    {"sw_scan_sse2_128_16",          sw_scan_sse2_128_16,        8},
    {"sw_scan_sse2_128_8",           sw_scan_sse2_128_8,        16},
    {"sw_diag_sse2_128_32",          sw_diag_sse2_128_32,        4},
    {"sw_diag_sse2_128_16",          sw_diag_sse2_128_16,        8},
    {"sw_diag_sse2_128_8",           sw_diag_sse2_128_8,        16},
    {"sw_striped_sse2_128_32",       sw_striped_sse2_128_32,     4},
    {"sw_striped_sse2_128_16",       sw_striped_sse2_128_16,     8},
    {"sw_striped_sse2_128_8",        sw_striped_sse2_128_8,     16},
#endif
#if HAVE_SSE41
    {"sw_scan_sse41_128_32",         sw_scan_sse41_128_32,       4},
    {"sw_scan_sse41_128_16",         sw_scan_sse41_128_16,       8},
    {"sw_scan_sse41_128_8",          sw_scan_sse41_128_8,       16},
    {"sw_diag_sse41_128_32",         sw_diag_sse41_128_32,       4},
    {"sw_diag_sse41_128_16",         sw_diag_sse41_128_16,       8},
    {"sw_diag_sse41_128_8",          sw_diag_sse41_128_8,       16},
    {"sw_striped_sse41_128_32",      sw_striped_sse41_128_32,    4},
    {"sw_striped_sse41_128_16",      sw_striped_sse41_128_16,    8},
    {"sw_striped_sse41_128_8",       sw_striped_sse41_128_8,    16},
#endif
#if HAVE_AVX2
    {"sw_scan_avx2_256_32",          sw_scan_avx2_256_32,        8},
    {"sw_scan_avx2_256_16",          sw_scan_avx2_256_16,       16},
    {"sw_scan_avx2_256_8",           sw_scan_avx2_256_8,        32},
    {"sw_diag_avx2_256_32",          sw_diag_avx2_256_32,        8},
    {"sw_diag_avx2_256_16",          sw_diag_avx2_256_16,       16},
    {"sw_diag_avx2_256_8",           sw_diag_avx2_256_8,        32},
    {"sw_striped_avx2_256_32",       sw_striped_avx2_256_32,     8},
    {"sw_striped_avx2_256_16",       sw_striped_avx2_256_16,    16},
    {"sw_striped_avx2_256_8",        sw_striped_avx2_256_8,     32},
#endif
#if HAVE_KNC
    {"sw_scan_knc_512_32",           sw_scan_knc_512_32,        16},
    {"sw_diag_knc_512_32",           sw_diag_knc_512_32,        16},
    {"sw_striped_knc_512_32",        sw_striped_knc_512_32,     16},
#endif


    {"nw_stats",                     nw_stats,                      1},
    {"nw_stats_scan",                nw_stats_scan,                 1},
#if HAVE_SSE2
    {"nw_stats_scan_sse2_128_32",    nw_stats_scan_sse2_128_32,     4},
    {"nw_stats_scan_sse2_128_16",    nw_stats_scan_sse2_128_16,     8},
    {"nw_stats_scan_sse2_128_8",     nw_stats_scan_sse2_128_8,     16},
    {"nw_stats_diag_sse2_128_32",    nw_stats_diag_sse2_128_32,     4},
    {"nw_stats_diag_sse2_128_16",    nw_stats_diag_sse2_128_16,     8},
    {"nw_stats_diag_sse2_128_8",     nw_stats_diag_sse2_128_8,     16},
    {"nw_stats_striped_sse2_128_32", nw_stats_striped_sse2_128_32,  4},
    {"nw_stats_striped_sse2_128_16", nw_stats_striped_sse2_128_16,  8},
    {"nw_stats_striped_sse2_128_8",  nw_stats_striped_sse2_128_8,  16},
#endif
#if HAVE_SSE41
    {"nw_stats_scan_sse41_128_32",   nw_stats_scan_sse41_128_32,    4},
    {"nw_stats_scan_sse41_128_16",   nw_stats_scan_sse41_128_16,    8},
    {"nw_stats_scan_sse41_128_8",    nw_stats_scan_sse41_128_8,    16},
    {"nw_stats_diag_sse41_128_32",   nw_stats_diag_sse41_128_32,    4},
    {"nw_stats_diag_sse41_128_16",   nw_stats_diag_sse41_128_16,    8},
    {"nw_stats_diag_sse41_128_8",    nw_stats_diag_sse41_128_8,    16},
    {"nw_stats_striped_sse41_128_32",nw_stats_striped_sse41_128_32, 4},
    {"nw_stats_striped_sse41_128_16",nw_stats_striped_sse41_128_16, 8},
    {"nw_stats_striped_sse41_128_8", nw_stats_striped_sse41_128_8, 16},
#endif
#if HAVE_AVX2
    {"nw_stats_scan_avx2_256_32",    nw_stats_scan_avx2_256_32,     8},
    {"nw_stats_scan_avx2_256_16",    nw_stats_scan_avx2_256_16,    16},
    {"nw_stats_scan_avx2_256_8",     nw_stats_scan_avx2_256_8,     32},
    {"nw_stats_diag_avx2_256_32",    nw_stats_diag_avx2_256_32,     8},
    {"nw_stats_diag_avx2_256_16",    nw_stats_diag_avx2_256_16,    16},
    {"nw_stats_diag_avx2_256_8",     nw_stats_diag_avx2_256_8,     32},
    {"nw_stats_striped_avx2_256_32", nw_stats_striped_avx2_256_32,  8},
    {"nw_stats_striped_avx2_256_16", nw_stats_striped_avx2_256_16, 16},
    {"nw_stats_striped_avx2_256_8",  nw_stats_striped_avx2_256_8,  32},
#endif
#if HAVE_KNC
    {"nw_stats_scan_knc_512_32",     nw_stats_scan_knc_512_32,     16},
    {"nw_stats_diag_knc_512_32",     nw_stats_diag_knc_512_32,     16},
    {"nw_stats_striped_knc_512_32",  nw_stats_striped_knc_512_32,  16},
#endif

    {"sg_stats",                     sg_stats,                      1},
    {"sg_stats_scan",                sg_stats_scan,                 1},
#if HAVE_SSE2
    {"sg_stats_scan_sse2_128_32",    sg_stats_scan_sse2_128_32,     4},
    {"sg_stats_scan_sse2_128_16",    sg_stats_scan_sse2_128_16,     8},
    {"sg_stats_scan_sse2_128_8",     sg_stats_scan_sse2_128_8,     16},
    {"sg_stats_diag_sse2_128_32",    sg_stats_diag_sse2_128_32,     4},
    {"sg_stats_diag_sse2_128_16",    sg_stats_diag_sse2_128_16,     8},
    {"sg_stats_diag_sse2_128_8",     sg_stats_diag_sse2_128_8,     16},
    {"sg_stats_striped_sse2_128_32", sg_stats_striped_sse2_128_32,  4},
    {"sg_stats_striped_sse2_128_16", sg_stats_striped_sse2_128_16,  8},
    {"sg_stats_striped_sse2_128_8",  sg_stats_striped_sse2_128_8,  16},
#endif
#if HAVE_SSE41
    {"sg_stats_scan_sse41_128_32",   sg_stats_scan_sse41_128_32,    4},
    {"sg_stats_scan_sse41_128_16",   sg_stats_scan_sse41_128_16,    8},
    {"sg_stats_scan_sse41_128_8",    sg_stats_scan_sse41_128_8,    16},
    {"sg_stats_diag_sse41_128_32",   sg_stats_diag_sse41_128_32,    4},
    {"sg_stats_diag_sse41_128_16",   sg_stats_diag_sse41_128_16,    8},
    {"sg_stats_diag_sse41_128_8",    sg_stats_diag_sse41_128_8,    16},
    {"sg_stats_striped_sse41_128_32",sg_stats_striped_sse41_128_32, 4},
    {"sg_stats_striped_sse41_128_16",sg_stats_striped_sse41_128_16, 8},
    {"sg_stats_striped_sse41_128_8", sg_stats_striped_sse41_128_8, 16},
#endif
#if HAVE_AVX2
    {"sg_stats_scan_avx2_256_32",    sg_stats_scan_avx2_256_32,     8},
    {"sg_stats_scan_avx2_256_16",    sg_stats_scan_avx2_256_16,    16},
    {"sg_stats_scan_avx2_256_8",     sg_stats_scan_avx2_256_8,     32},
    {"sg_stats_diag_avx2_256_32",    sg_stats_diag_avx2_256_32,     8},
    {"sg_stats_diag_avx2_256_16",    sg_stats_diag_avx2_256_16,    16},
    {"sg_stats_diag_avx2_256_8",     sg_stats_diag_avx2_256_8,     32},
    {"sg_stats_striped_avx2_256_32", sg_stats_striped_avx2_256_32,  8},
    {"sg_stats_striped_avx2_256_16", sg_stats_striped_avx2_256_16, 16},
    {"sg_stats_striped_avx2_256_8",  sg_stats_striped_avx2_256_8,  32},
#endif
#if HAVE_KNC
    {"sg_stats_scan_knc_512_32",     sg_stats_scan_knc_512_32,     16},
    {"sg_stats_diag_knc_512_32",     sg_stats_diag_knc_512_32,     16},
    {"sg_stats_striped_knc_512_32",  sg_stats_striped_knc_512_32,  16},
#endif

    {"sw_stats",                     sw_stats,                      1},
    {"sw_stats_scan",                sw_stats_scan,                 1},
#if HAVE_SSE2
    {"sw_stats_scan_sse2_128_32",    sw_stats_scan_sse2_128_32,     4},
    {"sw_stats_scan_sse2_128_16",    sw_stats_scan_sse2_128_16,     8},
    {"sw_stats_scan_sse2_128_8",     sw_stats_scan_sse2_128_8,     16},
    {"sw_stats_diag_sse2_128_32",    sw_stats_diag_sse2_128_32,     4},
    {"sw_stats_diag_sse2_128_16",    sw_stats_diag_sse2_128_16,     8},
    {"sw_stats_diag_sse2_128_8",     sw_stats_diag_sse2_128_8,     16},
    {"sw_stats_striped_sse2_128_32", sw_stats_striped_sse2_128_32,  4},
    {"sw_stats_striped_sse2_128_16", sw_stats_striped_sse2_128_16,  8},
    {"sw_stats_striped_sse2_128_8",  sw_stats_striped_sse2_128_8,  16},
#endif
#if HAVE_SSE41
    {"sw_stats_scan_sse41_128_32",   sw_stats_scan_sse41_128_32,    4},
    {"sw_stats_scan_sse41_128_16",   sw_stats_scan_sse41_128_16,    8},
    {"sw_stats_scan_sse41_128_8",    sw_stats_scan_sse41_128_8,    16},
    {"sw_stats_diag_sse41_128_32",   sw_stats_diag_sse41_128_32,    4},
    {"sw_stats_diag_sse41_128_16",   sw_stats_diag_sse41_128_16,    8},
    {"sw_stats_diag_sse41_128_8",    sw_stats_diag_sse41_128_8,    16},
    {"sw_stats_striped_sse41_128_32",sw_stats_striped_sse41_128_32, 4},
    {"sw_stats_striped_sse41_128_16",sw_stats_striped_sse41_128_16, 8},
    {"sw_stats_striped_sse41_128_8", sw_stats_striped_sse41_128_8, 16},
#endif
#if HAVE_AVX2
    {"sw_stats_scan_avx2_256_32",    sw_stats_scan_avx2_256_32,     8},
    {"sw_stats_scan_avx2_256_16",    sw_stats_scan_avx2_256_16,    16},
    {"sw_stats_scan_avx2_256_8",     sw_stats_scan_avx2_256_8,     32},
    {"sw_stats_diag_avx2_256_32",    sw_stats_diag_avx2_256_32,     8},
    {"sw_stats_diag_avx2_256_16",    sw_stats_diag_avx2_256_16,    16},
    {"sw_stats_diag_avx2_256_8",     sw_stats_diag_avx2_256_8,     32},
    {"sw_stats_striped_avx2_256_32", sw_stats_striped_avx2_256_32,  8},
    {"sw_stats_striped_avx2_256_16", sw_stats_striped_avx2_256_16, 16},
    {"sw_stats_striped_avx2_256_8",  sw_stats_striped_avx2_256_8,  32},
#endif
#if HAVE_KNC
    {"sw_stats_scan_knc_512_32",     sw_stats_scan_knc_512_32,     16},
    {"sw_stats_diag_knc_512_32",     sw_stats_diag_knc_512_32,     16},
    {"sw_stats_striped_knc_512_32",  sw_stats_striped_knc_512_32,  16},
#endif
    {"NULL", NULL, 0}
};

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
