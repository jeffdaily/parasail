#!/usr/bin/python

# Generate the various function dispatcher implemementations.

import os

def codegen():
    txt = """/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include "parasail.h"
#include "parasail/cpuid.h"

/* forward declare the dispatcher functions */
"""
    for stats in ["", "stats_"]:
        for width in ["64", "32", "16", "8", "sat"]:
            prefix = "parasail_profile_create_%s%s"%(stats, width)
            txt += "parasail_pcreator_t %s_dispatcher;\n" % prefix

    txt += """
/* declare and initialize the pointer to the dispatcher function */
"""
    for stats in ["", "stats_"]:
        for width in ["64", "32", "16", "8", "sat"]:
            prefix = "parasail_profile_create_%s%s"%(stats, width)
            txt += "parasail_pcreator_t * %s_pointer = %s_dispatcher;\n"%(prefix, prefix)

    txt += """
/* dispatcher function implementations */
"""
    for stats in ["", "stats_"]:
        for width in ["64", "32", "16", "8", "sat"]:
            prefix = "parasail_profile_create_%s%s"%(stats, width)
            params = {
                "PREFIX": prefix,
                "WIDTH": width
            }
            txt += """
parasail_profile_t* %(PREFIX)s_dispatcher(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        %(PREFIX)s_pointer = parasail_profile_create_avx_256_%(WIDTH)s;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        %(PREFIX)s_pointer = parasail_profile_create_sse_128_%(WIDTH)s;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        %(PREFIX)s_pointer = parasail_profile_create_sse_128_%(WIDTH)s;
    }
    else
#endif
#if HAVE_ALTIVEC
    if (parasail_can_use_altivec()) {
        %(PREFIX)s_pointer = parasail_profile_create_altivec_128_%(WIDTH)s;
    }
    else
#endif
#if HAVE_NEON
    if (parasail_can_use_neon()) {
        %(PREFIX)s_pointer = parasail_profile_create_neon_128_%(WIDTH)s;
    }
    else
#endif
    {
        /* no fallback; caller must check for non-NULL profile */
        return NULL;
    }
    return %(PREFIX)s_pointer(s1, s1Len, matrix);
}
""" % params

    txt += """
/* implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct impl */
"""
    for stats in ["", "stats_"]:
        for width in ["64", "32", "16", "8", "sat"]:
            prefix = "parasail_profile_create_%s%s"%(stats, width)
            params = {
                    "PREFIX": prefix,
                    "WIDTH": width
            }
            txt += """
parasail_profile_t* %(PREFIX)s(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    return %(PREFIX)s_pointer(s1, s1Len, matrix);
}
""" % params

    return txt

output_dir = "generated/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

output_filename = "%sdispatch_profile.c" % output_dir
writer = open(output_filename, "w")
writer.write(codegen())
writer.write("\n")
writer.close()
