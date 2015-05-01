#!/usr/bin/python

# Generate the various function dispatcher implemementations.

import os

def codegen(alg):
    txt = """/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include "parasail.h"
#include "parasail/cpuid.h"

/* forward declare the dispatcher functions */
"""
    for table in ["", "_table"]:
        for stats in ["", "_stats"]:
            for par in ["scan", "striped", "diag"]:
                for width in [64, 32, 16, 8]:
                    prefix = "parasail_%s%s%s_%s_%d"%(
                        alg, stats, table, par, width)
                    txt += "parasail_function_t %s_dispatcher;\n" % prefix

    txt += """
/* declare and initialize the pointer to the dispatcher function */
"""
    for table in ["", "_table"]:
        for stats in ["", "_stats"]:
            for par in ["scan", "striped", "diag"]:
                for width in [64, 32, 16, 8]:
                    prefix = "parasail_%s%s%s_%s_%d"%(
                        alg, stats, table, par, width)
                    txt += "parasail_function_t * %s_pointer = %s_dispatcher;\n"%(
                            prefix, prefix)
    txt += """
/* dispatcher function implementations */
"""
    for table in ["", "_table"]:
        for stats in ["", "_stats"]:
            for par in ["scan", "striped", "diag"]:
                for width in [64, 32, 16, 8]:
                    prefix = "parasail_%s%s%s_%s_%d"%(
                        alg, stats, table, par, width)
                    prefix2 = "parasail_%s%s%s_%s"%(
                        alg, stats, table, par)
                    base = alg
                    if par == "scan":
                        base += "_scan"
                    params = {
                            "ALG": alg,
                            "BASE": base,
                            "PREFIX": prefix,
                            "PREFIX2": prefix2,
                            "TABLE": table,
                            "STATS": stats,
                            "PAR": par,
                            "WIDTH": width
                    }
                    txt += """
parasail_result_t* %(PREFIX)s_dispatcher(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
#if HAVE_KNC
    %(PREFIX)s_pointer = %(PREFIX2)s_knc_512_32;
#else
#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        %(PREFIX)s_pointer = %(PREFIX2)s_avx2_256_%(WIDTH)s;
    }
    else
#endif
#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        %(PREFIX)s_pointer = %(PREFIX2)s_sse41_128_%(WIDTH)s;
    }
    else
#endif
#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        %(PREFIX)s_pointer = %(PREFIX2)s_sse2_128_%(WIDTH)s;
    }
    else
#endif
#endif
    {
        %(PREFIX)s_pointer = parasail_%(BASE)s;
    }
    return %(PREFIX)s_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}
""" % params
    txt += """
/* implementation which simply calls the pointer,
 * first time it's the dispatcher, otherwise it's correct impl */
"""
    for table in ["", "_table"]:
        for stats in ["", "_stats"]:
            for par in ["scan", "striped", "diag"]:
                for width in [64, 32, 16, 8]:
                    prefix = "parasail_%s%s%s_%s_%d"%(
                        alg, stats, table, par, width)
                    params = {
                            "ALG": alg,
                            "PREFIX": prefix,
                            "TABLE": table,
                            "STATS": stats,
                            "PAR": par,
                            "WIDTH": width
                    }
                    txt += """
parasail_result_t* %(PREFIX)s(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    return %(PREFIX)s_pointer(s1, s1Len, s2, s2Len, open, gap, matrix);
}
""" % params
    return txt

output_dir = "generated/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for alg in ["nw", "sg", "sw"]:
    output_filename = "%s%s_dispatch.c" % (output_dir, alg)
    writer = open(output_filename, "w")
    writer.write(codegen(alg))
    writer.write("\n")
    writer.close()
