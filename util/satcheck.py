#!/usr/bin/python

# Generate the various saturation checking and recalculating
# implementations.

import os

def codegen(alg):
    txt = """/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "parasail.h"

"""
    for table in ["", "_table", "_rowcol"]:
        for stats in ["", "_stats"]:
            for par in ["scan", "striped", "diag"]:
                prefix = "parasail_%s%s%s_%s"%(alg, stats, table, par)
                params = {"PREFIX":prefix}
                txt += """
parasail_result_t* %(PREFIX)s_sat(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t *matrix)
{
    parasail_result_t * result = NULL;
    
    result = %(PREFIX)s_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        parasail_result_free(result);
        result = %(PREFIX)s_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = %(PREFIX)s_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }

    return result;
}
""" % params

    for table in ["", "_table", "_rowcol"]:
        for stats in [""]:
            for par in ["scan_profile", "striped_profile"]:
                prefix = "parasail_%s%s%s_%s"%(alg, stats, table, par)
                params = {"PREFIX":prefix}
                txt += """
parasail_result_t* %(PREFIX)s_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;
    
    result = %(PREFIX)s_8(profile, s2, s2Len, open, gap);
    if (result->saturated) {
        parasail_result_free(result);
        parasail_profile_t *copy = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
        (void)memcpy(copy, profile, sizeof(parasail_profile_t));
        copy->profile = copy->profile2;
        result = %(PREFIX)s_16(copy, s2, s2Len, open, gap);
        free(copy);
    }

    return result;
}
""" % params

    for table in ["", "_table", "_rowcol"]:
        for stats in ["_stats"]:
            for par in ["scan_profile", "striped_profile"]:
                prefix = "parasail_%s%s%s_%s"%(alg, stats, table, par)
                params = {"PREFIX":prefix}
                txt += """
parasail_result_t* %(PREFIX)s_sat(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    parasail_result_t * result = NULL;
    
    result = %(PREFIX)s_8(profile, s2, s2Len, open, gap);
    if (result->saturated) {
        parasail_result_free(result);
        parasail_profile_t *copy = (parasail_profile_t*)malloc(sizeof(parasail_profile_t));
        (void)memcpy(copy, profile, sizeof(parasail_profile_t));
        copy->profile = copy->profile2;
        copy->profile_m = copy->profile2_m;
        copy->profile_s = copy->profile2_s;
        result = %(PREFIX)s_16(copy, s2, s2Len, open, gap);
        free(copy);
    }

    return result;
}
""" % params

    return txt

output_dir = "generated/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for alg in ["nw", "sg", "sw"]:
    output_filename = "%s%s_sat.c" % (output_dir, alg)
    writer = open(output_filename, "w")
    writer.write(codegen(alg))
    writer.write("\n")
    writer.close()

