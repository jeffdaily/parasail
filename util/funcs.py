#!/usr/bin/env python

print """/**
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

#include "parasail.h"

typedef struct func {
    parasail_function_t * pointer;
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

"""


def print_fmt(*args):
    fmt = '{%-36s %-38s %5s %10s %-8s %6s %5s %3s %1s %1s %1s},'
    new_args = [arg for arg in args]
    new_args[0] = '%s,'   % new_args[0]
    new_args[1] = '"%s",' % new_args[1]
    new_args[2] = '"%s",' % new_args[2]
    new_args[3] = '"%s",' % new_args[3]
    new_args[4] = '"%s",' % new_args[4]
    new_args[5] = '"%s",' % new_args[5]
    new_args[6] = '"%s",' % new_args[6]
    new_args[7] = '%d,'   % new_args[7]
    new_args[8] = '%d,'   % new_args[8]
    new_args[9] = '%d,'   % new_args[9]
    new_args[10]= '%d'    % new_args[10]
    print fmt % tuple(new_args)

def print_null():
    fmt = '{%s, "%s", "%s", "%s", "%s", "%s", "%s", %d, %d, %d, %d},'
    print fmt[:-1] % ("NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", 0, 0, 0, 0)

isa_to_bits = {
    "sse2"  : 128,
    "sse41" : 128,
    "avx2"  : 256,
    "knc"   : 512,
}

print "func_t functions[] = {"

for table in ["", "_table"]:
    for stats in ["", "_stats"]:
        for alg in ["nw", "sg", "sw"]:
            is_table = 0
            if table:
                is_table = 1
            is_stats = 0
            if stats:
                is_stats = 1
            pre = alg+stats+table
            print_fmt(pre,         pre,         alg+stats, "orig", "NA", "32", "32", 1, is_table, is_stats, 1)
            print_fmt(pre+"_scan", pre+"_scan", alg+stats, "scan", "NA", "32", "32", 1, is_table, is_stats, 0)
            for isa in ["sse2", "sse41", "avx2"]:
                print "#if HAVE_%s" % isa.upper()
                bits = isa_to_bits[isa]
                for par in ["scan", "striped", "diag"]:
                    widths = [64, 32, 16, 8]
                    # temporary hack until stats codegen catches up
                    if stats:
                        widths = [32, 16, 8]
                    for width in widths:
                        name = "%s_%s_%s_%s_%s" % (pre, par, isa, bits, width)
                        print_fmt(name, name, alg+stats, par, isa, bits, width, bits/width, is_table, is_stats, 0)
                print "#endif"
            for isa in ["knc"]:
                print "#if HAVE_%s" % isa.upper()
                bits = isa_to_bits[isa]
                for par in ["scan", "striped", "diag"]:
                    for width in [32]:
                        name = "%s_%s_%s_%s_%s" % (pre, par, isa, bits, width)
                        print_fmt(name, name, alg+stats, par, isa, bits, width, bits/width, is_table, is_stats, 0)
                print "#endif"
print_fmt("sw_blocked_sse41_128_32", "sw_blocked_sse41_128_32", "sw", "blocked", "sse41", "128", "32", 4, 0, 0, 0)
print_fmt("sw_blocked_sse41_128_16", "sw_blocked_sse41_128_16", "sw", "blocked", "sse41", "128", "16", 8, 0, 0, 0)
print_fmt("sw_table_blocked_sse41_128_32", "sw_blocked_sse41_128_32", "sw", "blocked", "sse41", "128", "32", 4, 1, 0, 0)
print_fmt("sw_table_blocked_sse41_128_16", "sw_blocked_sse41_128_16", "sw", "blocked", "sse41", "128", "16", 8, 1, 0, 0)
print_null()

print "};"

for table in ["", "_table"]:
    for stats in ["", "_stats"]:
        for alg in ["nw", "sg", "sw"]:
            is_table = 0
            if table:
                is_table = 1
            is_stats = 0
            if stats:
                is_stats = 1
            pre = alg+stats+table
            for isa in ["sse2", "sse41", "avx2"]:
                print "#if HAVE_%s" % isa.upper()
                print "func_t %s_%s_functions[] = {" % (pre, isa)
                print_fmt(pre,         pre,         alg+stats, "orig", "NA", "32", "32", 1, is_table, is_stats, 1)
                print_fmt(pre+"_scan", pre+"_scan", alg+stats, "scan", "NA", "32", "32", 1, is_table, is_stats, 0)
                bits = isa_to_bits[isa]
                for par in ["scan", "striped", "diag"]:
                    widths = [64, 32, 16, 8]
                    # temporary hack until stats codegen catches up
                    if stats:
                        widths = [32, 16, 8]
                    for width in widths:
                        name = "%s_%s_%s_%s_%s" % (pre, par, isa, bits, width)
                        print_fmt(name, name, alg+stats, par, isa, bits, width, bits/width, is_table, is_stats, 0)
                print_null()
                print "};"
                print 'funcs_t %s_%s = {"%s_%s", %s_%s_functions};' % ((pre, isa)*3)
                print "#endif"
            for isa in ["knc"]:
                print "#if HAVE_%s" % isa.upper()
                print "func_t %s_%s_functions[] = {" % (pre, isa)
                print_fmt(pre,         pre,         alg+stats, "orig", "NA", "32", "32", 1, is_table, is_stats, 1)
                print_fmt(pre+"_scan", pre+"_scan", alg+stats, "scan", "NA", "32", "32", 1, is_table, is_stats, 0)
                bits = isa_to_bits[isa]
                for par in ["scan", "striped", "diag"]:
                    for width in [32]:
                        name = "%s_%s_%s_%s_%s" % (pre, par, isa, bits, width)
                        print_fmt(name, name, alg+stats, par, isa, bits, width, bits/width, is_table, is_stats, 0)
                print_null()
                print "};"
                print 'funcs_t %s_%s = {"%s_%s", %s_%s_functions};' % ((pre, isa)*3)
                print "#endif"

print """
parasail_function_t * lookup_function(const char *funcname)
{
    parasail_function_t * function = NULL;

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

#endif /* _PARASAIL_FUNCTION_TYPE_H_ */"""
