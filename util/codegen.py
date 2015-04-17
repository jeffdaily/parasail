#!/usr/bin/python
#
# This file contains a C template for striped. The template is filled
# in as appropriate for various instruction sets and integer precisions.
#
# author jeff.daily@pnnl.gov
#
# Copyright (c) 2014 Battelle Memorial Institute.
#
# All rights reserved. No warranty, explicit or implicit, provided.


import copy
import os
import sys

from isa import sse2
from isa import sse41
from isa import avx2

keys = sse2.keys()

# gather templates
template_dir = "templates/"

template_filenames = [
"nw_diag.c",
"nw_scan.c",
"nw_striped.c",
"sg_diag.c",
"sg_scan.c",
"sg_striped.c",
"sw_diag.c",
"sw_scan.c",
"sw_striped.c",

"nw_stats_diag.c",
"sg_stats_diag.c",
"sw_stats_diag.c",
"nw_stats_striped.c",
"sg_stats_striped.c",
"sw_stats_striped.c",
]

special_templates = [
"sg_diag_8.c",
"sw_diag_8.c",
"sw_stats_diag_8.c",
]


output_dir = "generated/"
if not os.path.exists(output_dir):
        os.makedirs(output_dir)


def generate_printer(params):
    text = ""
    if "striped" in params["NAME"] or "scan" in params["NAME"]:
        for lane in range(params["LANES"]):
            params["LANE"] = lane
            if params["LANES"] / 10:
                text += "    array[(%(LANE)2d*seglen+t)*dlen + d] = (%(INT)s)%(VEXTRACT)s(vH, %(LANE)2d);\n" % params
            else:
                text += "    array[(%(LANE)s*seglen+t)*dlen + d] = (%(INT)s)%(VEXTRACT)s(vH, %(LANE)s);\n" % params
        del params["LANE"]
    elif "diag" in params["NAME"]:
        for lane in range(params["LANES"]):
            params["LANE"] = lane
            params["LANE_END"] = params["LANES"]-lane-1
            text += """
    if (0 <= i+%(LANE)s && i+%(LANE)s < s1Len && 0 <= j-%(LANE)s && j-%(LANE)s < s2Len) {
        array[(i+%(LANE)s)*s2Len + (j-%(LANE)s)] = (%(INT)s)%(VEXTRACT)s(vWscore, %(LANE_END)s);
    }\n"""[1:] % params
    else:
        print "bad printer name"
        sys.exit(1)
    params["PRINTER"] = text[:-1] # remove last newline
    return params


def generate_saturation_check_old(params):
    width = params["WIDTH"]
    if width == 8:
        
        params["SATURATION_CHECK_INIT"] = """
    %(VTYPE)s vSaturationCheck = %(VSET0)s();
    %(VTYPE)s vNegLimit = %(VSET1)s(INT8_MIN);
    %(VTYPE)s vPosLimit = %(VSET1)s(INT8_MAX);""".strip() % params

        params["SATURATION_CHECK_MID"] = """
            /* check for saturation */
            {
                vSaturationCheck = %(VOR)s(vSaturationCheck,
                        %(VOR)s(
                            %(VCMPEQ)s(vH, vNegLimit),
                            %(VCMPEQ)s(vH, vPosLimit)));
            }""".strip() % params
            
        params["SATURATION_CHECK_FINAL"] = """
    if (%(VMOVEMASK)s(vSaturationCheck)) {
        result->saturated = 1;
        score = INT8_MAX;
    }""".strip() % params

        params["STATS_SATURATION_CHECK_INIT"] = """
    %(VTYPE)s vSaturationCheck = %(VSET0)s();
    %(VTYPE)s vNegLimit = %(VSET1)s(INT8_MIN);
    %(VTYPE)s vPosLimit = %(VSET1)s(INT8_MAX);""".strip() % params

        params["STATS_SATURATION_CHECK_MID"] = """
            /* check for saturation */
            {
                vSaturationCheck = %(VOR)s(vSaturationCheck,
                        %(VOR)s(
                            %(VCMPEQ)s(vH, vNegLimit),
                            %(VCMPEQ)s(vH, vPosLimit)));
            }""".strip() % params

        params["STATS_SATURATION_CHECK_FINAL"] = """
    if (%(VMOVEMASK)s(vSaturationCheck)) {
        result->saturated = 1;
        score = INT8_MAX;
    }""".strip() % params

        params["NEG_INF"] = "INT8_MIN"
        params["VADD"] = params["VADDSx8"]
        params["VSUB"] = params["VSUBSx8"]
    else:
        params["SATURATION_CHECK_INIT"] = ""
        params["SATURATION_CHECK_MID"] = ""
        params["SATURATION_CHECK_FINAL"] = ""
        params["STATS_SATURATION_CHECK_INIT"] = ""
        params["STATS_SATURATION_CHECK_MID"] = ""
        params["STATS_SATURATION_CHECK_FINAL"] = ""
    return params


def generate_saturation_check(params):
    width = params["WIDTH"]
    if width == 8:
        params["SATURATION_CHECK_INIT"] = """
    %(VTYPE)s vNegLimit = %(VSET1)s(INT8_MIN);
    %(VTYPE)s vPosLimit = %(VSET1)s(INT8_MAX);
    %(VTYPE)s vSaturationCheckMin = vPosLimit;
    %(VTYPE)s vSaturationCheckMax = vNegLimit;""".strip() % params
        if "diag" in params["NAME"]:
            params["SATURATION_CHECK_MID"] = """
            /* check for saturation */
            {
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWscore);
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vWscore);
            }""".strip() % params
        else:
            params["SATURATION_CHECK_MID"] = """
            /* check for saturation */
            {
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vH);
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vH);
            }""".strip() % params

        params["SATURATION_CHECK_FINAL"] = """
    if (%(VMOVEMASK)s(%(VOR)s(
            %(VCMPEQ)s(vSaturationCheckMin, vNegLimit),
            %(VCMPEQ)s(vSaturationCheckMax, vPosLimit)))) {
        result->saturated = 1;
        score = INT8_MAX;
    }""".strip() % params

        params["STATS_SATURATION_CHECK_INIT"] = """
    %(VTYPE)s vNegLimit = %(VSET1)s(INT8_MIN);
    %(VTYPE)s vPosLimit = %(VSET1)s(INT8_MAX);
    %(VTYPE)s vSaturationCheckMin = vPosLimit;
    %(VTYPE)s vSaturationCheckMax = vNegLimit;""".strip() % params
        if "diag" in params["NAME"]:
            params["STATS_SATURATION_CHECK_MID"] = """
            /* check for saturation */
            {
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWscore);
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vWscore);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWmatch);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWsimilar);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vWlength);
            }""".strip() % params
        else:
            params["STATS_SATURATION_CHECK_MID"] = """
            /* check for saturation */
            {
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vH);
                vSaturationCheckMin = %(VMIN)s(vSaturationCheckMin, vH);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHM);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHS);
                vSaturationCheckMax = %(VMAX)s(vSaturationCheckMax, vHL);
            }""".strip() % params

        params["STATS_SATURATION_CHECK_FINAL"] = """
    if (%(VMOVEMASK)s(%(VOR)s(
            %(VCMPEQ)s(vSaturationCheckMin, vNegLimit),
            %(VCMPEQ)s(vSaturationCheckMax, vPosLimit)))) {
        result->saturated = 1;
        score = INT8_MAX;
        matches = 0;
        similar = 0;
        length = 0;
    }""".strip() % params

        params["NEG_INF"] = "INT8_MIN"
        params["VADD"] = params["VADDSx8"]
        params["VSUB"] = params["VSUBSx8"]
        for p in ["VMAX", "VMIN"]:
            if (params[p].endswith("_rpl")
                    and params[p] not in params["FIXES"]):
                params["FIXES"] += params[params[p]]
    else:
        params["SATURATION_CHECK_INIT"] = ""
        params["SATURATION_CHECK_MID"] = ""
        params["SATURATION_CHECK_FINAL"] = ""
        params["STATS_SATURATION_CHECK_INIT"] = ""
        params["STATS_SATURATION_CHECK_MID"] = ""
        params["STATS_SATURATION_CHECK_FINAL"] = ""
    return params


def generated_params_diag(params):
    lanes = params["LANES"]
    params["DIAG_I"] = ",".join(["%d"%i for i in range(lanes)])
    params["DIAG_ILO"] = ",".join(["%d"%i for i in range(lanes/2,lanes)])
    params["DIAG_IHI"] = ",".join(["%d"%i for i in range(lanes/2)])
    params["DIAG_J"] = ",".join(["%d"%-i for i in range(lanes)])
    params["DIAG_JLO"] = ",".join(["%d"%-i for i in range(lanes/2,lanes)])
    params["DIAG_JHI"] = ",".join(["%d"%-i for i in range(lanes/2)])
    params["DIAG_IBoundary"] = "            ".join(
            ["-open-%d*gap,\n"%(i)
                for i in range(lanes)])[:-2]
    params["DIAG_VS1"] = "                ".join(
            ["s1[i+%d],\n"%i
                for i in range(lanes)])[:-2]
    params["DIAG_MATROW_DECL"] = "        ".join(
            ["const int * const restrict matrow%d = matrix[s1[i+%d]];\n"%(i,i)
                for i in range(lanes)])[:-1]
    params["DIAG_MATROW_USE"] = "                    ".join(
            ["matrow%d[s2[j-%d]],\n"%(i,i)
                for i in range(lanes)])[:-2]
    return params


def generated_params_striped(params):
    return params


def generated_params_scan(params):
    params["SCAN_INSERT_MASK"] = "1"+",0"*(params["LANES"]-1)
    return params


def generated_params(params):
    # some params are generated from given params
    bits = params["BITS"]
    width = params["WIDTH"]
    params["INDEX"] = "int32_t"
    params["ALIGNMENT"] = bits/8
    params["BYTES"] = width/8
    params["LANES"] = bits/width
    params["LAST_POS"] = params["LANES"]-1
    params["INT"] = "int%(WIDTH)s_t" % params
    params["NEG_INF"] = "(INT%(WIDTH)s_MIN/(%(INT)s)(2))" % params
    if "diag" in params["NAME"]:
        params = generated_params_diag(params)
    elif "striped" in params["NAME"]:
        params = generated_params_striped(params)
    elif "scan" in params["NAME"]:
        params = generated_params_scan(params)
    # select appropriate vector functions for given width
    suffix = "x%s" % width
    for key in keys:
        if key.endswith(suffix):
            new_key = key.split('x')[0]
            params[new_key] = params[key]
    fixes = ""
    for param in params:
        if param in template and str(params[param]).endswith("_rpl"):
            fixes += params[params[param]]
    params["FIXES"] = fixes
    params = generate_printer(params)
    params = generate_saturation_check(params)
    return params


for template_filename in template_filenames:
    template = open(template_dir+template_filename).read()
    for width in [64,32,16,8]:
        for isa in [sse2,sse41,avx2]:
            params = copy.deepcopy(isa)
            params["WIDTH"] = width
            prefix = template_filename[:-2]
            parts = prefix.split('_')
            table_prefix = ""
            if len(parts) == 2:
                table_prefix = "%s_table_%s" % (parts[0], parts[1])
            if len(parts) == 3:
                table_prefix = "%s_%s_table_%s" % (parts[0], parts[1], parts[2])
            function_name = "%s_%s%s_%s_%s" % (prefix,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_table_name = "%s_%s%s_%s_%s" % (table_prefix,
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            params["NAME"] = function_name
            params["NAME_TABLE"] = function_table_name
            params = generated_params(params)
            output_filename = "%s%s.c" % (output_dir, function_name)
            result = template % params
            writer = open(output_filename, "w")
            writer.write(template % params)
            writer.write("\n")
            writer.close()


# some templates have specializations for certain bit widths, e.g., 8
for template_filename in special_templates:
    template = open(template_dir+template_filename).read()
    prefix = template_filename[:-2]
    parts = prefix.split('_')
    width = int(parts[-1])
    parts = parts[:-1]
    prefix = "_".join(parts)
    table_prefix = ""
    if len(parts) == 2:
        table_prefix = "%s_table_%s" % (parts[0], parts[1])
    if len(parts) == 3:
        table_prefix = "%s_%s_table_%s" % (parts[0], parts[1], parts[2])
    for isa in [sse2,sse41,avx2]:
        params = copy.deepcopy(isa)
        params["WIDTH"] = width
        function_name = "%s_%s%s_%s_%s" % (prefix,
                isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
        function_table_name = "%s_%s%s_%s_%s" % (table_prefix,
                isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
        params["NAME"] = function_name
        params["NAME_TABLE"] = function_table_name
        params = generated_params(params)
        output_filename = "%s%s.c" % (output_dir, function_name)
        result = template % params
        writer = open(output_filename, "w")
        writer.write(template % params)
        writer.write("\n")
        writer.close()
