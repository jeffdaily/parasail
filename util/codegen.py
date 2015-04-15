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
template_filenames = [name for name in os.listdir(template_dir)
                      if name[-2:] == ".c"]

output_dir = "generated/"
if not os.path.exists(output_dir):
        os.makedirs(output_dir)


def generate_printer(params):
    text = ""
    if "striped" in params["NAME"]:
        for lane in range(params["LANES"]):
            new_params = params
            params["LANE"] = lane
            if params["LANES"] / 10:
                text += "    array[(%(LANE)2d*seglen+t)*dlen + d] = (%(INT)s)%(VEXTRACT)s(vH, %(LANE)2d);\n" % params
            else:
                text += "    array[(%(LANE)s*seglen+t)*dlen + d] = (%(INT)s)%(VEXTRACT)s(vH, %(LANE)s);\n" % params
        del params["LANE"]
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
        params["NEG_INF"] = "INT8_MIN"
        params["VADD"] = params["VADDSx8"]
        params["VSUB"] = params["VSUBSx8"]
    else:
        params["SATURATION_CHECK_INIT"] = ""
        params["SATURATION_CHECK_MID"] = ""
        params["SATURATION_CHECK_FINAL"] = ""
    return params


def generate_saturation_check(params):
    width = params["WIDTH"]
    if width == 8:
        params["SATURATION_CHECK_INIT"] = """
    %(VTYPE)s vNegLimit = %(VSET1)s(INT8_MIN);
    %(VTYPE)s vPosLimit = %(VSET1)s(INT8_MAX);
    %(VTYPE)s vSaturationCheckMin = vPosLimit;
    %(VTYPE)s vSaturationCheckMax = vNegLimit;""".strip() % params
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
            function_name = "%s_%s%s_%s_%s" % (template_filename[:-2],
                    isa["ISA"], isa["ISA_VERSION"], isa["BITS"], width)
            function_table_name = "%s_table_%s" % (
                    function_name[:2], function_name[3:])
            params["NAME"] = function_name
            params["NAME_TABLE"] = function_table_name
            params = generated_params(params)
            output_filename = "%s%s.c" % (output_dir, function_name)
            result = template % params
            writer = open(output_filename, "w")
            writer.write(template % params)
            writer.write("\n")
            writer.close()
