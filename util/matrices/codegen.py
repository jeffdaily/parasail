#!/usr/bin/python

# This file extracts the raw substitution matrix contents located within
# the same directory, e.g., BLOSUM62, PAM250, and constructs files
# appropriate for inclusion within a C program.

import glob
import os

header = """/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 *
 * This file was converted to C code from the raw file found at
 * ftp://ftp.cbi.pku.edu.cn/pub/software/blast/matrices/%(FILENAME)s, the
 * Center for Bioinformatics, Peking University, China.
 */
#ifndef _PARASAIL_%(FILENAME)s_H_
#define _PARASAIL_%(FILENAME)s_H_

#include "parasail.h"
#include "%(BASE)s_map.h"

"""

footer = "#endif /* _PARASAIL_%s_H_ */"

output_dir = "generated/"
if not os.path.exists(output_dir):
        os.makedirs(output_dir)

filenames = []
filenames.extend(glob.glob("BLOSUM*"))
filenames.extend(glob.glob("PAM*"))
filenames = sorted(filenames)

def get_base(name):
    base = ""
    if "BLOSUM" in name:
        base = "blosum"
    elif "PAM" in name:
        base = "pam"
    assert base
    return base

def generate_mapper(name, line):
    base = get_base(name)
    text = """/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#ifndef _PARASAIL_%s_MAP_H_
#define _PARASAIL_%s_MAP_H_

static const int parasail_%s_map[256] = {
""" % (name, name, name.lower())
    parts = line.split()
    limit = len(parts) - 1
    lookup = {}
    for index,part in enumerate(parts):
        lookup[part.upper()] = index
        lookup[part.lower()] = index
    text += "    "
    for i in range(256):
        if i % 16 == 0 and i > 0:
            text += "\n    "
        c = str(chr(i))
        if c in lookup:
            text += "%3s," % lookup[c]
        else:
            text += "%3d," % limit
    text += """
};

#endif /* _PARASAIL_%s_MAP_H_ */
""" % name
    return text

for filename in filenames:
    base = get_base(filename)
    filename_lower = filename.lower()
    output_filename = output_dir + filename_lower + ".h"
    writer = open(output_filename, "w")
    writer.write(header % {"FILENAME":filename, "BASE":base})
    the_lines = []
    for line in open(filename):
        line = line.strip()
        if line[0] == '#':
            writer.write("/* %s */\n" % line)
        else:
            the_lines.append(line)
    count = len(the_lines)-1
    # write out the mat[] form
    writer.write("\n")
    writer.write("static const int8_t parasail_%s_[] = {\n" % filename_lower)
    writer.write("/*     " + ("%4s"*24) % tuple(the_lines[0].split()) + " */\n")
    text = ""
    for line in the_lines[1:]:
        parts = line.split()
        text += "/* %s */ " % parts[0]
        text += ("%3s,"*24) % tuple(parts[1:])
        text += "\n"
    writer.write(text[:-2])
    writer.write("\n};\n")
    # write out the mat[][] form
    writer.write("\n")
    writer.write("static const int parasail_%s__[%d][%d] = {\n" % (
        filename_lower, count, count))
    writer.write("/*     " + ("%4s"*24) % tuple(the_lines[0].split()) + " */\n")
    text = ""
    for line in the_lines[1:]:
        parts = line.split()
        text += "/* %s */{" % parts[0]
        text += ((("%3s,"*24)[:-1]+"},") % tuple(parts[1:]))
        text += "\n"
    writer.write(text[:-2])
    writer.write("\n};\n")
    writer.write("""
static const parasail_matrix_t parasail_%s = {
    "%s",
    parasail_%s_,
    parasail_%s__,
    parasail_%s_map,
    %d
};

""" % (filename_lower, filename_lower, filename_lower, filename_lower,
    base, count))
    writer.write(footer % filename)
    writer.write("\n")
    writer.close()
    # write out the mapper file
    if "BLOSUM" in filename:
        output_filename = output_dir + "blosum_map.h"
        writer = open(output_filename, "w")
        writer.write(generate_mapper("BLOSUM", the_lines[0]))
        writer.close()
    elif "PAM" in filename:
        output_filename = output_dir + "pam_map.h"
        writer = open(output_filename, "w")
        writer.write(generate_mapper("PAM", the_lines[0]))
        writer.close()

# write out the matrix lookup function
text = """/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#ifndef _PARASAIL_MATRIX_LOOKUP_H_
#define _PARASAIL_MATRIX_LOOKUP_H_

%(HEADERS)s
#include "parasail/matrices/blosum_map.h"
#include "parasail/matrices/pam_map.h"

const parasail_matrix_t * parasail_matrices[] = {
%(MATRICES)s
    NULL
};

const parasail_matrix_t* parasail_matrix_lookup(const char *matrixname)
{
    const parasail_matrix_t *matrix = NULL;

    if (matrixname) {
        int index = 0;
        const parasail_matrix_t *current = parasail_matrices[index++];
        while (current) {
            if (0 == strcmp(matrixname, current->name)) {
                matrix = current;
                break;
            }
            current = parasail_matrices[index++];
        }
    }

    return matrix;
}

#endif /* _PARASAIL_MATRIX_LOOKUP_H_ */
"""

headers = ""
matrices = ""
for name in filenames:
    base = get_base(name)
    headers += '#include "parasail/matrices/%s.h"\n' % name.lower()
    matrices += '    &parasail_%s,\n' % name.lower()
output_filename = output_dir + "matrix_lookup.h"
writer = open(output_filename, "w")
writer.write(text % {"MATRICES":matrices[:-1], "HEADERS":headers[:-1]})
writer.write("\n")
writer.close()
