#!/usr/bin/python

# Generate the various memory implementations.

import copy
import os

sse = {
        "ISA_HDRNAME" : "sse",
        "ISA_SRCNAME" : "sse",
        "HEADERS": "#include <emmintrin.h>",
        "VTYPE" : "__m128i",
        "VWIDTH" : "128",
        "VSTORE" : "_mm_store_si128",
        "ALIGN" : "16",
        "DECL" : (lambda bit : "__m128i_%s_t"%bit),
        "ACCESS" : (lambda bit : ".v"),
        "ASSIGN" : ".m"
}

avx2 = {
        "ISA_HDRNAME" : "avx",
        "ISA_SRCNAME" : "avx2",
        "HEADERS": "#include <immintrin.h>",
        "VTYPE" : "__m256i",
        "VWIDTH" : "256",
        "VSTORE" : "_mm256_store_si256",
        "ALIGN" : "32",
        "DECL" : (lambda bit : "__m256i_%s_t"%bit),
        "ACCESS" : (lambda bit : ".v"),
        "ASSIGN" : ".m"
}

altivec = {
        "ISA_HDRNAME" : "altivec",
        "ISA_SRCNAME" : "altivec",
        "HEADERS": "#include <stdint.h>\n#include <stdlib.h>",
        "VTYPE" : "vec128i",
        "VWIDTH" : "128",
        "VSTORE" : "_mm_store_si128",
        "ALIGN" : "16",
        "DECL" : (lambda bit : "vec128i_%s_t"%bit),
        "ACCESS" : (lambda bit : ".v"),
        "ASSIGN" : ".m"
}

neon = {
        "ISA_HDRNAME" : "neon",
        "ISA_SRCNAME" : "neon",
        "HEADERS": "#include <stdint.h>\n#include <stdlib.h>",
        "VTYPE" : "simde__m128i",
        "VWIDTH" : "128",
        "VSTORE" : "simde_mm_store_si128",
        "ALIGN" : "16",
        "DECL" : (lambda bit : "simde__m128i"),
        "ACCESS" : (lambda bit : ".i%s"%bit),
        "ASSIGN" : ""
}

def codegen(isa):
    txt = """/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

%(HEADERS)s

#include "parasail/memory.h"
#include "parasail/internal_%(ISA_HDRNAME)s.h"

%(VTYPE)s * parasail_memalign_%(VTYPE)s(size_t alignment, size_t size)
{
    return (%(VTYPE)s *) parasail_memalign(alignment, size*sizeof(%(VTYPE)s));
}

void parasail_memset_%(VTYPE)s(%(VTYPE)s *b, %(VTYPE)s c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        %(VSTORE)s(&b[i], c);
    }
}
""" % isa

    width = int(isa["VWIDTH"])
    for bitsize in [8,16,32,64]:
        params = copy.deepcopy(isa)
        params["BITSIZE"] = bitsize
        params["SEGWIDTH"] = width/bitsize
        params["DYNTYPE"] = isa["DECL"](bitsize)
        params["DYNACCESS"] = isa["ACCESS"](bitsize)
        txt += """
parasail_profile_t * parasail_profile_create_%(ISA_HDRNAME)s_%(VWIDTH)s_%(BITSIZE)d(
        const char * const restrict s1, const int _s1Len,
        const parasail_matrix_t *matrix)
{
    int s1Len = 0;
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    int32_t n = 0; /* number of amino acids in table */
    const int32_t segWidth = %(SEGWIDTH)d; /* number of values in vector unit */
    int32_t segLen = 0;
    %(VTYPE)s* restrict vProfile = NULL;
    int32_t index = 0;
    parasail_profile_t *profile = NULL;

    PARASAIL_CHECK_NULL(matrix);
    /* s1 is ignored for pssm, required for square */
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        PARASAIL_CHECK_NULL(s1);
    }

    s1Len = matrix->type == PARASAIL_MATRIX_TYPE_SQUARE ? _s1Len : matrix->length;
    n = matrix->size;
    segLen = (s1Len + segWidth - 1) / segWidth;
    vProfile = parasail_memalign_%(VTYPE)s(%(ALIGN)s, n * segLen);
    if (!vProfile) return NULL;
    profile = parasail_profile_new(s1, s1Len, matrix);
    if (!profile) return NULL;

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            %(DYNTYPE)s t;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
                    t%(DYNACCESS)s[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                }
                else {
                    t%(DYNACCESS)s[segNum] = j >= s1Len ? 0 : matrix->matrix[n*j+matrix->mapper[(unsigned char)matrix->alphabet[k]]];
                }
                j += segLen;
            }
            %(VSTORE)s(&vProfile[index], t%(ASSIGN)s);
            ++index;
        }
    }

    profile->profile%(BITSIZE)d.score = vProfile;
    profile->free = &parasail_free_%(VTYPE)s;
    return profile;
}
""" % params

    txt += """
parasail_profile_t* parasail_profile_create_%(ISA_HDRNAME)s_%(VWIDTH)s_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile8 = parasail_profile_create_%(ISA_HDRNAME)s_%(VWIDTH)s_8(s1, s1Len, matrix);
    parasail_profile_t *profile16 = parasail_profile_create_%(ISA_HDRNAME)s_%(VWIDTH)s_16(s1, s1Len, matrix);
    parasail_profile_t *profile32 = parasail_profile_create_%(ISA_HDRNAME)s_%(VWIDTH)s_32(s1, s1Len, matrix);
    profile8->profile16 = profile16->profile16;
    profile8->profile32 = profile32->profile32;
    free(profile16);
    free(profile32);

    return profile8;
}
""" % isa

    for bitsize in [8,16,32,64]:
        params = copy.deepcopy(isa)
        params["BITSIZE"] = bitsize
        params["SEGWIDTH"] = width/bitsize
        params["DYNTYPE"] = isa["DECL"](bitsize)
        params["DYNACCESS"] = isa["ACCESS"](bitsize)
        txt += """
parasail_profile_t * parasail_profile_create_stats_%(ISA_HDRNAME)s_%(VWIDTH)s_%(BITSIZE)d(
        const char * const restrict s1, const int _s1Len,
        const parasail_matrix_t *matrix)
{
    int s1Len = 0;
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t segNum = 0;
    int32_t n = 0;
    const int32_t segWidth = %(SEGWIDTH)d; /* number of values in vector unit */
    int32_t segLen = 0;
    %(VTYPE)s* restrict vProfile = NULL;
    %(VTYPE)s* restrict vProfileM = NULL;
    %(VTYPE)s* restrict vProfileS = NULL;
    int32_t index = 0;
    parasail_profile_t *profile = NULL;

    PARASAIL_CHECK_NULL(matrix);
    /* s1 is required for both pssm and square */
    PARASAIL_CHECK_NULL(s1);

    s1Len = matrix->type == PARASAIL_MATRIX_TYPE_SQUARE ? _s1Len : matrix->length;
    n = matrix->size; /* number of amino acids in table */
    segLen = (s1Len + segWidth - 1) / segWidth;
    vProfile = parasail_memalign_%(VTYPE)s(%(ALIGN)s, n * segLen);
    if (!vProfile) return NULL;
    vProfileM = parasail_memalign_%(VTYPE)s(%(ALIGN)s, n * segLen);
    if (!vProfileM) return NULL;
    vProfileS = parasail_memalign_%(VTYPE)s(%(ALIGN)s, n * segLen);
    if (!vProfileS) return NULL;
    profile = parasail_profile_new(s1, s1Len, matrix);
    if (!profile) return NULL;

    for (k=0; k<n; ++k) {
        for (i=0; i<segLen; ++i) {
            %(DYNTYPE)s p;
            %(DYNTYPE)s m;
            %(DYNTYPE)s s;
            j = i;
            for (segNum=0; segNum<segWidth; ++segNum) {
                if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
                    p%(DYNACCESS)s[segNum] = j >= s1Len ? 0 : matrix->matrix[n*k+matrix->mapper[(unsigned char)s1[j]]];
                }
                else {
                    p%(DYNACCESS)s[segNum] = j >= s1Len ? 0 : matrix->matrix[n*j+matrix->mapper[(unsigned char)matrix->alphabet[k]]];
                }
                m%(DYNACCESS)s[segNum] = j >= s1Len ? 0 : (k == matrix->mapper[(unsigned char)s1[j]]);
                s%(DYNACCESS)s[segNum] = p%(DYNACCESS)s[segNum] > 0;
                j += segLen;
            }
            %(VSTORE)s(&vProfile[index], p%(ASSIGN)s);
            %(VSTORE)s(&vProfileM[index], m%(ASSIGN)s);
            %(VSTORE)s(&vProfileS[index], s%(ASSIGN)s);
            ++index;
        }
    }

    profile->profile%(BITSIZE)d.score = vProfile;
    profile->profile%(BITSIZE)d.matches = vProfileM;
    profile->profile%(BITSIZE)d.similar = vProfileS;
    profile->free = &parasail_free_%(VTYPE)s;
    return profile;
}
""" % params

    txt += """
parasail_profile_t* parasail_profile_create_stats_%(ISA_HDRNAME)s_%(VWIDTH)s_sat(
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile8 = parasail_profile_create_stats_%(ISA_HDRNAME)s_%(VWIDTH)s_8(s1, s1Len, matrix);
    parasail_profile_t *profile16 = parasail_profile_create_stats_%(ISA_HDRNAME)s_%(VWIDTH)s_16(s1, s1Len, matrix);
    parasail_profile_t *profile32 = parasail_profile_create_stats_%(ISA_HDRNAME)s_%(VWIDTH)s_32(s1, s1Len, matrix);
    profile8->profile16 = profile16->profile16;
    profile8->profile32 = profile32->profile32;
    free(profile16);
    free(profile32);

    return profile8;
}

void parasail_free_%(VTYPE)s(void *ptr)
{
    parasail_free((%(VTYPE)s*)ptr);
}
""" % isa
    return txt

output_dir = "generated/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for isa in [sse,avx2,altivec,neon]:
    output_filename = "%smemory_%s.c" % (output_dir, isa["ISA_SRCNAME"])
    writer = open(output_filename, "w")
    writer.write(codegen(isa))
    writer.write("\n")
    writer.close()
