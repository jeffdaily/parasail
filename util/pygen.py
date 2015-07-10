#!/usr/bin/env python

# Creates the parasail.pyx file used for the python bindings.

import glob
import os

print """# cython: c_string_type=str, c_string_encoding=ascii

from libc.string cimport const_char
import numpy as np
cimport numpy as np
from cpython cimport PyObject, Py_INCREF

np.import_array()

cdef extern from "parasail.h":
    ctypedef struct parasail_result_t:
        int saturated
        int score
        int matches
        int similar
        int length
        int * score_table
        int * matches_table
        int * similar_table
        int * length_table
        int * score_row
        int * matches_row
        int * similar_row
        int * length_row
        int * score_col
        int * matches_col
        int * similar_col
        int * length_col

    ctypedef struct parasail_matrix_t:
        const char * name
        int * matrix
        #int * mapper
        int size
        int need_free

    ctypedef struct parasail_profile_t:
        const char *s1
        int s1Len
        const parasail_matrix_t *matrix
        #parasail_profile_data_t * profile8;
        #parasail_profile_data_t * profile16;
        #parasail_profile_data_t * profile32;
        #parasail_profile_data_t * profile64;
        #void (*free)(void * profile)

    void parasail_result_free(parasail_result_t * result)

    parasail_matrix_t* parasail_matrix_create(
            const char * alphabet, int match, int mismatch)

    void parasail_matrix_free(parasail_matrix_t *matrix)

    parasail_profile_t * parasail_profile_create_8(
            const char * s1, const int s1Len,
            const parasail_matrix_t *matrix);

    parasail_profile_t * parasail_profile_create_16(
            const char * s1, const int s1Len,
            const parasail_matrix_t *matrix);

    parasail_profile_t * parasail_profile_create_32(
            const char * s1, const int s1Len,
            const parasail_matrix_t *matrix);

    parasail_profile_t * parasail_profile_create_64(
            const char * s1, const int s1Len,
            const parasail_matrix_t *matrix);

    parasail_profile_t * parasail_profile_create_sat(
            const char * s1, const int s1Len,
            const parasail_matrix_t *matrix);

    parasail_profile_t * parasail_profile_create_stats_8(
            const char * s1, const int s1Len,
            const parasail_matrix_t *matrix);

    parasail_profile_t * parasail_profile_create_stats_16(
            const char * s1, const int s1Len,
            const parasail_matrix_t *matrix);

    parasail_profile_t * parasail_profile_create_stats_32(
            const char * s1, const int s1Len,
            const parasail_matrix_t *matrix);

    parasail_profile_t * parasail_profile_create_stats_64(
            const char * s1, const int s1Len,
            const parasail_matrix_t *matrix);

    parasail_profile_t * parasail_profile_create_stats_sat(
            const char * s1, const int s1Len,
            const parasail_matrix_t *matrix);

    void parasail_profile_free(parasail_profile_t *profile)
"""

# serial reference implementations (3x2x3 = 18 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol"]
for a in alg:
    for s in stats:
        for t in table:
            print ""
            print " "*4+"parasail_result_t* parasail_"+a+s+t+'('
            print " "*8+"const_char * s1, const int s1Len,"
            print " "*8+"const_char * s2, const int s2Len,"
            print " "*8+"const int open, const int gap,"
            print " "*8+"const parasail_matrix_t* matrix);"

# serial scan reference implementations (3x2x3 = 18 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol"]
for a in alg:
    for s in stats:
        for t in table:
            print ""
            print " "*4+"parasail_result_t* parasail_"+a+s+t+'_scan('
            print " "*8+"const_char * s1, const int s1Len,"
            print " "*8+"const_char * s2, const int s2Len,"
            print " "*8+"const int open, const int gap,"
            print " "*8+"const parasail_matrix_t* matrix);"

# vectorized implementations (3x2x3x3x4 = 216 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol"]
par = ["_scan", "_striped", "_diag"]
width = ["_64","_32","_16","_8","_sat"]
for a in alg:
    for s in stats:
        for t in table:
            for p in par:
                for w in width:
                    print ""
                    print " "*4+"parasail_result_t* parasail_"+a+s+t+p+w+'('
                    print " "*8+"const_char * s1, const int s1Len,"
                    print " "*8+"const_char * s2, const int s2Len,"
                    print " "*8+"const int open, const int gap,"
                    print " "*8+"const parasail_matrix_t* matrix);"

# vectorized profile implementations (3x2x3x2x4 = 144 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol"]
par = ["_scan_profile", "_striped_profile"]
width = ["_64","_32","_16","_8","_sat"]
for a in alg:
    for s in stats:
        for t in table:
            for p in par:
                for w in width:
                    print ""
                    print " "*4+"parasail_result_t* parasail_"+a+s+t+p+w+'('
                    print " "*8+"const parasail_profile_t* profile,"
                    print " "*8+"const_char * s2, const int s2Len,"
                    print " "*8+"const int open, const int gap);"

# vectorized implementations of blocked (just a couple)
alg = ["sw"]
stats = [""]
table = ["", "_table", "_rowcol"]
par = ["_blocked"]
width = ["_32","_16"]
for a in alg:
    for s in stats:
        for t in table:
            for p in par:
                for w in width:
                    print ""
                    print " "*4+"parasail_result_t* parasail_"+a+s+t+p+w+'('
                    print " "*8+"const_char * s1, const int s1Len,"
                    print " "*8+"const_char * s2, const int s2Len,"
                    print " "*8+"const int open, const int gap,"
                    print " "*8+"const parasail_matrix_t* matrix);"


# now for the actual implementations and not just the header info

print """

cdef np.ndarray to_ndarray(object base, int * data, np.npy_intp length):
    cdef np.ndarray[np.int32_t, ndim=1] ndarray
    ndarray = np.PyArray_SimpleNewFromData(1, &length, np.NPY_INT32, data)
    Py_INCREF(base)
    ndarray.base = <PyObject*>base
    return ndarray

cdef class Result:
    cdef parasail_result_t *_c_object
    cdef int s1Len
    cdef int s2Len
    def __dealloc__(self):
        if self._c_object is not NULL:
            parasail_result_free(self._c_object)
    @staticmethod
    cdef create(parasail_result_t *c_object, s1Len, s2Len):
        p = Result()
        p._c_object = c_object
        p.s1Len = s1Len
        p.s2Len = s2Len
        return p
    property saturated:
        def __get__(self): return self._c_object.saturated
    property score:
        def __get__(self): return self._c_object.score
    property matches:
        def __get__(self): return self._c_object.matches
    property similar:
        def __get__(self): return self._c_object.similar
    property length:
        def __get__(self): return self._c_object.length
    cdef get_table(self, int * table):
        if not table:
            return None
        else:
            array = to_ndarray(self, table, self.s1Len * self.s2Len)
            return array.reshape((self.s1Len, self.s2Len))
    property score_table:
        def __get__(self):
            return self.get_table(self._c_object.score_table)
    property matches_table:
        def __get__(self):
            return self.get_table(self._c_object.matches_table)
    property similar_table:
        def __get__(self):
            return self.get_table(self._c_object.similar_table)
    property length_table:
        def __get__(self):
            return self.get_table(self._c_object.length_table)
    cdef get_row(self, int * row):
        if not row:
            return None
        else:
            return to_ndarray(self, row, self.s2Len)
    property score_row:
        def __get__(self):
            return self.get_row(self._c_object.score_row)
    property matches_row:
        def __get__(self):
            return self.get_row(self._c_object.matches_row)
    property similar_row:
        def __get__(self):
            return self.get_row(self._c_object.similar_row)
    property length_row:
        def __get__(self):
            return self.get_row(self._c_object.length_row)
    cdef get_col(self, int * col):
        if not col:
            return None
        else:
            return to_ndarray(self, col, self.s1Len)
    property score_col:
        def __get__(self):
            return self.get_col(self._c_object.score_col)
    property matches_col:
        def __get__(self):
            return self.get_col(self._c_object.matches_col)
    property similar_col:
        def __get__(self):
            return self.get_col(self._c_object.similar_col)
    property length_col:
        def __get__(self):
            return self.get_col(self._c_object.length_col)
    def __str__(self):
        return str(self.score)
    def __int__(self):
        return self.score
    def __long__(self):
        return self.score
"""

print """
cdef class Matrix:
    cdef const parasail_matrix_t *_c_object
    def __dealloc__(self):
        if self._c_object is not NULL:
            if self._c_object.need_free:
                parasail_matrix_free(<parasail_matrix_t*>self._c_object)
    @staticmethod
    cdef create(const parasail_matrix_t *c_object):
        p = Matrix()
        p._c_object = c_object
        return p
    property name:
        def __get__(self): return self._c_object.name
    property size:
        def __get__(self): return self._c_object.size
    property matrix:
        def __get__(self):
            size = self.size
            array = to_ndarray(self, self._c_object.matrix, size*size)
            return array.reshape((size, size))
    def __str__(self):
        return str(self.matrix)

def matrix_create(const char* alphabet, int match, int mismatch):
    cdef parasail_matrix_t* matrix = parasail_matrix_create(
            alphabet, match, mismatch)
    return Matrix.create(matrix)

"""

prefix = "matrices/"
filenames = []
filenames.extend(glob.glob(prefix+"BLOSUM*"))
filenames.extend(glob.glob(prefix+"PAM*"))
names = [f[len(prefix):].lower() for f in filenames]

for name in names:
    print 'cdef extern from "parasail/matrices/%s.h":' % name
    print ' '*4+"const parasail_matrix_t parasail_"+name
for name in names:
    print "%s = Matrix.create(&parasail_%s)" % (name,name)

print """
cdef class Profile:
    cdef parasail_profile_t *_c_object
    cdef Matrix matrix
    def __init__(self, const_char * s1, Matrix matrix not None, bits):
        if bits not in [8,16,32,64,"sat"]:
            raise ValueError('bits must be one of [8,16,32,64,"sat"]')
        self.create(s1, matrix, bits)
    def __dealloc__(self):
        if self._c_object is not NULL:
            parasail_profile_free(self._c_object)
    cdef create(self, const_char * s1, Matrix matrix, bits):
        cdef size_t t1 = len(s1)
        cdef int l1 = <int>t1
        if 8 == bits:
            self._c_object = parasail_profile_create_8(s1, l1, matrix._c_object)
        elif 16 == bits:
            self._c_object = parasail_profile_create_16(s1, l1, matrix._c_object)
        elif 32 == bits:
            self._c_object = parasail_profile_create_32(s1, l1, matrix._c_object)
        elif 64 == bits:
            self._c_object = parasail_profile_create_64(s1, l1, matrix._c_object)
        elif "sat" == bits:
            self._c_object = parasail_profile_create_sat(s1, l1, matrix._c_object)
        self.matrix = matrix
    property s1:
        def __get__(self): return self._c_object.s1
    property matrix:
        def __get__(self): return self.matrix

"""

# serial reference implementations (3x2x3x2 = 36 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol"]
scan = ["", "_scan"]
for a in alg:
    for s in stats:
        for t in table:
            for x in scan:
                prefix = a+s+t+x
                print ""
                print "def "+prefix+'('
                print " "*8+"const_char * s1,"
                print " "*8+"const_char * s2,"
                print " "*8+"const int open, const int gap,"
                print " "*8+"Matrix matrix not None):"
                print " "*4+"cdef size_t t1 = len(s1)"
                print " "*4+"cdef size_t t2 = len(s2)"
                print " "*4+"cdef int l1 = <int>t1"
                print " "*4+"cdef int l2 = <int>t2"
                print " "*4+"assert l1 == t1"
                print " "*4+"assert l2 == t2"
                print " "*4+"return Result.create(parasail_"+prefix+"(s1,l1,s2,l2,open,gap,matrix._c_object), l1, l2)"

# vectorized implementations (3x2x3x3x4 = 216 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol"]
par = ["_scan", "_striped", "_diag"]
width = ["_64","_32","_16","_8","_sat"]
for a in alg:
    for s in stats:
        for t in table:
            for p in par:
                for w in width:
                    prefix = a+s+t+p+w
                    print ""
                    print "def "+prefix+'('
                    print " "*8+"const_char * s1,"
                    print " "*8+"const_char * s2,"
                    print " "*8+"const int open, const int gap,"
                    print " "*8+"Matrix matrix not None):"
                    print " "*4+"cdef size_t t1 = len(s1)"
                    print " "*4+"cdef size_t t2 = len(s2)"
                    print " "*4+"cdef int l1 = <int>t1"
                    print " "*4+"cdef int l2 = <int>t2"
                    print " "*4+"assert l1 == t1"
                    print " "*4+"assert l2 == t2"
                    print " "*4+"return Result.create(parasail_"+prefix+"(s1,l1,s2,l2,open,gap,matrix._c_object), l1, l2)"

# vectorized profile implementations (3x2x3x3x4 = 216 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol"]
par = ["_scan_profile", "_striped_profile"]
width = ["_64","_32","_16","_8","_sat"]
for a in alg:
    for s in stats:
        for t in table:
            for p in par:
                for w in width:
                    prefix = a+s+t+p+w
                    print ""
                    print "def "+prefix+'('
                    print " "*8+"Profile profile not None,"
                    print " "*8+"const_char * s2,"
                    print " "*8+"const int open, const int gap):"
                    print " "*4+"cdef size_t t2 = len(s2)"
                    print " "*4+"cdef int l1 = profile._c_object.s1Len"
                    print " "*4+"cdef int l2 = <int>t2"
                    print " "*4+"assert l2 == t2"
                    print " "*4+"return Result.create(parasail_"+prefix+"(profile._c_object,s2,l2,open,gap), l1, l2)"

