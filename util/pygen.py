#!/usr/bin/env python

# Creates the parasail.pyx file used for the python bindings.

import glob
import os

print """# cython: c_string_type=str, c_string_encoding=ascii

from libc.string cimport const_char

cdef extern from "parasail.h":
    ctypedef struct parasail_result_t:
        int saturated
        int score
        int matches
        int similar
        int length
        #int * score_table
        #int * matches_table
        #int * similar_table
        #int * length_table

    cdef struct parasail_matrix:
        const char * name
        #int * matrix
        #int * mapper
        int size
        int need_free
    ctypedef parasail_matrix parasail_matrix_t

    void parasail_result_free(parasail_result_t * result)

    parasail_matrix_t* parasail_matrix_create(
            const char * alphabet, int match, int mismatch)

    void parasail_matrix_free(parasail_matrix_t *matrix)
"""

# serial reference implementations (3x2x2x2 = 24 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table"]
for a in alg:
    for s in stats:
        for t in table:
            print ""
            print " "*4+"parasail_result_t* parasail_"+a+s+t+'('
            print " "*8+"const_char * s1, const int s1Len,"
            print " "*8+"const_char * s2, const int s2Len,"
            print " "*8+"const int open, const int gap,"
            print " "*8+"const parasail_matrix_t* matrix);"

# serial scan reference implementations (3x2x2x2 = 24 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table"]
for a in alg:
    for s in stats:
        for t in table:
            print ""
            print " "*4+"parasail_result_t* parasail_"+a+s+t+'_scan('
            print " "*8+"const_char * s1, const int s1Len,"
            print " "*8+"const_char * s2, const int s2Len,"
            print " "*8+"const int open, const int gap,"
            print " "*8+"const parasail_matrix_t* matrix);"

# vectorized implementations (3x2x2x3x13 = 468 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table"]
par = ["_scan", "_striped", "_diag"]
width = ["_64","_32","_16","_8"]

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

# vectorized implementations of blocked (just a couple)
alg = ["sw"]
stats = [""]
table = ["", "_table"]
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
cdef class Result:
    cdef parasail_result_t *_c_object
    def __dealloc__(self):
        if self._c_object is not NULL:
            parasail_result_free(self._c_object)
    @staticmethod
    cdef create(parasail_result_t *c_object):
        p = Result()
        p._c_object = c_object
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
                parasail_matrix_free(self._c_object)
    @staticmethod
    cdef create(const parasail_matrix_t *c_object):
        p = Matrix()
        p._c_object = c_object
        return p
    property name:
        def __get__(self): return self._c_object.name
    property size:
        def __get__(self): return self._c_object.size

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

# serial reference implementations (3x2x2x2 = 24 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table"]
for a in alg:
    for s in stats:
        for t in table:
            print ""
            print "def "+a+s+t+'('
            print " "*8+"const_char * s1,"
            print " "*8+"const_char * s2,"
            print " "*8+"const int open, const int gap,"
            print " "*8+"Matrix matrix not None):"
            print " "*4+"return Result.create(parasail_"+a+s+t+"(s1,len(s1),s2,len(s2),open,gap,matrix._c_object))"

