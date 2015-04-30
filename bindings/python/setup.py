#!/usr/bin/python
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import os

PARASAIL_PREFIX = "PARASAIL_PREFIX"
PARASAIL_INCLUDE = "PARASAIL_INCLUDE"
PARASAIL_LIB = "PARASAIL_LIB"

PREFIX = None
INCLUDE = None
LIB = None

if PARASAIL_PREFIX in os.environ:
    PREFIX = os.environ[PARASAIL_PREFIX]
    INCLUDE = PREFIX + "/include"
    LIB = PREFIX + "/lib"
else:
    INCLUDE = os.environ.get(PARASAIL_INCLUDE, "/usr/local/include")
    LIB = os.environ.get(PARASAIL_LIB, "/usr/local/lib")

ext_modules=[
    Extension("parasail",
        sources=["parasail.pyx"],
        libraries=["parasail"],
        include_dirs=[INCLUDE],
        library_dirs=[LIB],
    )
]

setup(
    ext_modules = cythonize(ext_modules)
)
