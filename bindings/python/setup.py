#!/usr/bin/python
from distutils.core import setup
from distutils.extension import Extension
from distutils.command.build_ext import build_ext
from distutils.command.config import config
from Cython.Build import cythonize

import os

# numpy is required -- attempt import
try:
    import numpy
except ImportError:
    print "numpy is required"
    raise

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

class my_build_ext(build_ext):
    def run(self):
        self.run_command("config")
        config_cmd = self.get_finalized_command("config")
        if self.define is None:
            self.define = [('restrict', config_cmd.restrict)]
        else:
            self.define += [('restrict', config_cmd.restrict)]
        build_ext.run(self)

# modelled after numpy inline test
class my_config(config):
    def initialize_options(self):
        config.initialize_options(self)
        self.restrict = None

    def finalize_options(self):
        config.finalize_options(self)
        if self.restrict is None:
            self.restrict = ""

    def check_restrict(self):
        """Return the restrict identifier (may be empty)."""
        self._check_compiler()
        body = """
#include <string.h>
void* nostatic_func (void * const %(restrict)s dst, const void * const %(restrict)s src)
{
    return memcpy(dst, src, 1);
}
"""

        for kw in ['restrict', '__restrict', '__restrict__', '_Restrict']:
            st = self.try_compile(body % {'restrict': kw}, None, None)
            if st:
                return kw

        return ''

    def run(self):
        config.run(self)
        self.restrict = self.check_restrict()

ext_modules=[
        Extension("parasail",
            sources=["parasail.pyx"],
            libraries=["parasail"],
            include_dirs=[INCLUDE, numpy.get_include()],
            library_dirs=[LIB]
            )
        ]

setup(
        name = "parasail",
        description = "Pairwise Sequence Alignment Library",
        author = "Jeff Daily",
        author_email = "jeff.daily@pnnl.gov",
        url = "https://github.com/jeffdaily/parasail",
        cmdclass={'config': my_config, "build_ext": my_build_ext},
        ext_modules = cythonize(ext_modules)
        )
