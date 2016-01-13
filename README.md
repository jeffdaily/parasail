# parasail: Pairwise Sequence Alignment Library

Author: Jeff Daily (jeff.daily@pnnl.gov)

## Table of Contents

  * [Introduction](#introduction)
  * [A Note About Instruction Sets and CPU Dispatching](#a-note-about-instruction-sets-and-cpu-dispatching)
  * [Installing](#installing)
  * [C Interface Example](#c-interface-example)
    * [Standard Function Naming Convention](#standard-function-naming-convention)
    * [Profile Function Naming Convention](#profile-function-naming-convention)
    * [Substitution Matrices](#substitution-matrices)
    * [Function Lookup](#function-lookup)
  * [Language Bindings](#language-bindings)
    * [C/C++](#cc)
    * [Python](#python)
  * [Example Applications](#example-applications)
  * [Windows Build](#windows-build)
  * [Citing parasail](#citing-parasail)
  * [License: Battelle BSD-style](#license-battelle-bsd-style)

## Introduction

[back to top](#parasail-pairwise-sequence-alignment-library)

parasail is a SIMD C (C99) library containing implementations of the Smith-Waterman (local), Needleman-Wunsch (global), and semi-global pairwise sequence alignment algorithms.  Here, semi-global means insertions before the start or after the end of either the query or target sequence are not penalized.  parasail implements most known algorithms for vectorized pairwise sequence alignment, including diagonal [Wozniak, 1997], blocked [Rognes and Seeberg, 2000], striped [Farrar, 2007], and prefix scan [Daily, 2015].  Therefore, parasail is a reference implementation for these algorithms in addition to providing an implementation of the best-performing algorithm(s) to date on today's most advanced CPUs.

parasail implements the above algorithms currently in two variants, 1) returning the alignment score and ending locations, and 2) additionally returning alignment statistics (number of exact matches, number of similarities, and alignment length).  The two variants exist because parasail is intended to be high-performing; caculating additional statistics is not free. Select the appropriate implementation for your needs.

Note: When any of the algorithms open a gap, only the gap open penalty alone is applied.

## A Note About Instruction Sets and CPU Dispatching

[back to top](#parasail-pairwise-sequence-alignment-library)

parasail supports the SSE2, SSE4.1, AVX2, and KNC (Xeon Phi) instruction sets.  In many cases, your compiler can compile source code for an instruction set which is not supported by your host CPU.  The code is still compiled, however, parasail uses CPU dispatching at runtime to correctly select the appropriate implementation for the highest level of instruction set supported.  This allows parasail to be compiled and distributed by a maintainer for the best available system while still allowing the distribution to run with a lesser CPU.

Rarely, configure will report that your compiler supports an instruction set, e.g., AVX2, when in fact it does not.  This is a bug; please file an issue and I will amend the configure script, if needed.  You likely won't discover this until you attempt to compile all of the sources.  If you run into issues (or if for some other reason you wish to disable one of the instruction sets), and after reporting the issue, you can disable the offending code as a temporary work-around using one or more of the following configure flags.
```bash
configure SSE2_FLAGS=choke SSE41_FLAGS=choke AVX2_CFLAGS=choke
```

## Installing

[back to top](#parasail-pairwise-sequence-alignment-library)

parasail follows the typical configure, make, make install steps of other GNU autotools-based installations.  There are no external dependencies.  There is an optional CMake build available, however the GNU autotools-based installation is the preferred method.

By default, running "make install" will install parasail into /usr/local. You will find the parasail.h header in /usr/local/include and the parasail library, e.g., libparasail.so, in /usr/local/lib. If you specify a different prefix during configure, for example `configure --prefix=/some/other/path`, then look within the include and lib directories there for the parasail.h header and
libparasail.so library, respectively.

Don't forget to link your application to the parasail library.  For example, `gcc foo.c -I/where/you/installed/include -L/where/you/installed/lib -lparasail`.  Otherwise, you'll see errors such as `undefined reference to 'parasail_sw'`.

## C Interface Example

[back to top](#parasail-pairwise-sequence-alignment-library)

All parasail functions have identical function signatures with respect to return types and parameters.

```C
parasail_result_t* the_parasail_function_name(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const parasail_matrix_t* matrix);
```

With respect to traditional database query use, s1 is the query sequence and s2 is the database sequence.  For the functions returning the DP table or last row and column, the dimensions of the DP table are s1Len x s2Len where s1Len is the number of rows and s2Len is the number of columns (in C row-major order).

The return type is a C struct and will be populated based on the function called; no single function will fill all members.

```C
typedef struct parasail_result {
    int saturated;  /* for the 8-bit functions, whether score overflowed and should be discarded */
    int score;      /* alignment score */
    int matches;    /* number of exactly matching characters in the alignment */
    int similar;    /* number of similar characters (positive substitutions) in the alignment */
    int length;     /* length of the alignment */
    int end_query;  /* end position of query sequence */
    int end_ref;    /* end position of reference sequence */
    int * restrict score_table;     /* DP table of scores */
    int * restrict matches_table;   /* DP table of exact match counts */
    int * restrict similar_table;   /* DP table of similar substitution counts */
    int * restrict length_table;    /* DP table of lengths */
    int * restrict score_row;       /* last row of DP table of scores */
    int * restrict matches_row;     /* last row of DP table of exact match counts */
    int * restrict similar_row;     /* last row of DP table of similar substitution counts */
    int * restrict length_row;      /* last row of DP table of lengths */
    int * restrict score_col;       /* last col of DP table of scores */
    int * restrict matches_col;     /* last col of DP table of exact match counts */
    int * restrict similar_col;     /* last col of DP table of similar substitution counts */
    int * restrict length_col;      /* last col of DP table of lengths */
} parasail_result_t;
```

You must free the returned parasail result using `void parasail_result_free(parasail_result_t *result)`.

### Standard Function Naming Convention

[back to top](#parasail-pairwise-sequence-alignment-library)

There are over 1,000 functions within the parasail library.  To make it easier to find the function you're looking for, the function names follow a naming convention.  The following will use set notation {} to indicate a selection must be made and brackets [] to indicate an optional part of the name.

`parasail_ {nw,sg,sw}_ [stats_] [{table,rowcol}_] {striped,scan,diag,blocked}_ [{sse2_128,sse4_128,avx2_256,knc_512}_] {8,16,32,64}`

Here is a breakdown of each section of the name:
  1. parasail_ -- prefix a.k.a. namespace
  2. {nw,sg,sw}_ -- the class of alignment; global, semi-global, or local, respectively
  3. [stats_] -- optionally if the additional statistics are requested
  4. [{table,rowcol}_] -- optionally if the DP table or last row and column of DP table should be returned
  5. {striped,scan,diag,blocked} -- the vectorized approach; striped is always a good choice
  6. [{sse2_128,sse4_128,avx2_256,knc_512}_] -- optionally the instruction set and vector width
  7. {8,16,32,64,sat} -- the integer width of the solution, a.k.a. the vector element widths; knc only supports _32; 16 is often a good choice; 'sat' is short for 'saturation check' -- the 8-bit solution is attempted and if the score overflows (saturates), the 16-bit solution is then attempted. In some cases this is faster than simply running the 16-bit solution.

For example:

- `parasail_nw_stats_striped_sse41_128_16` would use Needleman-Wunsch, with alignment statistics, using striped vectors for sse41 16-bit integers.
- `parasail_sg` would use semi-global, no alignment statistics, no vectorized code (i.e. serial).
- `parasail_sw_scan_avx2_256_8` would use Smith-Waterman, no alignment statistics, using prefix scan vectors for avx2 8-bit integers.
- `parasail_nw_stats_striped_16` would use Needleman-Wunsch, with alignment statistics, using striped vectors, dispatching to the best CPU, and 16-bit integers.
- `parasail_sw_scan_8` would use Smith-Waterman, no alignment statistics, using prefix scan vectors, dispatching to the best CPU, for 8-bit integers.
- `parasail_sg_rowcol_striped_16` would use semi-global, no alignment statistics, output the last row and column of the DP table, using striped vectors, dispatching to the best CPU, for 16-bit integers.

Note: The blocked vector implementation only exist for sse41 16-bit and 32-bit integer elements.

Note: The dispatcher for the KNC instruction set will always dispatch to the 32-bit integer element implementation since it is the only one supported on that platform.

### Profile Function Naming Convention

[back to top](#parasail-pairwise-sequence-alignment-library)

There is a special subset of functions that mimic the behavior of the [SSW library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library). For the striped and scan vector implementations *only*, a query profile can be created and reused for subsequent alignments. This can speed up applications such as database search noticeably.  To create a profile, use

```C
parasail_profile_t* parasail_profile_create[_stats][_{sse_128,avx_256,knc_512}]_{8,16,32,64,sat} (
        const char * const restrict s1, const int s1Len,
        const parasail_matrix_t* matrix);
```

Similar to the standard function naming convention, there are some variants for profile creation.  If you intend to use a stats-enabled alignment function, use a corresponding stats-enabled profile creation function.  If you intend to use CPU dispatching for your alignment function, use a corresponding profile creation function, e.g., `parasail_profile_create_16`; no need to specify the vector instruction set.  The 'sat' profile creation function is compatible with all integer width profile-based functions because it will contain profiles for each of the implemented widths -- just note it will of course use more memory instead of selecting a specific integer width.

You must not forget to free the profile(s) when you are finished.  There is only one function used to free the profile memory and it will handle all profiles created by the various profile creation functions.

```C
void parasail_profile_free(parasail_profile_t *profile);
```

The profile data structure is part of parasail's public interface, though you really shouldn't access its members.  Most of the attributes are opaque blocks of memory that hold the vaious profiles.  Occasionally, however, it may be useful to refer back to the parameters that were used during profile creation, namely s1, s1Len, and the substitution matrix -- these attributes can be accessed though they should be treated as read-only.

```C
typedef struct parasail_profile {
    const char *s1;
    int s1Len;
    const parasail_matrix_t *matrix;
    struct parasail_profile_data profile8;
    struct parasail_profile_data profile16;
    struct parasail_profile_data profile32;
    struct parasail_profile_data profile64;
    void (*free)(void * profile);
    int stop;
} parasail_profile_t;
```

### Substitution Matrices

[back to top](#parasail-pairwise-sequence-alignment-library)

parasails bundles a number of substitution matrices including PAM and BLOSUM.  To use them, include the appropriate header, or look them up by name (useful for command-line parsing). For example

```C
#include "parasail.h"
#include "parasail/matrices/blosum62.h"
#include "parasail/matrix_lookup.h"

int main(int argc, char **argv) {
        const char s1 = "asdf";
        const char s2 = "asdf";
        int s1Len = (int)strlen(s1);
        int s2Len = (int)strlen(s2);
        parasail_result_t *result = NULL;
        const parasail_matrix_t *matrix = NULL;
        
        /* note the address-of operator '&' */
        result = parasail_sw(s1, s1Len, s2, s2Len, -11, -1, &parasail_blosum62);
        parasail_result_free(result);
        
        matrix = parasail_matrix_lookup("pam100");
        result = parasail_sw(s1, s1Len, s2, s2Len, -11, -1, matrix);
        parasail_result_free(result);
}
```

### Function Lookup

[back to top](#parasail-pairwise-sequence-alignment-library)

Typically used for command-line or user-input parsing, the parasail functions can be accessed using their string name.  For example

```C
#include "parasail.h"
#include "parasail/matrices/blosum62.h"

int main(int argc, char **argv) {
        const char s1 = "asdf";
        const char s2 = "asdf";
        int s1Len = (int)strlen(s1);
        int s2Len = (int)strlen(s2);
        parasail_result_t *result = NULL;
        const parasail_matrix_t *matrix = NULL;
        parasail_function_t *function = NULL;
        
        function = parasail_lookup_function(argv[1]);
        result = function(s1, s1Len, s2, s2Len, -11, -1, &parasail_blosum62);
        parasail_result_free(result);
        
        /* 'parasail_' prefix is optional */
        function = parasail_lookup_function("nw_striped_32");
        result = function(s1, s1Len, s2, s2Len, -11, -1, &parasail_blosum62);
        parasail_result_free(result);
}
```

## Language Bindings

[back to top](#parasail-pairwise-sequence-alignment-library)

### C/C++

[back to top](#parasail-pairwise-sequence-alignment-library)

C is the native API for parasail.  C++ is supported directly because the parasial.h header uses the common C++ include guards (#ifdef __cplusplus) to extern "C" all of the functions.  Once you have installed parasail, #include "parasail.h" into your sources.

### Python

[back to top](#parasail-pairwise-sequence-alignment-library)

Once you have installed parasail into --prefix=$PREFIX, you can also compile the Python bindings.  Don't forget to add $PREFIX/lib to your LD_LIBRARY_PATH, if needed.  The Python bindings are in the <parasail>/bindings/python directory.  To build, run:

```
PARASAIL_PREFIX=$PREFIX python setup.py build
```

This will correctly setup the necessary CPPFLAGS, LDFLAGS, and LIBS variables during the build.  Becuase the parasail.h header uses C99 keywords, e.g., restrict, the setup.py process will test your C compiler for the correct use of restrict, automatically.

The Python interface only includes bindings for the dispatching functions, not the low-level instruction set-specific function calls.  The Python interface also includes wrappers for the various PAM and BLOSUM matrices included in the distribution.

Example:

```python
import parasail
result = parasail.sw_scan_16("asdf", "asdf", -11, -1, parasail.blosum62)
result = parasail.sw_stats_striped_8("asdf", "asdf", -11, -1, parasail.pam100)
```

## Windows Build

[back to top](#parasail-pairwise-sequence-alignment-library)

A parasail.dll was successfully created using a combination of Cygwin and its mingw64 package.

```bash
cd parasail
mkdir bld # it's cleaner to keep build and source directories separated
cd bld
../configure --prefix=`pwd` --host=x86_64-w64-mingw32 --build=x86_64-pc-cygwin
make
make install # installs DLL into `pwd`/bin
ldd bin/parasail.dll
        ntdll.dll => /cygdrive/c/Windows/SYSTEM32/ntdll.dll (0x77910000)
        kernel32.dll => /cygdrive/c/Windows/system32/kernel32.dll (0x777f0000)
        KERNELBASE.dll => /cygdrive/c/Windows/system32/KERNELBASE.dll (0x7fefd790000)
        msvcrt.dll => /cygdrive/c/Windows/system32/msvcrt.dll (0x7fefdbc0000)
        libwinpthread-1.dll => /usr/x86_64-w64-mingw32/sys-root/mingw/bin/libwinpthread-1.dll (0x64940000)
        USER32.dll => /cygdrive/c/Windows/system32/USER32.dll (0x776f0000)
        GDI32.dll => /cygdrive/c/Windows/system32/GDI32.dll (0x7fefdc60000)
        LPK.dll => /cygdrive/c/Windows/system32/LPK.dll (0x7fefde30000)
        USP10.dll => /cygdrive/c/Windows/system32/USP10.dll (0x7fefde40000)
        IMM32.DLL => /cygdrive/c/Windows/system32/IMM32.DLL (0x7feff370000)
        MSCTF.dll => /cygdrive/c/Windows/system32/MSCTF.dll (0x7feff000000)
        KATRK64.DLL => /cygdrive/c/Windows/KATRK64.DLL (0x180000000)
        WTSAPI32.dll => /cygdrive/c/Windows/system32/WTSAPI32.dll (0x7fefc920000)
```
Note that parasail.dll depends on a `libwinpthread-1.dll` which is located in the sys-root of the corresponding mingw installation.  See the `ldd` output above for an example.  If it weren't for the `parasail_time(void)` function, parasail would have zero external dependencies on Windows.

## Example Applications

[back to top](#parasail-pairwise-sequence-alignment-library)

In addition to the parasail library, there is one binary that is also compiled and installed.
See [the apps README](apps/README.md) for more details.

## Citing parasail

[back to top](#parasail-pairwise-sequence-alignment-library)

If needed, please cite the following paper.

Daily, Jeff. 2015. "Scalable Parallel Methods for Analyzing Metagenomic Data at Extreme Scale". PhD dissertation, Washington State University.  http://www.pnnl.gov/main/publications/external/technical_reports/PNNL-24266.pdf

## License: Battelle BSD-style

[back to top](#parasail-pairwise-sequence-alignment-library)

Copyright (c) 2015, Battelle Memorial Institute

1.  Battelle Memorial Institute (hereinafter Battelle) hereby grants
    permission to any person or entity lawfully obtaining a copy of this
    software and associated documentation files (hereinafter “the
    Software”) to redistribute and use the Software in source and binary
    forms, with or without modification.  Such person or entity may use,
    copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and may permit others to do so, subject to
    the following conditions:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimers.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.

    - Other than as used herein, neither the name Battelle Memorial
      Institute or Battelle may be used in any form whatsoever without
      the express written consent of Battelle.

    - Redistributions of the software in any form, and publications
      based on work performed using the software should include the
      following citation as a reference:

    Daily, Jeff. 2015. "Scalable Parallel Methods for Analyzing
    Metagenomic Data at Extreme Scale". PhD dissertation, Washington
    State University.  http://www.pnnl.gov/main/publications/external/technical_reports/PNNL-24266.pdf

2.  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL BATTELLE
    OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
    USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
    OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
    SUCH DAMAGE.

