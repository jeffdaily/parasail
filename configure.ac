#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
###############################################################################
# Init autoconf
###############################################################################
AC_PREREQ([2.68])
AC_INIT([Pairwise Sequence Alignment Library],
        [m4_esyscmd([./version.sh])],
        [jeff.daily@pnnl.gov],
        [parasail])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([parasail.h])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

# This lets us use PACKAGE_VERSION in Makefiles
AC_SUBST([PACKAGE_VERSION])

# Libtool -version-info (ABI version)
#
# Don't change this unless you know exactly what you're doing and have read and
# understand:
# http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
#
# Changes:
#
# parasail version 1.0.0: 1:0:0 (ABI version 1)
# parasail version 1.0.1: 1:1:0 (ABI version 1)
# parasail version 1.0.2: 1:2:0 (ABI version 1)
# parasail version 1.0.3: 1:3:0 (ABI version 1)
# parasail version 1.1.0: 1:4:0 (ABI version 1)
# parasail version 1.1.1: 1:5:0 (ABI version 1)
# parasail version 1.1.2: 1:5:0 (ABI version 1) no change
# parasail version 1.2.0: 2:0:0 (ABI version 2)
# parasail version 1.3.0: 2:0:0 (ABI version 2) no change
# parasail version 1.3.1: 2:0:0 (ABI version 2) no change
# parasail version 2.0.0: 3:0:0 (ABI version 3)
# parasail version 2.0.1: 3:1:0 (ABI version 3)
# parasail version 2.0.2: 3:2:0 (ABI version 3)
# parasail version 2.0.3: 3:2:0 (ABI version 3) no change
# parasail version 2.0.4: 3:3:0 (ABI version 3)
#
# libparasail -version-info current:revision:age
LTVER="3:3:0"
AC_SUBST([LTVER])

###############################################################################
# Init automake
###############################################################################
AM_INIT_AUTOMAKE([color-tests foreign parallel-tests silent-rules subdir-objects])

# Don't emit "rebuild rules" for configure, Makefile.ins, etc.
AM_MAINTAINER_MODE

###############################################################################
# Platform identity
###############################################################################
AC_CANONICAL_HOST

is_power_isa=0
AS_CASE(["$host_cpu"],
        [*powerpc*|*ppc*], [is_power_isa=1])
AM_CONDITIONAL([IS_POWER_ISA], [test "x$is_power_isa" = x1])

###############################################################################
# C compiler
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([C compiler])
AC_MSG_NOTICE

AC_PROG_CC

# Checks for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT
AC_CHECK_SIZEOF([int])
AC_SYS_LARGEFILE
AC_TYPE_INT8_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AX_OPENMP
AC_SUBST([OPENMP_CFLAGS])
AM_CONDITIONAL([HAVE_OPENMP],
               [test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" != "xunknown"])

AC_CACHE_CHECK([for SSE2 support], [pt_cv_sse2],
[pt_cv_sse2=no; pt_cv_sse2_val=0
for flag in default -msse2 -march=core2
do 
     save_CFLAGS="$CFLAGS"
     AS_IF([test "x$flag" != xdefault], [CFLAGS="$CFLAGS $flag"])
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <emmintrin.h>
int foo() {
    __m128i vOne = _mm_set1_epi16(1);
    __m128i result = _mm_add_epi16(vOne,vOne);
    return _mm_extract_epi16(result, 0);
}
int main(void) { return foo(); }
]])],
        [pt_cv_sse2="$flag"; pt_cv_sse2_val=1])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"
     AS_IF([test "x$pt_cv_sse2" != xno], [break])
done])
AC_DEFINE_UNQUOTED([HAVE_SSE2], [$pt_cv_sse2_val],
    [define to 1 if the C compiler supports SSE2])
AM_CONDITIONAL([HAVE_SSE2], [test "x$pt_cv_sse2_val" = x1])
AS_IF([test "x$pt_cv_sse2" == xdefault], [SSE2_CFLAGS=],
      [test "x$pt_cv_sse2" != xno],      [SSE2_CFLAGS="$pt_cv_sse2"],
      [SSE2_CFLAGS="not supported"])
AC_SUBST([SSE2_CFLAGS])

AC_CACHE_CHECK([for SSE2 _mm_set1_epi64x],
    [pt_cv_sse2_mm_set1_epi64x],
    [save_CFLAGS="$CFLAGS"
     CFLAGS="$SSE2_CFLAGS $CFLAGS"
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
#include <emmintrin.h>
__m128i foo() {
    __m128i vOne = _mm_set1_epi64x(1);
    return vOne;
}
]])],
        [pt_cv_sse2_mm_set1_epi64x=yes;
         pt_cv_sse2_mm_set1_epi64x_val=1],
        [pt_cv_sse2_mm_set1_epi64x=no;
         pt_cv_sse2_mm_set1_epi64x_val=0])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"])
AC_DEFINE_UNQUOTED([HAVE_SSE2_MM_SET1_EPI64X],
    [$pt_cv_sse2_mm_set1_epi64x_val],
    [define to 1 if the C compiler supports SSE2 _mm_set1_epi64x])

AC_CACHE_CHECK([for SSE2 _mm_set_epi64x],
    [pt_cv_sse2_mm_set_epi64x],
    [save_CFLAGS="$CFLAGS"
     CFLAGS="$SSE2_CFLAGS $CFLAGS"
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
#include <emmintrin.h>
__m128i foo() {
    __m128i vOne = _mm_set_epi64x(1,1);
    return vOne;
}
]])],
        [pt_cv_sse2_mm_set_epi64x=yes;
         pt_cv_sse2_mm_set_epi64x_val=1],
        [pt_cv_sse2_mm_set_epi64x=no;
         pt_cv_sse2_mm_set_epi64x_val=0])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"])
AC_DEFINE_UNQUOTED([HAVE_SSE2_MM_SET_EPI64X],
    [$pt_cv_sse2_mm_set_epi64x_val],
    [define to 1 if the C compiler supports SSE2 _mm_set_epi64x])

AC_CACHE_CHECK([for SSE4.1 support], [pt_cv_sse41],
[pt_cv_sse41=no; pt_cv_sse41_val=0
for flag in default -msse4.1 -msse4 -march=corei7
do
     save_CFLAGS="$CFLAGS"
     AS_IF([test "x$flag" != xdefault], [CFLAGS="$CFLAGS $flag"])
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <smmintrin.h>
int foo() {
    __m128i vOne = _mm_set1_epi8(1);
    __m128i result =  _mm_max_epi8(vOne,vOne);
    return _mm_extract_epi8(result, 0);
}
int main(void) { return foo(); }
]])],
        [pt_cv_sse41="$flag"; pt_cv_sse41_val=1])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"
     AS_IF([test "x$pt_cv_sse41" != xno], [break])
done])
AC_DEFINE_UNQUOTED([HAVE_SSE41], [$pt_cv_sse41_val],
    [define to 1 if the C compiler supports SSE4.1])
AM_CONDITIONAL([HAVE_SSE41], [test "x$pt_cv_sse41_val" = x1])
AS_IF([test "x$pt_cv_sse41" == xdefault], [SSE41_CFLAGS=],
      [test "x$pt_cv_sse41" != xno],      [SSE41_CFLAGS="$pt_cv_sse41"],
      [SSE41_CFLAGS="not supported"])
AC_SUBST([SSE41_CFLAGS])

AM_CONDITIONAL([HAVE_SSE], [test "x$pt_cv_sse2_val" = x1 || test "x$pt_cv_sse41_val" = x1])

AC_CACHE_CHECK([for SSE4.1 _mm_insert_epi64],
    [pt_cv_sse41_mm_insert_epi64],
    [save_CFLAGS="$CFLAGS"
     CFLAGS="$SSE41_CFLAGS $CFLAGS"
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
#include <smmintrin.h>
__m128i foo() {
    __m128i vOne = _mm_set1_epi8(1);
    return _mm_insert_epi64(vOne,INT64_MIN,0);
}
]])],
        [pt_cv_sse41_mm_insert_epi64=yes;
         pt_cv_sse41_mm_insert_epi64_val=1],
        [pt_cv_sse41_mm_insert_epi64=no;
         pt_cv_sse41_mm_insert_epi64_val=0])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"])
AC_DEFINE_UNQUOTED([HAVE_SSE41_MM_INSERT_EPI64],
    [$pt_cv_sse41_mm_insert_epi64_val],
    [define to 1 if the C compiler supports SSE4.1 _mm_insert_epi64])

AC_CACHE_CHECK([for SSE4.1 _mm_extract_epi64],
    [pt_cv_sse41_mm_extract_epi64],
    [save_CFLAGS="$CFLAGS"
     CFLAGS="$SSE41_CFLAGS $CFLAGS"
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
#include <immintrin.h>
int64_t foo() {
    __m128i vZero = _mm_set1_epi8(0);
    return _mm_extract_epi64(vZero,0);
}
int main(void) { return foo(); }
]])],
        [pt_cv_sse41_mm_extract_epi64=yes;
         pt_cv_sse41_mm_extract_epi64_val=1],
        [pt_cv_sse41_mm_extract_epi64=no;
         pt_cv_sse41_mm_extract_epi64_val=0])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"])
AC_DEFINE_UNQUOTED([HAVE_SSE41_MM_EXTRACT_EPI64],
    [$pt_cv_sse41_mm_extract_epi64_val],
    [define to 1 if the C compiler supports SSE4.1 _mm_extract_epi64])

AC_CACHE_CHECK([for AVX2 support], [pt_cv_avx2],
[pt_cv_avx2=no; pt_cv_avx2_val=0
for flag in default -mavx2 -march=core-avx2
do
     save_CFLAGS="$CFLAGS"
     AS_IF([test "x$flag" != xdefault], [CFLAGS="$CFLAGS $flag"])
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <immintrin.h>
void parasail_memset___m256i(__m256i *b, __m256i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm256_store_si256(&b[i], c);
    }
}

int foo() {
    __m256i vOne = _mm256_set1_epi8(1);
    __m256i result =  _mm256_add_epi8(vOne,vOne);
    return _mm_extract_epi16(_mm256_extracti128_si256(result,0),0);
}
int main(void) { return foo(); }
]])],
        [pt_cv_avx2="$flag"; pt_cv_avx2_val=1])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"
     AS_IF([test "x$pt_cv_avx2" != xno], [break])
done])
AC_DEFINE_UNQUOTED([HAVE_AVX2], [$pt_cv_avx2_val],
    [define to 1 if the C compiler supports AVX2])
AM_CONDITIONAL([HAVE_AVX2], [test "x$pt_cv_avx2_val" = x1])
AS_IF([test "x$pt_cv_avx2" == xdefault], [AVX2_CFLAGS=],
      [test "x$pt_cv_avx2" != xno],      [AVX2_CFLAGS="$pt_cv_avx2"],
      [AVX2_CFLAGS="not supported"])
AC_SUBST([AVX2_CFLAGS])

AC_CACHE_CHECK([for AVX2 _mm256_set1_epi64x],
    [pt_cv_avx2_mm256_set1_epi64x],
    [save_CFLAGS="$CFLAGS"
     CFLAGS="$AVX2_CFLAGS $CFLAGS"
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
#include <immintrin.h>
__m256i foo() {
    __m256i vOne = _mm256_set1_epi64x(1);
    return vOne;
}
]])],
        [pt_cv_avx2_mm256_set1_epi64x=yes;
         pt_cv_avx2_mm256_set1_epi64x_val=1],
        [pt_cv_avx2_mm256_set1_epi64x=no;
         pt_cv_avx2_mm256_set1_epi64x_val=0])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"])
AC_DEFINE_UNQUOTED([HAVE_AVX2_MM256_SET1_EPI64X],
    [$pt_cv_avx2_mm256_set1_epi64x_val],
    [define to 1 if the C compiler supports AVX2 _mm256_set1_epi64x])

AC_CACHE_CHECK([for AVX2 _mm256_set_epi64x],
    [pt_cv_avx2_mm256_set_epi64x],
    [save_CFLAGS="$CFLAGS"
     CFLAGS="$AVX2_CFLAGS $CFLAGS"
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
#include <immintrin.h>
__m256i foo() {
    __m256i vOne = _mm256_set_epi64x(1,1,1,1);
    return vOne;
}
]])],
        [pt_cv_avx2_mm256_set_epi64x=yes;
         pt_cv_avx2_mm256_set_epi64x_val=1],
        [pt_cv_avx2_mm256_set_epi64x=no;
         pt_cv_avx2_mm256_set_epi64x_val=0])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"])
AC_DEFINE_UNQUOTED([HAVE_AVX2_MM256_SET_EPI64X],
    [$pt_cv_avx2_mm256_set_epi64x_val],
    [define to 1 if the C compiler supports AVX2 _mm256_set_epi64x])

AC_CACHE_CHECK([for AVX2 _mm256_insert_epi64],
    [pt_cv_avx2_mm256_insert_epi64],
    [save_CFLAGS="$CFLAGS"
     CFLAGS="$AVX2_CFLAGS $CFLAGS"
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
#include <immintrin.h>
__m256i foo() {
    __m256i vOne = _mm256_set1_epi8(1);
    return _mm256_insert_epi64(vOne,INT64_MIN,0);
}
]])],
        [pt_cv_avx2_mm256_insert_epi64=yes;
         pt_cv_avx2_mm256_insert_epi64_val=1],
        [pt_cv_avx2_mm256_insert_epi64=no;
         pt_cv_avx2_mm256_insert_epi64_val=0])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"])
AC_DEFINE_UNQUOTED([HAVE_AVX2_MM256_INSERT_EPI64],
    [$pt_cv_avx2_mm256_insert_epi64_val],
    [define to 1 if the C compiler supports AVX2 _mm256_insert_epi64])

AC_CACHE_CHECK([for AVX2 _mm256_insert_epi32],
    [pt_cv_avx2_mm256_insert_epi32],
    [save_CFLAGS="$CFLAGS"
     CFLAGS="$AVX2_CFLAGS $CFLAGS"
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
#include <immintrin.h>
__m256i foo() {
    __m256i vOne = _mm256_set1_epi8(1);
    return _mm256_insert_epi32(vOne,INT32_MIN,0);
}
]])],
        [pt_cv_avx2_mm256_insert_epi32=yes;
         pt_cv_avx2_mm256_insert_epi32_val=1],
        [pt_cv_avx2_mm256_insert_epi32=no;
         pt_cv_avx2_mm256_insert_epi32_val=0])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"])
AC_DEFINE_UNQUOTED([HAVE_AVX2_MM256_INSERT_EPI32],
    [$pt_cv_avx2_mm256_insert_epi32_val],
    [define to 1 if the C compiler supports AVX2 _mm256_insert_epi32])

AC_CACHE_CHECK([for AVX2 _mm256_insert_epi16],
    [pt_cv_avx2_mm256_insert_epi16],
    [save_CFLAGS="$CFLAGS"
     CFLAGS="$AVX2_CFLAGS $CFLAGS"
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
#include <immintrin.h>
__m256i foo() {
    __m256i vZero = _mm256_set1_epi8(0);
    return _mm256_insert_epi16(vZero,INT16_MIN,0);
}
]])],
        [pt_cv_avx2_mm256_insert_epi16=yes;
         pt_cv_avx2_mm256_insert_epi16_val=1],
        [pt_cv_avx2_mm256_insert_epi16=no;
         pt_cv_avx2_mm256_insert_epi16_val=0])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"])
AC_DEFINE_UNQUOTED([HAVE_AVX2_MM256_INSERT_EPI16],
    [$pt_cv_avx2_mm256_insert_epi16_val],
    [define to 1 if the C compiler supports AVX2 _mm256_insert_epi16])

AC_CACHE_CHECK([for AVX2 _mm256_insert_epi8],
    [pt_cv_avx2_mm256_insert_epi8],
    [save_CFLAGS="$CFLAGS"
     CFLAGS="$AVX2_CFLAGS $CFLAGS"
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
#include <immintrin.h>
__m256i foo() {
    __m256i vZero = _mm256_set1_epi8(0);
    return _mm256_insert_epi8(vZero,INT8_MIN,0);
}
]])],
        [pt_cv_avx2_mm256_insert_epi8=yes;
         pt_cv_avx2_mm256_insert_epi8_val=1],
        [pt_cv_avx2_mm256_insert_epi8=no;
         pt_cv_avx2_mm256_insert_epi8_val=0])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"])
AC_DEFINE_UNQUOTED([HAVE_AVX2_MM256_INSERT_EPI8],
    [$pt_cv_avx2_mm256_insert_epi8_val],
    [define to 1 if the C compiler supports AVX2 _mm256_insert_epi8])

AC_CACHE_CHECK([for AVX2 _mm256_extract_epi64],
    [pt_cv_avx2_mm256_extract_epi64],
    [save_CFLAGS="$CFLAGS"
     CFLAGS="$AVX2_CFLAGS $CFLAGS"
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
#include <immintrin.h>
int64_t foo() {
    __m256i vZero = _mm256_set1_epi8(0);
    return _mm256_extract_epi64(vZero,0);
}
int main(void) { return foo(); }
]])],
        [pt_cv_avx2_mm256_extract_epi64=yes;
         pt_cv_avx2_mm256_extract_epi64_val=1],
        [pt_cv_avx2_mm256_extract_epi64=no;
         pt_cv_avx2_mm256_extract_epi64_val=0])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"])
AC_DEFINE_UNQUOTED([HAVE_AVX2_MM256_EXTRACT_EPI64],
    [$pt_cv_avx2_mm256_extract_epi64_val],
    [define to 1 if the C compiler supports AVX2 _mm256_extract_epi64])

AC_CACHE_CHECK([for AVX2 _mm256_extract_epi32],
    [pt_cv_avx2_mm256_extract_epi32],
    [save_CFLAGS="$CFLAGS"
     CFLAGS="$AVX2_CFLAGS $CFLAGS"
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
#include <immintrin.h>
int32_t foo() {
    __m256i vZero = _mm256_set1_epi8(0);
    return _mm256_extract_epi32(vZero,0);
}
int main(void) { return foo(); }
]])],
        [pt_cv_avx2_mm256_extract_epi32=yes;
         pt_cv_avx2_mm256_extract_epi32_val=1],
        [pt_cv_avx2_mm256_extract_epi32=no;
         pt_cv_avx2_mm256_extract_epi32_val=0])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"])
AC_DEFINE_UNQUOTED([HAVE_AVX2_MM256_EXTRACT_EPI32],
    [$pt_cv_avx2_mm256_extract_epi32_val],
    [define to 1 if the C compiler supports AVX2 _mm256_extract_epi32])

AC_CACHE_CHECK([for AVX2 _mm256_extract_epi16],
    [pt_cv_avx2_mm256_extract_epi16],
    [save_CFLAGS="$CFLAGS"
     CFLAGS="$AVX2_CFLAGS $CFLAGS"
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
#include <immintrin.h>
int16_t foo() {
    __m256i vZero = _mm256_set1_epi8(0);
    return _mm256_extract_epi16(vZero,0);
}
int main(void) { return foo(); }
]])],
        [pt_cv_avx2_mm256_extract_epi16=yes;
         pt_cv_avx2_mm256_extract_epi16_val=1],
        [pt_cv_avx2_mm256_extract_epi16=no;
         pt_cv_avx2_mm256_extract_epi16_val=0])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"])
AC_DEFINE_UNQUOTED([HAVE_AVX2_MM256_EXTRACT_EPI16],
    [$pt_cv_avx2_mm256_extract_epi16_val],
    [define to 1 if the C compiler supports AVX2 _mm256_extract_epi16])

AC_CACHE_CHECK([for AVX2 _mm256_extract_epi8],
    [pt_cv_avx2_mm256_extract_epi8],
    [save_CFLAGS="$CFLAGS"
     CFLAGS="$AVX2_CFLAGS $CFLAGS"
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
#include <immintrin.h>
int8_t foo() {
    __m256i vZero = _mm256_set1_epi8(0);
    return _mm256_extract_epi8(vZero,0);
}
int main(void) { return foo(); }
]])],
        [pt_cv_avx2_mm256_extract_epi8=yes;
         pt_cv_avx2_mm256_extract_epi8_val=1],
        [pt_cv_avx2_mm256_extract_epi8=no;
         pt_cv_avx2_mm256_extract_epi8_val=0])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"])
AC_DEFINE_UNQUOTED([HAVE_AVX2_MM256_EXTRACT_EPI8],
    [$pt_cv_avx2_mm256_extract_epi8_val],
    [define to 1 if the C compiler supports AVX2 _mm256_extract_epi8])

AC_CACHE_CHECK([for KNC],
    [pt_cv_knc],
    [AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#ifdef __MIC__
#include <immintrin.h>
#else
#error CHOKE
#endif
]])],
        [pt_cv_knc=yes; pt_cv_knc_val=1],
        [pt_cv_knc=no; pt_cv_knc_val=0])])
AC_DEFINE_UNQUOTED([HAVE_KNC], [$pt_cv_knc_val],
    [define to 1 if the C compiler supports KNC])
AM_CONDITIONAL([HAVE_KNC], [test "x$pt_cv_knc" = xyes])

AC_CACHE_CHECK([for AVX512F support], [pt_cv_avx512f],
[pt_cv_avx512f=no; pt_cv_avx512f_val=0
for flag in default -mavx512f -xMIC-AVX512 -march=knl
do
     save_CFLAGS="$CFLAGS"
     AS_IF([test "x$flag" != xdefault], [CFLAGS="$CFLAGS $flag"])
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <immintrin.h>
void parasail_memset___m512i(__m512i *b, __m512i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm512_store_si512(&b[i], c);
    }
}

int foo() {
    __m512i vOne = _mm512_set1_epi32(1);
    __m512i result =  _mm512_add_epi32(vOne,vOne);
    return 0;
}
int main(void) { return foo(); }
]])],
        [pt_cv_avx512f="$flag"; pt_cv_avx512f_val=1])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"
     AS_IF([test "x$pt_cv_avx512f" != xno], [break])
done])
AC_DEFINE_UNQUOTED([HAVE_AVX512F], [$pt_cv_avx512f_val],
    [define to 1 if the C compiler supports AVX512F])
AM_CONDITIONAL([HAVE_AVX512F], [test "x$pt_cv_avx512f_val" = x1])
AS_IF([test "x$pt_cv_avx512f" == xdefault], [AVX512F_CFLAGS=],
      [test "x$pt_cv_avx512f" != xno],      [AVX512F_CFLAGS="$pt_cv_avx512f"],
      [AVX512F_CFLAGS="not supported"])
AC_SUBST([AVX512F_CFLAGS])

AC_CACHE_CHECK([for AVX512BW support], [pt_cv_avx512bw],
[pt_cv_avx512bw=no; pt_cv_avx512bw_val=0
for flag in default -mavx512bw -xCORE-AVX512 -march=skylake-avx512
do
     save_CFLAGS="$CFLAGS"
     AS_IF([test "x$flag" != xdefault], [CFLAGS="$CFLAGS $flag"])
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <immintrin.h>
void parasail_memset___m512i(__m512i *b, __m512i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm512_store_si512(&b[i], c);
    }
}

int foo() {
    __m512i vOne = _mm512_set1_epi8(1);
    __m512i result =  _mm512_add_epi8(vOne,vOne);
    return 0;
}
int main(void) { return foo(); }
]])],
        [pt_cv_avx512bw="$flag"; pt_cv_avx512bw_val=1])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"
     AS_IF([test "x$pt_cv_avx512bw" != xno], [break])
done])
AC_DEFINE_UNQUOTED([HAVE_AVX512BW], [$pt_cv_avx512bw_val],
    [define to 1 if the C compiler supports AVX512BW])
AM_CONDITIONAL([HAVE_AVX512bw], [test "x$pt_cv_avx512bw_val" = x1])
AS_IF([test "x$pt_cv_avx512bw" == xdefault], [AVX512BW_CFLAGS=],
      [test "x$pt_cv_avx512bw" != xno],      [AVX512BW_CFLAGS="$pt_cv_avx512bw"],
      [AVX512BW_CFLAGS="not supported"])
AC_SUBST([AVX512BW_CFLAGS])

AC_CACHE_CHECK([for AVX512VBMI support], [pt_cv_avx512vbmi],
[pt_cv_avx512vbmi=no; pt_cv_avx512vbmi_val=0
for flag in default -mavx512vbmi -xCORE-AVX512 -march=skylake-avx512
do
     save_CFLAGS="$CFLAGS"
     AS_IF([test "x$flag" != xdefault], [CFLAGS="$CFLAGS $flag"])
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <immintrin.h>
void parasail_memset___m512i(__m512i *b, __m512i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm512_store_si512(&b[i], c);
    }
}

int foo() {
    __m512i vOne = _mm512_set1_epi8(1);
    __m512i result =  _mm512_permutex2var_epi8(vOne,vOne,vOne);
    return 0;
}
int main(void) { return foo(); }
]])],
        [pt_cv_avx512vbmi="$flag"; pt_cv_avx512vbmi_val=1])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"
     AS_IF([test "x$pt_cv_avx512vbmi" != xno], [break])
done])
AC_DEFINE_UNQUOTED([HAVE_AVX512VBMI], [$pt_cv_avx512vbmi_val],
    [define to 1 if the C compiler supports AVX512VBMI])
AM_CONDITIONAL([HAVE_AVX512vbmi], [test "x$pt_cv_avx512vbmi_val" = x1])
AS_IF([test "x$pt_cv_avx512vbmi" == xdefault], [AVX512VBMI_CFLAGS=],
      [test "x$pt_cv_avx512vbmi" != xno],      [AVX512VBMI_CFLAGS="$pt_cv_avx512vbmi"],
      [AVX512VBMI_CFLAGS="not supported"])
AC_SUBST([AVX512VBMI_CFLAGS])

AC_CACHE_CHECK([for xgetbv],
    [pt_cv_xgetbv],
    [AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <stdint.h>
int check_xcr0_ymm()
{
    uint32_t xcr0;
#if defined(_MSC_VER)
    xcr0 = (uint32_t)_xgetbv(0);  /* min VS2010 SP1 compiler is required */
#else
    __asm__ ("xgetbv" : "=a" (xcr0) : "c" (0) : "%edx" );
#endif
    return ((xcr0 & 6) == 6); /* checking if xmm and ymm state are enabled in XCR0 */
}
]])],
        [pt_cv_xgetbv=yes; pt_cv_xgetbv_val=1],
        [pt_cv_xgetbv=no; pt_cv_xgetbv_val=0])])
AC_DEFINE_UNQUOTED([HAVE_XGETBV], [$pt_cv_xgetbv_val],
    [define to 1 if the C compiler supports xgetbv])
AM_CONDITIONAL([HAVE_XGETBV], [test "x$pt_cv_xgetbv" = xyes])

AC_CACHE_CHECK([for AltiVec support], [pt_cv_altivec],
[pt_cv_altivec=no; pt_cv_altivec_val=0
for flag in default -maltivec -faltivec
do 
     save_CFLAGS="$CFLAGS"
     AS_IF([test "x$flag" != xdefault], [CFLAGS="$CFLAGS $flag"])
     psl_save_c_werror_flag=$ac_c_werror_flag
     ac_c_werror_flag=yes
     AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <altivec.h>
vector signed int vec_2sComp (vector signed int x)
{
    vector signed int one = (vector signed int) (1);
    x = vec_nor (x, x);
    x = vec_add (x, one);
    return x;
}
int main(void)
{
    vector signed int two = (vector signed int) (2);
    vector signed int result = vec_2sComp(two);
    return vec_extract(two,0);
}
]])],
        [pt_cv_altivec="$flag"; pt_cv_altivec_val=1])
     ac_c_werror_flag=$psl_save_c_werror_flag
     CFLAGS="$save_CFLAGS"
     AS_IF([test "x$pt_cv_altivec" != xno], [break])
done])
AC_DEFINE_UNQUOTED([HAVE_ALTIVEC], [$pt_cv_altivec_val],
    [define to 1 if the C compiler supports ALTIVEC])
AM_CONDITIONAL([HAVE_ALTIVEC], [test "x$pt_cv_altivec_val" = x1])
AS_IF([test "x$pt_cv_altivec" == xdefault], [ALTIVEC_CFLAGS=],
      [test "x$pt_cv_altivec" != xno],      [ALTIVEC_CFLAGS="$pt_cv_altivec"],
      [ALTIVEC_CFLAGS="not supported"])
AC_SUBST([ALTIVEC_CFLAGS])

# Checks for library functions.
AC_CHECK_FUNCS([_aligned_malloc])
AC_CHECK_FUNCS([aligned_alloc])
AC_CHECK_FUNCS([memalign])
AC_CHECK_FUNCS([memset])
AC_CHECK_FUNCS([posix_memalign])

# clock_gettime is in librt on linux, libpthread on windows
parasail_save_LIBS="$LIBS"
AC_SEARCH_LIBS([clock_gettime], [rt pthread])
LIBS="$parasail_save_LIBS"
AS_CASE(["$ac_cv_search_clock_gettime"],
        [*none*], [],
        [*no*], [],
        [CLOCK_LIBS="$ac_cv_search_clock_gettime"])
AC_SUBST([CLOCK_LIBS])

parasail_save_LIBS="$LIBS"
AC_SEARCH_LIBS([sqrt], [m])
LIBS="$parasail_save_LIBS"
AS_CASE(["$ac_cv_search_sqrt"],
        [*none*], [],
        [*no*], [],
        [MATH_LIBS="$ac_cv_search_sqrt"])
AC_SUBST([MATH_LIBS])

PKG_CHECK_MODULES([Z], [zlib], [
    zlib_okay_val=1
    ], [
    parasail_save_LIBS="$LIBS"
    AC_CHECK_HEADERS([zlib.h])
    AC_SEARCH_LIBS([gzopen], [z])
    LIBS="$parasail_save_LIBS"
    AS_CASE(["$ac_cv_search_gzopen"],
            [*none*], [],
            [*no*], [],
            [Z_LIBS="$ac_cv_search_gzopen"])
    AC_SUBST([Z_LIBS])
    AS_IF([test "x$ac_cv_search_gzopen" != xno && test "x$ac_cv_header_zlib_h" = xyes],
          [zlib_okay_val=1],
          [zlib_okay_val=0])
    ])

AM_CONDITIONAL([HAVE_ZLIB], [test "x$zlib_okay_val" = x1])
AC_DEFINE_UNQUOTED([HAVE_ZLIB], [$zlib_okay_val], [define to 1 if zlib is available])

# For host type checks
AC_CANONICAL_HOST

# OS-specific tests
psl_on_mingw32="no"
psl_dso_visibility="yes"
case "${host_os}" in
    *mingw32*)
        psl_on_mingw32="yes"
        psl_dso_visibility="no"
        if test "x$enable_static" = "xyes"; then
            AC_MSG_ERROR([Building static libraries is not supported under MinGW32])
        fi
        ;;
    *cygwin*)
        if test "x$enable_static" = "xyes"; then
            AC_MSG_ERROR([Building static libraries is not supported under Cygwin])
        fi
        ;;
esac

AM_CONDITIONAL([ON_MINGW], [test "x$psl_on_mingw32" = "xyes"])

#
# Check if the compiler supports -fvisibility=hidden flag.
# MinGW32 uses __declspec.
#
AC_PROG_GREP
AC_MSG_CHECKING([whether C compiler supports visibility flag])
if test "x$psl_dso_visibility" = "xyes"; then
    psl_save_c_werror_flag=$ac_c_werror_flag
    for flag in -fvisibility=hidden -xldscope=hidden; do
        AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM()],
            # This hack exist for ICC, which outputs unknown options as remarks
            # Remarks are not turned into errors even with -Werror on
            [if ($GREP 'ignoring unknown' conftest.err ||
                 $GREP 'not supported' conftest.err) >/dev/null 2>&1; then
                psl_dso_visibility="no"
            else
                psl_dso_visibility="$flag"
            fi],
            [psl_dso_visibility="no"])
        if test "x$psl_dso_visibility" != "xno"; then
            break
        fi
    done
    ac_c_werror_flag=$psl_save_c_werror_flag
fi
AC_MSG_RESULT([$psl_dso_visibility])
if test "x$psl_dso_visibility" = "xno"; then
    psl_dso_visibility=
fi
AC_SUBST([VISIBILITY_CFLAGS], [$psl_dso_visibility])

###############################################################################
# C++ compiler
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([C++ compiler])
AC_MSG_NOTICE
AC_LANG_PUSH([C++])

AC_PROG_CXX

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AX_OPENMP
AC_SUBST([OPENMP_CXXFLAGS])

AC_CACHE_CHECK([for variadic macros],
    [pt_cv_variadic],
    [psl_save_cxx_werror_flag=$ac_cxx_werror_flag
     ac_cxx_werror_flag=yes
     AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <wchar.h>
#include <stdio.h>
#define debug(...) fprintf (stderr, __VA_ARGS__)
#define showlist(...) puts (#__VA_ARGS__)
#define report(test,...) ((test) ? puts (#test) : printf (__VA_ARGS__))
static void
test_varargs_macros (void)
{
  int x = 1234;
  int y = 5678;
  debug ("Flag");
  debug ("X = %d\n", x);
  showlist (The first, second, and third items.);
  report (x>y, "x is %d but y is %d", x, y);
}
int main(void) { test_varargs_macros(); return 0; }
]])],
        [pt_cv_variadic=yes; pt_cv_variadic_val=1],
        [pt_cv_variadic=no; pt_cv_variadic_val=0])
     ac_cxx_werror_flag=$psl_save_cxx_werror_flag])
AC_DEFINE_UNQUOTED([HAVE_VARIADIC_MACROS], [$pt_cv_variadic_val],
    [define to 1 if the CXX compiler supports variadic macros])

# Checks for library functions.

AC_LANG_POP([C++])

###############################################################################
# Libtool setup -- no compiler/linker tests after this
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([Libtool setup])
AC_MSG_NOTICE

LT_INIT([win32-dll])

###############################################################################
# pkg-config
###############################################################################

AX_REQUIRE_DEFINED([PKG_INSTALLDIR])
PKG_INSTALLDIR

AX_CREATE_PKGCONFIG_INFO([parasail-1.pc], [], [-lparasail], [Pairwise Sequence Alignment Library], [], [])

###############################################################################
# The End
###############################################################################
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

# Report on what we found.
AC_MSG_NOTICE([])
AC_MSG_NOTICE([**************************************************************])
AC_MSG_NOTICE([ $PACKAGE_NAME])
AC_MSG_NOTICE([ configured as follows:])
AC_MSG_NOTICE([**************************************************************])
AC_MSG_NOTICE([])
AC_MSG_NOTICE([                    CC=$CC])
AC_MSG_NOTICE([                CFLAGS=$CFLAGS])
AC_MSG_NOTICE([         OPENMP_CFLAGS=$OPENMP_CFLAGS])
AC_MSG_NOTICE([           SSE2_CFLAGS=$SSE2_CFLAGS])
AC_MSG_NOTICE([          SSE41_CFLAGS=$SSE41_CFLAGS])
AC_MSG_NOTICE([           AVX2_CFLAGS=$AVX2_CFLAGS])
AC_MSG_NOTICE([        AVX512F_CFLAGS=$AVX512F_CFLAGS])
AC_MSG_NOTICE([       AVX512BW_CFLAGS=$AVX512BW_CFLAGS])
AC_MSG_NOTICE([     AVX512VBMI_CFLAGS=$AVX512VBMI_CFLAGS])
AC_MSG_NOTICE([        ALTIVEC_CFLAGS=$ALTIVEC_CFLAGS])
AC_MSG_NOTICE([                   CXX=$CXX])
AC_MSG_NOTICE([              CXXFLAGS=$CXXFLAGS])
AC_MSG_NOTICE([       OPENMP_CXXFLAGS=$OPENMP_CXXFLAGS])
AC_MSG_NOTICE([                   CPP=$CPP])
AC_MSG_NOTICE([              CPPFLAGS=$CPPFLAGS])
AC_MSG_NOTICE([               LDFLAGS=$LDFLAGS])
AC_MSG_NOTICE([                  LIBS=$LIBS])
AC_MSG_NOTICE([            CLOCK_LIBS=$CLOCK_LIBS])
AC_MSG_NOTICE([             MATH_LIBS=$MATH_LIBS])
AC_MSG_NOTICE([              Z_CFLAGS=$Z_CFLAGS])
AC_MSG_NOTICE([                Z_LIBS=$Z_LIBS])
