#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
###############################################################################
# Init autoconf
###############################################################################
AC_PREREQ([2.68])
AC_INIT([Pairwise Sequence Alignment Library],
        [0.1], [jeff.daily@pnnl.gov], parasail)
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/align.h])

###############################################################################
# Init automake
###############################################################################
AM_INIT_AUTOMAKE([color-tests foreign parallel-tests silent-rules subdir-objects])

###############################################################################
# C compiler
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([C compiler])
AC_MSG_NOTICE

# Checks for header files.
AC_CHECK_HEADERS([limits.h])
AC_CHECK_HEADERS([malloc.h])
AC_CHECK_HEADERS([sched.h])
AC_CHECK_HEADERS([stddef.h])
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_INT8_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AX_OPENMP
AC_SUBST([OPENMP_CFLAGS])

AC_CACHE_CHECK([for SSE2],
    [pt_cv_sse2],
    [AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <emmintrin.h>
__m128i foo() {
    __m128i vZero = _mm_set1_epi16(0);
    return _mm_add_epi16(vZero,vZero);
}
]])],
        [pt_cv_sse2=yes; pt_cv_sse2_val=1],
        [pt_cv_sse2=no]; pt_cv_sse2_val=0)])
AC_DEFINE_UNQUOTED([HAVE_SSE2], [$pt_cv_sse2_val],
    [define to 1 if the C compiler supports SSE2])
AM_CONDITIONAL([HAVE_SSE2], [test x$pt_cv_sse2 = xyes])

AC_CACHE_CHECK([for SSE4.1],
    [pt_cv_sse41],
    [AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <smmintrin.h>
__m128i foo() {
    __m128i vZero = _mm_set1_epi8(0);
    return _mm_max_epi16(vZero,vZero);
}
]])],
        [pt_cv_sse41=yes; pt_cv_sse41_val=1],
        [pt_cv_sse41=no]; pt_cv_sse41_val=0)])
AC_DEFINE_UNQUOTED([HAVE_SSE41], [$pt_cv_sse41_val],
    [define to 1 if the C compiler supports SSE4.1])
AM_CONDITIONAL([HAVE_SSE41], [test x$pt_cv_sse41 = xyes])

AC_CACHE_CHECK([for AVX2],
    [pt_cv_avx2],
    [AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <immintrin.h>
__m256i foo() {
    __m256i vZero = _mm256_set1_epi8(0);
    return _mm256_add_epi8(vZero,vZero);
}
]])],
        [pt_cv_avx2=yes; pt_cv_avx2_val=1],
        [pt_cv_avx2=no]; pt_cv_avx2_val=0)])
AC_DEFINE_UNQUOTED([HAVE_AVX2], [$pt_cv_avx2_val],
    [define to 1 if the C compiler supports AVX2])
AM_CONDITIONAL([HAVE_AVX2], [test x$pt_cv_avx2 = xyes])

AC_CACHE_CHECK([for KNC],
    [pt_cv_knc],
    [AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include <zmmintrin.h>
__mmask16 foo() {
    __mmask16 vMask;
    return _mm512_kand(vMask,vMask);
}
]])],
        [pt_cv_knc=yes; pt_cv_knc_val=1],
        [pt_cv_knc=no]; pt_cv_knc_val=0)])
AC_DEFINE_UNQUOTED([HAVE_KNC], [$pt_cv_knc_val],
    [define to 1 if the C compiler supports KNC])
AM_CONDITIONAL([HAVE_KNC], [test x$pt_cv_knc = xyes])

# Checks for library functions.
AC_CHECK_FUNCS([memset])

parasail_save_LIBS="$LIBS"
AC_SEARCH_LIBS([clock_gettime], [rt])
LIBS="$parasail_save_LIBS"
AS_CASE(["$ac_cv_search_clock_gettime"],
        [*none*], [],
        [*no*], [],
        [CLOCK_LIBS="$ac_cv_search_clock_gettime"])
AC_SUBST([CLOCK_LIBS])

###############################################################################
# Libtool setup -- no compiler/linker tests after this
###############################################################################
AC_MSG_NOTICE
AC_MSG_NOTICE([Libtool setup])
AC_MSG_NOTICE

LT_INIT([disable-shared])

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
AC_MSG_NOTICE([                   CPP=$CPP])
AC_MSG_NOTICE([              CPPFLAGS=$CPPFLAGS])
AC_MSG_NOTICE([               LDFLAGS=$LDFLAGS])
AC_MSG_NOTICE([                  LIBS=$LIBS])
