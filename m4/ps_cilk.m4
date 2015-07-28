# SYNOPSIS
#
#   PSL_CILK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile programs that use Intel Cilk.
#
#   On success, it sets the CILK_CFLAGS/CILK_CXXFLAGS/CILK_F77FLAGS
#   output variable to the flag (e.g. -fcilkplus) used both to compile
#   *and* link Cilk programs in the current language.
#   It may also set CILK_LIBS, if needed.
#
#   NOTE: You are assumed to not only compile your program with these flags,
#   but also link it with them as well.
#
#   If you want to compile everything with Cilk, you should set:
#
#     CFLAGS="$CFLAGS $CILK_CFLAGS"
#     #OR#  CXXFLAGS="$CXXFLAGS $CILK_CXXFLAGS"
#     #OR#  FFLAGS="$FFLAGS $CILK_FFLAGS"
#
#   (depending on the selected language).
#
#   The user can override the default choice by setting the corresponding
#   environment variable (e.g. CILK_CFLAGS).
#
#   ACTION-IF-FOUND is a list of shell commands to run if an Cilk flag is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_CILK.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 10

AC_DEFUN([PSL_CILK], [
AC_PREREQ([2.69]) dnl for _AC_LANG_PREFIX

AC_CACHE_CHECK([for Cilk flag of _AC_LANG compiler], psl_cv_[]_AC_LANG_ABBREV[]_cilk, [save[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
psl_cv_[]_AC_LANG_ABBREV[]_cilk=unknown
# Flags to try:  -fcilkplus (gcc), none
psl_cilk_flags="-fcilkplus none"
if test "x$CILK_[]_AC_LANG_PREFIX[]FLAGS" != x; then
  psl_cilk_flags="$CILK_[]_AC_LANG_PREFIX[]FLAGS $psl_cilk_flags"
fi
for psl_cilk_flag in $psl_cilk_flags; do
  case $psl_cilk_flag in
    none) []_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[] ;;
    *) []_AC_LANG_PREFIX[]FLAGS="$save[]_AC_LANG_PREFIX[]FLAGS $psl_cilk_flag" ;;
  esac
  try_again=0
  AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <stdio.h>
#include <cilk/cilk.h>

static void hello(){
    int i=0;
    for(i=0;i<1000000;i++)
        printf("");
    printf("Hello ");
}

static void world(){
    int i=0;
    for(i=0;i<1000000;i++)
        printf("");
    printf("world! ");
}

int main(){
    cilk_spawn hello();
    cilk_spawn world();
    printf("Done! ");
} 
]])],[psl_cv_[]_AC_LANG_ABBREV[]_cilk=$psl_cilk_flag; break],[try_again=1])
  AS_IF([test "x$try_again" = x1],[
  save_LIBS="$LIBS"
  LIBS="$LIBS -lcilkrts"
  AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <stdio.h>
#include <cilk/cilk.h>

static void hello(){
    int i=0;
    for(i=0;i<1000000;i++)
        printf("");
    printf("Hello ");
}

static void world(){
    int i=0;
    for(i=0;i<1000000;i++)
        printf("");
    printf("world! ");
}

int main(){
    cilk_spawn hello();
    cilk_spawn world();
    printf("Done! ");
    cilk_for(long i=0;i<1000000;i++)
        printf("");
} 
]])],[psl_cv_[]_AC_LANG_ABBREV[]_cilk=$psl_cilk_flag; CILK_LIBS=-lcilkrts],[])
  LIBS="$save_LIBS"
  if test "x$psl_cv_[]_AC_LANG_ABBREV[]_cilk" != "xunknown"; then
    break
  fi
  ])
done
[]_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[]FLAGS
])
if test "x$psl_cv_[]_AC_LANG_ABBREV[]_cilk" = "xunknown"; then
  m4_default([$2],:)
else
  if test "x$psl_cv_[]_AC_LANG_ABBREV[]_cilk" != "xnone"; then
    CILK_[]_AC_LANG_PREFIX[]FLAGS=$psl_cv_[]_AC_LANG_ABBREV[]_cilk
  fi
  m4_default([$1], [AC_DEFINE(HAVE_CILK,1,[Define if Cilk is enabled])])
fi
])dnl PSL_CILK
