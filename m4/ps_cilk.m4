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

AC_DEFUN([PSL_CILK_SOURCE], [AC_LANG_CONFTEST([AC_LANG_SOURCE(
[[#include <stdio.h>
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
}]])])])

AC_DEFUN([PSL_CILK], [
AC_PREREQ([2.69]) dnl for _AC_LANG_PREFIX

AC_CACHE_CHECK([for Cilk flag of _AC_LANG compiler], psl_cv_[]_AC_LANG_ABBREV[]_cilk, [
save[]_AC_LANG_ABBREV[]_werror_flag=$ac_[]_AC_LANG_ABBREV[]_werror_flag
ac_[]_AC_LANG_ABBREV[]_werror_flag=yes
save[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
psl_cv_[]_AC_LANG_ABBREV[]_cilk=unknown
# Flags to try:  -fcilkplus (gcc), none
psl_cilk_flags="-fcilkplus none"
AS_IF([test "x$CILK_[]_AC_LANG_PREFIX[]FLAGS" != x],
  [psl_cilk_flags="$CILK_[]_AC_LANG_PREFIX[]FLAGS $psl_cilk_flags"])
for psl_cilk_flag in $psl_cilk_flags; do
  AS_CASE([$psl_cilk_flag],
    [none], [[]_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[]],
            [[]_AC_LANG_PREFIX[]FLAGS="$save[]_AC_LANG_PREFIX[]FLAGS $psl_cilk_flag"])
  try_again=0
  PSL_CILK_SOURCE
  AC_LINK_IFELSE([],
    [psl_cv_[]_AC_LANG_ABBREV[]_cilk=$psl_cilk_flag; break],
    [try_again=1])
  AS_IF([test "x$try_again" = x1],[
    save_LIBS="$LIBS"
    LIBS="$LIBS -lcilkrts"
    PSL_CILK_SOURCE
    AC_LINK_IFELSE([],
      [psl_cv_[]_AC_LANG_ABBREV[]_cilk=$psl_cilk_flag; CILK_LIBS=-lcilkrts],
      [])
    LIBS="$save_LIBS"
    AS_IF([test "x$psl_cv_[]_AC_LANG_ABBREV[]_cilk" != "xunknown"], [break])])
done
[]_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[]FLAGS
ac_[]_AC_LANG_ABBREV[]_werror_flag=$save[]_AC_LANG_ABBREV[]_werror_flag
])
AS_IF([test "x$psl_cv_[]_AC_LANG_ABBREV[]_cilk" = "xunknown"],
  [m4_default([$2],:)],
  [AS_IF([test "x$psl_cv_[]_AC_LANG_ABBREV[]_cilk" != "xnone"],
    [CILK_[]_AC_LANG_PREFIX[]FLAGS=$psl_cv_[]_AC_LANG_ABBREV[]_cilk])
   m4_default([$1], [AC_DEFINE(HAVE_CILK,1,[Define if Cilk is enabled])])])
])dnl PSL_CILK

