# PS_COMPILER_VENDOR
# ------------------
# other compilers may define __GNUC__, like intel or pathscale
# that's why we put it later in the for loop, like we do with our compiler
# checks since GCC is so pervasive
AC_DEFUN([PS_COMPILER_VENDOR], [
AS_VAR_PUSHDEF([ps_cv_compiler_vendor],
               [ps_cv_[]_AC_LANG_ABBREV[]_compiler_vendor])
AC_CACHE_CHECK([for _AC_LANG compiler vendor], [ps_cv_compiler_vendor], [
ps_save_ac_ext="$ac_ext"
AC_LANG_CASE([Fortran],    [ac_ext=F])
AC_LANG_CASE([Fortran 77], [ac_ext=F])
ps_cv_compiler_vendor=unknown
ps_cpp_vendor_symbols=
for vendor in intel ibm pathscale amd cray clang gnu sun hp dec borland comeau kai lcc metrowerks sgi microsoft watcom portland fujitsu
do
AS_CASE([$vendor],
[amd],       [ps_cpp_vendor_symbols="defined(__OPEN64__)"],
[borland],   [ps_cpp_vendor_symbols="defined(__BORLANDC__) || defined(__TURBOC__)"],
[clang],     [ps_cpp_vendor_symbols="defined(__clang__)"],
[comeau],    [ps_cpp_vendor_symbols="defined(__COMO__)"],
[cray],      [ps_cpp_vendor_symbols="defined(_CRAYC) || defined(_ADDR64)"],
[dec],       [ps_cpp_vendor_symbols="defined(__DECC) || defined(__DECCXX) || defined(__DECC_VER) || defined(__DECCXX_VER)"],
[fujitsu],   [ps_cpp_vendor_symbols="defined(__fcc__) || defined(__fcc_version__) || defined(_FCC_VER) || defined(__FCC_VER_)"],
[gnu],       [ps_cpp_vendor_symbols="defined(__GNUC__)"],
[hp],        [ps_cpp_vendor_symbols="defined(__HP_cc) || defined(__HP_aCC)"],
[ibm],       [ps_cpp_vendor_symbols="defined(__xlc__) || defined(__xlC__) || defined(__IBMC__) || defined(__IBMCPP__)"],
[intel],     [ps_cpp_vendor_symbols="defined(__ICC) || defined(__ECC) || defined(__INTEL_COMPILER)"],
[kai],       [ps_cpp_vendor_symbols="defined(__KCC)"],
[lcc],       [ps_cpp_vendor_symbols="defined(__LCC__)"],
[metrowerks],[ps_cpp_vendor_symbols="defined(__MWERKS__)"],
[microsoft], [ps_cpp_vendor_symbols="defined(_MSC_VER)"],
[pathscale], [ps_cpp_vendor_symbols="defined(__PATHCC__) || defined(__PATHSCALE__)"],
[portland],  [ps_cpp_vendor_symbols="defined(__PGI)"],
[sgi],       [ps_cpp_vendor_symbols="defined(__sgi) || defined(sgi)"],
[sun],       [ps_cpp_vendor_symbols="defined(__SUNPRO_C) || defined(__SUNPRO_CC)"],
[watcom],    [ps_cpp_vendor_symbols="defined(__WATCOMC__)"])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
#if !($ps_cpp_vendor_symbols)
    chokeonthis
#endif
])], [ps_cv_compiler_vendor=$vendor; break])
done
ps_cpp_vendor_symbols=
ac_ext="$ps_save_ac_ext"
])
AS_VAR_POPDEF([ps_cv_compiler_vendor])
])dnl
