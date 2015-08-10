#.rst:
# FindAVX2
# --------
#
# Finds AVX2 support
#
# This module can be used to detect AVX2 support in a C/C++ compiler.  If
# the compiler supports AVX2, the flags required to compile with
# AVX2 support are returned in variables for the different languages.
# The variables may be empty if the compiler does not need a special
# flag to support AVX2.
#
# The following variables are set:
#
# ::
#
#    AVX2_C_FLAGS - flags to add to the C compiler for AVX2 support
#    AVX2_CXX_FLAGS - flags to add to the CXX compiler for AVX2 support
#    AVX2_FOUND - true if AVX2 is detected
#
#=============================================================================

set(_AVX2_REQUIRED_VARS)
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set(CMAKE_REQUIRED_QUIET ${AVX2_FIND_QUIETLY})

function(_AVX2_FLAG_CANDIDATES LANG)
  set(AVX2_FLAG_CANDIDATES
    #Empty, if compiler automatically accepts AVX2
    " "
    #GNU, Intel
    "-march=core-avx2"
    #clang
    "-mavx2"
  )

  set(AVX2_${LANG}_FLAG_CANDIDATES "${AVX2_FLAG_CANDIDATES}" PARENT_SCOPE)
endfunction()

# sample AVX2 source code to test
set(AVX2_C_TEST_SOURCE
"
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
    __m256i result2 =  _mm256_insert_epi64(vOne,0,0);
    parasail_memset___m256i(&result2, result, 1);
    return _mm_extract_epi16(_mm256_extracti128_si256(result,0),0);
}
int main(void) { return foo(); }
")

# check c compiler
if(CMAKE_C_COMPILER_LOADED)
  # if these are set then do not try to find them again,
  # by avoiding any try_compiles for the flags
  if(AVX2_C_FLAGS)
    unset(AVX2_C_FLAG_CANDIDATES)
  else()
    _AVX2_FLAG_CANDIDATES("C")
    include(CheckCSourceCompiles)
  endif()

  foreach(FLAG IN LISTS AVX2_C_FLAG_CANDIDATES)
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(AVX2_FLAG_DETECTED CACHE)
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Try AVX2 C flag = [${FLAG}]")
    endif()
    check_c_source_compiles("${AVX2_C_TEST_SOURCE}" AVX2_FLAG_DETECTED)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    if(AVX2_FLAG_DETECTED)
      set(AVX2_C_FLAGS_INTERNAL "${FLAG}")
      break()
    endif()
  endforeach()

  set(AVX2_C_FLAGS "${AVX2_C_FLAGS_INTERNAL}"
    CACHE STRING "C compiler flags for AVX2 intrinsics")

  list(APPEND _AVX2_REQUIRED_VARS AVX2_C_FLAGS)
  unset(AVX2_C_FLAG_CANDIDATES)
endif()

# check cxx compiler
if(CMAKE_CXX_COMPILER_LOADED)
  # if these are set then do not try to find them again,
  # by avoiding any try_compiles for the flags
  if(AVX2_CXX_FLAGS)
    unset(AVX2_CXX_FLAG_CANDIDATES)
  else()
    _AVX2_FLAG_CANDIDATES("CXX")
    include(CheckCXXSourceCompiles)

    # use the same source for CXX as C for now
    set(AVX2_CXX_TEST_SOURCE ${AVX2_C_TEST_SOURCE})
  endif()

  foreach(FLAG IN LISTS AVX2_CXX_FLAG_CANDIDATES)
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(AVX2_FLAG_DETECTED CACHE)
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Try AVX2 CXX flag = [${FLAG}]")
    endif()
    check_cxx_source_compiles("${AVX2_CXX_TEST_SOURCE}" AVX2_FLAG_DETECTED)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    if(AVX2_FLAG_DETECTED)
      set(AVX2_CXX_FLAGS_INTERNAL "${FLAG}")
      break()
    endif()
  endforeach()

  set(AVX2_CXX_FLAGS "${AVX2_CXX_FLAGS_INTERNAL}"
    CACHE STRING "C++ compiler flags for AVX2 parallization")

  list(APPEND _AVX2_REQUIRED_VARS AVX2_CXX_FLAGS)
  unset(AVX2_CXX_FLAG_CANDIDATES)
  unset(AVX2_CXX_TEST_SOURCE)
endif()

set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if(_AVX2_REQUIRED_VARS)
  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args(AVX2
                                    REQUIRED_VARS ${_AVX2_REQUIRED_VARS})

  mark_as_advanced(${_AVX2_REQUIRED_VARS})

  unset(_AVX2_REQUIRED_VARS)
else()
  message(SEND_ERROR "FindAVX2 requires C or CXX language to be enabled")
endif()
