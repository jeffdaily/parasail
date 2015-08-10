#.rst:
# FindSSE41
# ---------
#
# Finds SSE41 support
#
# This module can be used to detect SSE41 support in a C/C++ compiler.  If
# the compiler supports SSE41, the flags required to compile with
# SSE41 support are returned in variables for the different languages.
# The variables may be empty if the compiler does not need a special
# flag to support SSE41.
#
# The following variables are set:
#
# ::
#
#    SSE41_C_FLAGS - flags to add to the C compiler for SSE41 support
#    SSE41_CXX_FLAGS - flags to add to the CXX compiler for SSE41 support
#    SSE41_FOUND - true if SSE41 is detected
#
#=============================================================================

set(_SSE41_REQUIRED_VARS)
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set(CMAKE_REQUIRED_QUIET ${SSE41_FIND_QUIETLY})

function(_SSE41_FLAG_CANDIDATES LANG)
  set(SSE41_FLAG_CANDIDATES
    #Empty, if compiler automatically accepts SSE41
    " "
    #GNU, Intel
    "-march=corei7"
    #clang
    "-msse4"
  )

  set(SSE41_${LANG}_FLAG_CANDIDATES "${SSE41_FLAG_CANDIDATES}" PARENT_SCOPE)
endfunction()

# sample SSE41 source code to test
set(SSE41_C_TEST_SOURCE
"
#include <smmintrin.h>
int foo() {
    __m128i vOne = _mm_set1_epi8(1);
    __m128i result =  _mm_max_epi8(vOne,vOne);
    __m128i result2 =  _mm_insert_epi64(result,0,0);
    return _mm_extract_epi8(result2, 0);
}
int main(void) { return foo(); }
")

# check c compiler
if(CMAKE_C_COMPILER_LOADED)
  # if these are set then do not try to find them again,
  # by avoiding any try_compiles for the flags
  if(SSE41_C_FLAGS)
    unset(SSE41_C_FLAG_CANDIDATES)
  else()
    _SSE41_FLAG_CANDIDATES("C")
    include(CheckCSourceCompiles)
  endif()

  foreach(FLAG IN LISTS SSE41_C_FLAG_CANDIDATES)
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(SSE41_FLAG_DETECTED CACHE)
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Try SSE41 C flag = [${FLAG}]")
    endif()
    check_c_source_compiles("${SSE41_C_TEST_SOURCE}" SSE41_FLAG_DETECTED)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    if(SSE41_FLAG_DETECTED)
      set(SSE41_C_FLAGS_INTERNAL "${FLAG}")
      break()
    endif()
  endforeach()

  set(SSE41_C_FLAGS "${SSE41_C_FLAGS_INTERNAL}"
    CACHE STRING "C compiler flags for SSE41 intrinsics")

  list(APPEND _SSE41_REQUIRED_VARS SSE41_C_FLAGS)
  unset(SSE41_C_FLAG_CANDIDATES)
endif()

# check cxx compiler
if(CMAKE_CXX_COMPILER_LOADED)
  # if these are set then do not try to find them again,
  # by avoiding any try_compiles for the flags
  if(SSE41_CXX_FLAGS)
    unset(SSE41_CXX_FLAG_CANDIDATES)
  else()
    _SSE41_FLAG_CANDIDATES("CXX")
    include(CheckCXXSourceCompiles)

    # use the same source for CXX as C for now
    set(SSE41_CXX_TEST_SOURCE ${SSE41_C_TEST_SOURCE})
  endif()

  foreach(FLAG IN LISTS SSE41_CXX_FLAG_CANDIDATES)
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(SSE41_FLAG_DETECTED CACHE)
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Try SSE41 CXX flag = [${FLAG}]")
    endif()
    check_cxx_source_compiles("${SSE41_CXX_TEST_SOURCE}" SSE41_FLAG_DETECTED)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    if(SSE41_FLAG_DETECTED)
      set(SSE41_CXX_FLAGS_INTERNAL "${FLAG}")
      break()
    endif()
  endforeach()

  set(SSE41_CXX_FLAGS "${SSE41_CXX_FLAGS_INTERNAL}"
    CACHE STRING "C++ compiler flags for SSE41 parallization")

  list(APPEND _SSE41_REQUIRED_VARS SSE41_CXX_FLAGS)
  unset(SSE41_CXX_FLAG_CANDIDATES)
  unset(SSE41_CXX_TEST_SOURCE)
endif()

set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if(_SSE41_REQUIRED_VARS)
  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args(SSE41
                                    REQUIRED_VARS ${_SSE41_REQUIRED_VARS})

  mark_as_advanced(${_SSE41_REQUIRED_VARS})

  unset(_SSE41_REQUIRED_VARS)
else()
  message(SEND_ERROR "FindSSE41 requires C or CXX language to be enabled")
endif()
