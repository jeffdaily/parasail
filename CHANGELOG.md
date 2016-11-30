# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

This project follows the [Gitflow Workflow model](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow).

## [Unreleased]
The Unreleased section will be empty for tagged releases. Unreleased functionality appears in the develop branch.

## [1.1.1] - 2016-11-30

### Fixed
- libparasail now correctly links when pow() not in system C library

### Merged Pull Requests
- Allow injection via cmake submodule [\#27] ([armintoepfer])


## [1.1] - 2016-08-12

### Changed
- Stats functions are now affine, not linear.
- Semi-global and global alignments now use a more negative value to
  represent negative infinity instead of half the value of the smallest
  representable integer for the given bit width.
- end_query and end_ref reported for all routines.

### Fixed
- Stats functions are now affine, not linear.

### Closed Issues
- provide Java bindings [\#22]
- stats functions should be affine, not linear [\#10]
- parasail results off by one error [\#4]

## [1.0.3] - 2016-03-25

### Changed
- Added TravisCI support for autotools Linux and OSX builds.
- Added AppVeyor support for CMake Windows builds.
- PARASAIL_API and PARASAIL_LOCAL removed from all parasail functions.
- CMake build 
  - Added BUILD_SHARED_LIBS option.
  - Added parasail.def for MSVC DLL creation.
  - Set CMAKE_POSITION_INDEPENDENT_CODE to ON if BUILD_SHARED_LIBS is ON.
  - /arch:AVX is the correct flag for MSVC, not /arch:AVX2.

### Fixed
- parasail_free() was not being used to free ISA-specific sequence profiles. Caused MSVC 64-bit library to crash.
- CMake shared library build was basically not functional on any platform. It now works.

## [1.0.2] - 2016-03-17

### Changed
- 32-bit builds replace missing functionailty.
  - SSE2 _mm_set1_epi64x, _mm_set_epi64x
  - AVX2 _mm256_set1_epi64x, _mm256_set_epi64x

### Removed
- Python bindings and pygen.py generateor were removed. Now a stand-alone project [parasail-python].

### Fixed
- Multi-arch build for OSX now correctly detects SSE4.1 and AVX2 after fixing [\#20].

### Closed Issues
- -O3 optimization causes incorrect results on OSX clang for _mm256_blendv_epi8 [\#21]
- epi64 instructions not available on 32-bit platforms [\#20]
- python ctypes interface instead of cython [\#19] **wontfix** -- moved to [parasail-python] project.
- create python wheel for pip install [\#18] **wontfix** -- moved to [parasail-python] project.
- Adding example for python [\#18] **wontfix** -- moved to [parasail-python] project.

## [1.0.1] - 2016-03-01

### Changed
- Many improvements and bug fixes to the CMake build.
  - Needed to bump CMAKE_MINIMUM_REQUIRED to VERSION 3.1 to fix static linking.
  - Visual Studio, OSX, and Linux have been verified to work.
- Windows platform natively supported.
- If an instruction set, e.g., AVX2 is not detected, then the functions are stubbed out and return NULL and set errno to ENOSYS.
- restrict keyword is conditionally preprocessed away if it's not supported by the compiler (e.g., C++, C89).  parasail internally still uses restrict if there is a suitable extension (e.g., __restrict) but this change allows greater flexibility for external libraries and applications.
- parasail_aligner application now uses long instead of int for indexing. This supports larger input datasets.

### Fixed
- Changed C++ style comments to C style to support MSVC build.
- Corrected mixed declarations and code to support MSVC build.
- Fixed various warnings from gcc -Wall -Wextra, clang, icc. MSVC build still produces many warnings.

### Closed Issues
- incorrect default SSE41_CFLAGS for gcc 4.4.7 [\#17]
- test_isa should also report what the compiler supported [\#15]
- update README et al. for new citation [\#13]
- Adding flag to disable/enable binaries in CMakeLists.txt [\#9]
- Profile thread safety? [\#8]
- Can't get parasail\_aligner to use \> 1 thread [\#7]
- Missing \#include \<string.h\> in tests [\#6]
- Documentation [\#5]
- AVX2: no such instruction [\#1](https://github.com/jeffdaily/parasail/issues/1)

## [1.0.0] - 2015-09-16
First stable, production-ready version of parasail.

[parasail-python]: https://github.com/jeffdaily/parasail-python

[Unreleased]: https://github.com/jeffdaily/parasail/compare/v1.1.1...develop
[1.1.1]: https://github.com/jeffdaily/parasail/compare/v1.1...v1.1.1
[1.1]:   https://github.com/jeffdaily/parasail/compare/v1.0.3...v1.1
[1.0.3]: https://github.com/jeffdaily/parasail/compare/v1.0.2...v1.0.3
[1.0.2]: https://github.com/jeffdaily/parasail/compare/v1.0.1...v1.0.2
[1.0.1]: https://github.com/jeffdaily/parasail/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/jeffdaily/parasail/releases/tag/v1.0.0

[\#28]: https://github.com/jeffdaily/parasail/issues/28
[\#27]: https://github.com/jeffdaily/parasail/pull/27
[\#26]: https://github.com/jeffdaily/parasail/issues/26
[\#25]: https://github.com/jeffdaily/parasail/pull/25
[\#24]: https://github.com/jeffdaily/parasail/issues/24
[\#23]: https://github.com/jeffdaily/parasail/issues/23
[\#22]: https://github.com/jeffdaily/parasail/issues/22
[\#21]: https://github.com/jeffdaily/parasail/issues/21
[\#20]: https://github.com/jeffdaily/parasail/issues/20
[\#19]: https://github.com/jeffdaily/parasail/issues/19
[\#18]: https://github.com/jeffdaily/parasail/issues/18
[\#17]: https://github.com/jeffdaily/parasail/issues/17
[\#16]: https://github.com/jeffdaily/parasail/issues/16
[\#15]: https://github.com/jeffdaily/parasail/issues/15
[\#14]: https://github.com/jeffdaily/parasail/issues/14
[\#13]: https://github.com/jeffdaily/parasail/issues/13
[\#12]: https://github.com/jeffdaily/parasail/issues/12
[\#11]: https://github.com/jeffdaily/parasail/issues/11
[\#10]: https://github.com/jeffdaily/parasail/issues/10
[\#9]: https://github.com/jeffdaily/parasail/issues/9
[\#8]: https://github.com/jeffdaily/parasail/issues/8
[\#7]: https://github.com/jeffdaily/parasail/issues/7
[\#6]: https://github.com/jeffdaily/parasail/issues/6
[\#5]: https://github.com/jeffdaily/parasail/issues/5
[\#4]: https://github.com/jeffdaily/parasail/issues/4
[\#3]: https://github.com/jeffdaily/parasail/issues/3
[\#2]: https://github.com/jeffdaily/parasail/issues/2
[\#1]: https://github.com/jeffdaily/parasail/issues/1

[armintoepfer]: https://github.com/armintoepfer