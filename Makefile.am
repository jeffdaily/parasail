ACLOCAL_AMFLAGS = -I m4

include_HEADERS = src/parasail.h
lib_LTLIBRARIES = libparasail.la
noinst_LTLIBRARIES =
check_PROGRAMS =
AM_CFLAGS =
AM_CPPFLAGS = -I$(top_srcdir)/src
AM_LDFLAGS =
LDADD = libparasail.la

noinst_LTLIBRARIES += libparasail_memory_sse2.la
noinst_LTLIBRARIES += libparasail_memory_sse41.la
noinst_LTLIBRARIES += libparasail_memory_avx2.la

noinst_LTLIBRARIES += libparasail_novec.la
noinst_LTLIBRARIES += libparasail_novec_table.la
noinst_LTLIBRARIES += libparasail_sse2.la
noinst_LTLIBRARIES += libparasail_sse2_table.la
noinst_LTLIBRARIES += libparasail_sse41.la
noinst_LTLIBRARIES += libparasail_sse41_table.la
noinst_LTLIBRARIES += libparasail_avx2.la
noinst_LTLIBRARIES += libparasail_avx2_table.la
noinst_LTLIBRARIES += libparasail_knc.la
noinst_LTLIBRARIES += libparasail_knc_table.la

libparasail_la_SOURCES =
libparasail_la_SOURCES += src/parasail_cpuid.h
libparasail_la_SOURCES += src/parasail_internal.h
libparasail_la_SOURCES += src/parasail_internal_sse.h
libparasail_la_SOURCES += src/parasail_internal_avx.h
libparasail_la_SOURCES += src/parasail_internal_knc.h
libparasail_la_SOURCES += src/cpuid.c
libparasail_la_SOURCES += src/memory.c

MEM_SRC_SSE2 =
MEM_SRC_AVX2 =
MEM_SRC_KNC =
SRC_NOVEC =
SRC_SSE2 =
SRC_SSE41 =
SRC_AVX2 =
SRC_KNC =

#################
# vector memory
#################

if HAVE_SSE2
MEM_SRC_SSE2 += src/memory_sse2.c
endif

if HAVE_AVX2
MEM_SRC_AVX2 += src/memory_avx2.c
endif

if HAVE_KNC
MEM_SRC_KNC += src/memory_knc.c
endif

#################
# serial methods
#################

SRC_NOVEC += src/nw.c
SRC_NOVEC += src/sg.c
SRC_NOVEC += src/sw.c

SRC_NOVEC += src/nw_scan.c
SRC_NOVEC += src/sg_scan.c
SRC_NOVEC += src/sw_scan.c

SRC_NOVEC += src/nw_stats.c
SRC_NOVEC += src/sg_stats.c
SRC_NOVEC += src/sw_stats.c

SRC_NOVEC += src/nw_stats_scan.c
SRC_NOVEC += src/sg_stats_scan.c
SRC_NOVEC += src/sw_stats_scan.c

###################
# parallel methods
###################

# parallel scan methods

if HAVE_SSE2
SRC_SSE2 += src/nw_scan_sse2_128_32.c
SRC_SSE2 += src/sg_scan_sse2_128_32.c
SRC_SSE2 += src/sw_scan_sse2_128_32.c

SRC_SSE2 += src/nw_scan_sse2_128_16.c
SRC_SSE2 += src/sg_scan_sse2_128_16.c
SRC_SSE2 += src/sw_scan_sse2_128_16.c

SRC_SSE2 += src/nw_scan_sse2_128_8.c
SRC_SSE2 += src/sg_scan_sse2_128_8.c
SRC_SSE2 += src/sw_scan_sse2_128_8.c
endif

if HAVE_SSE41
SRC_SSE41 += src/nw_scan_sse41_128_32.c
SRC_SSE41 += src/sg_scan_sse41_128_32.c
SRC_SSE41 += src/sw_scan_sse41_128_32.c

SRC_SSE41 += src/nw_scan_sse41_128_16.c
SRC_SSE41 += src/sg_scan_sse41_128_16.c
SRC_SSE41 += src/sw_scan_sse41_128_16.c

SRC_SSE41 += src/nw_scan_sse41_128_8.c
SRC_SSE41 += src/sg_scan_sse41_128_8.c
SRC_SSE41 += src/sw_scan_sse41_128_8.c
endif

if HAVE_AVX2
SRC_AVX2 += src/nw_scan_avx2_256_32.c
SRC_AVX2 += src/sg_scan_avx2_256_32.c
SRC_AVX2 += src/sw_scan_avx2_256_32.c

SRC_AVX2 += src/nw_scan_avx2_256_16.c
SRC_AVX2 += src/sg_scan_avx2_256_16.c
SRC_AVX2 += src/sw_scan_avx2_256_16.c

SRC_AVX2 += src/nw_scan_avx2_256_8.c
SRC_AVX2 += src/sg_scan_avx2_256_8.c
SRC_AVX2 += src/sw_scan_avx2_256_8.c
endif

if HAVE_KNC
SRC_KNC += src/nw_scan_knc_512_32.c
SRC_KNC += src/sg_scan_knc_512_32.c
SRC_KNC += src/sw_scan_knc_512_32.c
endif

# parallel diag methods

if HAVE_SSE2
SRC_SSE2 += src/nw_diag_sse2_128_32.c
SRC_SSE2 += src/sg_diag_sse2_128_32.c
SRC_SSE2 += src/sw_diag_sse2_128_32.c

SRC_SSE2 += src/nw_diag_sse2_128_16.c
SRC_SSE2 += src/sg_diag_sse2_128_16.c
SRC_SSE2 += src/sw_diag_sse2_128_16.c

SRC_SSE2 += src/nw_diag_sse2_128_8.c
SRC_SSE2 += src/sg_diag_sse2_128_8.c
SRC_SSE2 += src/sw_diag_sse2_128_8.c
endif

if HAVE_SSE41
SRC_SSE41 += src/nw_diag_sse41_128_32.c
SRC_SSE41 += src/sg_diag_sse41_128_32.c
SRC_SSE41 += src/sw_diag_sse41_128_32.c

SRC_SSE41 += src/nw_diag_sse41_128_16.c
SRC_SSE41 += src/sg_diag_sse41_128_16.c
SRC_SSE41 += src/sw_diag_sse41_128_16.c

SRC_SSE41 += src/nw_diag_sse41_128_8.c
SRC_SSE41 += src/sg_diag_sse41_128_8.c
SRC_SSE41 += src/sw_diag_sse41_128_8.c
endif

if HAVE_AVX2
SRC_AVX2 += src/nw_diag_avx2_256_32.c
SRC_AVX2 += src/sg_diag_avx2_256_32.c
SRC_AVX2 += src/sw_diag_avx2_256_32.c

SRC_AVX2 += src/nw_diag_avx2_256_16.c
SRC_AVX2 += src/sg_diag_avx2_256_16.c
SRC_AVX2 += src/sw_diag_avx2_256_16.c

#SRC_AVX2 += src/nw_diag_avx2_256_8.c
#SRC_AVX2 += src/sg_diag_avx2_256_8.c
#SRC_AVX2 += src/sw_diag_avx2_256_8.c
endif

# parallel striped methods

if HAVE_SSE2
SRC_SSE2 += src/nw_striped_sse2_128_32.c
SRC_SSE2 += src/sg_striped_sse2_128_32.c
SRC_SSE2 += src/sw_striped_sse2_128_32.c

SRC_SSE2 += src/nw_striped_sse2_128_16.c
SRC_SSE2 += src/sg_striped_sse2_128_16.c
SRC_SSE2 += src/sw_striped_sse2_128_16.c

SRC_SSE2 += src/nw_striped_sse2_128_8.c
SRC_SSE2 += src/sg_striped_sse2_128_8.c
SRC_SSE2 += src/sw_striped_sse2_128_8.c
endif

if HAVE_SSE41
SRC_SSE41 += src/nw_striped_sse41_128_32.c
SRC_SSE41 += src/sg_striped_sse41_128_32.c
SRC_SSE41 += src/sw_striped_sse41_128_32.c

SRC_SSE41 += src/nw_striped_sse41_128_16.c
SRC_SSE41 += src/sg_striped_sse41_128_16.c
SRC_SSE41 += src/sw_striped_sse41_128_16.c

SRC_SSE41 += src/nw_striped_sse41_128_8.c
SRC_SSE41 += src/sg_striped_sse41_128_8.c
SRC_SSE41 += src/sw_striped_sse41_128_8.c
endif

if HAVE_AVX2
SRC_AVX2 += src/nw_striped_avx2_256_32.c
SRC_AVX2 += src/sg_striped_avx2_256_32.c
SRC_AVX2 += src/sw_striped_avx2_256_32.c

SRC_AVX2 += src/nw_striped_avx2_256_16.c
SRC_AVX2 += src/sg_striped_avx2_256_16.c
SRC_AVX2 += src/sw_striped_avx2_256_16.c

#SRC_AVX2 += src/nw_striped_avx2_256_8.c
#SRC_AVX2 += src/sg_striped_avx2_256_8.c
#SRC_AVX2 += src/sw_striped_avx2_256_8.c
endif

#########################
# parallel stats methods
#########################

# parallel scan methods
#
if HAVE_SSE2
SRC_SSE2 += src/nw_stats_scan_sse2_128_32.c
SRC_SSE2 += src/sg_stats_scan_sse2_128_32.c
SRC_SSE2 += src/sw_stats_scan_sse2_128_32.c

SRC_SSE2 += src/nw_stats_scan_sse2_128_16.c
SRC_SSE2 += src/sg_stats_scan_sse2_128_16.c
SRC_SSE2 += src/sw_stats_scan_sse2_128_16.c

SRC_SSE2 += src/nw_stats_scan_sse2_128_8.c
SRC_SSE2 += src/sg_stats_scan_sse2_128_8.c
SRC_SSE2 += src/sw_stats_scan_sse2_128_8.c
endif

if HAVE_SSE41
SRC_SSE41 += src/nw_stats_scan_sse41_128_32.c
SRC_SSE41 += src/sg_stats_scan_sse41_128_32.c
SRC_SSE41 += src/sw_stats_scan_sse41_128_32.c

SRC_SSE41 += src/nw_stats_scan_sse41_128_16.c
SRC_SSE41 += src/sg_stats_scan_sse41_128_16.c
SRC_SSE41 += src/sw_stats_scan_sse41_128_16.c

SRC_SSE41 += src/nw_stats_scan_sse41_128_8.c
SRC_SSE41 += src/sg_stats_scan_sse41_128_8.c
SRC_SSE41 += src/sw_stats_scan_sse41_128_8.c
endif

if HAVE_AVX2
SRC_AVX2 += src/nw_stats_scan_avx2_256_32.c
SRC_AVX2 += src/sg_stats_scan_avx2_256_32.c
SRC_AVX2 += src/sw_stats_scan_avx2_256_32.c

SRC_AVX2 += src/nw_stats_scan_avx2_256_16.c
SRC_AVX2 += src/sg_stats_scan_avx2_256_16.c
SRC_AVX2 += src/sw_stats_scan_avx2_256_16.c

SRC_AVX2 += src/nw_stats_scan_avx2_256_8.c
SRC_AVX2 += src/sg_stats_scan_avx2_256_8.c
SRC_AVX2 += src/sw_stats_scan_avx2_256_8.c
endif

# parallel diag methods

if HAVE_SSE2
SRC_SSE2 += src/nw_stats_diag_sse2_128_32.c
SRC_SSE2 += src/sg_stats_diag_sse2_128_32.c
SRC_SSE2 += src/sw_stats_diag_sse2_128_32.c

SRC_SSE2 += src/nw_stats_diag_sse2_128_16.c
SRC_SSE2 += src/sg_stats_diag_sse2_128_16.c
SRC_SSE2 += src/sw_stats_diag_sse2_128_16.c

SRC_SSE2 += src/nw_stats_diag_sse2_128_8.c
SRC_SSE2 += src/sg_stats_diag_sse2_128_8.c
SRC_SSE2 += src/sw_stats_diag_sse2_128_8.c
endif

if HAVE_SSE41
SRC_SSE41 += src/nw_stats_diag_sse41_128_32.c
SRC_SSE41 += src/sg_stats_diag_sse41_128_32.c
SRC_SSE41 += src/sw_stats_diag_sse41_128_32.c

SRC_SSE41 += src/nw_stats_diag_sse41_128_16.c
SRC_SSE41 += src/sg_stats_diag_sse41_128_16.c
SRC_SSE41 += src/sw_stats_diag_sse41_128_16.c

SRC_SSE41 += src/nw_stats_diag_sse41_128_8.c
SRC_SSE41 += src/sg_stats_diag_sse41_128_8.c
SRC_SSE41 += src/sw_stats_diag_sse41_128_8.c
endif

if HAVE_AVX2
#SRC_AVX2 += src/nw_stats_diag_avx2_256_32.c
#SRC_AVX2 += src/sg_stats_diag_avx2_256_32.c
#SRC_AVX2 += src/sw_stats_diag_avx2_256_32.c

#SRC_AVX2 += src/nw_stats_diag_avx2_256_16.c
#SRC_AVX2 += src/sg_stats_diag_avx2_256_16.c
#SRC_AVX2 += src/sw_stats_diag_avx2_256_16.c

#SRC_AVX2 += src/nw_stats_diag_avx2_256_8.c
#SRC_AVX2 += src/sg_stats_diag_avx2_256_8.c
#SRC_AVX2 += src/sw_stats_diag_avx2_256_8.c
endif

# parallel striped methods

if HAVE_SSE2
SRC_SSE2 += src/nw_stats_striped_sse2_128_32.c
SRC_SSE2 += src/sg_stats_striped_sse2_128_32.c
SRC_SSE2 += src/sw_stats_striped_sse2_128_32.c

SRC_SSE2 += src/nw_stats_striped_sse2_128_16.c
SRC_SSE2 += src/sg_stats_striped_sse2_128_16.c
SRC_SSE2 += src/sw_stats_striped_sse2_128_16.c

SRC_SSE2 += src/nw_stats_striped_sse2_128_8.c
SRC_SSE2 += src/sg_stats_striped_sse2_128_8.c
SRC_SSE2 += src/sw_stats_striped_sse2_128_8.c
endif

if HAVE_SSE41
SRC_SSE41 += src/nw_stats_striped_sse41_128_32.c
SRC_SSE41 += src/sg_stats_striped_sse41_128_32.c
SRC_SSE41 += src/sw_stats_striped_sse41_128_32.c

SRC_SSE41 += src/nw_stats_striped_sse41_128_16.c
SRC_SSE41 += src/sg_stats_striped_sse41_128_16.c
SRC_SSE41 += src/sw_stats_striped_sse41_128_16.c

SRC_SSE41 += src/nw_stats_striped_sse41_128_8.c
SRC_SSE41 += src/sg_stats_striped_sse41_128_8.c
SRC_SSE41 += src/sw_stats_striped_sse41_128_8.c
endif

if HAVE_AVX2
#SRC_AVX2 += src/nw_stats_striped_avx2_256_32.c
#SRC_AVX2 += src/sg_stats_striped_avx2_256_32.c
#SRC_AVX2 += src/sw_stats_striped_avx2_256_32.c

#SRC_AVX2 += src/nw_stats_striped_avx2_256_16.c
#SRC_AVX2 += src/sg_stats_striped_avx2_256_16.c
#SRC_AVX2 += src/sw_stats_striped_avx2_256_16.c

#SRC_AVX2 += src/nw_stats_striped_avx2_256_8.c
#SRC_AVX2 += src/sg_stats_striped_avx2_256_8.c
#SRC_AVX2 += src/sw_stats_striped_avx2_256_8.c
endif

libparasail_novec_la_SOURCES = $(SRC_NOVEC)
libparasail_sse2_la_SOURCES = $(SRC_SSE2)
libparasail_sse41_la_SOURCES = $(SRC_SSE41)
libparasail_avx2_la_SOURCES = $(SRC_AVX2)
libparasail_knc_la_SOURCES = $(SRC_KNC)

libparasail_novec_la_CFLAGS = $(AM_CFLAGS)
libparasail_sse2_la_CFLAGS = $(AM_CFLAGS) $(SSE2_CFLAGS)
libparasail_sse41_la_CFLAGS = $(AM_CFLAGS) $(SSE41_CFLAGS)
libparasail_avx2_la_CFLAGS = $(AM_CFLAGS) $(AVX2_CFLAGS)
libparasail_knc_la_CFLAGS = $(AM_CFLAGS)

libparasail_novec_table_la_SOURCES = $(SRC_NOVEC)
libparasail_sse2_table_la_SOURCES = $(SRC_SSE2)
libparasail_sse41_table_la_SOURCES = $(SRC_SSE41)
libparasail_avx2_table_la_SOURCES = $(SRC_AVX2)
libparasail_knc_table_la_SOURCES = $(SRC_KNC)

libparasail_novec_table_la_CFLAGS = $(AM_CFLAGS)
libparasail_sse2_table_la_CFLAGS = $(AM_CFLAGS) $(SSE2_CFLAGS)
libparasail_sse41_table_la_CFLAGS = $(AM_CFLAGS) $(SSE41_CFLAGS)
libparasail_avx2_table_la_CFLAGS = $(AM_CFLAGS) $(AVX2_CFLAGS)
libparasail_knc_table_la_CFLAGS = $(AM_CFLAGS)

libparasail_novec_table_la_CPPFLAGS = $(AM_CPPFLAGS) -DPARASAIL_TABLE=1
libparasail_sse2_table_la_CPPFLAGS = $(AM_CPPFLAGS) -DPARASAIL_TABLE=1
libparasail_sse41_table_la_CPPFLAGS = $(AM_CPPFLAGS) -DPARASAIL_TABLE=1
libparasail_avx2_table_la_CPPFLAGS = $(AM_CPPFLAGS) -DPARASAIL_TABLE=1
libparasail_knc_table_la_CPPFLAGS = $(AM_CPPFLAGS) -DPARASAIL_TABLE=1

libparasail_memory_sse2_la_SOURCES = $(MEM_SRC_SSE2)
libparasail_memory_sse41_la_SOURCES = $(MEM_SRC_SSE41)
libparasail_memory_avx2_la_SOURCES = $(MEM_SRC_AVX2)

libparasail_memory_sse2_la_CFLAGS = $(AM_CFLAGS) $(SSE2_CFLAGS)
libparasail_memory_sse41_la_CFLAGS = $(AM_CFLAGS) $(SSE41_CFLAGS)
libparasail_memory_avx2_la_CFLAGS = $(AM_CFLAGS) $(AVX2_CFLAGS)

libparasail_la_LIBADD =
libparasail_la_LIBADD += libparasail_memory_sse2.la
libparasail_la_LIBADD += libparasail_memory_sse41.la
libparasail_la_LIBADD += libparasail_memory_avx2.la
libparasail_la_LIBADD += libparasail_novec.la
libparasail_la_LIBADD += libparasail_novec_table.la
libparasail_la_LIBADD += libparasail_sse2.la
libparasail_la_LIBADD += libparasail_sse2_table.la
libparasail_la_LIBADD += libparasail_sse41.la
libparasail_la_LIBADD += libparasail_sse41_table.la
libparasail_la_LIBADD += libparasail_avx2.la
libparasail_la_LIBADD += libparasail_avx2_table.la
libparasail_la_LIBADD += libparasail_knc.la
libparasail_la_LIBADD += libparasail_knc_table.la

check_PROGRAMS += tests/test_align
check_PROGRAMS += tests/test_isa
if HAVE_OPENMP
check_PROGRAMS += tests/test_openmp
endif

tests_test_align_SOURCES        = tests/test_align.c
tests_test_align_LDADD          = $(LDADD) $(CLOCK_LIBS) $(MATH_LIBS)

tests_test_isa_SOURCES          = tests/test_isa.c

tests_test_openmp_SOURCES       = tests/test_openmp.c
tests_test_openmp_CFLAGS        = $(AM_CFLAGS) $(OPENMP_CFLAGS)
tests_test_openmp_LDFLAGS       = $(AM_LDFLAGS) $(OPENMP_CFLAGS)
tests_test_openmp_LDADD         = $(LDADD) $(CLOCK_LIBS) $(MATH_LIBS)

