ACLOCAL_AMFLAGS = -I m4

include_HEADERS = src/parasail.h
lib_LTLIBRARIES = libparasail.la
noinst_LTLIBRARIES = libparasail_table.la
check_PROGRAMS =
AM_CFLAGS =
AM_CPPFLAGS = -I$(top_srcdir)/src
AM_LDFLAGS =
LDADD = libparasail.la

libparasail_la_SOURCES =
libparasail_la_SOURCES += src/parasail_internal.h
libparasail_la_SOURCES += src/parasail_internal_sse.h
libparasail_la_SOURCES += src/parasail_internal_avx.h
libparasail_la_SOURCES += src/memory.c

ALIGNMENT_SRC = 

#################
# serial methods
#################

ALIGNMENT_SRC += src/nw.c
ALIGNMENT_SRC += src/sg.c
ALIGNMENT_SRC += src/sw.c

ALIGNMENT_SRC += src/nw_scan.c
ALIGNMENT_SRC += src/sg_scan.c
ALIGNMENT_SRC += src/sw_scan.c

ALIGNMENT_SRC += src/nw_stats.c
ALIGNMENT_SRC += src/sg_stats.c
ALIGNMENT_SRC += src/sw_stats.c

ALIGNMENT_SRC += src/nw_stats_scan.c
ALIGNMENT_SRC += src/sg_stats_scan.c
ALIGNMENT_SRC += src/sw_stats_scan.c

###################
# parallel methods
###################

# parallel scan methods
#
if HAVE_SSE2
ALIGNMENT_SRC += src/nw_scan_sse2_128_32.c
ALIGNMENT_SRC += src/sg_scan_sse2_128_32.c
ALIGNMENT_SRC += src/sw_scan_sse2_128_32.c

ALIGNMENT_SRC += src/nw_scan_sse2_128_16.c
ALIGNMENT_SRC += src/sg_scan_sse2_128_16.c
ALIGNMENT_SRC += src/sw_scan_sse2_128_16.c

ALIGNMENT_SRC += src/nw_scan_sse2_128_8.c
ALIGNMENT_SRC += src/sg_scan_sse2_128_8.c
ALIGNMENT_SRC += src/sw_scan_sse2_128_8.c
endif

if HAVE_SSE41
ALIGNMENT_SRC += src/nw_scan_sse41_128_32.c
ALIGNMENT_SRC += src/sg_scan_sse41_128_32.c
ALIGNMENT_SRC += src/sw_scan_sse41_128_32.c

ALIGNMENT_SRC += src/nw_scan_sse41_128_16.c
ALIGNMENT_SRC += src/sg_scan_sse41_128_16.c
ALIGNMENT_SRC += src/sw_scan_sse41_128_16.c

ALIGNMENT_SRC += src/nw_scan_sse41_128_8.c
ALIGNMENT_SRC += src/sg_scan_sse41_128_8.c
ALIGNMENT_SRC += src/sw_scan_sse41_128_8.c
endif

if HAVE_AVX2
ALIGNMENT_SRC += src/nw_scan_avx2_256_32.c
#ALIGNMENT_SRC += src/sg_scan_avx2_256_32.c
#ALIGNMENT_SRC += src/sw_scan_avx2_256_32.c

#ALIGNMENT_SRC += src/nw_scan_avx2_256_16.c
#ALIGNMENT_SRC += src/sg_scan_avx2_256_16.c
#ALIGNMENT_SRC += src/sw_scan_avx2_256_16.c

#ALIGNMENT_SRC += src/nw_scan_avx2_256_8.c
#ALIGNMENT_SRC += src/sg_scan_avx2_256_8.c
#ALIGNMENT_SRC += src/sw_scan_avx2_256_8.c
endif

# parallel diag methods

if HAVE_SSE2
ALIGNMENT_SRC += src/nw_diag_sse2_128_32.c
ALIGNMENT_SRC += src/sg_diag_sse2_128_32.c
ALIGNMENT_SRC += src/sw_diag_sse2_128_32.c

ALIGNMENT_SRC += src/nw_diag_sse2_128_16.c
ALIGNMENT_SRC += src/sg_diag_sse2_128_16.c
ALIGNMENT_SRC += src/sw_diag_sse2_128_16.c

ALIGNMENT_SRC += src/nw_diag_sse2_128_8.c
ALIGNMENT_SRC += src/sg_diag_sse2_128_8.c
ALIGNMENT_SRC += src/sw_diag_sse2_128_8.c
endif

if HAVE_SSE41
ALIGNMENT_SRC += src/nw_diag_sse41_128_32.c
ALIGNMENT_SRC += src/sg_diag_sse41_128_32.c
ALIGNMENT_SRC += src/sw_diag_sse41_128_32.c

ALIGNMENT_SRC += src/nw_diag_sse41_128_16.c
ALIGNMENT_SRC += src/sg_diag_sse41_128_16.c
ALIGNMENT_SRC += src/sw_diag_sse41_128_16.c

ALIGNMENT_SRC += src/nw_diag_sse41_128_8.c
ALIGNMENT_SRC += src/sg_diag_sse41_128_8.c
ALIGNMENT_SRC += src/sw_diag_sse41_128_8.c
endif

# parallel striped methods

if HAVE_SSE2
ALIGNMENT_SRC += src/nw_striped_sse2_128_32.c
ALIGNMENT_SRC += src/sg_striped_sse2_128_32.c
ALIGNMENT_SRC += src/sw_striped_sse2_128_32.c

ALIGNMENT_SRC += src/nw_striped_sse2_128_16.c
ALIGNMENT_SRC += src/sg_striped_sse2_128_16.c
ALIGNMENT_SRC += src/sw_striped_sse2_128_16.c

ALIGNMENT_SRC += src/nw_striped_sse2_128_8.c
ALIGNMENT_SRC += src/sg_striped_sse2_128_8.c
ALIGNMENT_SRC += src/sw_striped_sse2_128_8.c
endif

if HAVE_SSE41
ALIGNMENT_SRC += src/nw_striped_sse41_128_32.c
ALIGNMENT_SRC += src/sg_striped_sse41_128_32.c
ALIGNMENT_SRC += src/sw_striped_sse41_128_32.c

ALIGNMENT_SRC += src/nw_striped_sse41_128_16.c
ALIGNMENT_SRC += src/sg_striped_sse41_128_16.c
ALIGNMENT_SRC += src/sw_striped_sse41_128_16.c

ALIGNMENT_SRC += src/nw_striped_sse41_128_8.c
ALIGNMENT_SRC += src/sg_striped_sse41_128_8.c
ALIGNMENT_SRC += src/sw_striped_sse41_128_8.c
endif

#########################
# parallel stats methods
#########################

# parallel scan methods
#
if HAVE_SSE2
ALIGNMENT_SRC += src/nw_stats_scan_sse2_128_32.c
#ALIGNMENT_SRC += src/sg_stats_scan_sse2_128_32.c
#ALIGNMENT_SRC += src/sw_stats_scan_sse2_128_32.c

ALIGNMENT_SRC += src/nw_stats_scan_sse2_128_16.c
ALIGNMENT_SRC += src/sg_stats_scan_sse2_128_16.c
#ALIGNMENT_SRC += src/sw_stats_scan_sse2_128_16.c

ALIGNMENT_SRC += src/nw_stats_scan_sse2_128_8.c
#ALIGNMENT_SRC += src/sg_stats_scan_sse2_128_8.c
#ALIGNMENT_SRC += src/sw_stats_scan_sse2_128_8.c
endif

if HAVE_SSE41
ALIGNMENT_SRC += src/nw_stats_scan_sse41_128_32.c
#ALIGNMENT_SRC += src/sg_stats_scan_sse41_128_32.c
#ALIGNMENT_SRC += src/sw_stats_scan_sse41_128_32.c

ALIGNMENT_SRC += src/nw_stats_scan_sse41_128_16.c
#ALIGNMENT_SRC += src/sg_stats_scan_sse41_128_16.c
#ALIGNMENT_SRC += src/sw_stats_scan_sse41_128_16.c

ALIGNMENT_SRC += src/nw_stats_scan_sse41_128_8.c
#ALIGNMENT_SRC += src/sg_stats_scan_sse41_128_8.c
#ALIGNMENT_SRC += src/sw_stats_scan_sse41_128_8.c
endif

libparasail_la_SOURCES += $(ALIGNMENT_SRC)

libparasail_table_la_SOURCES = $(ALIGNMENT_SRC)
libparasail_table_la_CPPFLAGS = $(AM_CPPFLAGS) -DPARASAIL_TABLE

libparasail_la_LIBADD = libparasail_table.la

check_PROGRAMS += tests/test_align
check_PROGRAMS += tests/test_align_stats
check_PROGRAMS += tests/test_align_table
if HAVE_OPENMP
check_PROGRAMS += tests/test_openmp
endif

tests_test_align_SOURCES        = tests/test_align.c
tests_test_align_stats_SOURCES  = tests/test_align_stats.c
tests_test_align_table_SOURCES  = tests/test_align_table.c
tests_test_openmp_SOURCES       = tests/test_openmp.c
tests_test_openmp_CFLAGS        = $(AM_CFLAGS) $(OPENMP_CFLAGS)
tests_test_openmp_LDFLAGS       = $(AM_LDFLAGS) $(OPENMP_CFLAGS)
tests_test_openmp_LDADD         = $(LDADD) $(CLOCK_LIBS) $(MATH_LIBS)

