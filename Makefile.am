ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES =
check_PROGRAMS =
AM_CFLAGS =
AM_CPPFLAGS =
AM_LDFLAGS =
LDADD = libparasail.la

AM_CPPFLAGS += -I$(top_srcdir)/src

lib_LTLIBRARIES += libparasail.la
libparasail_la_SOURCES =
libparasail_la_SOURCES += src/align.h
libparasail_la_SOURCES += src/nw.c
libparasail_la_SOURCES += src/nw_scan.c
libparasail_la_SOURCES += src/nw_scan_row.c
libparasail_la_SOURCES += src/nw_stats.c
libparasail_la_SOURCES += src/nw_stats_scan.c
libparasail_la_SOURCES += src/sg.c
libparasail_la_SOURCES += src/sg_scan.c
libparasail_la_SOURCES += src/sg_stats.c
libparasail_la_SOURCES += src/sg_stats_scan.c
libparasail_la_SOURCES += src/sw.c
libparasail_la_SOURCES += src/sw_scan.c
libparasail_la_SOURCES += src/sw_stats.c
libparasail_la_SOURCES += src/sw_stats_scan.c
if HAVE_SSE2
libparasail_la_SOURCES += src/align_scan_128_16.h
libparasail_la_SOURCES += src/align_striped_128_16.h
libparasail_la_SOURCES += src/align_wozniak_128_16.h
libparasail_la_SOURCES += src/nw_scan_128_16.c
libparasail_la_SOURCES += src/nw_striped_128_16.c
libparasail_la_SOURCES += src/nw_wozniak_128_16.c
libparasail_la_SOURCES += src/nw_stats_scan_128_16.c
libparasail_la_SOURCES += src/nw_stats_striped_128_16.c
libparasail_la_SOURCES += src/nw_stats_wozniak_128_16.c
libparasail_la_SOURCES += src/sg_scan_128_16.c
libparasail_la_SOURCES += src/sg_striped_128_16.c
libparasail_la_SOURCES += src/sg_wozniak_128_16.c
libparasail_la_SOURCES += src/sg_stats_scan_128_16.c
libparasail_la_SOURCES += src/sg_stats_striped_128_16.c
libparasail_la_SOURCES += src/sg_stats_wozniak_128_16.c
libparasail_la_SOURCES += src/sw_scan_128_16.c
libparasail_la_SOURCES += src/sw_striped_128_16.c
libparasail_la_SOURCES += src/sw_wozniak_128_16.c
libparasail_la_SOURCES += src/sw_stats_scan_128_16.c
libparasail_la_SOURCES += src/sw_stats_striped_128_16.c
libparasail_la_SOURCES += src/sw_stats_wozniak_128_16.c
endif
if HAVE_SSE41
libparasail_la_SOURCES += src/align_scan_128_8.h
libparasail_la_SOURCES += src/align_striped_128_8.h
libparasail_la_SOURCES += src/align_wozniak_128_8.h
libparasail_la_SOURCES += src/nw_scan_128_8.c
libparasail_la_SOURCES += src/nw_striped_128_8.c
libparasail_la_SOURCES += src/nw_wozniak_128_8.c
libparasail_la_SOURCES += src/nw_stats_scan_128_8.c
libparasail_la_SOURCES += src/nw_stats_striped_128_8.c
libparasail_la_SOURCES += src/nw_stats_wozniak_128_8.c
libparasail_la_SOURCES += src/sg_scan_128_8.c
libparasail_la_SOURCES += src/sg_striped_128_8.c
libparasail_la_SOURCES += src/sg_wozniak_128_8.c
libparasail_la_SOURCES += src/sg_stats_scan_128_8.c
libparasail_la_SOURCES += src/sg_stats_wozniak_128_8.c
libparasail_la_SOURCES += src/sw_scan_128_8.c
libparasail_la_SOURCES += src/sw_stats_scan_128_8.c
endif
if HAVE_KNC
libparasail_la_SOURCES += src/nw_scan_512_32.c
endif

check_PROGRAMS += tests/test_align
check_PROGRAMS += tests/test_align_stats
#check_PROGRAMS += tests/test_align_debug

tests_test_align_SOURCES       = tests/test_align.c
tests_test_align_stats_SOURCES = tests/test_align_stats.c
#tests_test_align_debug_SOURCES = tests/test_align_debug.c

