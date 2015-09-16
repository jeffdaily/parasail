ACLOCAL_AMFLAGS = -I m4

include_HEADERS = parasail.h
nobase_include_HEADERS =
lib_LTLIBRARIES = libparasail.la
noinst_LTLIBRARIES =
bin_PROGRAMS =
check_PROGRAMS =
check_LTLIBRARIES = 
AM_CFLAGS =
AM_CPPFLAGS =
AM_LDFLAGS =
LDADD = libparasail.la
EXTRA_DIST = version.sh

AM_CFLAGS += $(VISIBILITY_FLAGS)

if ON_MINGW
libparasail_la_LDFLAGS = $(AM_LDFLAGS) -no-undefined -avoid-version -version-info @LTVER@
else
libparasail_la_LDFLAGS = $(AM_LDFLAGS) -version-info @LTVER@
endif

nobase_include_HEADERS += parasail/cpuid.h
nobase_include_HEADERS += parasail/io.h
nobase_include_HEADERS += parasail/function_lookup.h
nobase_include_HEADERS += parasail/matrix_lookup.h

nobase_include_HEADERS += parasail/matrices/blosum100.h
nobase_include_HEADERS += parasail/matrices/blosum30.h
nobase_include_HEADERS += parasail/matrices/blosum35.h
nobase_include_HEADERS += parasail/matrices/blosum40.h
nobase_include_HEADERS += parasail/matrices/blosum45.h
nobase_include_HEADERS += parasail/matrices/blosum50.h
nobase_include_HEADERS += parasail/matrices/blosum55.h
nobase_include_HEADERS += parasail/matrices/blosum60.h
nobase_include_HEADERS += parasail/matrices/blosum62.h
nobase_include_HEADERS += parasail/matrices/blosum65.h
nobase_include_HEADERS += parasail/matrices/blosum70.h
nobase_include_HEADERS += parasail/matrices/blosum75.h
nobase_include_HEADERS += parasail/matrices/blosum80.h
nobase_include_HEADERS += parasail/matrices/blosum85.h
nobase_include_HEADERS += parasail/matrices/blosum90.h
nobase_include_HEADERS += parasail/matrices/blosum_map.h
nobase_include_HEADERS += parasail/matrices/pam10.h
nobase_include_HEADERS += parasail/matrices/pam100.h
nobase_include_HEADERS += parasail/matrices/pam110.h
nobase_include_HEADERS += parasail/matrices/pam120.h
nobase_include_HEADERS += parasail/matrices/pam130.h
nobase_include_HEADERS += parasail/matrices/pam140.h
nobase_include_HEADERS += parasail/matrices/pam150.h
nobase_include_HEADERS += parasail/matrices/pam160.h
nobase_include_HEADERS += parasail/matrices/pam170.h
nobase_include_HEADERS += parasail/matrices/pam180.h
nobase_include_HEADERS += parasail/matrices/pam190.h
nobase_include_HEADERS += parasail/matrices/pam20.h
nobase_include_HEADERS += parasail/matrices/pam200.h
nobase_include_HEADERS += parasail/matrices/pam210.h
nobase_include_HEADERS += parasail/matrices/pam220.h
nobase_include_HEADERS += parasail/matrices/pam230.h
nobase_include_HEADERS += parasail/matrices/pam240.h
nobase_include_HEADERS += parasail/matrices/pam250.h
nobase_include_HEADERS += parasail/matrices/pam260.h
nobase_include_HEADERS += parasail/matrices/pam270.h
nobase_include_HEADERS += parasail/matrices/pam280.h
nobase_include_HEADERS += parasail/matrices/pam290.h
nobase_include_HEADERS += parasail/matrices/pam30.h
nobase_include_HEADERS += parasail/matrices/pam300.h
nobase_include_HEADERS += parasail/matrices/pam310.h
nobase_include_HEADERS += parasail/matrices/pam320.h
nobase_include_HEADERS += parasail/matrices/pam330.h
nobase_include_HEADERS += parasail/matrices/pam340.h
nobase_include_HEADERS += parasail/matrices/pam350.h
nobase_include_HEADERS += parasail/matrices/pam360.h
nobase_include_HEADERS += parasail/matrices/pam370.h
nobase_include_HEADERS += parasail/matrices/pam380.h
nobase_include_HEADERS += parasail/matrices/pam390.h
nobase_include_HEADERS += parasail/matrices/pam40.h
nobase_include_HEADERS += parasail/matrices/pam400.h
nobase_include_HEADERS += parasail/matrices/pam410.h
nobase_include_HEADERS += parasail/matrices/pam420.h
nobase_include_HEADERS += parasail/matrices/pam430.h
nobase_include_HEADERS += parasail/matrices/pam440.h
nobase_include_HEADERS += parasail/matrices/pam450.h
nobase_include_HEADERS += parasail/matrices/pam460.h
nobase_include_HEADERS += parasail/matrices/pam470.h
nobase_include_HEADERS += parasail/matrices/pam480.h
nobase_include_HEADERS += parasail/matrices/pam490.h
nobase_include_HEADERS += parasail/matrices/pam50.h
nobase_include_HEADERS += parasail/matrices/pam500.h
nobase_include_HEADERS += parasail/matrices/pam60.h
nobase_include_HEADERS += parasail/matrices/pam70.h
nobase_include_HEADERS += parasail/matrices/pam80.h
nobase_include_HEADERS += parasail/matrices/pam90.h
nobase_include_HEADERS += parasail/matrices/pam_map.h

AM_CPPFLAGS += -I$(top_srcdir)/contrib
AM_CPPFLAGS += -I$(top_srcdir)/contrib/sais-lite-lcp

noinst_LTLIBRARIES += libparasail_novec.la
noinst_LTLIBRARIES += libparasail_novec_table.la
noinst_LTLIBRARIES += libparasail_novec_rowcol.la
if HAVE_SSE
noinst_LTLIBRARIES += libparasail_sse_memory.la
endif
if HAVE_SSE2
noinst_LTLIBRARIES += libparasail_sse2.la
noinst_LTLIBRARIES += libparasail_sse2_table.la
noinst_LTLIBRARIES += libparasail_sse2_rowcol.la
endif
if HAVE_SSE41
noinst_LTLIBRARIES += libparasail_sse41.la
noinst_LTLIBRARIES += libparasail_sse41_table.la
noinst_LTLIBRARIES += libparasail_sse41_rowcol.la
endif
if HAVE_AVX2
noinst_LTLIBRARIES += libparasail_avx2_memory.la
noinst_LTLIBRARIES += libparasail_avx2.la
noinst_LTLIBRARIES += libparasail_avx2_table.la
noinst_LTLIBRARIES += libparasail_avx2_rowcol.la
endif
if HAVE_KNC
noinst_LTLIBRARIES += libparasail_knc_memory.la
noinst_LTLIBRARIES += libparasail_knc.la
noinst_LTLIBRARIES += libparasail_knc_table.la
noinst_LTLIBRARIES += libparasail_knc_rowcol.la
endif

SRC_CORE =
MEM_SRC_SSE =
MEM_SRC_AVX2 =
MEM_SRC_KNC =
SRC_NOVEC =
SRC_SSE2 =
SRC_SSE41 =
SRC_AVX2 =
SRC_KNC =

#################
# core sources
#################

SRC_CORE += src/cpuid.c
SRC_CORE += src/function_lookup.c
SRC_CORE += src/io.c
SRC_CORE += src/matrix_lookup.c
SRC_CORE += src/memory.c
SRC_CORE += src/time.c
SRC_CORE += src/nw_dispatch.c
SRC_CORE += src/sg_dispatch.c
SRC_CORE += src/sw_dispatch.c
SRC_CORE += src/dispatch_profile.c
SRC_CORE += src/satcheck.c
SRC_CORE += parasail/memory.h
SRC_CORE += parasail/stats.h

#################
# vector memory
#################

MEM_SRC_SSE += src/memory_sse.c
MEM_SRC_SSE += parasail/internal_sse.h

MEM_SRC_AVX2 += src/memory_avx2.c
MEM_SRC_AVX2 += parasail/internal_avx.h

MEM_SRC_KNC += src/memory_knc.c
MEM_SRC_KNC += parasail/internal_knc.h

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

SRC_SSE2 += src/nw_scan_sse2_128_64.c
SRC_SSE2 += src/sg_scan_sse2_128_64.c
SRC_SSE2 += src/sw_scan_sse2_128_64.c

SRC_SSE2 += src/nw_scan_sse2_128_32.c
SRC_SSE2 += src/sg_scan_sse2_128_32.c
SRC_SSE2 += src/sw_scan_sse2_128_32.c

SRC_SSE2 += src/nw_scan_sse2_128_16.c
SRC_SSE2 += src/sg_scan_sse2_128_16.c
SRC_SSE2 += src/sw_scan_sse2_128_16.c

SRC_SSE2 += src/nw_scan_sse2_128_8.c
SRC_SSE2 += src/sg_scan_sse2_128_8.c
SRC_SSE2 += src/sw_scan_sse2_128_8.c

SRC_SSE41 += src/nw_scan_sse41_128_64.c
SRC_SSE41 += src/sg_scan_sse41_128_64.c
SRC_SSE41 += src/sw_scan_sse41_128_64.c

SRC_SSE41 += src/nw_scan_sse41_128_32.c
SRC_SSE41 += src/sg_scan_sse41_128_32.c
SRC_SSE41 += src/sw_scan_sse41_128_32.c

SRC_SSE41 += src/nw_scan_sse41_128_16.c
SRC_SSE41 += src/sg_scan_sse41_128_16.c
SRC_SSE41 += src/sw_scan_sse41_128_16.c

SRC_SSE41 += src/nw_scan_sse41_128_8.c
SRC_SSE41 += src/sg_scan_sse41_128_8.c
SRC_SSE41 += src/sw_scan_sse41_128_8.c

SRC_AVX2 += src/nw_scan_avx2_256_64.c
SRC_AVX2 += src/sg_scan_avx2_256_64.c
SRC_AVX2 += src/sw_scan_avx2_256_64.c

SRC_AVX2 += src/nw_scan_avx2_256_32.c
SRC_AVX2 += src/sg_scan_avx2_256_32.c
SRC_AVX2 += src/sw_scan_avx2_256_32.c

SRC_AVX2 += src/nw_scan_avx2_256_16.c
SRC_AVX2 += src/sg_scan_avx2_256_16.c
SRC_AVX2 += src/sw_scan_avx2_256_16.c

SRC_AVX2 += src/nw_scan_avx2_256_8.c
SRC_AVX2 += src/sg_scan_avx2_256_8.c
SRC_AVX2 += src/sw_scan_avx2_256_8.c

SRC_KNC += src/nw_scan_knc_512_32.c
SRC_KNC += src/sg_scan_knc_512_32.c
SRC_KNC += src/sw_scan_knc_512_32.c

# parallel diag methods

SRC_SSE2 += src/nw_diag_sse2_128_64.c
SRC_SSE2 += src/sg_diag_sse2_128_64.c
SRC_SSE2 += src/sw_diag_sse2_128_64.c

SRC_SSE2 += src/nw_diag_sse2_128_32.c
SRC_SSE2 += src/sg_diag_sse2_128_32.c
SRC_SSE2 += src/sw_diag_sse2_128_32.c

SRC_SSE2 += src/nw_diag_sse2_128_16.c
SRC_SSE2 += src/sg_diag_sse2_128_16.c
SRC_SSE2 += src/sw_diag_sse2_128_16.c

SRC_SSE2 += src/nw_diag_sse2_128_8.c
SRC_SSE2 += src/sg_diag_sse2_128_8.c
SRC_SSE2 += src/sw_diag_sse2_128_8.c

SRC_SSE41 += src/nw_diag_sse41_128_64.c
SRC_SSE41 += src/sg_diag_sse41_128_64.c
SRC_SSE41 += src/sw_diag_sse41_128_64.c

SRC_SSE41 += src/nw_diag_sse41_128_32.c
SRC_SSE41 += src/sg_diag_sse41_128_32.c
SRC_SSE41 += src/sw_diag_sse41_128_32.c

SRC_SSE41 += src/nw_diag_sse41_128_16.c
SRC_SSE41 += src/sg_diag_sse41_128_16.c
SRC_SSE41 += src/sw_diag_sse41_128_16.c

SRC_SSE41 += src/nw_diag_sse41_128_8.c
SRC_SSE41 += src/sg_diag_sse41_128_8.c
SRC_SSE41 += src/sw_diag_sse41_128_8.c

SRC_AVX2 += src/nw_diag_avx2_256_64.c
SRC_AVX2 += src/sg_diag_avx2_256_64.c
SRC_AVX2 += src/sw_diag_avx2_256_64.c

SRC_AVX2 += src/nw_diag_avx2_256_32.c
SRC_AVX2 += src/sg_diag_avx2_256_32.c
SRC_AVX2 += src/sw_diag_avx2_256_32.c

SRC_AVX2 += src/nw_diag_avx2_256_16.c
SRC_AVX2 += src/sg_diag_avx2_256_16.c
SRC_AVX2 += src/sw_diag_avx2_256_16.c

SRC_AVX2 += src/nw_diag_avx2_256_8.c
SRC_AVX2 += src/sg_diag_avx2_256_8.c
SRC_AVX2 += src/sw_diag_avx2_256_8.c

SRC_KNC += src/nw_diag_knc_512_32.c
SRC_KNC += src/sg_diag_knc_512_32.c
SRC_KNC += src/sw_diag_knc_512_32.c

# parallel striped methods

SRC_SSE2 += src/nw_striped_sse2_128_64.c
SRC_SSE2 += src/sg_striped_sse2_128_64.c
SRC_SSE2 += src/sw_striped_sse2_128_64.c

SRC_SSE2 += src/nw_striped_sse2_128_32.c
SRC_SSE2 += src/sg_striped_sse2_128_32.c
SRC_SSE2 += src/sw_striped_sse2_128_32.c

SRC_SSE2 += src/nw_striped_sse2_128_16.c
SRC_SSE2 += src/sg_striped_sse2_128_16.c
SRC_SSE2 += src/sw_striped_sse2_128_16.c

SRC_SSE2 += src/nw_striped_sse2_128_8.c
SRC_SSE2 += src/sg_striped_sse2_128_8.c
SRC_SSE2 += src/sw_striped_sse2_128_8.c

SRC_SSE41 += src/nw_striped_sse41_128_64.c
SRC_SSE41 += src/sg_striped_sse41_128_64.c
SRC_SSE41 += src/sw_striped_sse41_128_64.c

SRC_SSE41 += src/nw_striped_sse41_128_32.c
SRC_SSE41 += src/sg_striped_sse41_128_32.c
SRC_SSE41 += src/sw_striped_sse41_128_32.c

SRC_SSE41 += src/nw_striped_sse41_128_16.c
SRC_SSE41 += src/sg_striped_sse41_128_16.c
SRC_SSE41 += src/sw_striped_sse41_128_16.c

SRC_SSE41 += src/nw_striped_sse41_128_8.c
SRC_SSE41 += src/sg_striped_sse41_128_8.c
SRC_SSE41 += src/sw_striped_sse41_128_8.c

SRC_AVX2 += src/nw_striped_avx2_256_64.c
SRC_AVX2 += src/sg_striped_avx2_256_64.c
SRC_AVX2 += src/sw_striped_avx2_256_64.c

SRC_AVX2 += src/nw_striped_avx2_256_32.c
SRC_AVX2 += src/sg_striped_avx2_256_32.c
SRC_AVX2 += src/sw_striped_avx2_256_32.c

SRC_AVX2 += src/nw_striped_avx2_256_16.c
SRC_AVX2 += src/sg_striped_avx2_256_16.c
SRC_AVX2 += src/sw_striped_avx2_256_16.c

SRC_AVX2 += src/nw_striped_avx2_256_8.c
SRC_AVX2 += src/sg_striped_avx2_256_8.c
SRC_AVX2 += src/sw_striped_avx2_256_8.c

SRC_KNC += src/nw_striped_knc_512_32.c
SRC_KNC += src/sg_striped_knc_512_32.c
SRC_KNC += src/sw_striped_knc_512_32.c

# parallel blocked methods

#SRC_SSE2 += src/nw_blocked_sse2_128_32.c
#SRC_SSE2 += src/sg_blocked_sse2_128_32.c
#SRC_SSE2 += src/sw_blocked_sse2_128_32.c
#
#SRC_SSE2 += src/nw_blocked_sse2_128_16.c
#SRC_SSE2 += src/sg_blocked_sse2_128_16.c
#SRC_SSE2 += src/sw_blocked_sse2_128_16.c
#
#SRC_SSE2 += src/nw_blocked_sse2_128_8.c
#SRC_SSE2 += src/sg_blocked_sse2_128_8.c
#SRC_SSE2 += src/sw_blocked_sse2_128_8.c
#
#SRC_SSE41 += src/nw_blocked_sse41_128_32.c
#SRC_SSE41 += src/sg_blocked_sse41_128_32.c
SRC_SSE41 += src/sw_blocked_sse41_128_32.c
#
#SRC_SSE41 += src/nw_blocked_sse41_128_16.c
#SRC_SSE41 += src/sg_blocked_sse41_128_16.c
SRC_SSE41 += src/sw_blocked_sse41_128_16.c
#
#SRC_SSE41 += src/nw_blocked_sse41_128_8.c
#SRC_SSE41 += src/sg_blocked_sse41_128_8.c
#SRC_SSE41 += src/sw_blocked_sse41_128_8.c
#
#SRC_AVX2 += src/nw_blocked_avx2_256_32.c
#SRC_AVX2 += src/sg_blocked_avx2_256_32.c
#SRC_AVX2 += src/sw_blocked_avx2_256_32.c
#
#SRC_AVX2 += src/nw_blocked_avx2_256_16.c
#SRC_AVX2 += src/sg_blocked_avx2_256_16.c
#SRC_AVX2 += src/sw_blocked_avx2_256_16.c
#
#SRC_AVX2 += src/nw_blocked_avx2_256_8.c
#SRC_AVX2 += src/sg_blocked_avx2_256_8.c
#SRC_AVX2 += src/sw_blocked_avx2_256_8.c
#
#SRC_KNC += src/nw_blocked_knc_512_32.c
#SRC_KNC += src/sg_blocked_knc_512_32.c
#SRC_KNC += src/sw_blocked_knc_512_32.c

#########################
# parallel stats methods
#########################

# parallel scan methods

SRC_SSE2 += src/nw_stats_scan_sse2_128_64.c
SRC_SSE2 += src/sg_stats_scan_sse2_128_64.c
SRC_SSE2 += src/sw_stats_scan_sse2_128_64.c

SRC_SSE2 += src/nw_stats_scan_sse2_128_32.c
SRC_SSE2 += src/sg_stats_scan_sse2_128_32.c
SRC_SSE2 += src/sw_stats_scan_sse2_128_32.c

SRC_SSE2 += src/nw_stats_scan_sse2_128_16.c
SRC_SSE2 += src/sg_stats_scan_sse2_128_16.c
SRC_SSE2 += src/sw_stats_scan_sse2_128_16.c

SRC_SSE2 += src/nw_stats_scan_sse2_128_8.c
SRC_SSE2 += src/sg_stats_scan_sse2_128_8.c
SRC_SSE2 += src/sw_stats_scan_sse2_128_8.c

SRC_SSE41 += src/nw_stats_scan_sse41_128_64.c
SRC_SSE41 += src/sg_stats_scan_sse41_128_64.c
SRC_SSE41 += src/sw_stats_scan_sse41_128_64.c

SRC_SSE41 += src/nw_stats_scan_sse41_128_32.c
SRC_SSE41 += src/sg_stats_scan_sse41_128_32.c
SRC_SSE41 += src/sw_stats_scan_sse41_128_32.c

SRC_SSE41 += src/nw_stats_scan_sse41_128_16.c
SRC_SSE41 += src/sg_stats_scan_sse41_128_16.c
SRC_SSE41 += src/sw_stats_scan_sse41_128_16.c

SRC_SSE41 += src/nw_stats_scan_sse41_128_8.c
SRC_SSE41 += src/sg_stats_scan_sse41_128_8.c
SRC_SSE41 += src/sw_stats_scan_sse41_128_8.c

SRC_AVX2 += src/nw_stats_scan_avx2_256_64.c
SRC_AVX2 += src/sg_stats_scan_avx2_256_64.c
SRC_AVX2 += src/sw_stats_scan_avx2_256_64.c

SRC_AVX2 += src/nw_stats_scan_avx2_256_32.c
SRC_AVX2 += src/sg_stats_scan_avx2_256_32.c
SRC_AVX2 += src/sw_stats_scan_avx2_256_32.c

SRC_AVX2 += src/nw_stats_scan_avx2_256_16.c
SRC_AVX2 += src/sg_stats_scan_avx2_256_16.c
SRC_AVX2 += src/sw_stats_scan_avx2_256_16.c

SRC_AVX2 += src/nw_stats_scan_avx2_256_8.c
SRC_AVX2 += src/sg_stats_scan_avx2_256_8.c
SRC_AVX2 += src/sw_stats_scan_avx2_256_8.c

SRC_KNC += src/nw_stats_scan_knc_512_32.c
SRC_KNC += src/sg_stats_scan_knc_512_32.c
SRC_KNC += src/sw_stats_scan_knc_512_32.c

# parallel diag methods

SRC_SSE2 += src/nw_stats_diag_sse2_128_64.c
SRC_SSE2 += src/sg_stats_diag_sse2_128_64.c
SRC_SSE2 += src/sw_stats_diag_sse2_128_64.c

SRC_SSE2 += src/nw_stats_diag_sse2_128_32.c
SRC_SSE2 += src/sg_stats_diag_sse2_128_32.c
SRC_SSE2 += src/sw_stats_diag_sse2_128_32.c

SRC_SSE2 += src/nw_stats_diag_sse2_128_16.c
SRC_SSE2 += src/sg_stats_diag_sse2_128_16.c
SRC_SSE2 += src/sw_stats_diag_sse2_128_16.c

SRC_SSE2 += src/nw_stats_diag_sse2_128_8.c
SRC_SSE2 += src/sg_stats_diag_sse2_128_8.c
SRC_SSE2 += src/sw_stats_diag_sse2_128_8.c

SRC_SSE41 += src/nw_stats_diag_sse41_128_64.c
SRC_SSE41 += src/sg_stats_diag_sse41_128_64.c
SRC_SSE41 += src/sw_stats_diag_sse41_128_64.c

SRC_SSE41 += src/nw_stats_diag_sse41_128_32.c
SRC_SSE41 += src/sg_stats_diag_sse41_128_32.c
SRC_SSE41 += src/sw_stats_diag_sse41_128_32.c

SRC_SSE41 += src/nw_stats_diag_sse41_128_16.c
SRC_SSE41 += src/sg_stats_diag_sse41_128_16.c
SRC_SSE41 += src/sw_stats_diag_sse41_128_16.c

SRC_SSE41 += src/nw_stats_diag_sse41_128_8.c
SRC_SSE41 += src/sg_stats_diag_sse41_128_8.c
SRC_SSE41 += src/sw_stats_diag_sse41_128_8.c

SRC_AVX2 += src/nw_stats_diag_avx2_256_64.c
SRC_AVX2 += src/sg_stats_diag_avx2_256_64.c
SRC_AVX2 += src/sw_stats_diag_avx2_256_64.c

SRC_AVX2 += src/nw_stats_diag_avx2_256_32.c
SRC_AVX2 += src/sg_stats_diag_avx2_256_32.c
SRC_AVX2 += src/sw_stats_diag_avx2_256_32.c

SRC_AVX2 += src/nw_stats_diag_avx2_256_16.c
SRC_AVX2 += src/sg_stats_diag_avx2_256_16.c
SRC_AVX2 += src/sw_stats_diag_avx2_256_16.c

SRC_AVX2 += src/nw_stats_diag_avx2_256_8.c
SRC_AVX2 += src/sg_stats_diag_avx2_256_8.c
SRC_AVX2 += src/sw_stats_diag_avx2_256_8.c

SRC_KNC += src/nw_stats_diag_knc_512_32.c
SRC_KNC += src/sg_stats_diag_knc_512_32.c
SRC_KNC += src/sw_stats_diag_knc_512_32.c

# parallel striped methods

SRC_SSE2 += src/nw_stats_striped_sse2_128_64.c
SRC_SSE2 += src/sg_stats_striped_sse2_128_64.c
SRC_SSE2 += src/sw_stats_striped_sse2_128_64.c

SRC_SSE2 += src/nw_stats_striped_sse2_128_32.c
SRC_SSE2 += src/sg_stats_striped_sse2_128_32.c
SRC_SSE2 += src/sw_stats_striped_sse2_128_32.c

SRC_SSE2 += src/nw_stats_striped_sse2_128_16.c
SRC_SSE2 += src/sg_stats_striped_sse2_128_16.c
SRC_SSE2 += src/sw_stats_striped_sse2_128_16.c

SRC_SSE2 += src/nw_stats_striped_sse2_128_8.c
SRC_SSE2 += src/sg_stats_striped_sse2_128_8.c
SRC_SSE2 += src/sw_stats_striped_sse2_128_8.c

SRC_SSE41 += src/nw_stats_striped_sse41_128_64.c
SRC_SSE41 += src/sg_stats_striped_sse41_128_64.c
SRC_SSE41 += src/sw_stats_striped_sse41_128_64.c

SRC_SSE41 += src/nw_stats_striped_sse41_128_32.c
SRC_SSE41 += src/sg_stats_striped_sse41_128_32.c
SRC_SSE41 += src/sw_stats_striped_sse41_128_32.c

SRC_SSE41 += src/nw_stats_striped_sse41_128_16.c
SRC_SSE41 += src/sg_stats_striped_sse41_128_16.c
SRC_SSE41 += src/sw_stats_striped_sse41_128_16.c

SRC_SSE41 += src/nw_stats_striped_sse41_128_8.c
SRC_SSE41 += src/sg_stats_striped_sse41_128_8.c
SRC_SSE41 += src/sw_stats_striped_sse41_128_8.c

SRC_AVX2 += src/nw_stats_striped_avx2_256_64.c
SRC_AVX2 += src/sg_stats_striped_avx2_256_64.c
SRC_AVX2 += src/sw_stats_striped_avx2_256_64.c

SRC_AVX2 += src/nw_stats_striped_avx2_256_32.c
SRC_AVX2 += src/sg_stats_striped_avx2_256_32.c
SRC_AVX2 += src/sw_stats_striped_avx2_256_32.c

SRC_AVX2 += src/nw_stats_striped_avx2_256_16.c
SRC_AVX2 += src/sg_stats_striped_avx2_256_16.c
SRC_AVX2 += src/sw_stats_striped_avx2_256_16.c

SRC_AVX2 += src/nw_stats_striped_avx2_256_8.c
SRC_AVX2 += src/sg_stats_striped_avx2_256_8.c
SRC_AVX2 += src/sw_stats_striped_avx2_256_8.c

SRC_KNC += src/nw_stats_striped_knc_512_32.c
SRC_KNC += src/sg_stats_striped_knc_512_32.c
SRC_KNC += src/sw_stats_striped_knc_512_32.c

#########################
# library magic
#########################

libparasail_la_SOURCES = $(SRC_CORE)

libparasail_novec_la_SOURCES = $(SRC_NOVEC)
libparasail_sse2_la_SOURCES  = $(SRC_SSE2)
libparasail_sse41_la_SOURCES = $(SRC_SSE41)
libparasail_avx2_la_SOURCES  = $(SRC_AVX2)
libparasail_knc_la_SOURCES   = $(SRC_KNC)

libparasail_novec_la_CFLAGS = $(AM_CFLAGS)
libparasail_sse2_la_CFLAGS  = $(AM_CFLAGS) $(SSE2_CFLAGS)
libparasail_sse41_la_CFLAGS = $(AM_CFLAGS) $(SSE41_CFLAGS)
libparasail_avx2_la_CFLAGS  = $(AM_CFLAGS) $(AVX2_CFLAGS)
libparasail_knc_la_CFLAGS   = $(AM_CFLAGS)

libparasail_novec_table_la_SOURCES = $(SRC_NOVEC)
libparasail_sse2_table_la_SOURCES  = $(SRC_SSE2)
libparasail_sse41_table_la_SOURCES = $(SRC_SSE41)
libparasail_avx2_table_la_SOURCES  = $(SRC_AVX2)
libparasail_knc_table_la_SOURCES   = $(SRC_KNC)

libparasail_novec_table_la_CFLAGS = $(AM_CFLAGS)
libparasail_sse2_table_la_CFLAGS  = $(AM_CFLAGS) $(SSE2_CFLAGS)
libparasail_sse41_table_la_CFLAGS = $(AM_CFLAGS) $(SSE41_CFLAGS)
libparasail_avx2_table_la_CFLAGS  = $(AM_CFLAGS) $(AVX2_CFLAGS)
libparasail_knc_table_la_CFLAGS   = $(AM_CFLAGS)

libparasail_novec_table_la_CPPFLAGS = $(AM_CPPFLAGS) -DPARASAIL_TABLE=1
libparasail_sse2_table_la_CPPFLAGS  = $(AM_CPPFLAGS) -DPARASAIL_TABLE=1
libparasail_sse41_table_la_CPPFLAGS = $(AM_CPPFLAGS) -DPARASAIL_TABLE=1
libparasail_avx2_table_la_CPPFLAGS  = $(AM_CPPFLAGS) -DPARASAIL_TABLE=1
libparasail_knc_table_la_CPPFLAGS   = $(AM_CPPFLAGS) -DPARASAIL_TABLE=1

libparasail_novec_rowcol_la_SOURCES = $(SRC_NOVEC)
libparasail_sse2_rowcol_la_SOURCES  = $(SRC_SSE2)
libparasail_sse41_rowcol_la_SOURCES = $(SRC_SSE41)
libparasail_avx2_rowcol_la_SOURCES  = $(SRC_AVX2)
libparasail_knc_rowcol_la_SOURCES   = $(SRC_KNC)

libparasail_novec_rowcol_la_CFLAGS = $(AM_CFLAGS)
libparasail_sse2_rowcol_la_CFLAGS  = $(AM_CFLAGS) $(SSE2_CFLAGS)
libparasail_sse41_rowcol_la_CFLAGS = $(AM_CFLAGS) $(SSE41_CFLAGS)
libparasail_avx2_rowcol_la_CFLAGS  = $(AM_CFLAGS) $(AVX2_CFLAGS)
libparasail_knc_rowcol_la_CFLAGS   = $(AM_CFLAGS)

libparasail_novec_rowcol_la_CPPFLAGS = $(AM_CPPFLAGS) -DPARASAIL_ROWCOL=1
libparasail_sse2_rowcol_la_CPPFLAGS  = $(AM_CPPFLAGS) -DPARASAIL_ROWCOL=1
libparasail_sse41_rowcol_la_CPPFLAGS = $(AM_CPPFLAGS) -DPARASAIL_ROWCOL=1
libparasail_avx2_rowcol_la_CPPFLAGS  = $(AM_CPPFLAGS) -DPARASAIL_ROWCOL=1
libparasail_knc_rowcol_la_CPPFLAGS   = $(AM_CPPFLAGS) -DPARASAIL_ROWCOL=1

libparasail_sse_memory_la_SOURCES  = $(MEM_SRC_SSE)
libparasail_avx2_memory_la_SOURCES = $(MEM_SRC_AVX2)
libparasail_knc_memory_la_SOURCES  = $(MEM_SRC_KNC)

if HAVE_SSE2
libparasail_sse_memory_la_CFLAGS  = $(AM_CFLAGS) $(SSE2_CFLAGS)
else
if HAVE_SSE41
libparasail_sse_memory_la_CFLAGS  = $(AM_CFLAGS) $(SSE41_CFLAGS)
endif
endif
libparasail_avx2_memory_la_CFLAGS = $(AM_CFLAGS) $(AVX2_CFLAGS)
libparasail_knc_memory_la_CFLAGS  = $(AM_CFLAGS)

libparasail_la_LIBADD =
libparasail_la_LIBADD += libparasail_novec.la
libparasail_la_LIBADD += libparasail_novec_table.la
libparasail_la_LIBADD += libparasail_novec_rowcol.la
if HAVE_SSE
libparasail_la_LIBADD += libparasail_sse_memory.la
endif
if HAVE_SSE2
libparasail_la_LIBADD += libparasail_sse2.la
libparasail_la_LIBADD += libparasail_sse2_table.la
libparasail_la_LIBADD += libparasail_sse2_rowcol.la
endif
if HAVE_SSE41
libparasail_la_LIBADD += libparasail_sse41.la
libparasail_la_LIBADD += libparasail_sse41_table.la
libparasail_la_LIBADD += libparasail_sse41_rowcol.la
endif
if HAVE_AVX2
libparasail_la_LIBADD += libparasail_avx2_memory.la
libparasail_la_LIBADD += libparasail_avx2.la
libparasail_la_LIBADD += libparasail_avx2_table.la
libparasail_la_LIBADD += libparasail_avx2_rowcol.la
endif
if HAVE_KNC
libparasail_la_LIBADD += libparasail_knc_memory.la
libparasail_la_LIBADD += libparasail_knc.la
libparasail_la_LIBADD += libparasail_knc_table.la
libparasail_la_LIBADD += libparasail_knc_rowcol.la
endif
libparasail_la_LIBADD += $(CLOCK_LIBS)

if HAVE_SSE2
check_LTLIBRARIES += tests/libssw.la
tests_libssw_la_SOURCES =
tests_libssw_la_SOURCES += contrib/ssw.c
tests_libssw_la_SOURCES += contrib/ssw.h
tests_libssw_la_CFLAGS  = $(AM_CFLAGS) $(SSE2_CFLAGS)
endif

check_LTLIBRARIES += tests/libsais.la
tests_libsais_la_SOURCES =
tests_libsais_la_SOURCES += contrib/sais-lite-lcp/sais.c
tests_libsais_la_SOURCES += contrib/sais-lite-lcp/sais.h

bin_PROGRAMS += apps/parasail_aligner
bin_PROGRAMS += apps/parasail_stats
if HAVE_CILK
bin_PROGRAMS += apps/parasail_aligner_cilk
endif

check_PROGRAMS += tests/test_align
check_PROGRAMS += tests/test_gcups
check_PROGRAMS += tests/test_isa
check_PROGRAMS += tests/test_openmp
check_PROGRAMS += tests/test_query
check_PROGRAMS += tests/test_scatter
check_PROGRAMS += tests/test_sw_dispatch
check_PROGRAMS += tests/test_verify
check_PROGRAMS += tests/test_verify_tables
check_PROGRAMS += tests/test_verify_rowcols
if HAVE_SSE2
check_PROGRAMS += tests/test_ssw
endif

#LDADD += $(CLOCK_LIBS)
LDADD += $(MATH_LIBS)
LDADD += $(Z_LIBS)

apps_parasail_aligner_SOURCES = apps/parasail_aligner.cpp
apps_parasail_aligner_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
apps_parasail_aligner_LDFLAGS  = $(AM_LDFLAGS) $(OPENMP_CXXFLAGS)
apps_parasail_aligner_LDADD =
apps_parasail_aligner_LDADD += $(LDADD)
apps_parasail_aligner_LDADD += tests/libsais.la

apps_parasail_aligner_cilk_SOURCES = apps/parasail_aligner.cpp
apps_parasail_aligner_cilk_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_CILK=1
apps_parasail_aligner_cilk_CXXFLAGS = $(AM_CXXFLAGS) $(CILK_CXXFLAGS)
apps_parasail_aligner_cilk_LDFLAGS  = $(AM_LDFLAGS) $(CILK_CXXFLAGS)
apps_parasail_aligner_cilk_LDADD =
apps_parasail_aligner_cilk_LDADD += $(LDADD)
apps_parasail_aligner_cilk_LDADD += $(CILK_LIBS)
apps_parasail_aligner_cilk_LDADD += tests/libsais.la

apps_parasail_stats_SOURCES = apps/parasail_stats.c

tests_test_align_SOURCES = tests/test_align.c

tests_test_gcups_SOURCES = tests/test_gcups.c

tests_test_isa_SOURCES = tests/test_isa.c

tests_test_openmp_SOURCES = tests/test_openmp.c
tests_test_openmp_CFLAGS  = $(AM_CFLAGS) $(OPENMP_CFLAGS)
tests_test_openmp_LDFLAGS = $(AM_LDFLAGS) $(OPENMP_CFLAGS)
tests_test_openmp_LDADD   =
tests_test_openmp_LDADD   += $(LDADD)
if HAVE_SSE2
tests_test_openmp_LDADD   += tests/libssw.la
endif

tests_test_query_SOURCES = tests/test_query.c
tests_test_query_CFLAGS  = $(AM_CFLAGS) $(OPENMP_CFLAGS)
tests_test_query_LDFLAGS = $(AM_LDFLAGS) $(OPENMP_CFLAGS)

tests_test_scatter_SOURCES = tests/test_scatter.c
tests_test_scatter_CFLAGS  = $(AM_CFLAGS) $(OPENMP_CFLAGS)
tests_test_scatter_LDFLAGS = $(AM_LDFLAGS) $(OPENMP_CFLAGS)
tests_test_scatter_LDADD   =
tests_test_scatter_LDADD   += $(LDADD)
if HAVE_SSE2
tests_test_scatter_LDADD   += tests/libssw.la
endif

tests_test_ssw_SOURCES  = tests/test_ssw.cpp
tests_test_ssw_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
tests_test_ssw_LDFLAGS  = $(AM_LDFLAGS) $(OPENMP_CXXFLAGS)
tests_test_ssw_LDADD    =
tests_test_ssw_LDADD    += $(LDADD)
tests_test_ssw_LDADD    += tests/libsais.la
tests_test_ssw_LDADD    += tests/libssw.la

tests_test_sw_dispatch_SOURCES = tests/test_sw_dispatch.c

tests_test_verify_SOURCES = tests/test_verify.c
tests_test_verify_CFLAGS  = $(AM_CFLAGS) $(OPENMP_CFLAGS)
tests_test_verify_LDFLAGS = $(AM_LDFLAGS) $(OPENMP_CFLAGS)

tests_test_verify_tables_SOURCES = tests/test_verify_tables.c
tests_test_verify_tables_CFLAGS  = $(AM_CFLAGS) $(OPENMP_CFLAGS)
tests_test_verify_tables_LDFLAGS = $(AM_LDFLAGS) $(OPENMP_CFLAGS)

tests_test_verify_rowcols_SOURCES = tests/test_verify_rowcols.c
tests_test_verify_rowcols_CFLAGS  = $(AM_CFLAGS) $(OPENMP_CFLAGS)
tests_test_verify_rowcols_LDFLAGS = $(AM_LDFLAGS) $(OPENMP_CFLAGS)

EXTRA_DIST += DISCLAIMER
EXTRA_DIST += README.md
EXTRA_DIST += bindings/python/parasail.pyx
EXTRA_DIST += bindings/python/setup.py
EXTRA_DIST += bindings/python/test.py
EXTRA_DIST += contrib/cpuid.c
EXTRA_DIST += contrib/cpuid_gen4.c
EXTRA_DIST += contrib/kseq.h
EXTRA_DIST += contrib/sais-lite-lcp/COPYING
EXTRA_DIST += contrib/sais-lite-lcp/Makefile
EXTRA_DIST += contrib/sais-lite-lcp/README
EXTRA_DIST += contrib/sais-lite-lcp/sais.c
EXTRA_DIST += contrib/sais-lite-lcp/sais.h
EXTRA_DIST += contrib/sais-lite-lcp/sais.hxx
EXTRA_DIST += contrib/sais-lite-lcp/suftest.c
EXTRA_DIST += contrib/sais-lite-lcp/test.c
EXTRA_DIST += CMakeLists.txt
EXTRA_DIST += cmake/FindAVX2.cmake
EXTRA_DIST += cmake/FindSSE2.cmake
EXTRA_DIST += cmake/FindSSE41.cmake
EXTRA_DIST += cmake/Findxgetbv.cmake
EXTRA_DIST += cmake/config.h.in
EXTRA_DIST += util/codegen.py
EXTRA_DIST += util/diff_all.sh
EXTRA_DIST += util/dispatcher.py
EXTRA_DIST += util/func_group_tables.py
EXTRA_DIST += util/func_groups.py
EXTRA_DIST += util/funcs.py
EXTRA_DIST += util/gap_tester.py
EXTRA_DIST += util/isa.py
EXTRA_DIST += util/isa.pyc
EXTRA_DIST += util/names.py
EXTRA_DIST += util/pygen.py
EXTRA_DIST += util/tester.py
EXTRA_DIST += util/matrices/BLOSUM100
EXTRA_DIST += util/matrices/BLOSUM30
EXTRA_DIST += util/matrices/BLOSUM35
EXTRA_DIST += util/matrices/BLOSUM40
EXTRA_DIST += util/matrices/BLOSUM45
EXTRA_DIST += util/matrices/BLOSUM50
EXTRA_DIST += util/matrices/BLOSUM55
EXTRA_DIST += util/matrices/BLOSUM60
EXTRA_DIST += util/matrices/BLOSUM62
EXTRA_DIST += util/matrices/BLOSUM65
EXTRA_DIST += util/matrices/BLOSUM70
EXTRA_DIST += util/matrices/BLOSUM75
EXTRA_DIST += util/matrices/BLOSUM80
EXTRA_DIST += util/matrices/BLOSUM85
EXTRA_DIST += util/matrices/BLOSUM90
EXTRA_DIST += util/matrices/DAYHOFF
EXTRA_DIST += util/matrices/GONNET
EXTRA_DIST += util/matrices/PAM10
EXTRA_DIST += util/matrices/PAM100
EXTRA_DIST += util/matrices/PAM110
EXTRA_DIST += util/matrices/PAM120
EXTRA_DIST += util/matrices/PAM130
EXTRA_DIST += util/matrices/PAM140
EXTRA_DIST += util/matrices/PAM150
EXTRA_DIST += util/matrices/PAM160
EXTRA_DIST += util/matrices/PAM170
EXTRA_DIST += util/matrices/PAM180
EXTRA_DIST += util/matrices/PAM190
EXTRA_DIST += util/matrices/PAM20
EXTRA_DIST += util/matrices/PAM200
EXTRA_DIST += util/matrices/PAM210
EXTRA_DIST += util/matrices/PAM220
EXTRA_DIST += util/matrices/PAM230
EXTRA_DIST += util/matrices/PAM240
EXTRA_DIST += util/matrices/PAM250
EXTRA_DIST += util/matrices/PAM260
EXTRA_DIST += util/matrices/PAM270
EXTRA_DIST += util/matrices/PAM280
EXTRA_DIST += util/matrices/PAM290
EXTRA_DIST += util/matrices/PAM30
EXTRA_DIST += util/matrices/PAM300
EXTRA_DIST += util/matrices/PAM310
EXTRA_DIST += util/matrices/PAM320
EXTRA_DIST += util/matrices/PAM330
EXTRA_DIST += util/matrices/PAM340
EXTRA_DIST += util/matrices/PAM350
EXTRA_DIST += util/matrices/PAM360
EXTRA_DIST += util/matrices/PAM370
EXTRA_DIST += util/matrices/PAM380
EXTRA_DIST += util/matrices/PAM390
EXTRA_DIST += util/matrices/PAM40
EXTRA_DIST += util/matrices/PAM400
EXTRA_DIST += util/matrices/PAM410
EXTRA_DIST += util/matrices/PAM420
EXTRA_DIST += util/matrices/PAM430
EXTRA_DIST += util/matrices/PAM440
EXTRA_DIST += util/matrices/PAM450
EXTRA_DIST += util/matrices/PAM460
EXTRA_DIST += util/matrices/PAM470
EXTRA_DIST += util/matrices/PAM480
EXTRA_DIST += util/matrices/PAM490
EXTRA_DIST += util/matrices/PAM50
EXTRA_DIST += util/matrices/PAM500
EXTRA_DIST += util/matrices/PAM60
EXTRA_DIST += util/matrices/PAM70
EXTRA_DIST += util/matrices/PAM80
EXTRA_DIST += util/matrices/PAM90
EXTRA_DIST += util/matrices/README
EXTRA_DIST += util/matrices/codegen.py
EXTRA_DIST += util/templates/nw_diag.c
EXTRA_DIST += util/templates/nw_scan.c
EXTRA_DIST += util/templates/nw_stats_diag.c
EXTRA_DIST += util/templates/nw_stats_scan.c
EXTRA_DIST += util/templates/nw_stats_striped.c
EXTRA_DIST += util/templates/nw_striped.c
EXTRA_DIST += util/templates/sg_diag.c
EXTRA_DIST += util/templates/sg_diag_8.c
EXTRA_DIST += util/templates/sg_scan.c
EXTRA_DIST += util/templates/sg_stats_diag.c
EXTRA_DIST += util/templates/sg_stats_scan.c
EXTRA_DIST += util/templates/sg_stats_striped.c
EXTRA_DIST += util/templates/sg_striped.c
EXTRA_DIST += util/templates/sw_diag.c
EXTRA_DIST += util/templates/sw_diag_8.c
EXTRA_DIST += util/templates/sw_scan.c
EXTRA_DIST += util/templates/sw_stats_diag.c
EXTRA_DIST += util/templates/sw_stats_diag_8.c
EXTRA_DIST += util/templates/sw_stats_scan.c
EXTRA_DIST += util/templates/sw_stats_striped.c
EXTRA_DIST += util/templates/sw_stats_striped_bias.c
EXTRA_DIST += util/templates/sw_striped.c
EXTRA_DIST += util/templates/sw_striped_bias.c
