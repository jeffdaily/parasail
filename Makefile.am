ACLOCAL_AMFLAGS = -I m4

include_HEADERS = parasail.h
nobase_include_HEADERS =
lib_LTLIBRARIES = libparasail.la
noinst_LTLIBRARIES =
check_PROGRAMS =
check_LTLIBRARIES = 
AM_CFLAGS =
AM_CPPFLAGS =
AM_LDFLAGS =
LDADD = libparasail.la
libparasail_la_LDFLAGS = $(AM_LDFLAGS) -version-info @LTVER@

nobase_include_HEADERS += parasail/cpuid.h
nobase_include_HEADERS += parasail/function_lookup.h
nobase_include_HEADERS += parasail/matrix_lookup.h
nobase_include_HEADERS += parasail/memory.h

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
if HAVE_SSE
noinst_LTLIBRARIES += libparasail_sse_memory.la
endif
if HAVE_SSE2
noinst_LTLIBRARIES += libparasail_sse2.la
noinst_LTLIBRARIES += libparasail_sse2_table.la
endif
if HAVE_SSE41
noinst_LTLIBRARIES += libparasail_sse41.la
noinst_LTLIBRARIES += libparasail_sse41_table.la
endif
if HAVE_AVX2
noinst_LTLIBRARIES += libparasail_avx2_memory.la
noinst_LTLIBRARIES += libparasail_avx2.la
noinst_LTLIBRARIES += libparasail_avx2_table.la
endif
if HAVE_KNC
noinst_LTLIBRARIES += libparasail_knc_memory.la
noinst_LTLIBRARIES += libparasail_knc.la
noinst_LTLIBRARIES += libparasail_knc_table.la
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

SRC_CORE += src/parasail_cpuid.h
SRC_CORE += src/parasail_internal.h
SRC_CORE += src/parasail_internal_sse.h
SRC_CORE += src/parasail_internal_avx.h
SRC_CORE += src/parasail_internal_knc.h
SRC_CORE += src/cpuid.c
SRC_CORE += src/memory.c
SRC_CORE += src/nw_dispatch.c
SRC_CORE += src/sg_dispatch.c
SRC_CORE += src/sw_dispatch.c

#################
# vector memory
#################

MEM_SRC_SSE += src/memory_sse.c

MEM_SRC_AVX2 += src/memory_avx2.c

MEM_SRC_KNC += src/memory_knc.c

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

SRC_SSE2 += src/nw_stats_scan_sse2_128_32.c
SRC_SSE2 += src/sg_stats_scan_sse2_128_32.c
SRC_SSE2 += src/sw_stats_scan_sse2_128_32.c

SRC_SSE2 += src/nw_stats_scan_sse2_128_16.c
SRC_SSE2 += src/sg_stats_scan_sse2_128_16.c
SRC_SSE2 += src/sw_stats_scan_sse2_128_16.c

SRC_SSE2 += src/nw_stats_scan_sse2_128_8.c
SRC_SSE2 += src/sg_stats_scan_sse2_128_8.c
SRC_SSE2 += src/sw_stats_scan_sse2_128_8.c

SRC_SSE41 += src/nw_stats_scan_sse41_128_32.c
SRC_SSE41 += src/sg_stats_scan_sse41_128_32.c
SRC_SSE41 += src/sw_stats_scan_sse41_128_32.c

SRC_SSE41 += src/nw_stats_scan_sse41_128_16.c
SRC_SSE41 += src/sg_stats_scan_sse41_128_16.c
SRC_SSE41 += src/sw_stats_scan_sse41_128_16.c

SRC_SSE41 += src/nw_stats_scan_sse41_128_8.c
SRC_SSE41 += src/sg_stats_scan_sse41_128_8.c
SRC_SSE41 += src/sw_stats_scan_sse41_128_8.c

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

SRC_SSE2 += src/nw_stats_diag_sse2_128_32.c
SRC_SSE2 += src/sg_stats_diag_sse2_128_32.c
SRC_SSE2 += src/sw_stats_diag_sse2_128_32.c

SRC_SSE2 += src/nw_stats_diag_sse2_128_16.c
SRC_SSE2 += src/sg_stats_diag_sse2_128_16.c
SRC_SSE2 += src/sw_stats_diag_sse2_128_16.c

SRC_SSE2 += src/nw_stats_diag_sse2_128_8.c
SRC_SSE2 += src/sg_stats_diag_sse2_128_8.c
SRC_SSE2 += src/sw_stats_diag_sse2_128_8.c

SRC_SSE41 += src/nw_stats_diag_sse41_128_32.c
SRC_SSE41 += src/sg_stats_diag_sse41_128_32.c
SRC_SSE41 += src/sw_stats_diag_sse41_128_32.c

SRC_SSE41 += src/nw_stats_diag_sse41_128_16.c
SRC_SSE41 += src/sg_stats_diag_sse41_128_16.c
SRC_SSE41 += src/sw_stats_diag_sse41_128_16.c

SRC_SSE41 += src/nw_stats_diag_sse41_128_8.c
SRC_SSE41 += src/sg_stats_diag_sse41_128_8.c
SRC_SSE41 += src/sw_stats_diag_sse41_128_8.c

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

SRC_SSE2 += src/nw_stats_striped_sse2_128_32.c
SRC_SSE2 += src/sg_stats_striped_sse2_128_32.c
SRC_SSE2 += src/sw_stats_striped_sse2_128_32.c

SRC_SSE2 += src/nw_stats_striped_sse2_128_16.c
SRC_SSE2 += src/sg_stats_striped_sse2_128_16.c
SRC_SSE2 += src/sw_stats_striped_sse2_128_16.c

SRC_SSE2 += src/nw_stats_striped_sse2_128_8.c
SRC_SSE2 += src/sg_stats_striped_sse2_128_8.c
SRC_SSE2 += src/sw_stats_striped_sse2_128_8.c

SRC_SSE41 += src/nw_stats_striped_sse41_128_32.c
SRC_SSE41 += src/sg_stats_striped_sse41_128_32.c
SRC_SSE41 += src/sw_stats_striped_sse41_128_32.c

SRC_SSE41 += src/nw_stats_striped_sse41_128_16.c
SRC_SSE41 += src/sg_stats_striped_sse41_128_16.c
SRC_SSE41 += src/sw_stats_striped_sse41_128_16.c

SRC_SSE41 += src/nw_stats_striped_sse41_128_8.c
SRC_SSE41 += src/sg_stats_striped_sse41_128_8.c
SRC_SSE41 += src/sw_stats_striped_sse41_128_8.c

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
if HAVE_SSE
libparasail_la_LIBADD += libparasail_sse_memory.la
endif
if HAVE_SSE2
libparasail_la_LIBADD += libparasail_sse2.la
libparasail_la_LIBADD += libparasail_sse2_table.la
endif
if HAVE_SSE41
libparasail_la_LIBADD += libparasail_sse41.la
libparasail_la_LIBADD += libparasail_sse41_table.la
endif
if HAVE_AVX2
libparasail_la_LIBADD += libparasail_avx2_memory.la
libparasail_la_LIBADD += libparasail_avx2.la
libparasail_la_LIBADD += libparasail_avx2_table.la
endif
if HAVE_KNC
libparasail_la_LIBADD += libparasail_knc_memory.la
libparasail_la_LIBADD += libparasail_knc.la
libparasail_la_LIBADD += libparasail_knc_table.la
endif

if HAVE_SSE2
check_LTLIBRARIES += tests/libssw.la
tests_libssw_la_SOURCES =
tests_libssw_la_SOURCES += contrib/ssw.c
tests_libssw_la_SOURCES += contrib/ssw.h
endif

check_LTLIBRARIES += tests/libsais.la
tests_libsais_la_SOURCES =
tests_libsais_la_SOURCES += contrib/sais-lite-lcp/sais.c
tests_libsais_la_SOURCES += contrib/sais-lite-lcp/sais.h

check_PROGRAMS += tests/test_align
check_PROGRAMS += tests/test_gcups
check_PROGRAMS += tests/test_isa
check_PROGRAMS += tests/test_openmp
check_PROGRAMS += tests/test_openmp_sais
check_PROGRAMS += tests/test_query
check_PROGRAMS += tests/test_scatter
check_PROGRAMS += tests/test_sw_dispatch
check_PROGRAMS += tests/test_verify
check_PROGRAMS += tests/test_verify_tables

LDADD += $(CLOCK_LIBS)
LDADD += $(MATH_LIBS)
LDADD += $(Z_LIBS)

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

tests_test_openmp_sais_SOURCES  = tests/test_openmp_sais.cpp
tests_test_openmp_sais_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
tests_test_openmp_sais_LDFLAGS  = $(AM_LDFLAGS) $(OPENMP_CXXFLAGS)
tests_test_openmp_sais_LDADD    =
tests_test_openmp_sais_LDADD    += $(LDADD)
tests_test_openmp_sais_LDADD    += tests/libsais.la

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

tests_test_sw_dispatch_SOURCES = tests/test_sw_dispatch.c

tests_test_verify_SOURCES = tests/test_verify.c
tests_test_verify_CFLAGS  = $(AM_CFLAGS) $(OPENMP_CFLAGS)
tests_test_verify_LDFLAGS = $(AM_LDFLAGS) $(OPENMP_CFLAGS)

tests_test_verify_tables_SOURCES = tests/test_verify_tables.c
tests_test_verify_tables_CFLAGS  = $(AM_CFLAGS) $(OPENMP_CFLAGS)
tests_test_verify_tables_LDFLAGS = $(AM_LDFLAGS) $(OPENMP_CFLAGS)

