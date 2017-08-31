#!/usr/bin/env python

print """LIBRARY parasail
EXPORTS
; from parasail.h
    parasail_profile_free
    parasail_version
    parasail_result_free
    parasail_lookup_function
    parasail_lookup_pfunction
    parasail_lookup_pcreator
    parasail_lookup_function_info
    parasail_lookup_pfunction_info
    parasail_time
    parasail_matrix_lookup
    parasail_matrix_create
    parasail_matrix_from_file
    parasail_matrix_free
    parasail_matrix_copy
    parasail_matrix_set_value
    parasail_nw_banded
    parasail_traceback
    parasail_cigar
; from parasail/io.h
    parasail_open
    parasail_close
    parasail_is_fasta
    parasail_is_fastq
    parasail_stat
    parasail_stat_fasta
    parasail_stat_fastq
    parasail_read
    parasail_pack
    parasail_pack_fasta
    parasail_pack_fastq
    parasail_is_fasta_buffer
    parasail_is_fastq_buffer
    parasail_stat_buffer
    parasail_stat_fasta_buffer
    parasail_stat_fastq_buffer
    parasail_pack_buffer
    parasail_pack_fasta_buffer
    parasail_pack_fastq_buffer
; from parasail/cpuid.h
    parasail_can_use_avx512vbmi
    parasail_can_use_avx512bw
    parasail_can_use_avx512f
    parasail_can_use_avx2
    parasail_can_use_sse41
    parasail_can_use_sse2
; from parasail/memory.h (mostly internal functions)
    parasail_memalign
    parasail_memalign_int
    parasail_memalign_int8_t
    parasail_memalign_int16_t
    parasail_memalign_int32_t
    parasail_memalign_int64_t
    parasail_free
    parasail_memset
    parasail_memset_int
    parasail_memset_int8_t
    parasail_memset_int16_t
    parasail_memset_int32_t
    parasail_memset_int64_t
    parasail_result_new
    parasail_result_new_table1
    parasail_result_new_table3
    parasail_result_new_rowcol1
    parasail_result_new_rowcol3
    parasail_result_new_trace
    parasail_profile_new
    parasail_reverse
    parasail_striped_unwind
; from parasail.h, generated names"""

# serial reference implementations (3x2x3 = 18 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
for a in alg:
    for s in stats:
        for t in table:
            if 'stats' in s and 'trace' in t: continue
            print "    parasail_"+a+s+t

# serial scan reference implementations (3x2x3 = 18 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
for a in alg:
    for s in stats:
        for t in table:
            if 'stats' in s and 'trace' in t: continue
            print "    parasail_"+a+s+t+"_scan"

# vectorized implementations (3x2x3x3x13 = 702 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan", "_striped", "_diag"]
isa = [
    "_sse2_128_64", "_sse2_128_32", "_sse2_128_16", "_sse2_128_8", "_sse2_128_sat",
    "_sse41_128_64", "_sse41_128_32", "_sse41_128_16", "_sse41_128_8", "_sse41_128_sat",
    "_avx2_256_64", "_avx2_256_32", "_avx2_256_16", "_avx2_256_8", "_avx2_256_sat",
    "_knc_512_32"
    ]
for a in alg:
    for s in stats:
        for t in table:
            if 'stats' in s and 'trace' in t: continue
            for p in par:
                for i in isa:
                    print "    parasail_"+a+s+t+p+i

# vectorized profile implementations (3x2x3x2x13 = 468 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan_profile", "_striped_profile"]
isa = [
    "_sse2_128_64", "_sse2_128_32", "_sse2_128_16", "_sse2_128_8", "_sse2_128_sat",
    "_sse41_128_64", "_sse41_128_32", "_sse41_128_16", "_sse41_128_8", "_sse41_128_sat",
    "_avx2_256_64", "_avx2_256_32", "_avx2_256_16", "_avx2_256_8", "_avx2_256_sat",
    "_knc_512_32"
    ]
for a in alg:
    for s in stats:
        for t in table:
            if 'stats' in s and 'trace' in t: continue
            for p in par:
                for i in isa:
                    print "    parasail_"+a+s+t+p+i

# vectorized implementations of blocked (1x1x3x1x2 = 6 impl)
alg = ["sw"]
stats = [""]
table = ["", "_table", "_rowcol"]
par = ["_blocked"]
isa = ["_sse41_128_32", "_sse41_128_16"]
for a in alg:
    for s in stats:
        for t in table:
            for p in par:
                for i in isa:
                    print "    parasail_"+a+s+t+p+i

# dispatching implementations (3x2x3x3x4 = 216 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan", "_striped", "_diag"]
width = ["_64", "_32", "_16", "_8", "_sat"]
for a in alg:
    for s in stats:
        for t in table:
            if 'stats' in s and 'trace' in t: continue
            for p in par:
                for w in width:
                    print "    parasail_"+a+s+t+p+w

# dispatching profile implementations (3x2x3x2x4 = 144 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan_profile", "_striped_profile"]
width = ["_64", "_32", "_16", "_8", "_sat"]
for a in alg:
    for s in stats:
        for t in table:
            if 'stats' in s and 'trace' in t: continue
            for p in par:
                for w in width:
                    print "    parasail_"+a+s+t+p+w

# profile creation functions (2x13 = 26 impl)
stats = ["", "_stats"]
isa = [
    "_sse_128_64", "_sse_128_32", "_sse_128_16", "_sse_128_8", "_sse_128_sat",
    "_avx_256_64", "_avx_256_32", "_avx_256_16", "_avx_256_8", "_avx_256_sat",
    "_knc_512_32",
    "_64", "_32", "_16", "_8", "_sat"
    ]
for s in stats:
    for i in isa:
        print "    parasail_profile_create"+s+i

# dispatching saturation check implementations (3x2x3x3 = 54 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan", "_striped", "_diag"]
for a in alg:
    for s in stats:
        for t in table:
            if 'stats' in s and 'trace' in t: continue
            for p in par:
                print "    parasail_"+a+s+t+p+"_sat"

