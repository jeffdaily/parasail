#!/usr/bin/env python

# serial reference implementations
alg = ["nw", "sg", "sw", "sg_qb", "sg_qe", "sg_qx", "sg_db", "sg_de", "sg_dx", "sg_qb_de", "sg_qe_db"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            print ""
            print "extern parasail_result_t* parasail_"+a+s+t+'('
            print " "*8+"const char * const restrict s1, const int s1Len,"
            print " "*8+"const char * const restrict s2, const int s2Len,"
            print " "*8+"const int open, const int gap,"
            print " "*8+"const parasail_matrix_t* matrix);"

# serial reference implementations of internal sg function
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
for s in stats:
    for t in table:
        if 'trace' in t and 'stats' in s: continue
        print ""
        print "extern parasail_result_t* parasail_sg_flags"+s+t+'('
        print " "*8+"const char * const restrict s1, const int s1Len,"
        print " "*8+"const char * const restrict s2, const int s2Len,"
        print " "*8+"const int open, const int gap,"
        print " "*8+"const parasail_matrix_t* matrix,"
        print " "*8+"int s1_beg, int s1_end, int s2_beg, int s2_end);"

# serial scan reference implementations
alg = ["nw", "sg", "sw", "sg_qb", "sg_qe", "sg_qx", "sg_db", "sg_de", "sg_dx", "sg_qb_de", "sg_qe_db"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            print ""
            print "extern parasail_result_t* parasail_"+a+s+t+'_scan('
            print " "*8+"const char * const restrict s1, const int s1Len,"
            print " "*8+"const char * const restrict s2, const int s2Len,"
            print " "*8+"const int open, const int gap,"
            print " "*8+"const parasail_matrix_t* matrix);"

# serial scan reference implementations of internal sg function
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
for s in stats:
    for t in table:
        if 'trace' in t and 'stats' in s: continue
        print ""
        print "extern parasail_result_t* parasail_sg_flags"+s+t+'_scan('
        print " "*8+"const char * const restrict s1, const int s1Len,"
        print " "*8+"const char * const restrict s2, const int s2Len,"
        print " "*8+"const int open, const int gap,"
        print " "*8+"const parasail_matrix_t* matrix,"
        print " "*8+"int s1_beg, int s1_end, int s2_beg, int s2_end);"

# vectorized implementations
alg = ["nw", "sg", "sw", "sg_qb", "sg_qe", "sg_qx", "sg_db", "sg_de", "sg_dx", "sg_qb_de", "sg_qe_db"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan", "_striped", "_diag"]
isa = [
    "_sse2_128_64", "_sse2_128_32", "_sse2_128_16", "_sse2_128_8", "_sse2_128_sat",
    "_sse41_128_64", "_sse41_128_32", "_sse41_128_16", "_sse41_128_8", "_sse41_128_sat",
    "_avx2_256_64", "_avx2_256_32", "_avx2_256_16", "_avx2_256_8", "_avx2_256_sat",
    "_altivec_128_64", "_altivec_128_32", "_altivec_128_16", "_altivec_128_8", "_altivec_128_sat",
    "_neon_128_64", "_neon_128_32", "_neon_128_16", "_neon_128_8", "_neon_128_sat"
    ]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            for p in par:
                for i in isa:
                    print ""
                    print "extern parasail_result_t* parasail_"+a+s+t+p+i+'('
                    print " "*8+"const char * const restrict s1, const int s1Len,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap,"
                    print " "*8+"const parasail_matrix_t* matrix);"

# vectorized implementations of internal sg function
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan", "_striped", "_diag"]
isa = [
    "_sse2_128_64", "_sse2_128_32", "_sse2_128_16", "_sse2_128_8", "_sse2_128_sat",
    "_sse41_128_64", "_sse41_128_32", "_sse41_128_16", "_sse41_128_8", "_sse41_128_sat",
    "_avx2_256_64", "_avx2_256_32", "_avx2_256_16", "_avx2_256_8", "_avx2_256_sat",
    "_altivec_128_64", "_altivec_128_32", "_altivec_128_16", "_altivec_128_8", "_altivec_128_sat",
    "_neon_128_64", "_neon_128_32", "_neon_128_16", "_neon_128_8", "_neon_128_sat"
    ]
for s in stats:
    for t in table:
        if 'trace' in t and 'stats' in s: continue
        for p in par:
            for i in isa:
                print ""
                print "extern parasail_result_t* parasail_sg_flags"+s+t+p+i+'('
                print " "*8+"const char * const restrict s1, const int s1Len,"
                print " "*8+"const char * const restrict s2, const int s2Len,"
                print " "*8+"const int open, const int gap,"
                print " "*8+"const parasail_matrix_t* matrix,"
                print " "*8+"int s1_beg, int s1_end, int s2_beg, int s2_end);"

# vectorized profile implementations
alg = ["nw", "sg", "sw", "sg_qb", "sg_qe", "sg_qx", "sg_db", "sg_de", "sg_dx", "sg_qb_de", "sg_qe_db"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan_profile", "_striped_profile"]
isa = [
    "_sse2_128_64", "_sse2_128_32", "_sse2_128_16", "_sse2_128_8", "_sse2_128_sat",
    "_sse41_128_64", "_sse41_128_32", "_sse41_128_16", "_sse41_128_8", "_sse41_128_sat",
    "_avx2_256_64", "_avx2_256_32", "_avx2_256_16", "_avx2_256_8", "_avx2_256_sat",
    "_altivec_128_64", "_altivec_128_32", "_altivec_128_16", "_altivec_128_8", "_altivec_128_sat",
    "_neon_128_64", "_neon_128_32", "_neon_128_16", "_neon_128_8", "_neon_128_sat"
    ]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            for p in par:
                for i in isa:
                    print ""
                    print "extern parasail_result_t* parasail_"+a+s+t+p+i+'('
                    print " "*8+"const parasail_profile_t * const restrict profile,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap);"

# vectorized profile implementations of internal sg function
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan_profile", "_striped_profile"]
isa = [
    "_sse2_128_64", "_sse2_128_32", "_sse2_128_16", "_sse2_128_8", "_sse2_128_sat",
    "_sse41_128_64", "_sse41_128_32", "_sse41_128_16", "_sse41_128_8", "_sse41_128_sat",
    "_avx2_256_64", "_avx2_256_32", "_avx2_256_16", "_avx2_256_8", "_avx2_256_sat",
    "_altivec_128_64", "_altivec_128_32", "_altivec_128_16", "_altivec_128_8", "_altivec_128_sat",
    "_neon_128_64", "_neon_128_32", "_neon_128_16", "_neon_128_8", "_neon_128_sat"
    ]
for s in stats:
    for t in table:
        if 'trace' in t and 'stats' in s: continue
        for p in par:
            for i in isa:
                print ""
                print "extern parasail_result_t* parasail_sg_flags"+s+t+p+i+'('
                print " "*8+"const parasail_profile_t * const restrict profile,"
                print " "*8+"const char * const restrict s2, const int s2Len,"
                print " "*8+"const int open, const int gap,"
                print " "*8+"int s1_beg, int s1_end, int s2_beg, int s2_end);"

# vectorized implementations of blocked
alg = ["sw"]
stats = [""]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_blocked"]
isa = ["_sse41_128_32", "_sse41_128_16"]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t: continue
            for p in par:
                for i in isa:
                    print ""
                    print "extern parasail_result_t* parasail_"+a+s+t+p+i+'('
                    print " "*8+"const char * const restrict s1, const int s1Len,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap,"
                    print " "*8+"const parasail_matrix_t* matrix);"
# dispatching implementations
alg = ["nw", "sg", "sw", "sg_qb", "sg_qe", "sg_qx", "sg_db", "sg_de", "sg_dx", "sg_qb_de", "sg_qe_db"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan", "_striped", "_diag"]
width = ["_64", "_32", "_16", "_8", "_sat"]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            for p in par:
                for w in width:
                    print ""
                    print "extern parasail_result_t* parasail_"+a+s+t+p+w+'('
                    print " "*8+"const char * const restrict s1, const int s1Len,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap,"
                    print " "*8+"const parasail_matrix_t* matrix);"

# dispatching implementations of internal sg function
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan", "_striped", "_diag"]
width = ["_64", "_32", "_16", "_8", "_sat"]
for s in stats:
    for t in table:
        if 'trace' in t and 'stats' in s: continue
        for p in par:
            for w in width:
                print ""
                print "extern parasail_result_t* parasail_sg_flags"+s+t+p+w+'('
                print " "*8+"const char * const restrict s1, const int s1Len,"
                print " "*8+"const char * const restrict s2, const int s2Len,"
                print " "*8+"const int open, const int gap,"
                print " "*8+"const parasail_matrix_t* matrix,"
                print " "*8+"int s1_beg, int s1_end, int s2_beg, int s2_end);"

# dispatching profile implementations
alg = ["nw", "sg", "sw", "sg_qb", "sg_qe", "sg_qx", "sg_db", "sg_de", "sg_dx", "sg_qb_de", "sg_qe_db"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan_profile", "_striped_profile"]
width = ["_64", "_32", "_16", "_8", "_sat"]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            for p in par:
                for w in width:
                    print ""
                    print "extern parasail_result_t* parasail_"+a+s+t+p+w+'('
                    print " "*8+"const parasail_profile_t * const restrict profile,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap);"

# dispatching profile implementations of internal sg function
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan_profile", "_striped_profile"]
width = ["_64", "_32", "_16", "_8", "_sat"]
for s in stats:
    for t in table:
        if 'trace' in t and 'stats' in s: continue
        for p in par:
            for w in width:
                print ""
                print "extern parasail_result_t* parasail_sg_flags"+s+t+p+w+'('
                print " "*8+"const parasail_profile_t * const restrict profile,"
                print " "*8+"const char * const restrict s2, const int s2Len,"
                print " "*8+"const int open, const int gap,"
                print " "*8+"int s1_beg, int s1_end, int s2_beg, int s2_end);"

# profile creation functions
stats = ["", "_stats"]
isa = [
    "_sse_128_64", "_sse_128_32", "_sse_128_16", "_sse_128_8", "_sse_128_sat",
    "_avx_256_64", "_avx_256_32", "_avx_256_16", "_avx_256_8", "_avx_256_sat",
    "_altivec_128_64", "_altivec_128_32", "_altivec_128_16", "_altivec_128_8", "_altivec_128_sat",
    "_neon_128_64", "_neon_128_32", "_neon_128_16", "_neon_128_8", "_neon_128_sat",
    "_64", "_32", "_16", "_8", "_sat"
    ]
for s in stats:
    for i in isa:
        print ""
        print "extern parasail_profile_t* parasail_profile_create"+s+i+'('
        print " "*8+"const char * const restrict s1, const int s1Len,"
        print " "*8+"const parasail_matrix_t* matrix);"

# dispatching saturation check implementations
alg = ["nw", "sg", "sw", "sg_qb", "sg_qe", "sg_qx", "sg_db", "sg_de", "sg_dx", "sg_qb_de", "sg_qe_db"]
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan", "_striped", "_diag"]
for a in alg:
    for s in stats:
        for t in table:
            if 'trace' in t and 'stats' in s: continue
            for p in par:
                print ""
                print "extern parasail_result_t* parasail_"+a+s+t+p+'_sat('
                print " "*8+"const char * const restrict s1, const int s1Len,"
                print " "*8+"const char * const restrict s2, const int s2Len,"
                print " "*8+"const int open, const int gap,"
                print " "*8+"const parasail_matrix_t* matrix);"

# dispatching saturation check implementations of internal sg function
stats = ["", "_stats"]
table = ["", "_table", "_rowcol", "_trace"]
par = ["_scan", "_striped", "_diag"]
for s in stats:
    for t in table:
        if 'trace' in t and 'stats' in s: continue
        for p in par:
            print ""
            print "extern parasail_result_t* parasail_sg_flags"+s+t+p+'_sat('
            print " "*8+"const char * const restrict s1, const int s1Len,"
            print " "*8+"const char * const restrict s2, const int s2Len,"
            print " "*8+"const int open, const int gap,"
            print " "*8+"const parasail_matrix_t* matrix,"
            print " "*8+"int s1_beg, int s1_end, int s2_beg, int s2_end);"

