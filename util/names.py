#!/usr/bin/env python

# serial reference implementations (3x2x2x2 = 24 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table"]
for a in alg:
    for s in stats:
        for t in table:
            print ""
            print "extern PARASAIL_API"
            print "parasail_result_t* parasail_"+a+s+t+'('
            print " "*8+"const char * const restrict s1, const int s1Len,"
            print " "*8+"const char * const restrict s2, const int s2Len,"
            print " "*8+"const int open, const int gap,"
            print " "*8+"const parasail_matrix_t* matrix);"

# serial scan reference implementations (3x2x2x2 = 24 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table"]
for a in alg:
    for s in stats:
        for t in table:
            print ""
            print "extern PARASAIL_API"
            print "parasail_result_t* parasail_"+a+s+t+'_scan('
            print " "*8+"const char * const restrict s1, const int s1Len,"
            print " "*8+"const char * const restrict s2, const int s2Len,"
            print " "*8+"const int open, const int gap,"
            print " "*8+"const parasail_matrix_t* matrix);"

# vectorized implementations (3x2x2x3x13 = 468 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table"]
par = ["_scan", "_striped", "_diag"]
isa = [
    "_sse2_128_64", "_sse2_128_32", "_sse2_128_16", "_sse2_128_8",
    "_sse41_128_64", "_sse41_128_32", "_sse41_128_16", "_sse41_128_8",
    "_avx2_256_64", "_avx2_256_32", "_avx2_256_16", "_avx2_256_8",
    "_knc_512_32"
    ]

for a in alg:
    for s in stats:
        for t in table:
            for p in par:
                for i in isa:
                    print ""
                    print "extern PARASAIL_API"
                    print "parasail_result_t* parasail_"+a+s+t+p+i+'('
                    print " "*8+"const char * const restrict s1, const int s1Len,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap,"
                    print " "*8+"const parasail_matrix_t* matrix);"

# vectorized implementations of blocked (just a couple)
alg = ["sw"]
stats = [""]
table = ["", "_table"]
par = ["_blocked"]
isa = ["_sse41_128_32", "_sse41_128_16"]

for a in alg:
    for s in stats:
        for t in table:
            for p in par:
                for i in isa:
                    print ""
                    print "extern PARASAIL_API"
                    print "parasail_result_t* parasail_"+a+s+t+p+i+'('
                    print " "*8+"const char * const restrict s1, const int s1Len,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap,"
                    print " "*8+"const parasail_matrix_t* matrix);"

# dispatching implementations (3x2x3x4 = 72 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table"]
par = ["_scan", "_striped", "_diag"]
width = ["", "_64", "_32", "_16", "_8"]

for a in alg:
    for s in stats:
        for t in table:
            for p in par:
                for w in width:
                    print ""
                    print "extern PARASAIL_API"
                    print "parasail_result_t* parasail_"+a+s+t+p+w+'('
                    print " "*8+"const char * const restrict s1, const int s1Len,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap,"
                    print " "*8+"const parasail_matrix_t* matrix);"

