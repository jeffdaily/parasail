#!/usr/bin/env python

# serial reference implementations (3x2x2x2 = 24 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table"]
for a in alg:
    for s in stats:
        for t in table:
            print ""
            print "parasail_result_t* "+a+s+t+'('
            print " "*8+"const char * const restrict s1, const int s1Len,"
            print " "*8+"const char * const restrict s2, const int s2Len,"
            print " "*8+"const int open, const int gap, const int matrix[24][24]);"

# serial scan reference implementations (3x2x2x2 = 24 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table"]
for a in alg:
    for s in stats:
        for t in table:
            print ""
            print "parasail_result_t* "+a+s+t+'_scan('
            print " "*8+"const char * const restrict s1, const int s1Len,"
            print " "*8+"const char * const restrict s2, const int s2Len,"
            print " "*8+"const int open, const int gap, const int matrix[24][24]);"

# vectorized implementations (3x2x2x2x3x7 = 504 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table"]
par = ["_scan", "_striped", "_blocked", "_diag"]
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
                    print "parasail_result_t* "+a+s+t+p+i+'('
                    print " "*8+"const char * const restrict s1, const int s1Len,"
                    print " "*8+"const char * const restrict s2, const int s2Len,"
                    print " "*8+"const int open, const int gap, const int matrix[24][24]);"

# dispatching implementations (3x2x3x4 = 72 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
par = ["_scan", "_striped", "_diag"]
width = ["", "_64", "_32", "_16", "_8"]

for a in alg:
    for s in stats:
        for p in par:
            for w in width:
                print ""
                print "parasail_result_t* parasail_"+a+s+p+w+'('
                print " "*8+"const char * const restrict s1, const int s1Len,"
                print " "*8+"const char * const restrict s2, const int s2Len,"
                print " "*8+"const int open, const int gap, const int matrix[24][24]);"

