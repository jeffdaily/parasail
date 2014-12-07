#!/usr/bin/env python

# serial reference implementations (3x2x2x2 = 24 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table"]
ext = ["", "_ext"]
for a in alg:
    for s in stats:
        for t in table:
            for e in ext:
                print ""
                print "void "+a+s+t+e+'('
                print " "*8+"const char * const restrict s1, const int s1Len,"
                print " "*8+"const char * const restrict s2, const int s2Len,"
                print " "*8+"const int open, const int gap, const int matrix[24][24],"
                if e:
                    print " "*8+"parasail_result_t *result,"
                    print " "*8+"parasail_workspace_t *workspace);"
                else:
                    print " "*8+"parasail_result_t *result);"

# serial scan reference implementations (3x2x2x2 = 24 impl)
alg = ["nw_scan", "sg_scan", "sw_scan"]
stats = ["", "_stats"]
table = ["", "_table"]
ext = ["", "_ext"]
for a in alg:
    for s in stats:
        for t in table:
            for e in ext:
                print ""
                print "void "+a+s+t+e+'('
                print " "*8+"const char * const restrict s1, const int s1Len,"
                print " "*8+"const char * const restrict s2, const int s2Len,"
                print " "*8+"const int open, const int gap, const int matrix[24][24],"
                if e:
                    print " "*8+"parasail_result_t *result,"
                    print " "*8+"parasail_workspace_t *workspace);"
                else:
                    print " "*8+"parasail_result_t *result);"

# vectorized implementations (3x2x2x2x3x7 = 504 impl)
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
table = ["", "_table"]
ext = ["", "_ext"]
par = ["_scan", "_striped", "_diag"]
isa = [
    "_sse2_128_16",
    "_sse41_128_32", "_sse41_128_8",
    "_avx2_256_32", "_avx2_256_16", "_avx2_256_8",
    "_knc_512_32"
    ]

for a in alg:
    for s in stats:
        for t in table:
            for e in ext:
                for p in par:
                    for i in isa:
                        print ""
                        print "void "+a+s+t+e+p+i+'('
                        print " "*8+"const char * const restrict s1, const int s1Len,"
                        print " "*8+"const char * const restrict s2, const int s2Len,"
                        print " "*8+"const int open, const int gap, const int matrix[24][24],"
                        if e:
                            print " "*8+"parasail_result_t *result,"
                            print " "*8+"parasail_workspace_t *workspace);"
                        else:
                            print " "*8+"parasail_result_t *result);"

