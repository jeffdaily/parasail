#!/usr/bin/env python

import sys

# serial reference implementations
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
par = ["", "_scan"]
blosum = ["40","45","50","62","75","80","90"]
for a in alg:
    for s in stats:
        for p in par:
            for b in blosum:
                txt = "./test_openmp -a "+a+s+p+" -b blosum"+b+" -f "+sys.argv[1]
                print "echo '"+txt+"'"
                print txt

# vectorized implementations
alg = ["nw", "sg", "sw"]
stats = ["", "_stats"]
par = ["_scan", "_striped", "_diag"]
isa = ["_sse2_128_16", "_sse41_128_16", "_avx2_256_16"]
blosum = ["40","45","50","62","75","80","90"]
for a in alg:
    for s in stats:
        for p in par:
            for i in isa:
                for b in blosum:
                    txt = "./test_openmp -a "+a+s+p+i+" -b blosum"+b+" -f "+sys.argv[1]
                    print "echo '"+txt+"'"
                    print txt

