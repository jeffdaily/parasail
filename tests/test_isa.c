#include "config.h"

#include <stdio.h>

#include "parasail/cpuid.h"

#define UNUSED(expr) do { (void)(expr); } while (0)

int main(int argc, char **argv)
{
    UNUSED(argc);
    UNUSED(argv);

    if (parasail_can_use_sse2()) {
        printf("sse2 supported\n");
    }
    else {
        printf("sse2 NOT supported\n");
    }

    if (parasail_can_use_sse41()) {
        printf("sse41 supported\n");
    }
    else {
        printf("sse41 NOT supported\n");
    }

    if (parasail_can_use_avx2()) {
        printf("avx2 supported\n");
    }
    else {
        printf("avx2 NOT supported\n");
    }

#ifdef __MIC__
    printf("knc supported\n");
#else
    printf("knc NOT supported\n");
#endif

    return 0;
}

