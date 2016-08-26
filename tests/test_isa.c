#include "config.h"

#include <stdio.h>

#include "parasail/cpuid.h"

#define UNUSED(expr) do { (void)(expr); } while (0)

static inline char* yesno(int value)
{
    return value ? "yes" : "no ";
}


int main(int argc, char **argv)
{
    int cpu_sse2;
    int cpu_sse41;
    int cpu_avx2;
    int cpu_avx512;
    int cpu_knl;
    int cpu_knc;
    int cc_sse2;
    int cc_sse41;
    int cc_avx2;
    int cc_avx512;
    int cc_knl;
    int cc_knc;

    UNUSED(argc);
    UNUSED(argv);

    cpu_sse2 = parasail_can_use_sse2();
    cpu_sse41 = parasail_can_use_sse41();
    cpu_avx2 = parasail_can_use_avx2();
    cpu_avx512 = parasail_can_use_avx512();
    cpu_knl = parasail_can_use_knl();

#if HAVE_SSE2
    cc_sse2 = 1;
#else
    cc_sse2 = 0;
#endif

#if HAVE_SSE41
    cc_sse41 = 1;
#else
    cc_sse41 = 0;
#endif

#if HAVE_AVX2
    cc_avx2 = 1;
#else
    cc_avx2 = 0;
#endif

#if HAVE_AVX512
    cc_avx512 = 1;
#else
    cc_avx512 = 0;
#endif

#if HAVE_KNL
    cc_knl = 1;
#else
    cc_knl = 0;
#endif

#ifdef __MIC__
    cpu_knc = 1;
    cc_knc = 1;
#else
    cpu_knc = 0;
    cc_knc = 0;
#endif

    printf(" ISA    | Compiler | CPU  \n");
    printf("--------------------------\n");
    printf(" SSE2   |  %s     | %s\n", yesno(cc_sse2), yesno(cpu_sse2));
    printf(" SSE41  |  %s     | %s\n", yesno(cc_sse41), yesno(cpu_sse41));
    printf(" AVX2   |  %s     | %s\n", yesno(cc_avx2), yesno(cpu_avx2));
    printf(" AVX512 |  %s     | %s\n", yesno(cc_avx512), yesno(cpu_avx512));
    printf(" KNL    |  %s     | %s\n", yesno(cc_knl), yesno(cpu_knl));
    printf(" KNC    |  %s     | %s\n", yesno(cc_knc), yesno(cpu_knc));

    return 0;
}

