#include "config.h"

#include <stdio.h>

#include "parasail/internal_altivec.h"

static inline void print_8(vec128i v)
{
    vec128i_8_t a;
    a.m = v;
    printf( "%u" ",%u" ",%u" ",%u"
           ",%u" ",%u" ",%u" ",%u"
           ",%u" ",%u" ",%u" ",%u"
           ",%u" ",%u" ",%u" ",%u",
           a.v[0], a.v[1], a.v[2], a.v[3],
           a.v[4], a.v[5], a.v[6], a.v[7],
           a.v[8], a.v[9], a.v[10], a.v[11],
           a.v[12], a.v[13], a.v[14], a.v[15] );
}

int main(int argc, char **argv)
{
    vec128i vOne8 = (vec128i)vec_splats((signed char)1);

    printf("vOne=");print_8(vOne8);printf("\n");

    return 0;
}
