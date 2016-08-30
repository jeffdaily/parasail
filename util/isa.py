#!/usr/bin/python

sse2_header = """#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <emmintrin.h>
#endif"""

sse41_header = """#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <emmintrin.h>
#include <smmintrin.h>
#endif"""

sse2 = {
    "ISA"         : "sse",
    "ISA_VERSION" : "2",
    "HEADER"      : sse2_header,
    "BITS"        : 128,
    "VTYPE"       : "__m128i",
    "VADDx8"      : "_mm_add_epi8",
    "VADDx16"     : "_mm_add_epi16",
    "VADDx32"     : "_mm_add_epi32",
    "VADDx64"     : "_mm_add_epi64",
    "VADDSx8"     : "_mm_adds_epi8",
    "VADDSx16"    : "_mm_adds_epi16",
    "VBLEND"      : "_mm_blendv_epi8_rpl",
    "VCMPEQx8"    : "_mm_cmpeq_epi8",
    "VCMPEQx16"   : "_mm_cmpeq_epi16",
    "VCMPEQx32"   : "_mm_cmpeq_epi32",
    "VCMPEQx64"   : "_mm_cmpeq_epi64_rpl",
    "VCMPGTx8"    : "_mm_cmpgt_epi8",
    "VCMPGTx16"   : "_mm_cmpgt_epi16",
    "VCMPGTx32"   : "_mm_cmpgt_epi32",
    "VCMPGTx64"   : "_mm_cmpgt_epi64_rpl",
    "VCMPLTx8"    : "_mm_cmplt_epi8",
    "VCMPLTx16"   : "_mm_cmplt_epi16",
    "VCMPLTx32"   : "_mm_cmplt_epi32",
    "VCMPLTx64"   : "_mm_cmplt_epi64_rpl",
    "VEXTRACTx8"  : "_mm_extract_epi8_rpl",
    "VEXTRACTx16" : "_mm_extract_epi16",
    "VEXTRACTx32" : "_mm_extract_epi32_rpl",
    "VEXTRACTx64" : "_mm_extract_epi64_rpl",
    "VINSERTx8"   : "_mm_insert_epi8_rpl",
    "VINSERTx16"  : "_mm_insert_epi16",
    "VINSERTx32"  : "_mm_insert_epi32_rpl",
    "VINSERTx64"  : "_mm_insert_epi64_rpl",
    "VLOAD"       : "_mm_load_si128",
    "VHMAXx8"     : "_mm_hmax_epi8_rpl",
    "VHMAXx16"    : "_mm_hmax_epi16_rpl",
    "VHMAXx32"    : "_mm_hmax_epi32_rpl",
    "VHMAXx64"    : "_mm_hmax_epi64_rpl",
    "VMAXx8"      : "_mm_max_epi8_rpl",
    "VMAXx16"     : "_mm_max_epi16",
    "VMAXx32"     : "_mm_max_epi32_rpl",
    "VMAXx64"     : "_mm_max_epi64_rpl",
    "VMINx8"      : "_mm_min_epi8_rpl",
    "VMINx16"     : "_mm_min_epi16",
    "VMINx32"     : "_mm_min_epi32_rpl",
    "VMINx64"     : "_mm_min_epi64_rpl",
    "VMOVEMASK"   : "_mm_movemask_epi8",
    "VPACKS"      : "_mm_packs_epi16",
    "VUNPACKLO"   : "_mm_unpacklo_epi8",
    "VUNPACKHI"   : "_mm_unpackhi_epi8",
    "VAND"        : "_mm_and_si128",
    "VANDNOT"     : "_mm_andnot_si128",
    "VOR"         : "_mm_or_si128",
    "VSET0"       : "_mm_setzero_si128",
    "VSET1x8"     : "_mm_set1_epi8",
    "VSET1x16"    : "_mm_set1_epi16",
    "VSET1x32"    : "_mm_set1_epi32",
    "VSET1x64"    : "_mm_set1_epi64x_rpl",
    "VSETx8"      : "_mm_set_epi8",
    "VSETx16"     : "_mm_set_epi16",
    "VSETx32"     : "_mm_set_epi32",
    "VSETx64"     : "_mm_set_epi64x_rpl",
    "VRSHIFT"     : "_mm_srli_si128",
    "VSHIFT"      : "_mm_slli_si128",
    "VSTORE"      : "_mm_store_si128",
    "VSUBx8"      : "_mm_sub_epi8",
    "VSUBx16"     : "_mm_sub_epi16",
    "VSUBx32"     : "_mm_sub_epi32",
    "VSUBx64"     : "_mm_sub_epi64",
    "VSUBSx8"     : "_mm_subs_epi8",
    "VSUBSx16"    : "_mm_subs_epi16",
    "_mm_blendv_epi8_rpl" : """
static inline __m128i _mm_blendv_epi8_rpl(__m128i a, __m128i b, __m128i mask) {
    a = _mm_andnot_si128(mask, a);
    a = _mm_or_si128(a, _mm_and_si128(mask, b));
    return a;
}
""",
    "_mm_cmpeq_epi64_rpl" : """
static inline __m128i _mm_cmpeq_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]==B.v[0]) ? 0xFFFFFFFFFFFFFFFF : 0;
    A.v[1] = (A.v[1]==B.v[1]) ? 0xFFFFFFFFFFFFFFFF : 0;
    return A.m;
}
""",
    "_mm_cmpgt_epi64_rpl" : """
static inline __m128i _mm_cmpgt_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? 0xFFFFFFFFFFFFFFFF : 0;
    A.v[1] = (A.v[1]>B.v[1]) ? 0xFFFFFFFFFFFFFFFF : 0;
    return A.m;
}
""",
    "_mm_cmplt_epi64_rpl" : """
static inline __m128i _mm_cmplt_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]<B.v[0]) ? 0xFFFFFFFFFFFFFFFF : 0;
    A.v[1] = (A.v[1]<B.v[1]) ? 0xFFFFFFFFFFFFFFFF : 0;
    return A.m;
}
""",
    "_mm_extract_epi8_rpl" : """
static inline int8_t _mm_extract_epi8_rpl(__m128i a, const int imm) {
    __m128i_8_t A;
    A.m = a;
    return A.v[imm];
}
""",
    "_mm_extract_epi32_rpl" : """
static inline int32_t _mm_extract_epi32_rpl(__m128i a, const int imm) {
    __m128i_32_t A;
    A.m = a;
    return A.v[imm];
}
""",
    "_mm_extract_epi64_rpl" : """
static inline int64_t _mm_extract_epi64_rpl(__m128i a, const int imm) {
    __m128i_64_t A;
    A.m = a;
    return A.v[imm];
}
""",
    "_mm_insert_epi8_rpl" : """
static inline __m128i _mm_insert_epi8_rpl(__m128i a, int8_t i, const int imm) {
    __m128i_8_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
""",
    "_mm_insert_epi32_rpl" : """
static inline __m128i _mm_insert_epi32_rpl(__m128i a, int32_t i, const int imm) {
    __m128i_32_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
""",
    "_mm_insert_epi64_rpl" : """
static inline __m128i _mm_insert_epi64_rpl(__m128i a, int64_t i, const int imm) {
    __m128i_64_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
""",
    "_mm_hmax_epi8_rpl" : """
static inline int8_t _mm_hmax_epi8_rpl(__m128i a) {
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 8));
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 4));
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 2));
    a = _mm_max_epi8_rpl(a, _mm_srli_si128(a, 1));
    return _mm_extract_epi8_rpl(a, 0);
}
""",
    "_mm_hmax_epi16_rpl" : """
static inline int16_t _mm_hmax_epi16_rpl(__m128i a) {
    a = _mm_max_epi16(a, _mm_srli_si128(a, 8));
    a = _mm_max_epi16(a, _mm_srli_si128(a, 4));
    a = _mm_max_epi16(a, _mm_srli_si128(a, 2));
    return _mm_extract_epi16(a, 0);
}
""",
    "_mm_hmax_epi32_rpl" : """
static inline int32_t _mm_hmax_epi32_rpl(__m128i a) {
    a = _mm_max_epi32_rpl(a, _mm_srli_si128(a, 8));
    a = _mm_max_epi32_rpl(a, _mm_srli_si128(a, 4));
    return _mm_extract_epi32_rpl(a, 0);
}
""",
    "_mm_hmax_epi64_rpl" : """
static inline int64_t _mm_hmax_epi64_rpl(__m128i a) {
    a = _mm_max_epi64_rpl(a, _mm_srli_si128(a, 8));
    return _mm_extract_epi64_rpl(a, 0);
}
""",
    "_mm_max_epi8_rpl" : """
static inline __m128i _mm_max_epi8_rpl(__m128i a, __m128i b) {
    __m128i mask = _mm_cmpgt_epi8(a, b);
    a = _mm_and_si128(a, mask);
    b = _mm_andnot_si128(mask, b);
    return _mm_or_si128(a, b);
}
""",
    "_mm_max_epi32_rpl" : """
static inline __m128i _mm_max_epi32_rpl(__m128i a, __m128i b) {
    __m128i mask = _mm_cmpgt_epi32(a, b);
    a = _mm_and_si128(a, mask);
    b = _mm_andnot_si128(mask, b);
    return _mm_or_si128(a, b);
}
""",
    "_mm_max_epi64_rpl" : """
static inline __m128i _mm_max_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]>B.v[1]) ? A.v[1] : B.v[1];
    return A.m;
}
""",
    "_mm_min_epi8_rpl" : """
static inline __m128i _mm_min_epi8_rpl(__m128i a, __m128i b) {
    __m128i mask = _mm_cmpgt_epi8(b, a);
    a = _mm_and_si128(a, mask);
    b = _mm_andnot_si128(mask, b);
    return _mm_or_si128(a, b);
}
""",
    "_mm_min_epi32_rpl" : """
static inline __m128i _mm_min_epi32_rpl(__m128i a, __m128i b) {
    __m128i mask = _mm_cmpgt_epi32(b, a);
    a = _mm_and_si128(a, mask);
    b = _mm_andnot_si128(mask, b);
    return _mm_or_si128(a, b);
}
""",
    "_mm_min_epi64_rpl" : """
static inline __m128i _mm_min_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]<B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]<B.v[1]) ? A.v[1] : B.v[1];
    return A.m;
}
""",
    "_mm_set1_epi64x_rpl" : """
#if HAVE_SSE2_MM_SET1_EPI64X
#define _mm_set1_epi64x_rpl _mm_set1_epi64x
#else
static inline __m128i _mm_set1_epi64x_rpl(int64_t i) {
    __m128i_64_t A;
    A.v[0] = i;
    A.v[1] = i;
    return A.m;
}
#endif
""",
    "_mm_set_epi64x_rpl" : """
#if HAVE_SSE2_MM_SET_EPI64X
#define _mm_set_epi64x_rpl _mm_set_epi64x
#else
static inline __m128i _mm_set_epi64x_rpl(int64_t e1, int64_t e0) {
    __m128i_64_t A;
    A.v[0] = e0;
    A.v[1] = e1;
    return A.m;
}
#endif
""",
}


sse41 = {
    "ISA"         : "sse",
    "ISA_VERSION" : "41",
    "HEADER"      : sse41_header,
    "BITS"        : 128,
    "VTYPE"       : "__m128i",
    "VADDx8"      : "_mm_add_epi8",
    "VADDx16"     : "_mm_add_epi16",
    "VADDx32"     : "_mm_add_epi32",
    "VADDx64"     : "_mm_add_epi64",
    "VADDSx8"     : "_mm_adds_epi8",
    "VADDSx16"    : "_mm_adds_epi16",
    "VBLEND"      : "_mm_blendv_epi8",
    "VCMPEQx8"    : "_mm_cmpeq_epi8",
    "VCMPEQx16"   : "_mm_cmpeq_epi16",
    "VCMPEQx32"   : "_mm_cmpeq_epi32",
    "VCMPEQx64"   : "_mm_cmpeq_epi64",
    "VCMPGTx8"    : "_mm_cmpgt_epi8",
    "VCMPGTx16"   : "_mm_cmpgt_epi16",
    "VCMPGTx32"   : "_mm_cmpgt_epi32",
    "VCMPGTx64"   : "_mm_cmpgt_epi64_rpl",
    "VCMPLTx8"    : "_mm_cmplt_epi8",
    "VCMPLTx16"   : "_mm_cmplt_epi16",
    "VCMPLTx32"   : "_mm_cmplt_epi32",
    "VCMPLTx64"   : "_mm_cmplt_epi64_rpl",
    "VEXTRACTx8"  : "_mm_extract_epi8",
    "VEXTRACTx16" : "_mm_extract_epi16",
    "VEXTRACTx32" : "_mm_extract_epi32",
    "VEXTRACTx64" : "_mm_extract_epi64_rpl",
    "VINSERTx8"   : "_mm_insert_epi8",
    "VINSERTx16"  : "_mm_insert_epi16",
    "VINSERTx32"  : "_mm_insert_epi32",
    "VINSERTx64"  : "_mm_insert_epi64_rpl",
    "VLOAD"       : "_mm_load_si128",
    "VHMAXx8"     : "_mm_hmax_epi8_rpl",
    "VHMAXx16"    : "_mm_hmax_epi16_rpl",
    "VHMAXx32"    : "_mm_hmax_epi32_rpl",
    "VHMAXx64"    : "_mm_hmax_epi64_rpl",
    "VMAXx8"      : "_mm_max_epi8",
    "VMAXx16"     : "_mm_max_epi16",
    "VMAXx32"     : "_mm_max_epi32",
    "VMAXx64"     : "_mm_max_epi64_rpl",
    "VMINx8"      : "_mm_min_epi8",
    "VMINx16"     : "_mm_min_epi16",
    "VMINx32"     : "_mm_min_epi32",
    "VMINx64"     : "_mm_min_epi64_rpl",
    "VMOVEMASK"   : "_mm_movemask_epi8",
    "VPACKS"      : "_mm_packs_epi16",
    "VUNPACKLO"   : "_mm_unpacklo_epi8",
    "VUNPACKHI"   : "_mm_unpackhi_epi8",
    "VAND"        : "_mm_and_si128",
    "VANDNOT"     : "_mm_andnot_si128",
    "VOR"         : "_mm_or_si128",
    "VSET0"       : "_mm_setzero_si128",
    "VSET1x8"     : "_mm_set1_epi8",
    "VSET1x16"    : "_mm_set1_epi16",
    "VSET1x32"    : "_mm_set1_epi32",
    "VSET1x64"    : "_mm_set1_epi64x_rpl",
    "VSETx8"      : "_mm_set_epi8",
    "VSETx16"     : "_mm_set_epi16",
    "VSETx32"     : "_mm_set_epi32",
    "VSETx64"     : "_mm_set_epi64x_rpl",
    "VRSHIFT"     : "_mm_srli_si128",
    "VSHIFT"      : "_mm_slli_si128",
    "VSTORE"      : "_mm_store_si128",
    "VSUBx8"      : "_mm_sub_epi8",
    "VSUBx16"     : "_mm_sub_epi16",
    "VSUBx32"     : "_mm_sub_epi32",
    "VSUBx64"     : "_mm_sub_epi64",
    "VSUBSx8"     : "_mm_subs_epi8",
    "VSUBSx16"    : "_mm_subs_epi16",
    "_mm_cmpgt_epi64_rpl" : """
static inline __m128i _mm_cmpgt_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? 0xFFFFFFFFFFFFFFFF : 0;
    A.v[1] = (A.v[1]>B.v[1]) ? 0xFFFFFFFFFFFFFFFF : 0;
    return A.m;
}
""",
    "_mm_cmplt_epi64_rpl" : """
static inline __m128i _mm_cmplt_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]<B.v[0]) ? 0xFFFFFFFFFFFFFFFF : 0;
    A.v[1] = (A.v[1]<B.v[1]) ? 0xFFFFFFFFFFFFFFFF : 0;
    return A.m;
}
""",
    "_mm_hmax_epi8_rpl" : """
static inline int8_t _mm_hmax_epi8_rpl(__m128i a) {
    a = _mm_max_epi8(a, _mm_srli_si128(a, 8));
    a = _mm_max_epi8(a, _mm_srli_si128(a, 4));
    a = _mm_max_epi8(a, _mm_srli_si128(a, 2));
    a = _mm_max_epi8(a, _mm_srli_si128(a, 1));
    return _mm_extract_epi8(a, 0);
}
""",
    "_mm_hmax_epi16_rpl" : """
static inline int16_t _mm_hmax_epi16_rpl(__m128i a) {
    a = _mm_max_epi16(a, _mm_srli_si128(a, 8));
    a = _mm_max_epi16(a, _mm_srli_si128(a, 4));
    a = _mm_max_epi16(a, _mm_srli_si128(a, 2));
    return _mm_extract_epi16(a, 0);
}
""",
    "_mm_hmax_epi32_rpl" : """
static inline int32_t _mm_hmax_epi32_rpl(__m128i a) {
    a = _mm_max_epi32(a, _mm_srli_si128(a, 8));
    a = _mm_max_epi32(a, _mm_srli_si128(a, 4));
    return _mm_extract_epi32(a, 0);
}
""",
    "_mm_hmax_epi64_rpl" : """
static inline int64_t _mm_hmax_epi64_rpl(__m128i a) {
    a = _mm_max_epi64_rpl(a, _mm_srli_si128(a, 8));
    return _mm_extract_epi64_rpl(a, 0);
}
""",
    "_mm_max_epi64_rpl" : """
static inline __m128i _mm_max_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]>B.v[1]) ? A.v[1] : B.v[1];
    return A.m;
}
""",
    "_mm_min_epi64_rpl" : """
static inline __m128i _mm_min_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]<B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]<B.v[1]) ? A.v[1] : B.v[1];
    return A.m;
}
""",
    "_mm_extract_epi64_rpl" : """
#if HAVE_SSE41_MM_EXTRACT_EPI64
#define _mm_extract_epi64_rpl _mm_extract_epi64
#else
static inline int64_t _mm_extract_epi64_rpl(__m128i a, int imm) {
    __m128i_64_t A;
    A.m = a;
    return A.v[imm];
}
#endif
""",
    "_mm_insert_epi64_rpl" : """
static inline __m128i _mm_insert_epi64_rpl(__m128i a, int64_t i, const int imm) {
    __m128i_64_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
""",
    "_mm_set1_epi64x_rpl" : """
#if HAVE_SSE2_MM_SET1_EPI64X
#define _mm_set1_epi64x_rpl _mm_set1_epi64x
#else
static inline __m128i _mm_set1_epi64x_rpl(int64_t i) {
    __m128i_64_t A;
    A.v[0] = i;
    A.v[1] = i;
    return A.m;
}
#endif
""",
    "_mm_set_epi64x_rpl" : """
#if HAVE_SSE2_MM_SET_EPI64X
#define _mm_set_epi64x_rpl _mm_set_epi64x
#else
static inline __m128i _mm_set_epi64x_rpl(int64_t e1, int64_t e0) {
    __m128i_64_t A;
    A.v[0] = e0;
    A.v[1] = e1;
    return A.m;
}
#endif
""",
}


avx2 = {
    "ISA"         : "avx",
    "ISA_VERSION" : "2",
    "HEADER"      : "#include <immintrin.h>",
    "BITS"        : 256,
    "VTYPE"       : "__m256i",
    "VADDx8"      : "_mm256_add_epi8",
    "VADDx16"     : "_mm256_add_epi16",
    "VADDx32"     : "_mm256_add_epi32",
    "VADDx64"     : "_mm256_add_epi64",
    "VADDSx8"     : "_mm256_adds_epi8",
    "VADDSx16"    : "_mm256_adds_epi16",
    "VBLEND"      : "_mm256_blendv_epi8",
    "VCMPEQx8"    : "_mm256_cmpeq_epi8",
    "VCMPEQx16"   : "_mm256_cmpeq_epi16",
    "VCMPEQx32"   : "_mm256_cmpeq_epi32",
    "VCMPEQx64"   : "_mm256_cmpeq_epi64",
    "VCMPGTx8"    : "_mm256_cmpgt_epi8",
    "VCMPGTx16"   : "_mm256_cmpgt_epi16",
    "VCMPGTx32"   : "_mm256_cmpgt_epi32",
    "VCMPGTx64"   : "_mm256_cmpgt_epi64",
    "VCMPLTx8"    : "_mm256_cmplt_epi8_rpl",
    "VCMPLTx16"   : "_mm256_cmplt_epi16_rpl",
    "VCMPLTx32"   : "_mm256_cmplt_epi32_rpl",
    "VCMPLTx64"   : "_mm256_cmplt_epi64_rpl",
    "VEXTRACTx8"  : "_mm256_extract_epi8_rpl",
    "VEXTRACTx16" : "_mm256_extract_epi16_rpl",
    "VEXTRACTx32" : "_mm256_extract_epi32_rpl",
    "VEXTRACTx64" : "_mm256_extract_epi64_rpl",
    "VINSERTx8"   : "_mm256_insert_epi8_rpl",
    "VINSERTx16"  : "_mm256_insert_epi16_rpl",
    "VINSERTx32"  : "_mm256_insert_epi32_rpl",
    "VINSERTx64"  : "_mm256_insert_epi64_rpl",
    "VLOAD"       : "_mm256_load_si256",
    "VHMAXx8"     : "_mm256_hmax_epi8_rpl",
    "VHMAXx16"    : "_mm256_hmax_epi16_rpl",
    "VHMAXx32"    : "_mm256_hmax_epi32_rpl",
    "VHMAXx64"    : "_mm256_hmax_epi64_rpl",
    "VMAXx8"      : "_mm256_max_epi8",
    "VMAXx16"     : "_mm256_max_epi16",
    "VMAXx32"     : "_mm256_max_epi32",
    "VMAXx64"     : "_mm256_max_epi64_rpl",
    "VMINx8"      : "_mm256_min_epi8",
    "VMINx16"     : "_mm256_min_epi16",
    "VMINx32"     : "_mm256_min_epi32",
    "VMINx64"     : "_mm256_min_epi64_rpl",
    "VMOVEMASK"   : "_mm256_movemask_epi8",
    "VPACKS"      : "_mm256_packs_epi16_rpl",
    "VUNPACKLO"   : "_mm256_unpacklo_epi8_rpl",
    "VUNPACKHI"   : "_mm256_unpackhi_epi8_rpl",
    "VAND"        : "_mm256_and_si256",
    "VANDNOT"     : "_mm256_andnot_si256",
    "VOR"         : "_mm256_or_si256",
    "VSET0"       : "_mm256_setzero_si256",
    "VSET1x8"     : "_mm256_set1_epi8",
    "VSET1x16"    : "_mm256_set1_epi16",
    "VSET1x32"    : "_mm256_set1_epi32",
    "VSET1x64"    : "_mm256_set1_epi64x_rpl",
    "VSETx8"      : "_mm256_set_epi8",
    "VSETx16"     : "_mm256_set_epi16",
    "VSETx32"     : "_mm256_set_epi32",
    "VSETx64"     : "_mm256_set_epi64x_rpl",
    "VRSHIFT"     : "_mm256_srli_si256_rpl",
    "VSHIFT"      : "_mm256_slli_si256_rpl",
    "VSTORE"      : "_mm256_store_si256",
    "VSUBx8"      : "_mm256_sub_epi8",
    "VSUBx16"     : "_mm256_sub_epi16",
    "VSUBx32"     : "_mm256_sub_epi32",
    "VSUBx64"     : "_mm256_sub_epi64",
    "VSUBSx8"     : "_mm256_subs_epi8",
    "VSUBSx16"    : "_mm256_subs_epi16",
    "_mm256_cmplt_epi8_rpl" : """
#define _mm256_cmplt_epi8_rpl(a,b) _mm256_cmpgt_epi8(b,a)
""",
    "_mm256_cmplt_epi16_rpl" : """
#define _mm256_cmplt_epi16_rpl(a,b) _mm256_cmpgt_epi16(b,a)
""",
    "_mm256_cmplt_epi32_rpl" : """
#define _mm256_cmplt_epi32_rpl(a,b) _mm256_cmpgt_epi32(b,a)
""",
    "_mm256_cmplt_epi64_rpl" : """
#define _mm256_cmplt_epi64_rpl(a,b) _mm256_cmpgt_epi64(b,a)
""",
    "_mm256_hmax_epi8_rpl" : """
static inline int8_t _mm256_hmax_epi8_rpl(__m256i a) {
    a = _mm256_max_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 8));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 4));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 2));
    a = _mm256_max_epi8(a, _mm256_slli_si256(a, 1));
    return _mm256_extract_epi8_rpl(a, 31);
}
""",
    "_mm256_hmax_epi16_rpl" : """
static inline int16_t _mm256_hmax_epi16_rpl(__m256i a) {
    a = _mm256_max_epi16(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
    a = _mm256_max_epi16(a, _mm256_slli_si256(a, 8));
    a = _mm256_max_epi16(a, _mm256_slli_si256(a, 4));
    a = _mm256_max_epi16(a, _mm256_slli_si256(a, 2));
    return _mm256_extract_epi16_rpl(a, 15);
}
""",
    "_mm256_hmax_epi32_rpl" : """
static inline int32_t _mm256_hmax_epi32_rpl(__m256i a) {
    a = _mm256_max_epi32(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
    a = _mm256_max_epi32(a, _mm256_slli_si256(a, 8));
    a = _mm256_max_epi32(a, _mm256_slli_si256(a, 4));
    return _mm256_extract_epi32_rpl(a, 7);
}
""",
    "_mm256_hmax_epi64_rpl" : """
static inline int64_t _mm256_hmax_epi64_rpl(__m256i a) {
    a = _mm256_max_epi64_rpl(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,0,0)));
    a = _mm256_max_epi64_rpl(a, _mm256_slli_si256(a, 8));
    return _mm256_extract_epi64_rpl(a, 3);
}
""",
    "_mm256_max_epi64_rpl" : """
static inline __m256i _mm256_max_epi64_rpl(__m256i a, __m256i b) {
    __m256i_64_t A;
    __m256i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]>B.v[1]) ? A.v[1] : B.v[1];
    A.v[2] = (A.v[2]>B.v[2]) ? A.v[2] : B.v[2];
    A.v[3] = (A.v[3]>B.v[3]) ? A.v[3] : B.v[3];
    return A.m;
}
""",
    "_mm256_min_epi64_rpl" : """
static inline __m256i _mm256_min_epi64_rpl(__m256i a, __m256i b) {
    __m256i_64_t A;
    __m256i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]<B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]<B.v[1]) ? A.v[1] : B.v[1];
    A.v[2] = (A.v[2]<B.v[2]) ? A.v[2] : B.v[2];
    A.v[3] = (A.v[3]<B.v[3]) ? A.v[3] : B.v[3];
    return A.m;
}
""",
    "_mm256_packs_epi16_rpl" : """
static inline __m256i _mm256_packs_epi16_rpl(__m256i a, __m256i b) {
    return _mm256_permute4x64_epi64(
            _mm256_packs_epi16(a, b),
            _MM_SHUFFLE(3,1,2,0));
}
""",
    "_mm256_unpacklo_epi8_rpl" : """
static inline __m256i _mm256_unpacklo_epi8_rpl(__m256i a, __m256i b) {
    __m256i an = _mm256_permute4x64_epi64(a, _MM_SHUFFLE(1,1,0,0));
    __m256i bn = _mm256_permute4x64_epi64(b, _MM_SHUFFLE(1,1,0,0));
    return _mm256_unpacklo_epi8(an, bn);
}
""",
    "_mm256_unpackhi_epi8_rpl" : """
static inline __m256i _mm256_unpackhi_epi8_rpl(__m256i a, __m256i b) {
    __m256i an = _mm256_permute4x64_epi64(a, _MM_SHUFFLE(3,3,2,2));
    __m256i bn = _mm256_permute4x64_epi64(b, _MM_SHUFFLE(3,3,2,2));
    return _mm256_unpackhi_epi8(an, bn);
}
""",
    "_mm256_srli_si256_rpl" : """
#define _mm256_srli_si256_rpl(a,imm) _mm256_or_si256(_mm256_slli_si256(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(3,0,0,1)), 16-imm), _mm256_srli_si256(a, imm))
""",
    "_mm256_slli_si256_rpl" : """
#define _mm256_slli_si256_rpl(a,imm) _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)
""",
    "_mm256_insert_epi64_rpl" : """
#if HAVE_AVX2_MM256_INSERT_EPI64
#define _mm256_insert_epi64_rpl _mm256_insert_epi64
#else
static inline __m256i _mm256_insert_epi64_rpl(__m256i a, int64_t i, int imm) {
    __m256i_64_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif
""",
    "_mm256_insert_epi32_rpl" : """
#if HAVE_AVX2_MM256_INSERT_EPI32
#define _mm256_insert_epi32_rpl _mm256_insert_epi32
#else
static inline __m256i _mm256_insert_epi32_rpl(__m256i a, int32_t i, int imm) {
    __m256i_32_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif
""",
    "_mm256_insert_epi16_rpl" : """
#if HAVE_AVX2_MM256_INSERT_EPI16
#define _mm256_insert_epi16_rpl _mm256_insert_epi16
#else
static inline __m256i _mm256_insert_epi16_rpl(__m256i a, int16_t i, int imm) {
    __m256i_16_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif
""",
    "_mm256_insert_epi8_rpl" : """
#if HAVE_AVX2_MM256_INSERT_EPI8
#define _mm256_insert_epi8_rpl _mm256_insert_epi8
#else
static inline __m256i _mm256_insert_epi8_rpl(__m256i a, int8_t i, int imm) {
    __m256i_8_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif
""",
    "_mm256_extract_epi64_rpl" : """
#if HAVE_AVX2_MM256_EXTRACT_EPI64
#define _mm256_extract_epi64_rpl _mm256_extract_epi64
#else
static inline int64_t _mm256_extract_epi64_rpl(__m256i a, int imm) {
    __m256i_64_t A;
    A.m = a;
    return A.v[imm];
}
#endif
""",
    "_mm256_extract_epi32_rpl" : """
#if HAVE_AVX2_MM256_EXTRACT_EPI32
#define _mm256_extract_epi32_rpl _mm256_extract_epi32
#else
static inline int32_t _mm256_extract_epi32_rpl(__m256i a, int imm) {
    __m256i_32_t A;
    A.m = a;
    return A.v[imm];
}
#endif
""",
    "_mm256_extract_epi16_rpl" : """
#if HAVE_AVX2_MM256_EXTRACT_EPI16
#define _mm256_extract_epi16_rpl _mm256_extract_epi16
#else
static inline int16_t _mm256_extract_epi16_rpl(__m256i a, int imm) {
    __m256i_16_t A;
    A.m = a;
    return A.v[imm];
}
#endif
""",
    "_mm256_extract_epi8_rpl" : """
#if HAVE_AVX2_MM256_EXTRACT_EPI8
#define _mm256_extract_epi8_rpl _mm256_extract_epi8
#else
static inline int8_t _mm256_extract_epi8_rpl(__m256i a, int imm) {
    __m256i_8_t A;
    A.m = a;
    return A.v[imm];
}
#endif
""",
    "_mm256_set1_epi64x_rpl" : """
#if HAVE_AVX2_MM256_SET1_EPI64X
#define _mm256_set1_epi64x_rpl _mm256_set1_epi64x
#else
static inline __m256i _mm256_set1_epi64x_rpl(int64_t i) {
    __m256i_64_t A;
    A.v[0] = i;
    A.v[1] = i;
    A.v[2] = i;
    A.v[3] = i;
    return A.m;
}
#endif
""",
    "_mm256_set_epi64x_rpl" : """
#if HAVE_AVX2_MM256_SET_EPI64X
#define _mm256_set_epi64x_rpl _mm256_set_epi64x
#else
static inline __m256i _mm256_set_epi64x_rpl(int64_t e3, int64_t e2, int64_t e1, int64_t e0) {
    __m256i_64_t A;
    A.v[0] = e0;
    A.v[1] = e1;
    A.v[2] = e2;
    A.v[3] = e3;
    return A.m;
}
#endif
""",
}


avx512 = {
    "ISA"         : "avx",
    "ISA_VERSION" : "512",
    "HEADER"      : "#include <immintrin.h>",
    "BITS"        : 512,
    "VTYPE"       : "__m512i",
    "VADDx8"      : "_mm512_add_epi8",
    "VADDx16"     : "_mm512_add_epi16",
    "VADDx32"     : "_mm512_add_epi32",
    "VADDx64"     : "_mm512_add_epi64",
    "VADDSx8"     : "_mm512_adds_epi8",
    "VADDSx16"    : "_mm512_adds_epi16",
    "VBLENDx8"    : "_mm512_mask_blend_epi8",
    "VBLENDx16"   : "_mm512_mask_blend_epi16",
    "VBLENDx32"   : "_mm512_mask_blend_epi32",
    "VBLENDx64"   : "_mm512_mask_blend_epi64",
    "VCMPEQx8"    : "_mm512_cmpeq_epi8_mask",
    "VCMPEQx16"   : "_mm512_cmpeq_epi16_mask",
    "VCMPEQx32"   : "_mm512_cmpeq_epi32_mask",
    "VCMPEQx64"   : "_mm512_cmpeq_epi64_mask",
    "VCMPGTx8"    : "_mm512_cmpgt_epi8_mask",
    "VCMPGTx16"   : "_mm512_cmpgt_epi16_mask",
    "VCMPGTx32"   : "_mm512_cmpgt_epi32_mask",
    "VCMPGTx64"   : "_mm512_cmpgt_epi64_mask",
    "VCMPLTx8"    : "_mm512_cmplt_epi8_mask",
    "VCMPLTx16"   : "_mm512_cmplt_epi16_mask",
    "VCMPLTx32"   : "_mm512_cmplt_epi32_mask",
    "VCMPLTx64"   : "_mm512_cmplt_epi64_mask",
    "VEXTRACTx8"  : "_mm512_extract_epi8_rpl",
    "VEXTRACTx16" : "_mm512_extract_epi16_rpl",
    "VEXTRACTx32" : "_mm512_extract_epi32_rpl",
    "VEXTRACTx64" : "_mm512_extract_epi64_rpl",
    "VINSERTx8"   : "_mm512_insert_epi8_rpl",
    "VINSERTx16"  : "_mm512_insert_epi16_rpl",
    "VINSERTx32"  : "_mm512_insert_epi32_rpl",
    "VINSERTx64"  : "_mm512_insert_epi64_rpl",
    "VLOAD"       : "_mm512_load_si512",
    "VHMAXx8"     : "_mm512_hmax_epi8_rpl",
    "VHMAXx16"    : "_mm512_hmax_epi16_rpl",
    "VHMAXx32"    : "_mm512_hmax_epi32_rpl",
    "VHMAXx64"    : "_mm512_hmax_epi64_rpl",
    "VMAXx8"      : "_mm512_max_epi8",
    "VMAXx16"     : "_mm512_max_epi16",
    "VMAXx32"     : "_mm512_max_epi32",
    "VMAXx64"     : "_mm512_max_epi64",
    "VMINx8"      : "_mm512_min_epi8",
    "VMINx16"     : "_mm512_min_epi16",
    "VMINx32"     : "_mm512_min_epi32",
    "VMINx64"     : "_mm512_min_epi64",
    "VMOVEMASK"   : "(unsigned long long)",
    "VKANDx8"     : "_mm512_kand",
    "VKANDx16"    : "_mm512_kand",
    "VKANDx32"    : "_mm512_kand",
    "VKANDx64"    : "_mm512_kand",
    "VKANDNOTx8"  : "_mm512_kandn",
    "VKANDNOTx16" : "_mm512_kandn",
    "VKANDNOTx32" : "_mm512_kandn",
    "VKANDNOTx64" : "_mm512_kandn",
    "VKORx8"      : "_mm512_kor",
    "VKORx16"     : "_mm512_kor",
    "VKORx32"     : "_mm512_kor",
    "VKORx64"     : "_mm512_kor",
    "VSET0"       : "_mm512_setzero_si512",
    "VSET1x8"     : "_mm512_set1_epi8",
    "VSET1x16"    : "_mm512_set1_epi16",
    "VSET1x32"    : "_mm512_set1_epi32",
    "VSET1x64"    : "_mm512_set1_epi64",
    "VSETx8"      : "_mm512_set_epi8",
    "VSETx16"     : "_mm512_set_epi16",
    "VSETx32"     : "_mm512_set_epi32",
    "VSETx64"     : "_mm512_set_epi64",
    "VRSHIFT"     : "_mm512_srli_si512_rpl",
    "VSHIFT"      : "_mm512_slli_si512_rpl",
    "VSTORE"      : "_mm512_store_si512",
    "VSUBx8"      : "_mm512_sub_epi8",
    "VSUBx16"     : "_mm512_sub_epi16",
    "VSUBx32"     : "_mm512_sub_epi32",
    "VSUBx64"     : "_mm512_sub_epi64",
    "VSUBSx8"     : "_mm512_subs_epi8",
    "VSUBSx16"    : "_mm512_subs_epi16",
    "_mm512_hmax_epi8_rpl" : """
static inline int8_t _mm512_hmax_epi8_rpl(__m512i a) {
    /* TODO */
    assert(0);
    return 0;
}
""",
    "_mm512_hmax_epi16_rpl" : """
static inline int16_t _mm512_hmax_epi16_rpl(__m512i a) {
    /* TODO */
    assert(0);
    return 0;
}
""",
    "_mm512_hmax_epi32_rpl" : """
static inline int32_t _mm512_hmax_epi32_rpl(__m512i a) {
    __m512i idx1 = _mm512_set_epi32(7,6,5,4,3,2,1,0,15,14,13,12,11,10,9,8);
    __m512i idx2 = _mm512_set_epi32(14,15,12,13,10,11,8,9,6,7,4,5,2,3,0,1);
    /* first, shuffle across 256-bit lanes */
    a = _mm512_max_epi32(a, _mm512_permutexvar_epi32(idx1, a));
    /* shift within 256-bit lanes, 128-bit shuffle */
    a = _mm512_max_epi32(a, _mm512_permutex_epi64(a, _MM_SHUFFLE(1,0,3,2)));
    /* shift within 256-bit lanes, 64-bit shuffle */
    a = _mm512_max_epi32(a, _mm512_permutex_epi64(a, _MM_SHUFFLE(2,3,0,1)));
    /* shift within 256-bit lanes, 32-bit shuffle */
    a = _mm512_max_epi32(a, _mm512_permutexvar_epi32(idx2, a));
    /* result is in any element of either 256-bit lane */
    return _mm512_extract_epi32_rpl(a, 15);
}
""",
    "_mm512_hmax_epi64_rpl" : """
static inline int64_t _mm512_hmax_epi64_rpl(__m512i a) {
    __m512i idx;

    /* first, shuffle across 256-bit lanes */
    idx = _mm512_set_epi64(3,2,1,0,7,6,5,4);
    a = _mm512_max_epi64(a, _mm512_permutexvar_epi64(idx, a));
    /* shift within 256-bit lanes, 128-bit shuffle */
    a = _mm512_max_epi64(a, _mm512_permutex_epi64(a, _MM_SHUFFLE(1,0,3,2)));
    /* shift within 256-bit lanes, 64-bit shuffle */
    a = _mm512_max_epi64(a, _mm512_permutex_epi64(a, _MM_SHUFFLE(2,3,0,1)));
    /* result is in any element of either 256-bit lane */
    return _mm512_extract_epi64_rpl(a, 7);
}
""",
    "_mm512_packs_epi16_rpl" : """
static inline __m512i _mm512_packs_epi16_rpl(__m512i a, __m512i b) {
    return _mm512_permute4x64_epi64(
            _mm512_packs_epi16(a, b),
            _MM_SHUFFLE(3,1,2,0));
}
""",
    "_mm512_unpacklo_epi8_rpl" : """
static inline __m512i _mm512_unpacklo_epi8_rpl(__m512i a, __m512i b) {
    __m512i an = _mm512_permute4x64_epi64(a, _MM_SHUFFLE(1,1,0,0));
    __m512i bn = _mm512_permute4x64_epi64(b, _MM_SHUFFLE(1,1,0,0));
    return _mm512_unpacklo_epi8(an, bn);
}
""",
    "_mm512_unpackhi_epi8_rpl" : """
static inline __m512i _mm512_unpackhi_epi8_rpl(__m512i a, __m512i b) {
    __m512i an = _mm512_permute4x64_epi64(a, _MM_SHUFFLE(3,3,2,2));
    __m512i bn = _mm512_permute4x64_epi64(b, _MM_SHUFFLE(3,3,2,2));
    return _mm512_unpackhi_epi8(an, bn);
}
""",
    "_mm512_srli_si512_rpl" : """
#define _mm512_srli_si512_rpl(a,imm) _mm512_or_si512(_mm512_slli_si512(_mm512_permute2x128_si512(a, a, _MM_SHUFFLE(3,0,0,1)), 16-imm), _mm512_srli_si512(a, imm))
""",
    "_mm512_slli_si512_rpl" : """
#define _mm512_slli_si512_rpl(a,imm) _mm512_alignr_epi8(a, _mm512_permute2x128_si512(a, a, _MM_SHUFFLE(0,0,3,0)), 16-imm)
""",
    "_mm512_insert_epi64_rpl" : """
#if HAVE_AVX2_MM512_INSERT_EPI64
#define _mm512_insert_epi64_rpl _mm512_insert_epi64
#else
static inline __m512i _mm512_insert_epi64_rpl(__m512i a, int64_t i, int imm) {
    __m512i_64_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif
""",
    "_mm512_insert_epi32_rpl" : """
#if HAVE_AVX2_MM512_INSERT_EPI32
#define _mm512_insert_epi32_rpl _mm512_insert_epi32
#else
static inline __m512i _mm512_insert_epi32_rpl(__m512i a, int32_t i, int imm) {
    __m512i_32_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif
""",
    "_mm512_insert_epi16_rpl" : """
#if HAVE_AVX2_MM512_INSERT_EPI16
#define _mm512_insert_epi16_rpl _mm512_insert_epi16
#else
static inline __m512i _mm512_insert_epi16_rpl(__m512i a, int16_t i, int imm) {
    __m512i_16_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif
""",
    "_mm512_insert_epi8_rpl" : """
#if HAVE_AVX2_MM512_INSERT_EPI8
#define _mm512_insert_epi8_rpl _mm512_insert_epi8
#else
static inline __m512i _mm512_insert_epi8_rpl(__m512i a, int8_t i, int imm) {
    __m512i_8_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif
""",
    "_mm512_extract_epi64_rpl" : """
#if HAVE_AVX2_MM512_EXTRACT_EPI64
#define _mm512_extract_epi64_rpl _mm512_extract_epi64
#else
static inline int64_t _mm512_extract_epi64_rpl(__m512i a, int imm) {
    __m512i_64_t A;
    A.m = a;
    return A.v[imm];
}
#endif
""",
    "_mm512_extract_epi32_rpl" : """
#if HAVE_AVX2_MM512_EXTRACT_EPI32
#define _mm512_extract_epi32_rpl _mm512_extract_epi32
#else
static inline int32_t _mm512_extract_epi32_rpl(__m512i a, int imm) {
    __m512i_32_t A;
    A.m = a;
    return A.v[imm];
}
#endif
""",
    "_mm512_extract_epi16_rpl" : """
#if HAVE_AVX2_MM512_EXTRACT_EPI16
#define _mm512_extract_epi16_rpl _mm512_extract_epi16
#else
static inline int16_t _mm512_extract_epi16_rpl(__m512i a, int imm) {
    __m512i_16_t A;
    A.m = a;
    return A.v[imm];
}
#endif
""",
    "_mm512_extract_epi8_rpl" : """
#if HAVE_AVX2_MM512_EXTRACT_EPI8
#define _mm512_extract_epi8_rpl _mm512_extract_epi8
#else
static inline int8_t _mm512_extract_epi8_rpl(__m512i a, int imm) {
    __m512i_8_t A;
    A.m = a;
    return A.v[imm];
}
#endif
""",
}


isa = {
        "sse2"  : sse2,
        "sse41" : sse41,
        "avx2"  : avx2,
        "avx512"  : avx512,
}
