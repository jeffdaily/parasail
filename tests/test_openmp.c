#include "config.h"

#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#include <omp.h>

#include "kseq.h"
KSEQ_INIT(int, read)

#if HAVE_SSE2
#include "ssw.h"
#endif

#include "parasail.h"
#include "parasail_internal.h"
#include "blosum/blosum40.h"
#include "blosum/blosum45.h"
#include "blosum/blosum50.h"
#include "blosum/blosum62.h"
#include "blosum/blosum75.h"
#include "blosum/blosum80.h"
#include "blosum/blosum90.h"
#include "timer.h"
#include "timer_real.h"

typedef parasail_result_t* (*pf)(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24]);

#if HAVE_SSE2
parasail_result_t* parasail_ssw_(
        const char * const restrict s1, const int s1_len,
        const char * const restrict s2, const int s2_len,
        const int open, const int gap, const int matrix[24][24],
        int score_size)
{
    parasail_result_t *result = parasail_result_new();
    s_profile *profile = NULL;
    int8_t *s1_num = (int8_t*)malloc(sizeof(int8_t) * s1_len);
    int8_t *s2_num = (int8_t*)malloc(sizeof(int8_t) * s2_len);
    s_align *ssw_result = NULL;
    size_t m = 0;

    /* This table is used to transform amino acid letters into numbers. */
    static const int8_t table[128] = {
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
        14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
        23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
        14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
    };

    /* initialize score matrix */
    for (m = 0; m < s1_len; ++m) s1_num[m] = table[(int)s1[m]];
    for (m = 0; m < s2_len; ++m) s2_num[m] = table[(int)s2[m]];
    profile = ssw_init(s1_num, s1_len, blosum62__, 24, score_size);
    ssw_result = ssw_align(profile, s2_num, s2_len, -open, -gap, 2, 0, 0, s1_len/2);
    result->score = ssw_result->score1;
    result->saturated = ssw_result->saturated;
    align_destroy(ssw_result);
    init_destroy(profile);

    return result;
}

parasail_result_t* parasail_ssw(
        const char * const restrict s1, const int s1_len,
        const char * const restrict s2, const int s2_len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_ssw_(s1, s1_len, s2, s2_len, open, gap, matrix, 2);
}

parasail_result_t* parasail_ssw_16(
        const char * const restrict s1, const int s1_len,
        const char * const restrict s2, const int s2_len,
        const int open, const int gap, const int matrix[24][24])
{
    return parasail_ssw_(s1, s1_len, s2, s2_len, open, gap, matrix, 1);
}
#endif

parasail_result_t* parasail_sw(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const int matrix[24][24])
{
    int saturated = 0;
    parasail_result_t *result;

    result = parasail_sw_scan_8(s1, s1Len, s2, s2Len, open, gap, matrix);
    if (result->saturated) {
        saturated = 1;
        parasail_result_free(result);
        result = parasail_sw_scan_16(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    if (result->saturated) {
        parasail_result_free(result);
        result = parasail_sw_scan_32(s1, s1Len, s2, s2Len, open, gap, matrix);
    }
    result->saturated = saturated;

    return result;
}

typedef struct blosum {
    const char * name;
    const int (*blosum)[24];
} blosum_t;

blosum_t blosums[] = {
    {"blosum40",blosum40},
    {"blosum45",blosum45},
    {"blosum50",blosum50},
    {"blosum62",blosum62},
    {"blosum75",blosum75},
    {"blosum80",blosum80},
    {"blosum90",blosum90},
    {"NULL",NULL},
};

typedef struct func {
    const char * name;
    pf f;
} func_t;

func_t functions[] = {
    {"nw",                           nw,                        },
    {"nw_scan",                      nw_scan,                   },
#if HAVE_SSE2
    {"nw_scan_sse2_128_32",          nw_scan_sse2_128_32,       },
    {"nw_scan_sse2_128_16",          nw_scan_sse2_128_16,       },
    {"nw_scan_sse2_128_8",           nw_scan_sse2_128_8,        },
    {"nw_diag_sse2_128_32",          nw_diag_sse2_128_32,       },
    {"nw_diag_sse2_128_16",          nw_diag_sse2_128_16,       },
    {"nw_diag_sse2_128_8",           nw_diag_sse2_128_8,        },
    {"nw_striped_sse2_128_32",       nw_striped_sse2_128_32,    },
    {"nw_striped_sse2_128_16",       nw_striped_sse2_128_16,    },
    {"nw_striped_sse2_128_8",        nw_striped_sse2_128_8,     },
#endif
#if HAVE_SSE41
    {"nw_scan_sse41_128_32",         nw_scan_sse41_128_32,      },
    {"nw_scan_sse41_128_16",         nw_scan_sse41_128_16,      },
    {"nw_scan_sse41_128_8",          nw_scan_sse41_128_8,       },
    {"nw_diag_sse41_128_32",         nw_diag_sse41_128_32,      },
    {"nw_diag_sse41_128_16",         nw_diag_sse41_128_16,      },
    {"nw_diag_sse41_128_8",          nw_diag_sse41_128_8,       },
    {"nw_striped_sse41_128_32",      nw_striped_sse41_128_32,   },
    {"nw_striped_sse41_128_16",      nw_striped_sse41_128_16,   },
    {"nw_striped_sse41_128_8",       nw_striped_sse41_128_8,    },
#endif
#if HAVE_AVX2
    {"nw_scan_avx2_256_32",          nw_scan_avx2_256_32,       },
    {"nw_scan_avx2_256_16",          nw_scan_avx2_256_16,       },
    {"nw_scan_avx2_256_8",           nw_scan_avx2_256_8,        },
    {"nw_diag_avx2_256_32",          nw_diag_avx2_256_32,       },
    {"nw_diag_avx2_256_16",          nw_diag_avx2_256_16,       },
    {"nw_diag_avx2_256_8",           nw_diag_avx2_256_8,        },
    {"nw_striped_avx2_256_32",       nw_striped_avx2_256_32,    },
    {"nw_striped_avx2_256_16",       nw_striped_avx2_256_16,    },
    {"nw_striped_avx2_256_8",        nw_striped_avx2_256_8,     },
#endif
#if HAVE_KNC
    {"nw_scan_knc_512_32",           nw_scan_knc_512_32,        },
    {"nw_diag_knc_512_32",           nw_diag_knc_512_32,        },
    {"nw_striped_knc_512_32",        nw_striped_knc_512_32,     },
#endif

    {"sg",                           sg,                        },
    {"sg_scan",                      sg_scan,                   },
#if HAVE_SSE2
    {"sg_scan_sse2_128_32",          sg_scan_sse2_128_32,       },
    {"sg_scan_sse2_128_16",          sg_scan_sse2_128_16,       },
    {"sg_scan_sse2_128_8",           sg_scan_sse2_128_8,        },
    {"sg_diag_sse2_128_32",          sg_diag_sse2_128_32,       },
    {"sg_diag_sse2_128_16",          sg_diag_sse2_128_16,       },
    {"sg_diag_sse2_128_8",           sg_diag_sse2_128_8,        },
    {"sg_striped_sse2_128_32",       sg_striped_sse2_128_32,    },
    {"sg_striped_sse2_128_16",       sg_striped_sse2_128_16,    },
    {"sg_striped_sse2_128_8",        sg_striped_sse2_128_8,     },
#endif
#if HAVE_SSE41
    {"sg_scan_sse41_128_32",         sg_scan_sse41_128_32,      },
    {"sg_scan_sse41_128_16",         sg_scan_sse41_128_16,      },
    {"sg_scan_sse41_128_8",          sg_scan_sse41_128_8,       },
    {"sg_diag_sse41_128_32",         sg_diag_sse41_128_32,      },
    {"sg_diag_sse41_128_16",         sg_diag_sse41_128_16,      },
    {"sg_diag_sse41_128_8",          sg_diag_sse41_128_8,       },
    {"sg_striped_sse41_128_32",      sg_striped_sse41_128_32,   },
    {"sg_striped_sse41_128_16",      sg_striped_sse41_128_16,   },
    {"sg_striped_sse41_128_8",       sg_striped_sse41_128_8,    },
#endif
#if HAVE_AVX2
    {"sg_scan_avx2_256_32",          sg_scan_avx2_256_32,       },
    {"sg_scan_avx2_256_16",          sg_scan_avx2_256_16,       },
    {"sg_scan_avx2_256_8",           sg_scan_avx2_256_8,        },
    {"sg_diag_avx2_256_32",          sg_diag_avx2_256_32,       },
    {"sg_diag_avx2_256_16",          sg_diag_avx2_256_16,       },
    {"sg_diag_avx2_256_8",           sg_diag_avx2_256_8,        },
    {"sg_striped_avx2_256_32",       sg_striped_avx2_256_32,    },
    {"sg_striped_avx2_256_16",       sg_striped_avx2_256_16,    },
    {"sg_striped_avx2_256_8",        sg_striped_avx2_256_8,     },
#endif
#if HAVE_KNC
    {"sg_scan_knc_512_32",           sg_scan_knc_512_32,        },
    {"sg_diag_knc_512_32",           sg_diag_knc_512_32,        },
    {"sg_striped_knc_512_32",        sg_striped_knc_512_32,     },
#endif

    {"sw",                           sw,                        },
    {"sw_scan",                      sw_scan,                   },
#if HAVE_SSE2
    {"sw_scan_sse2_128_32",          sw_scan_sse2_128_32,       },
    {"sw_scan_sse2_128_16",          sw_scan_sse2_128_16,       },
    {"sw_scan_sse2_128_8",           sw_scan_sse2_128_8,        },
    {"sw_diag_sse2_128_32",          sw_diag_sse2_128_32,       },
    {"sw_diag_sse2_128_16",          sw_diag_sse2_128_16,       },
    {"sw_diag_sse2_128_8",           sw_diag_sse2_128_8,        },
    {"sw_striped_sse2_128_32",       sw_striped_sse2_128_32,    },
    {"sw_striped_sse2_128_16",       sw_striped_sse2_128_16,    },
    {"sw_striped_sse2_128_8",        sw_striped_sse2_128_8,     },
    {"ssw_16",                       parasail_ssw_16,           },
    {"ssw_8",                        parasail_ssw,              },
#endif
#if HAVE_SSE41
    {"sw_scan_sse41_128_32",         sw_scan_sse41_128_32,      },
    {"sw_scan_sse41_128_16",         sw_scan_sse41_128_16,      },
    {"sw_scan_sse41_128_8",          sw_scan_sse41_128_8,       },
    {"sw_diag_sse41_128_32",         sw_diag_sse41_128_32,      },
    {"sw_diag_sse41_128_16",         sw_diag_sse41_128_16,      },
    {"sw_diag_sse41_128_8",          sw_diag_sse41_128_8,       },
    {"sw_striped_sse41_128_32",      sw_striped_sse41_128_32,   },
    {"sw_striped_sse41_128_16",      sw_striped_sse41_128_16,   },
    {"sw_striped_sse41_128_8",       sw_striped_sse41_128_8,    },
#endif
#if HAVE_AVX2
    {"sw_scan_avx2_256_32",          sw_scan_avx2_256_32,       },
    {"sw_scan_avx2_256_16",          sw_scan_avx2_256_16,       },
    {"sw_scan_avx2_256_8",           sw_scan_avx2_256_8,        },
    {"sw_diag_avx2_256_32",          sw_diag_avx2_256_32,       },
    {"sw_diag_avx2_256_16",          sw_diag_avx2_256_16,       },
    {"sw_diag_avx2_256_8",           sw_diag_avx2_256_8,        },
    {"sw_striped_avx2_256_32",       sw_striped_avx2_256_32,    },
    {"sw_striped_avx2_256_16",       sw_striped_avx2_256_16,    },
    {"sw_striped_avx2_256_8",        sw_striped_avx2_256_8,     },
#endif
#if HAVE_KNC
    {"sw_scan_knc_512_32",           sw_scan_knc_512_32,        },
    {"sw_diag_knc_512_32",           sw_diag_knc_512_32,        },
    {"sw_striped_knc_512_32",        sw_striped_knc_512_32,     },
#endif


    {"nw_stats",                     nw_stats,                     },
    {"nw_stats_scan",                nw_stats_scan,                },
#if HAVE_SSE2
    {"nw_stats_scan_sse2_128_32",    nw_stats_scan_sse2_128_32,    },
    {"nw_stats_scan_sse2_128_16",    nw_stats_scan_sse2_128_16,    },
    {"nw_stats_scan_sse2_128_8",     nw_stats_scan_sse2_128_8,     },
    {"nw_stats_diag_sse2_128_32",    nw_stats_diag_sse2_128_32,    },
    {"nw_stats_diag_sse2_128_16",    nw_stats_diag_sse2_128_16,    },
    {"nw_stats_diag_sse2_128_8",     nw_stats_diag_sse2_128_8,     },
    {"nw_stats_striped_sse2_128_32", nw_stats_striped_sse2_128_32, },
    {"nw_stats_striped_sse2_128_16", nw_stats_striped_sse2_128_16, },
    {"nw_stats_striped_sse2_128_8",  nw_stats_striped_sse2_128_8,  },
#endif
#if HAVE_SSE41
    {"nw_stats_scan_sse41_128_32",   nw_stats_scan_sse41_128_32,   },
    {"nw_stats_scan_sse41_128_16",   nw_stats_scan_sse41_128_16,   },
    {"nw_stats_scan_sse41_128_8",    nw_stats_scan_sse41_128_8,    },
    {"nw_stats_diag_sse41_128_32",   nw_stats_diag_sse41_128_32,   },
    {"nw_stats_diag_sse41_128_16",   nw_stats_diag_sse41_128_16,   },
    {"nw_stats_diag_sse41_128_8",    nw_stats_diag_sse41_128_8,    },
    {"nw_stats_striped_sse41_128_32",nw_stats_striped_sse41_128_32,},
    {"nw_stats_striped_sse41_128_16",nw_stats_striped_sse41_128_16,},
    {"nw_stats_striped_sse41_128_8", nw_stats_striped_sse41_128_8, },
#endif
#if HAVE_AVX2
    {"nw_stats_scan_avx2_256_32",    nw_stats_scan_avx2_256_32,    },
    {"nw_stats_scan_avx2_256_16",    nw_stats_scan_avx2_256_16,    },
    {"nw_stats_scan_avx2_256_8",     nw_stats_scan_avx2_256_8,     },
    {"nw_stats_diag_avx2_256_32",    nw_stats_diag_avx2_256_32,    },
    {"nw_stats_diag_avx2_256_16",    nw_stats_diag_avx2_256_16,    },
    {"nw_stats_diag_avx2_256_8",     nw_stats_diag_avx2_256_8,     },
    {"nw_stats_striped_avx2_256_32", nw_stats_striped_avx2_256_32, },
    {"nw_stats_striped_avx2_256_16", nw_stats_striped_avx2_256_16, },
    {"nw_stats_striped_avx2_256_8",  nw_stats_striped_avx2_256_8,  },
#endif
#if HAVE_KNC
    {"nw_stats_scan_knc_512_32",     nw_stats_scan_knc_512_32,     },
    {"nw_stats_diag_knc_512_32",     nw_stats_diag_knc_512_32,     },
    {"nw_stats_striped_knc_512_32",  nw_stats_striped_knc_512_32,  },
#endif

    {"sg_stats",                     sg_stats,                     },
    {"sg_stats_scan",                sg_stats_scan,                },
#if HAVE_SSE2
    {"sg_stats_scan_sse2_128_32",    sg_stats_scan_sse2_128_32,    },
    {"sg_stats_scan_sse2_128_16",    sg_stats_scan_sse2_128_16,    },
    {"sg_stats_scan_sse2_128_8",     sg_stats_scan_sse2_128_8,     },
    {"sg_stats_diag_sse2_128_32",    sg_stats_diag_sse2_128_32,    },
    {"sg_stats_diag_sse2_128_16",    sg_stats_diag_sse2_128_16,    },
    {"sg_stats_diag_sse2_128_8",     sg_stats_diag_sse2_128_8,     },
    {"sg_stats_striped_sse2_128_32", sg_stats_striped_sse2_128_32, },
    {"sg_stats_striped_sse2_128_16", sg_stats_striped_sse2_128_16, },
    {"sg_stats_striped_sse2_128_8",  sg_stats_striped_sse2_128_8,  },
#endif
#if HAVE_SSE41
    {"sg_stats_scan_sse41_128_32",   sg_stats_scan_sse41_128_32,   },
    {"sg_stats_scan_sse41_128_16",   sg_stats_scan_sse41_128_16,   },
    {"sg_stats_scan_sse41_128_8",    sg_stats_scan_sse41_128_8,    },
    {"sg_stats_diag_sse41_128_32",   sg_stats_diag_sse41_128_32,   },
    {"sg_stats_diag_sse41_128_16",   sg_stats_diag_sse41_128_16,   },
    {"sg_stats_diag_sse41_128_8",    sg_stats_diag_sse41_128_8,    },
    {"sg_stats_striped_sse41_128_32",sg_stats_striped_sse41_128_32,},
    {"sg_stats_striped_sse41_128_16",sg_stats_striped_sse41_128_16,},
    {"sg_stats_striped_sse41_128_8", sg_stats_striped_sse41_128_8, },
#endif
#if HAVE_AVX2
    {"sg_stats_scan_avx2_256_32",    sg_stats_scan_avx2_256_32,    },
    {"sg_stats_scan_avx2_256_16",    sg_stats_scan_avx2_256_16,    },
    {"sg_stats_scan_avx2_256_8",     sg_stats_scan_avx2_256_8,     },
    {"sg_stats_diag_avx2_256_32",    sg_stats_diag_avx2_256_32,    },
    {"sg_stats_diag_avx2_256_16",    sg_stats_diag_avx2_256_16,    },
    {"sg_stats_diag_avx2_256_8",     sg_stats_diag_avx2_256_8,     },
    {"sg_stats_striped_avx2_256_32", sg_stats_striped_avx2_256_32, },
    {"sg_stats_striped_avx2_256_16", sg_stats_striped_avx2_256_16, },
    {"sg_stats_striped_avx2_256_8",  sg_stats_striped_avx2_256_8,  },
#endif
#if HAVE_KNC
    {"sg_stats_scan_knc_512_32",     sg_stats_scan_knc_512_32,     },
    {"sg_stats_diag_knc_512_32",     sg_stats_diag_knc_512_32,     },
    {"sg_stats_striped_knc_512_32",  sg_stats_striped_knc_512_32,  },
#endif

    {"sw_stats",                     sw_stats,                     },
    {"sw_stats_scan",                sw_stats_scan,                },
#if HAVE_SSE2
    {"sw_stats_scan_sse2_128_32",    sw_stats_scan_sse2_128_32,    },
    {"sw_stats_scan_sse2_128_16",    sw_stats_scan_sse2_128_16,    },
    {"sw_stats_scan_sse2_128_8",     sw_stats_scan_sse2_128_8,     },
    {"sw_stats_diag_sse2_128_32",    sw_stats_diag_sse2_128_32,    },
    {"sw_stats_diag_sse2_128_16",    sw_stats_diag_sse2_128_16,    },
    {"sw_stats_diag_sse2_128_8",     sw_stats_diag_sse2_128_8,     },
    {"sw_stats_striped_sse2_128_32", sw_stats_striped_sse2_128_32, },
    {"sw_stats_striped_sse2_128_16", sw_stats_striped_sse2_128_16, },
    {"sw_stats_striped_sse2_128_8",  sw_stats_striped_sse2_128_8,  },
#endif
#if HAVE_SSE41
    {"sw_stats_scan_sse41_128_32",   sw_stats_scan_sse41_128_32,   },
    {"sw_stats_scan_sse41_128_16",   sw_stats_scan_sse41_128_16,   },
    {"sw_stats_scan_sse41_128_8",    sw_stats_scan_sse41_128_8,    },
    {"sw_stats_diag_sse41_128_32",   sw_stats_diag_sse41_128_32,   },
    {"sw_stats_diag_sse41_128_16",   sw_stats_diag_sse41_128_16,   },
    {"sw_stats_diag_sse41_128_8",    sw_stats_diag_sse41_128_8,    },
    {"sw_stats_striped_sse41_128_32",sw_stats_striped_sse41_128_32,},
    {"sw_stats_striped_sse41_128_16",sw_stats_striped_sse41_128_16,},
    {"sw_stats_striped_sse41_128_8", sw_stats_striped_sse41_128_8, },
#endif
#if HAVE_AVX2
    {"sw_stats_scan_avx2_256_32",    sw_stats_scan_avx2_256_32,    },
    {"sw_stats_scan_avx2_256_16",    sw_stats_scan_avx2_256_16,    },
    {"sw_stats_scan_avx2_256_8",     sw_stats_scan_avx2_256_8,     },
    {"sw_stats_diag_avx2_256_32",    sw_stats_diag_avx2_256_32,    },
    {"sw_stats_diag_avx2_256_16",    sw_stats_diag_avx2_256_16,    },
    {"sw_stats_diag_avx2_256_8",     sw_stats_diag_avx2_256_8,     },
    {"sw_stats_striped_avx2_256_32", sw_stats_striped_avx2_256_32, },
    {"sw_stats_striped_avx2_256_16", sw_stats_striped_avx2_256_16, },
    {"sw_stats_striped_avx2_256_8",  sw_stats_striped_avx2_256_8,  },
#endif
#if HAVE_KNC
    {"sw_stats_scan_knc_512_32",     sw_stats_scan_knc_512_32,     },
    {"sw_stats_diag_knc_512_32",     sw_stats_diag_knc_512_32,     },
    {"sw_stats_striped_knc_512_32",  sw_stats_striped_knc_512_32,  },
#endif
    {"NULL", NULL}
};

static inline void parse_sequences(
        const char *filename, char ***strings_, size_t **sizes_, size_t *count_)
{
    FILE* fp;
    kseq_t *seq = NULL;
    int l = 0;
    char **strings = NULL;
    size_t *sizes = NULL;
    size_t count = 0;
    size_t memory = 1000;
    size_t i = 0;

    fp = fopen(filename, "r");
    if(fp == NULL) {
        perror("fopen");
        exit(1);
    }
    strings = malloc(sizeof(char*) * memory);
    sizes = malloc(sizeof(size_t) * memory);
    seq = kseq_init(fileno(fp));
    while ((l = kseq_read(seq)) >= 0) {
        strings[count] = strdup(seq->seq.s);
        if (NULL == strings[count]) {
            perror("strdup");
            exit(1);
        }
        sizes[count] = seq->seq.l;
        ++count;
        if (count >= memory) {
            char **new_strings = NULL;
            size_t *new_sizes = NULL;
            memory *= 2;
            new_strings = realloc(strings, sizeof(char*) * memory);
            if (NULL == new_strings) {
                perror("realloc");
                exit(1);
            }
            strings = new_strings;
            new_sizes = realloc(sizes, sizeof(size_t) * memory);
            if (NULL == new_sizes) {
                perror("realloc");
                exit(1);
            }
            sizes = new_sizes;
        }
    }
    kseq_destroy(seq);
    fclose(fp);

    *strings_ = strings;
    *sizes_ = sizes;
    *count_ = count;
}

static inline char* rand_string(size_t size)
{
    char *str = NULL;
    const char charset[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    if (size) {
        size_t n;
        --size;
        str = malloc(size + 1);
        for (n = 0; n < size; n++) {
            int key = rand() % (int) (sizeof charset - 1);
            str[n] = charset[key];
        }
        str[size] = '\0';
    }
    return str;
}

static inline unsigned long binomial_coefficient(unsigned long n, unsigned long k)
{
    /* from http://blog.plover.com/math/choose.html */
    unsigned long r = 1;
    unsigned long d;
    if (k > n) {
        return 0;
    }
    for (d = 1; d <= k; d++) {
        r *= n--;
        r /= d;
    }
    return r;
}

static inline void k_combination2(unsigned long pos, unsigned long *a, unsigned long *b)
{
    double s;
    double i = floor(sqrt(2.0 * pos)) - 1.0;
    if (i <= 1.0) {
        i = 1.0;
    }
    s = i * (i - 1.0) / 2.0;
    while (pos - s >= i) {
        s += i;
        i += 1;
    }
    *a = (unsigned long)(pos - s);
    *b = (unsigned long)(i);
}

int main(int argc, char **argv)
{
    int shortest = INT_MAX;
    int longest = 0;
    double timer_clock = 0.0;
    unsigned long i = 0;
    size_t seq_count = 10;
    size_t limit = 0;
    char **sequences = NULL;
    size_t *sizes = NULL;
    char *endptr = NULL;
    char *funcname = NULL;
    pf function = NULL;
    char *filename = NULL;
    int c = 0;
    char *blosumname = NULL;
    const int (*blosum)[24] = blosum62;
    int gap_open = 10;
    int gap_extend = 1;
    int N = 1;

    while ((c = getopt(argc, argv, "a:b:f:n:o:e:")) != -1) {
        switch (c) {
            case 'a':
                funcname = optarg;
                break;
            case 'b':
                blosumname = optarg;
                break;
            case 'f':
                filename = optarg;
                break;
            case 'n':
                errno = 0;
                seq_count = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case 'o':
                errno = 0;
                gap_open = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case 'e':
                errno = 0;
                gap_extend = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol");
                    exit(1);
                }
                break;
            case '?':
                if (optopt == 'f' || optopt == 'n') {
                    fprintf(stderr,
                            "Option -%c requires an argument.\n",
                            optopt);
                }
                else if (isprint(optopt)) {
                    fprintf(stderr, "Unknown option `-%c'.\n",
                            optopt);
                }
                else {
                    fprintf(stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                }
                exit(1);
            default:
                fprintf(stderr, "default case in getopt\n");
                exit(1);
        }
    }

    /* select the function */
    if (funcname) {
        int index = 0;
        func_t f;
        f = functions[index++];
        while (f.f) {
            if (0 == strcmp(funcname, f.name)) {
                function = f.f;
                printf("function: %s\n", funcname);
                break;
            }
            f = functions[index++];
        }
        if (NULL == function) {
            fprintf(stderr, "Specified function not found.\n");
            exit(1);
        }
    }
    else {
        fprintf(stderr, "No alignment function specified.\n");
        exit(1);
    }

    /* select the blosum matrix */
    if (blosumname) {
        int index = 0;
        blosum_t b;
        b = blosums[index++];
        while (b.blosum) {
            if (0 == strcmp(blosumname, b.name)) {
                blosum = b.blosum;
                printf("blosum: %s\n", blosumname);
                break;
            }
            b = blosums[index++];
        }
        if (NULL == b.blosum) {
            fprintf(stderr, "Specified blosum matrix not found.\n");
            fprintf(stderr, "Choices are {"
                    "blosum40,"
                    "blosum45,"
                    "blosum50,"
                    "blosum62,"
                    "blosum75,"
                    "blosum80,"
                    "blosum90}\n");
            exit(1);
        }
    }

    if (filename) {
        parse_sequences(filename, &sequences, &sizes, &seq_count);
    }
    else {
        /* generate 'seq_count' number of random strings */
        sequences = (char**)malloc(sizeof(char*)*seq_count);
        sizes = (size_t*)malloc(sizeof(size_t)*seq_count);
        for (i=0; i<seq_count; ++i) {
            sizes[i] = (rand()%32767)+10;
            shortest = sizes[i] < shortest ? sizes[i] : shortest;
            longest = sizes[i] > longest ? sizes[i] : longest;
            sequences[i] = rand_string(sizes[i]);
        }
        printf("done generating %lu random srings, "
                "shortest is %d, longest is %d\n",
                seq_count, shortest, longest);
    }

    limit = binomial_coefficient(seq_count, 2);
    printf("%lu choose 2 is %lu\n", seq_count, limit);

    timer_init();
    printf("%s timer\n", timer_name());

#pragma omp parallel
    {
#pragma omp single
        {
            N = omp_get_max_threads();
            printf("omp_get_max_threads()=%d\n", N);
        }
    }

    timer_clock = timer_real();
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        unsigned long a=0;
        unsigned long b=1;
#pragma omp for schedule(dynamic)
        for (i=0; i<limit; ++i) {
            parasail_result_t *result = NULL;
            k_combination2(i, &a, &b);
            result = function(sequences[a], sizes[a], sequences[b], sizes[b],
                    gap_open, gap_extend, blosum);
            parasail_result_free(result);
        }
    }
    timer_clock = timer_real() - timer_clock;
    printf("%s\t%s\t%d\t%f\n", funcname, blosumname, N, timer_clock);

    return 0;
}

