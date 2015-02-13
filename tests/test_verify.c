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

#include "kseq.h"
KSEQ_INIT(int, read)

#include "parasail.h"
#include "parasail_internal.h"
#include "parasail_cpuid.h"
#include "blosum/blosum40.h"
#include "blosum/blosum45.h"
#include "blosum/blosum50.h"
#include "blosum/blosum62.h"
#include "blosum/blosum75.h"
#include "blosum/blosum80.h"
#include "blosum/blosum90.h"

typedef parasail_result_t* (*pf)(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        const int matrix[24][24]);

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

typedef struct gap_score {
    int open;
    int extend;
} gap_score_t;

gap_score_t gap_scores[] = {
    {10,1},
    {10,2},
    {14,2},
    {40,2},
    {INT_MIN,INT_MIN}
};

typedef struct func {
    const char * name;
    pf f;
} func_t;

typedef struct funcs {
    const char * name;
    func_t *fs;
} funcs_t;

#if HAVE_SSE2
func_t nw_sse2_functions[] = {
    {"nw",                           nw,                        },
    {"nw_scan",                      nw_scan,                   },
    {"nw_scan_sse2_128_32",          nw_scan_sse2_128_32,       },
    {"nw_scan_sse2_128_16",          nw_scan_sse2_128_16,       },
    {"nw_scan_sse2_128_8",           nw_scan_sse2_128_8,        },
    {"nw_diag_sse2_128_32",          nw_diag_sse2_128_32,       },
    {"nw_diag_sse2_128_16",          nw_diag_sse2_128_16,       },
    {"nw_diag_sse2_128_8",           nw_diag_sse2_128_8,        },
    {"nw_striped_sse2_128_32",       nw_striped_sse2_128_32,    },
    {"nw_striped_sse2_128_16",       nw_striped_sse2_128_16,    },
    {"nw_striped_sse2_128_8",        nw_striped_sse2_128_8,     },
    {"NULL", NULL}
};
funcs_t nw_sse2 = {"nw_sse2", nw_sse2_functions};
#endif

#if HAVE_SSE41
func_t nw_sse41_functions[] = {
    {"nw",                           nw,                        },
    {"nw_scan",                      nw_scan,                   },
    {"nw_scan_sse41_128_32",         nw_scan_sse41_128_32,      },
    {"nw_scan_sse41_128_16",         nw_scan_sse41_128_16,      },
    {"nw_scan_sse41_128_8",          nw_scan_sse41_128_8,       },
    {"nw_diag_sse41_128_32",         nw_diag_sse41_128_32,      },
    {"nw_diag_sse41_128_16",         nw_diag_sse41_128_16,      },
    {"nw_diag_sse41_128_8",          nw_diag_sse41_128_8,       },
    {"nw_striped_sse41_128_32",      nw_striped_sse41_128_32,   },
    {"nw_striped_sse41_128_16",      nw_striped_sse41_128_16,   },
    {"nw_striped_sse41_128_8",       nw_striped_sse41_128_8,    },
    {"NULL", NULL}
};
funcs_t nw_sse41 = {"nw_sse41", nw_sse41_functions};
#endif

#if HAVE_AVX2
func_t nw_avx2_functions[] = {
    {"nw",                           nw,                        },
    {"nw_scan",                      nw_scan,                   },
    {"nw_scan_avx2_256_32",          nw_scan_avx2_256_32,       },
    {"nw_scan_avx2_256_16",          nw_scan_avx2_256_16,       },
    {"nw_scan_avx2_256_8",           nw_scan_avx2_256_8,        },
    {"nw_diag_avx2_256_32",          nw_diag_avx2_256_32,       },
    {"nw_diag_avx2_256_16",          nw_diag_avx2_256_16,       },
    {"nw_diag_avx2_256_8",           nw_diag_avx2_256_8,        },
    {"nw_striped_avx2_256_32",       nw_striped_avx2_256_32,    },
    {"nw_striped_avx2_256_16",       nw_striped_avx2_256_16,    },
    {"nw_striped_avx2_256_8",        nw_striped_avx2_256_8,     },
    {"NULL", NULL}
};
funcs_t nw_avx2 = {"nw_avx2", nw_avx2_functions};
#endif

#if HAVE_KNC
func_t nw_knc_functions[] = {
    {"nw",                           nw,                        },
    {"nw_scan",                      nw_scan,                   },
    {"nw_scan_knc_512_32",           nw_scan_knc_512_32,        },
    {"nw_diag_knc_512_32",           nw_diag_knc_512_32,        },
    {"nw_striped_knc_512_32",        nw_striped_knc_512_32,     },
    {"NULL", NULL}
};
funcs_t nw_knc = {"nw_knc", nw_knc_functions};
#endif

#if HAVE_SSE2
func_t sg_sse2_functions[] = {
    {"sg",                           sg,                        },
    {"sg_scan",                      sg_scan,                   },
    {"sg_scan_sse2_128_32",          sg_scan_sse2_128_32,       },
    {"sg_scan_sse2_128_16",          sg_scan_sse2_128_16,       },
    {"sg_scan_sse2_128_8",           sg_scan_sse2_128_8,        },
    {"sg_diag_sse2_128_32",          sg_diag_sse2_128_32,       },
    {"sg_diag_sse2_128_16",          sg_diag_sse2_128_16,       },
    {"sg_diag_sse2_128_8",           sg_diag_sse2_128_8,        },
    {"sg_striped_sse2_128_32",       sg_striped_sse2_128_32,    },
    {"sg_striped_sse2_128_16",       sg_striped_sse2_128_16,    },
    {"sg_striped_sse2_128_8",        sg_striped_sse2_128_8,     },
    {"NULL", NULL}
};
funcs_t sg_sse2 = {"sg_sse2", sg_sse2_functions};
#endif

#if HAVE_SSE41
func_t sg_sse41_functions[] = {
    {"sg",                           sg,                        },
    {"sg_scan",                      sg_scan,                   },
    {"sg_scan_sse41_128_32",         sg_scan_sse41_128_32,      },
    {"sg_scan_sse41_128_16",         sg_scan_sse41_128_16,      },
    {"sg_scan_sse41_128_8",          sg_scan_sse41_128_8,       },
    {"sg_diag_sse41_128_32",         sg_diag_sse41_128_32,      },
    {"sg_diag_sse41_128_16",         sg_diag_sse41_128_16,      },
    {"sg_diag_sse41_128_8",          sg_diag_sse41_128_8,       },
    {"sg_striped_sse41_128_32",      sg_striped_sse41_128_32,   },
    {"sg_striped_sse41_128_16",      sg_striped_sse41_128_16,   },
    {"sg_striped_sse41_128_8",       sg_striped_sse41_128_8,    },
    {"NULL", NULL}
};
funcs_t sg_sse41 = {"sg_sse41", sg_sse41_functions};
#endif

#if HAVE_AVX2
func_t sg_avx2_functions[] = {
    {"sg",                           sg,                        },
    {"sg_scan",                      sg_scan,                   },
    {"sg_scan_avx2_256_32",          sg_scan_avx2_256_32,       },
    {"sg_scan_avx2_256_16",          sg_scan_avx2_256_16,       },
    {"sg_scan_avx2_256_8",           sg_scan_avx2_256_8,        },
    {"sg_diag_avx2_256_32",          sg_diag_avx2_256_32,       },
    {"sg_diag_avx2_256_16",          sg_diag_avx2_256_16,       },
    {"sg_diag_avx2_256_8",           sg_diag_avx2_256_8,        },
    {"sg_striped_avx2_256_32",       sg_striped_avx2_256_32,    },
    {"sg_striped_avx2_256_16",       sg_striped_avx2_256_16,    },
    {"sg_striped_avx2_256_8",        sg_striped_avx2_256_8,     },
    {"NULL", NULL}
};
funcs_t sg_avx2 = {"sg_avx2", sg_avx2_functions};
#endif

#if HAVE_KNC
func_t sg_knc_functions[] = {
    {"sg",                           sg,                        },
    {"sg_scan",                      sg_scan,                   },
    {"sg_scan_knc_512_32",           sg_scan_knc_512_32,        },
    {"sg_diag_knc_512_32",           sg_diag_knc_512_32,        },
    {"sg_striped_knc_512_32",        sg_striped_knc_512_32,     },
    {"NULL", NULL}
};
funcs_t sg_avx2 = {"sg_avx2", sg_avx2_functions};
#endif

#if HAVE_SSE2
func_t sw_sse2_functions[] = {
    {"sw",                           sw,                        },
    {"sw_scan",                      sw_scan,                   },
    {"sw_scan_sse2_128_32",          sw_scan_sse2_128_32,       },
    {"sw_scan_sse2_128_16",          sw_scan_sse2_128_16,       },
    {"sw_scan_sse2_128_8",           sw_scan_sse2_128_8,        },
    {"sw_diag_sse2_128_32",          sw_diag_sse2_128_32,       },
    {"sw_diag_sse2_128_16",          sw_diag_sse2_128_16,       },
    {"sw_diag_sse2_128_8",           sw_diag_sse2_128_8,        },
    {"sw_striped_sse2_128_32",       sw_striped_sse2_128_32,    },
    {"sw_striped_sse2_128_16",       sw_striped_sse2_128_16,    },
    {"sw_striped_sse2_128_8",        sw_striped_sse2_128_8,     },
    {"NULL", NULL}
};
funcs_t sw_sse2 = {"sw_sse2", sw_sse2_functions};
#endif

#if HAVE_SSE41
func_t sw_sse41_functions[] = {
    {"sw",                           sw,                        },
    {"sw_scan",                      sw_scan,                   },
    {"sw_scan_sse41_128_32",         sw_scan_sse41_128_32,      },
    {"sw_scan_sse41_128_16",         sw_scan_sse41_128_16,      },
    {"sw_scan_sse41_128_8",          sw_scan_sse41_128_8,       },
    {"sw_diag_sse41_128_32",         sw_diag_sse41_128_32,      },
    {"sw_diag_sse41_128_16",         sw_diag_sse41_128_16,      },
    {"sw_diag_sse41_128_8",          sw_diag_sse41_128_8,       },
    {"sw_striped_sse41_128_32",      sw_striped_sse41_128_32,   },
    {"sw_striped_sse41_128_16",      sw_striped_sse41_128_16,   },
    {"sw_striped_sse41_128_8",       sw_striped_sse41_128_8,    },
    {"NULL", NULL}
};
funcs_t sw_sse41 = {"sw_sse41", sw_sse41_functions};
#endif

#if HAVE_AVX2
func_t sw_avx2_functions[] = {
    {"sw",                           sw,                        },
    {"sw_scan",                      sw_scan,                   },
    {"sw_scan_avx2_256_32",          sw_scan_avx2_256_32,       },
    {"sw_scan_avx2_256_16",          sw_scan_avx2_256_16,       },
    {"sw_scan_avx2_256_8",           sw_scan_avx2_256_8,        },
    {"sw_diag_avx2_256_32",          sw_diag_avx2_256_32,       },
    {"sw_diag_avx2_256_16",          sw_diag_avx2_256_16,       },
    {"sw_diag_avx2_256_8",           sw_diag_avx2_256_8,        },
    {"sw_striped_avx2_256_32",       sw_striped_avx2_256_32,    },
    {"sw_striped_avx2_256_16",       sw_striped_avx2_256_16,    },
    {"sw_striped_avx2_256_8",        sw_striped_avx2_256_8,     },
    {"NULL", NULL}
};
funcs_t sw_avx2 = {"sw_avx2", sw_avx2_functions};
#endif

#if HAVE_KNC
func_t sw_knc_functions[] = {
    {"sw",                           sw,                        },
    {"sw_scan",                      sw_scan,                   },
    {"sw_scan_knc_512_32",           sw_scan_knc_512_32,        },
    {"sw_diag_knc_512_32",           sw_diag_knc_512_32,        },
    {"sw_striped_knc_512_32",        sw_striped_knc_512_32,     },
    {"NULL", NULL}
};
funcs_t sw_knc = {"sw_knc", sw_knc_functions};
#endif

#if HAVE_SSE2
func_t nw_stats_sse2_functions[] = {
    {"nw_stats",                     nw_stats,                     },
    {"nw_stats_scan",                nw_stats_scan,                },
    {"nw_stats_scan_sse2_128_32",    nw_stats_scan_sse2_128_32,    },
    {"nw_stats_scan_sse2_128_16",    nw_stats_scan_sse2_128_16,    },
    {"nw_stats_scan_sse2_128_8",     nw_stats_scan_sse2_128_8,     },
    {"nw_stats_diag_sse2_128_32",    nw_stats_diag_sse2_128_32,    },
    {"nw_stats_diag_sse2_128_16",    nw_stats_diag_sse2_128_16,    },
    {"nw_stats_diag_sse2_128_8",     nw_stats_diag_sse2_128_8,     },
    {"nw_stats_striped_sse2_128_32", nw_stats_striped_sse2_128_32, },
    {"nw_stats_striped_sse2_128_16", nw_stats_striped_sse2_128_16, },
    {"nw_stats_striped_sse2_128_8",  nw_stats_striped_sse2_128_8,  },
    {"NULL", NULL}
};
funcs_t nw_stats_sse2 = {"nw_stats_sse2", nw_stats_sse2_functions};
#endif

#if HAVE_SSE41
func_t nw_stats_sse41_functions[] = {
    {"nw_stats",                     nw_stats,                     },
    {"nw_stats_scan",                nw_stats_scan,                },
    {"nw_stats_scan_sse41_128_32",   nw_stats_scan_sse41_128_32,   },
    {"nw_stats_scan_sse41_128_16",   nw_stats_scan_sse41_128_16,   },
    {"nw_stats_scan_sse41_128_8",    nw_stats_scan_sse41_128_8,    },
    {"nw_stats_diag_sse41_128_32",   nw_stats_diag_sse41_128_32,   },
    {"nw_stats_diag_sse41_128_16",   nw_stats_diag_sse41_128_16,   },
    {"nw_stats_diag_sse41_128_8",    nw_stats_diag_sse41_128_8,    },
    {"nw_stats_striped_sse41_128_32",nw_stats_striped_sse41_128_32,},
    {"nw_stats_striped_sse41_128_16",nw_stats_striped_sse41_128_16,},
    {"nw_stats_striped_sse41_128_8", nw_stats_striped_sse41_128_8, },
    {"NULL", NULL}
};
funcs_t nw_stats_sse41 = {"nw_stats_sse41", nw_stats_sse41_functions};
#endif

#if HAVE_AVX2
func_t nw_stats_avx2_functions[] = {
    {"nw_stats",                     nw_stats,                     },
    {"nw_stats_scan",                nw_stats_scan,                },
    {"nw_stats_scan_avx2_256_32",    nw_stats_scan_avx2_256_32,    },
    {"nw_stats_scan_avx2_256_16",    nw_stats_scan_avx2_256_16,    },
    {"nw_stats_scan_avx2_256_8",     nw_stats_scan_avx2_256_8,     },
    {"nw_stats_diag_avx2_256_32",    nw_stats_diag_avx2_256_32,    },
    {"nw_stats_diag_avx2_256_16",    nw_stats_diag_avx2_256_16,    },
    {"nw_stats_diag_avx2_256_8",     nw_stats_diag_avx2_256_8,     },
    {"nw_stats_striped_avx2_256_32", nw_stats_striped_avx2_256_32, },
    {"nw_stats_striped_avx2_256_16", nw_stats_striped_avx2_256_16, },
    {"nw_stats_striped_avx2_256_8",  nw_stats_striped_avx2_256_8,  },
    {"NULL", NULL}
};
funcs_t nw_stats_avx2 = {"nw_stats_avx2", nw_stats_avx2_functions};
#endif

#if HAVE_KNC
func_t nw_stats_knc_functions[] = {
    {"nw_stats",                     nw_stats,                     },
    {"nw_stats_scan",                nw_stats_scan,                },
    {"nw_stats_scan_knc_512_32",     nw_stats_scan_knc_512_32,     },
    {"nw_stats_diag_knc_512_32",     nw_stats_diag_knc_512_32,     },
    {"nw_stats_striped_knc_512_32",  nw_stats_striped_knc_512_32,  },
    {"NULL", NULL}
};
funcs_t nw_stats_knc = {"nw_stats_knc", nw_stats_knc_functions};
#endif

#if HAVE_SSE2
func_t sg_stats_sse2_functions[] = {
    {"sg_stats",                     sg_stats,                     },
    {"sg_stats_scan",                sg_stats_scan,                },
    {"sg_stats_scan_sse2_128_32",    sg_stats_scan_sse2_128_32,    },
    {"sg_stats_scan_sse2_128_16",    sg_stats_scan_sse2_128_16,    },
    {"sg_stats_scan_sse2_128_8",     sg_stats_scan_sse2_128_8,     },
    {"sg_stats_diag_sse2_128_32",    sg_stats_diag_sse2_128_32,    },
    {"sg_stats_diag_sse2_128_16",    sg_stats_diag_sse2_128_16,    },
    {"sg_stats_diag_sse2_128_8",     sg_stats_diag_sse2_128_8,     },
    {"sg_stats_striped_sse2_128_32", sg_stats_striped_sse2_128_32, },
    {"sg_stats_striped_sse2_128_16", sg_stats_striped_sse2_128_16, },
    {"sg_stats_striped_sse2_128_8",  sg_stats_striped_sse2_128_8,  },
    {"NULL", NULL}
};
funcs_t sg_stats_sse2 = {"sg_stats_sse2", sg_stats_sse2_functions};
#endif

#if HAVE_SSE41
func_t sg_stats_sse41_functions[] = {
    {"sg_stats",                     sg_stats,                     },
    {"sg_stats_scan",                sg_stats_scan,                },
    {"sg_stats_scan_sse41_128_32",   sg_stats_scan_sse41_128_32,   },
    {"sg_stats_scan_sse41_128_16",   sg_stats_scan_sse41_128_16,   },
    {"sg_stats_scan_sse41_128_8",    sg_stats_scan_sse41_128_8,    },
    {"sg_stats_diag_sse41_128_32",   sg_stats_diag_sse41_128_32,   },
    {"sg_stats_diag_sse41_128_16",   sg_stats_diag_sse41_128_16,   },
    {"sg_stats_diag_sse41_128_8",    sg_stats_diag_sse41_128_8,    },
    {"sg_stats_striped_sse41_128_32",sg_stats_striped_sse41_128_32,},
    {"sg_stats_striped_sse41_128_16",sg_stats_striped_sse41_128_16,},
    {"sg_stats_striped_sse41_128_8", sg_stats_striped_sse41_128_8, },
    {"NULL", NULL}
};
funcs_t sg_stats_sse41 = {"sg_stats_sse41", sg_stats_sse41_functions};
#endif

#if HAVE_AVX2
func_t sg_stats_avx2_functions[] = {
    {"sg_stats",                     sg_stats,                     },
    {"sg_stats_scan",                sg_stats_scan,                },
    {"sg_stats_scan_avx2_256_32",    sg_stats_scan_avx2_256_32,    },
    {"sg_stats_scan_avx2_256_16",    sg_stats_scan_avx2_256_16,    },
    {"sg_stats_scan_avx2_256_8",     sg_stats_scan_avx2_256_8,     },
    {"sg_stats_diag_avx2_256_32",    sg_stats_diag_avx2_256_32,    },
    {"sg_stats_diag_avx2_256_16",    sg_stats_diag_avx2_256_16,    },
    {"sg_stats_diag_avx2_256_8",     sg_stats_diag_avx2_256_8,     },
    {"sg_stats_striped_avx2_256_32", sg_stats_striped_avx2_256_32, },
    {"sg_stats_striped_avx2_256_16", sg_stats_striped_avx2_256_16, },
    {"sg_stats_striped_avx2_256_8",  sg_stats_striped_avx2_256_8,  },
    {"NULL", NULL}
};
funcs_t sg_stats_avx2 = {"sg_stats_avx2", sg_stats_avx2_functions};
#endif

#if HAVE_KNC
func_t sg_stats_knc_functions[] = {
    {"sg_stats",                     sg_stats,                     },
    {"sg_stats_scan",                sg_stats_scan,                },
    {"sg_stats_scan_knc_512_32",     sg_stats_scan_knc_512_32,     },
    {"sg_stats_diag_knc_512_32",     sg_stats_diag_knc_512_32,     },
    {"sg_stats_striped_knc_512_32",  sg_stats_striped_knc_512_32,  },
    {"NULL", NULL}
};
funcs_t sg_stats_knc = {"sg_stats_knc", sg_stats_knc_functions};
#endif

#if HAVE_SSE2
func_t sw_stats_sse2_functions[] = {
    {"sw_stats",                     sw_stats,                     },
    {"sw_stats_scan",                sw_stats_scan,                },
    {"sw_stats_scan_sse2_128_32",    sw_stats_scan_sse2_128_32,    },
    {"sw_stats_scan_sse2_128_16",    sw_stats_scan_sse2_128_16,    },
    {"sw_stats_scan_sse2_128_8",     sw_stats_scan_sse2_128_8,     },
    {"sw_stats_diag_sse2_128_32",    sw_stats_diag_sse2_128_32,    },
    {"sw_stats_diag_sse2_128_16",    sw_stats_diag_sse2_128_16,    },
    {"sw_stats_diag_sse2_128_8",     sw_stats_diag_sse2_128_8,     },
    {"sw_stats_striped_sse2_128_32", sw_stats_striped_sse2_128_32, },
    {"sw_stats_striped_sse2_128_16", sw_stats_striped_sse2_128_16, },
    {"sw_stats_striped_sse2_128_8",  sw_stats_striped_sse2_128_8,  },
    {"NULL", NULL}
};
funcs_t sw_stats_sse2 = {"sw_stats_sse2", sw_stats_sse2_functions};
#endif

#if HAVE_SSE41
func_t sw_stats_sse41_functions[] = {
    {"sw_stats",                     sw_stats,                     },
    {"sw_stats_scan",                sw_stats_scan,                },
    {"sw_stats_scan_sse41_128_32",   sw_stats_scan_sse41_128_32,   },
    {"sw_stats_scan_sse41_128_16",   sw_stats_scan_sse41_128_16,   },
    {"sw_stats_scan_sse41_128_8",    sw_stats_scan_sse41_128_8,    },
    {"sw_stats_diag_sse41_128_32",   sw_stats_diag_sse41_128_32,   },
    {"sw_stats_diag_sse41_128_16",   sw_stats_diag_sse41_128_16,   },
    {"sw_stats_diag_sse41_128_8",    sw_stats_diag_sse41_128_8,    },
    {"sw_stats_striped_sse41_128_32",sw_stats_striped_sse41_128_32,},
    {"sw_stats_striped_sse41_128_16",sw_stats_striped_sse41_128_16,},
    {"sw_stats_striped_sse41_128_8", sw_stats_striped_sse41_128_8, },
    {"NULL", NULL}
};
funcs_t sw_stats_sse41 = {"sw_stats_sse41", sw_stats_sse41_functions};
#endif

#if HAVE_AVX2
func_t sw_stats_avx2_functions[] = {
    {"sw_stats",                     sw_stats,                     },
    {"sw_stats_scan",                sw_stats_scan,                },
    {"sw_stats_scan_avx2_256_32",    sw_stats_scan_avx2_256_32,    },
    {"sw_stats_scan_avx2_256_16",    sw_stats_scan_avx2_256_16,    },
    {"sw_stats_scan_avx2_256_8",     sw_stats_scan_avx2_256_8,     },
    {"sw_stats_diag_avx2_256_32",    sw_stats_diag_avx2_256_32,    },
    {"sw_stats_diag_avx2_256_16",    sw_stats_diag_avx2_256_16,    },
    {"sw_stats_diag_avx2_256_8",     sw_stats_diag_avx2_256_8,     },
    {"sw_stats_striped_avx2_256_32", sw_stats_striped_avx2_256_32, },
    {"sw_stats_striped_avx2_256_16", sw_stats_striped_avx2_256_16, },
    {"sw_stats_striped_avx2_256_8",  sw_stats_striped_avx2_256_8,  },
    {"NULL", NULL}
};
funcs_t sw_stats_avx2 = {"sw_stats_avx2", sw_stats_avx2_functions};
#endif

#if HAVE_KNC
func_t sw_stats_knc_functions[] = {
    {"sw_stats",                     sw_stats,                     },
    {"sw_stats_scan",                sw_stats_scan,                },
    {"sw_stats_scan_knc_512_32",     sw_stats_scan_knc_512_32,     },
    {"sw_stats_diag_knc_512_32",     sw_stats_diag_knc_512_32,     },
    {"sw_stats_striped_knc_512_32",  sw_stats_striped_knc_512_32,  },
    {"NULL", NULL}
};
funcs_t sw_stats_knc = {"sw_stats_knc", sw_stats_knc_functions};
#endif

static inline void parse_sequences(
        const char *filename,
        char ***strings_,
        unsigned long **sizes_,
        unsigned long *count_)
{
    FILE* fp;
    kseq_t *seq = NULL;
    int l = 0;
    char **strings = NULL;
    unsigned long *sizes = NULL;
    unsigned long count = 0;
    unsigned long memory = 1000;
    unsigned long i = 0;

    fp = fopen(filename, "r");
    if(fp == NULL) {
        perror("fopen");
        exit(1);
    }
    strings = malloc(sizeof(char*) * memory);
    sizes = malloc(sizeof(unsigned long) * memory);
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
            unsigned long *new_sizes = NULL;
            memory *= 2;
            new_strings = realloc(strings, sizeof(char*) * memory);
            if (NULL == new_strings) {
                perror("realloc");
                exit(1);
            }
            strings = new_strings;
            new_sizes = realloc(sizes, sizeof(unsigned long) * memory);
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

static inline unsigned long binomial_coefficient(
        unsigned long n,
        unsigned long k)
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

static inline void k_combination2(
        unsigned long pos,
        unsigned long *a,
        unsigned long *b)
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

static inline void check_functions(
        funcs_t f,
        char **sequences,
        unsigned long *sizes,
        unsigned long count,
        unsigned long pair_limit)
{
    func_t *functions = f.fs;
    unsigned long blosum_index = 0;
    unsigned long gap_index = 0;
    unsigned long function_index = 0;
    unsigned long pair_index = 0;
    pf reference_function = NULL;

    printf("checking %s functions\n", f.name);
    for (blosum_index=0; NULL!=blosums[blosum_index].blosum; ++blosum_index) {
        //printf("\t%s\n", blosums[blosum_index].name);
        for (gap_index=0; INT_MIN!=gap_scores[gap_index].open; ++gap_index) {
            //printf("\t\topen=%d extend=%d\n",
            //        gap_scores[gap_index].open,
            //        gap_scores[gap_index].extend);
            reference_function = functions[0].f;
            for (function_index=1;
                    NULL!=functions[function_index].f;
                    ++function_index) {
                //printf("\t\t\t%s\n", functions[function_index].name);
#pragma omp parallel for
                for (pair_index=0; pair_index<pair_limit; ++pair_index) {
                    parasail_result_t *reference_result = NULL;
                    parasail_result_t *result = NULL;
                    unsigned long a = 0;
                    unsigned long b = 1;
                    int open = gap_scores[gap_index].open;
                    int extend = gap_scores[gap_index].extend;
                    k_combination2(pair_index, &a, &b);
                    //printf("\t\t\t\tpair=%lu (%lu,%lu)\n", pair_index, a, b);
                    reference_result = reference_function(
                            sequences[a], sizes[a],
                            sequences[b], sizes[b],
                            open, extend,
                            blosums[blosum_index].blosum);
                    result = functions[function_index].f(
                            sequences[a], sizes[a],
                            sequences[b], sizes[b],
                            open, extend,
                            blosums[blosum_index].blosum);
                    if (result->saturated) {
                        /* no point in comparing a result that saturated */
                        parasail_result_free(reference_result);
                        parasail_result_free(result);
                        continue;
                    }
                    if (reference_result->score != result->score) {
#pragma omp critical(printer)
                        {
                            printf("%s(%lu,%lu,%d,%d,%s) wrong score\n",
                                    functions[function_index].name,
                                    a, b, open, extend,
                                    blosums[blosum_index].name);
                        }
                    }
                    parasail_result_free(reference_result);
                    parasail_result_free(result);
                }
            }
        }
    }
}

int main(int argc, char **argv)
{
    unsigned long i = 0;
    unsigned long seq_count = 0;
    unsigned long limit = 0;
    char **sequences = NULL;
    unsigned long *sizes = NULL;
    char *endptr = NULL;
    char *funcname = NULL;
    pf function = NULL;
    char *filename = NULL;
    int c = 0;


    while ((c = getopt(argc, argv, "f:n:")) != -1) {
        switch (c) {
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

    if (filename) {
        parse_sequences(filename, &sequences, &sizes, &seq_count);
    }
    else {
        fprintf(stderr, "no filename specified\n");
        exit(1);
    }

    limit = binomial_coefficient(seq_count, 2);
    printf("%lu choose 2 is %lu\n", seq_count, limit);


#if HAVE_SSE2
    if (parasail_can_use_sse2()) {
        check_functions(nw_sse2, sequences, sizes, seq_count, limit);
        check_functions(sg_sse2, sequences, sizes, seq_count, limit);
        check_functions(sw_sse2, sequences, sizes, seq_count, limit);
        check_functions(nw_stats_sse2, sequences, sizes, seq_count, limit);
        check_functions(sg_stats_sse2, sequences, sizes, seq_count, limit);
        check_functions(sw_stats_sse2, sequences, sizes, seq_count, limit);
    }
#endif

#if HAVE_SSE41
    if (parasail_can_use_sse41()) {
        check_functions(nw_sse41, sequences, sizes, seq_count, limit);
        check_functions(sg_sse41, sequences, sizes, seq_count, limit);
        check_functions(sw_sse41, sequences, sizes, seq_count, limit);
        check_functions(nw_stats_sse41, sequences, sizes, seq_count, limit);
        check_functions(sg_stats_sse41, sequences, sizes, seq_count, limit);
        check_functions(sw_stats_sse41, sequences, sizes, seq_count, limit);
    }
#endif

#if HAVE_AVX2
    if (parasail_can_use_avx2()) {
        check_functions(nw_avx2, sequences, sizes, seq_count, limit);
        check_functions(sg_avx2, sequences, sizes, seq_count, limit);
        check_functions(sw_avx2, sequences, sizes, seq_count, limit);
        check_functions(nw_stats_avx2, sequences, sizes, seq_count, limit);
        check_functions(sg_stats_avx2, sequences, sizes, seq_count, limit);
        check_functions(sw_stats_avx2, sequences, sizes, seq_count, limit);
    }
#endif

    for (i=0; i<seq_count; ++i) {
        free(sequences[i]);
    }
    free(sequences);
    free(sizes);

    return 0;
}

