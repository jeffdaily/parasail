/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#ifndef _PARASAIL_CPUID_H_
#define _PARASAIL_CPUID_H_

#include "parasail.h"

#ifdef __cplusplus
extern "C" {
#endif

extern PARASAIL_API int parasail_can_use_avx2();
extern PARASAIL_API int parasail_can_use_sse41();
extern PARASAIL_API int parasail_can_use_sse2();

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_CPUID_H_ */

