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

#ifdef __cplusplus
extern "C" {
#endif

int parasail_can_use_avx2();
int parasail_can_use_sse41();
int parasail_can_use_sse2();

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_CPUID_H_ */

