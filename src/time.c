/**
 * @file timer_real.h
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Portable use of nanosecond precision realtime clock.
 */
#include "config.h"

#include "parasail.h"

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#elif defined(__FreeBSD__)
#include <time.h>
#else
#include <time.h>
#endif

#include <assert.h>

double parasail_time(void)
{
#ifdef __MACH__
    /* OS X does not have clock_gettime, use clock_get_time */
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    return (double)(mts.tv_sec) + (double)(mts.tv_nsec)/1000000000.0;
#elif defined(__FreeBSD__)
    struct timespec ts;
    /* Works on FreeBSD */
    long retval = clock_gettime(CLOCK_MONOTONIC, &ts);
    assert(0 == retval);
    return (double)(ts.tv_sec) + (double)(ts.tv_nsec)/1000000000.0;
#else
    struct timespec ts;
    /* Works on Linux */
    long retval = clock_gettime(CLOCK_REALTIME, &ts);
    assert(0 == retval);
    return (double)(ts.tv_sec) + (double)(ts.tv_nsec)/1000000000.0;
#endif
}

