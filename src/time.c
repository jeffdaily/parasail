/**
 * @file timer_real.h
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Portable use of nanosecond precision realtime clock.
 */
#include "config.h"

#include "parasail.h"

#if HAVE_GETSYSTEMTIMEASFILETIME && !defined(INT64_LITERAL_SUFFIX_UNKNOWN)
#include <windows.h>
#elif HAVE_CLOCK_GET_TIME
#include <mach/clock.h>
#include <mach/mach.h>
#elif HAVE_CLOCK_GETTIME_MONOTONIC || HAVE_CLOCK_GETTIME_REALTIME
/* undef if needed to avoid compiler warning */
#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif
#define _POSIX_C_SOURCE 199309L
#include <time.h>
#endif

double parasail_time(void)
{
#if HAVE_GETSYSTEMTIMEASFILETIME && !defined(INT64_LITERAL_SUFFIX_UNKNOWN)
    __int64 wintime;
    double sec;
    double nsec;
    GetSystemTimeAsFileTime((FILETIME*)&wintime);
#if defined(INT64_LITERAL_SUFFIX_I64)
    wintime -=116444736000000000i64;   /*1jan1601 to 1jan1970*/
    sec  = wintime / 10000000i64;      /*seconds*/
    nsec = wintime % 10000000i64 *100; /*nano-seconds*/
#elif defined(INT64_LITERAL_SUFFIX_LL)
    wintime -=116444736000000000LL;   /*1jan1601 to 1jan1970*/
    sec  = wintime / 10000000LL;      /*seconds*/
    nsec = wintime % 10000000LL *100; /*nano-seconds*/
#else
#error unknown int64 literal suffix
#endif
    return sec + nsec/1000000000.0;
#elif HAVE_CLOCK_GET_TIME
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    return (double)(mts.tv_sec) + (double)(mts.tv_nsec)/1000000000.0;
#elif HAVE_CLOCK_GETTIME_MONOTONIC
    struct timespec ts;
    /* Works on FreeBSD */
    long retval = clock_gettime(CLOCK_MONOTONIC, &ts);
    if (0 != retval) {
        fprintf(stderr, "%s: clock_gettime failed\n", __func__);
        return 0.0;
    }
    return (double)(ts.tv_sec) + (double)(ts.tv_nsec)/1000000000.0;
#elif HAVE_CLOCK_GETTIME_REALTIME
    struct timespec ts;
    /* Works on Linux */
    long retval = clock_gettime(CLOCK_REALTIME, &ts);
    if (0 != retval) {
        fprintf(stderr, "%s: clock_gettime failed\n", __func__);
        return 0.0;
    }
    return (double)(ts.tv_sec) + (double)(ts.tv_nsec)/1000000000.0;
#endif
    return 0.0; /* avoid warnings */
}

