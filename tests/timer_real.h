/**
 * @file timer_real.h
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Portable use of nanosecond precision realtime clock.
 */

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#elif defined(__FreeBSD__)
#include <time.h>
#else
#include <time.h>
#endif

inline static double timer_real()
{
    timespec ts;
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts.tv_sec = mts.tv_sec;
    ts.tv_nsec = mts.tv_nsec;
#elif defined(__FreeBSD__)
    long retval = clock_gettime(CLOCK_MONOTONIC, &ts); // Works on FreeBSD
    assert(0 == retval);
#else
    long retval = clock_gettime(CLOCK_REALTIME, &ts); // Works on Linux
    assert(0 == retval);
#endif
    return double(ts.tv_sec) + double(ts.tv_nsec)/1000000000.0;
}
