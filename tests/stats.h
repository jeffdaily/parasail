/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2014 Battelle Memorial Institute.
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#ifndef _STATS_H_
#define _STATS_H_

#include <math.h>

/** http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance */

typedef struct stats {
    unsigned long _n;
    double _mean;
    double _M2;
    double _sum;
    double _min;
    double _max;
} stats_t;

static inline void stats_clear(stats_t *stats) {
    stats->_n = 0UL;
    stats->_mean = 0.0;
    stats->_M2 = 0.0;
    stats->_sum = 0.0;
    stats->_min = 0.0;
    stats->_max = 0.0;
}

static inline stats_t* stats_new() {
    stats_t *stats = (stats_t*)malloc(sizeof(stats_t));
    stats_clear(stats);
    return stats;
}

static inline void stats_sample_value(stats_t *stats, const double x) {
    double delta = 0;

    /* extra stats */
    stats->_sum = stats->_sum + x;
    if (0UL == stats->_n) {
        stats->_min = x;
        stats->_max = x;
    }
    else {
        stats->_min = stats->_min < x ? stats->_min : x;
        stats->_max = stats->_max > x ? stats->_max : x;
    }

    stats->_n = stats->_n + 1UL;
    delta = x - stats->_mean;
    stats->_mean = stats->_mean + delta/stats->_n;
    stats->_M2 = stats->_M2 + delta * (x - stats->_mean);
}

static inline void stats_sample_stats(stats_t *stats, const stats_t * const B) {
    if (B->_n == 0) {
        return;
    }
    else if (stats->_n == 0) {
        stats->_n = B->_n;
        stats->_mean = B->_mean;
        stats->_M2 = B->_M2;
        stats->_sum = B->_sum;
        stats->_min = B->_min;
        stats->_max = B->_max;
    }
    else {
        double delta = B->_mean - stats->_mean;
        unsigned long X_n = stats->_n + B->_n;
        //double X_mean = stats->_mean + delta*(B->_n/X_n);
        double X_mean = (stats->_n*stats->_mean + B->_n*B->_mean) / X_n;
        double X_M2 = stats->_M2 + B->_M2 + delta*delta*stats->_n*B->_n/X_n;

        stats->_n = X_n;
        stats->_mean = X_mean;
        stats->_M2 = X_M2;
        stats->_sum += B->_sum;
        stats->_min = stats->_min < B->_min ? stats->_min : B->_min;
        stats->_max = stats->_max > B->_max ? stats->_max : B->_max;
    }
}

static inline double stats_variance(const stats_t * const stats) {
    return stats->_M2/(stats->_n-1);
}

static inline double stats_stddev(const stats_t * const stats) {
    return pow(stats_variance(stats),0.5);
}

#endif /* _STATS_H_ */

