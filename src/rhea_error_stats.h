/** RHEA_ERROR_STATS
 *
 * Statistics for numerical errors.
 */

#ifndef RHEA_ERROR_STATS_H
#define RHEA_ERROR_STATS_H

typedef int       (*rhea_error_stats_fn_t) (double *error_mean,
                                            double *error_stddev,
                                            void *data);

int                 rhea_error_stats_default_fn (double *error_mean,
                                                 double *error_stddev,
                                                 void *data);

#endif /* RHEA_ERROR_STATS_H */
