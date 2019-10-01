#include <rhea_error_stats.h>
#include <rhea_base.h>

int
rhea_error_stats_default_fn (double *error_mean, double *error_stddev,
                             void *data)
{
  *error_mean = NAN;
  *error_stddev = NAN;
  return 0;
}
