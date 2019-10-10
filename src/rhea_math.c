#include <rhea_math.h>
#include <rhea_base.h>

double
rhea_math_smin_logexp_nondim (const double x, const double y, const double p,
                              const double dim_est)
{
  double              minimum, maximum;

  if (x < y) {
    minimum = x;
    maximum = y;
  }
  else if (y < x) {
    minimum = y;
    maximum = x;
  }
  else { /* otherwise x=y */
    return x;
  }

  RHEA_ASSERT (p < 0.0);
  RHEA_ASSERT (0.0 < dim_est);
  RHEA_ASSERT ( isfinite (exp (p/dim_est*(maximum - minimum))) );

  /* return smooth minimum */
  return minimum + log (1.0 + exp (p/dim_est*(maximum - minimum))) * dim_est/p;
}

double
rhea_math_smin_logexp (const double x, const double y, const double p)
{
  return rhea_math_smin_logexp_nondim (x, y, p, SC_MAX (fabs (x), fabs (y)));
}

double
rhea_math_smax_logexp_nondim (const double x, const double y, const double p,
                              const double dim_est)
{
  double              minimum, maximum;

  if (x < y) {
    minimum = x;
    maximum = y;
  }
  else if (y < x) {
    minimum = y;
    maximum = x;
  }
  else { /* otherwise x=y */
    return x;
  }

  RHEA_ASSERT (0.0 < p);
  RHEA_ASSERT (0.0 < dim_est);
  RHEA_ASSERT ( isfinite (exp (p/dim_est*(minimum - maximum))) );

  /* return smooth maximum */
  return maximum + log (1.0 + exp (p/dim_est*(minimum - maximum))) * dim_est/p;
}

double
rhea_math_smax_logexp (const double x, const double y, const double p)
{
  return rhea_math_smax_logexp_nondim (x, y, p, SC_MAX (fabs (x), fabs (y)));
}

double
rhea_math_smin_gpm_nondim (const double x, const double y, const double p,
                           const double dim_est)
{
  RHEA_ASSERT (p < 0.0);
  RHEA_ASSERT (0.0 < dim_est);
  return dim_est * pow (pow (x/dim_est, p) + pow (y/dim_est, p), 1.0/p);
}

double
rhea_math_smin_gpm_dx_nondim (const double x, const double y, const double p,
                              const double dim_est)
{
  RHEA_ASSERT (p < 0.0);
  RHEA_ASSERT (0.0 < dim_est);
  return pow (x/dim_est, p - 1.0) *
         pow (pow (x/dim_est, p) + pow (y/dim_est, p), 1.0/p - 1.0);
}

double
rhea_math_smin_gpm (const double x, const double y, const double p)
{
  return rhea_math_smin_gpm_nondim (x, y, p, SC_MAX (fabs (x), fabs (y)));
}

double
rhea_math_smax_gpm_nondim (const double x, const double y, const double p,
                           const double dim_est)
{
  RHEA_ASSERT (0.0 < p);
  RHEA_ASSERT (0.0 < dim_est);
  return dim_est * pow (pow (x/dim_est, p) + pow (y/dim_est, p), 1.0/p);
}

double
rhea_math_smax_gpm (const double x, const double y, const double p)
{
  return rhea_math_smax_gpm_nondim (x, y, p, SC_MAX (fabs (x), fabs (y)));
}
